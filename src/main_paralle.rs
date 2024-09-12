use std::time::Instant;
use std::{
    borrow::{Borrow, Cow},
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Error, Read, Seek},
    path::Path,
    sync::{Arc, Mutex},
    thread,
};



use d4::{find_tracks, ssio::D4TrackReader, Chrom};
use d4_framefile::{Directory, OpenResult};
use regex::Regex;
use clap::{Arg, App};
use rayon::prelude::*;

use rand::seq::SliceRandom;
use rand::thread_rng;

#[derive(Debug)]
struct Range {
    length: u32,
    value: u16
}

fn parse_region_spec<T: Iterator<Item = String>>(
    regions: Option<T>,
    chrom_list: &[Chrom],
) -> std::io::Result<Vec<(usize, u32, u32)>> {
    let region_pattern = Regex::new(r"^(?P<CHR>[^:]+)((:(?P<FROM>\d+)-)?(?P<TO>\d+)?)?$").unwrap();
    let chr_map: HashMap<_, _> = chrom_list
        .iter()
        .enumerate()
        .map(|(a, b)| (b.name.to_string(), a))
        .collect();

    let ret = regions.map_or_else(
        || chrom_list.iter().enumerate().map(|(id, chrom)| (id, 0, chrom.size as u32)).collect(),
        |regions| regions.filter_map(|region_spec| {
            region_pattern.captures(&region_spec).and_then(|captures| {
                let chr = captures.name("CHR").unwrap().as_str();
                let start: u32 = captures.name("FROM").map_or(0, |x| x.as_str().parse().unwrap_or(0));
                let end: u32 = captures.name("TO").map_or_else(
                    || chr_map.get(chr).map_or(!0, |&id| chrom_list[id].size as u32),
                    |x| x.as_str().parse().unwrap_or(!0)
                );
                chr_map.get(chr).map(|&chr| (chr, start, end))
            }).or_else(|| {
                eprintln!("Warning: ignore chromosome {} which is not defined in d4 file", region_spec);
                None
            })
        }).collect()
    );

    Ok(ret)
}

fn add_one_interval(
    _chr: &str,
    left: u32,
    right: u32,
    values: &[i32],
    vect: &mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64,
    depthx1: &mut u64
) {
    let length = right - left;
    *all_length += length as u64;
    
    for &value in values {
        let value_u32 = value as u32;
        *totaldepth += (length * value_u32) as u64;
        if value_u32 >= 1 {
            *depthx1 += length as u64;
        }
        vect.push(Range { length, value: value as u16 });
    }
}

fn add_some_interval<R: Read + Seek>(
    inputs: &mut [D4TrackReader<R>],
    regions: &[(usize, u32, u32)],
    vect: &mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64,
    depthx1: &mut u64
) -> u32 {
    if inputs.is_empty() {
        return 1;
    }
    
    for &(cid, begin, end) in regions {
        let chrom = inputs[0].chrom_list()[cid].name.as_str().to_string();
        let mut values = vec![0; inputs.len()];
        let mut prev_values = vec![0; inputs.len()];
        let mut views: Vec<_> = inputs
            .iter_mut()
            .map(|x| x.get_view(&chrom, begin, end).unwrap())
            .collect();

        let mut value_changed = false;
        let mut last_pos = begin;

        for pos in begin..end {
            for (input_id, input) in views.iter_mut().enumerate() {
                if let Ok((reported_pos, value)) = input.next().unwrap() {
                    assert_eq!(reported_pos, pos);
                    if values[input_id] != value {
                        if !value_changed {
                            prev_values.clone_from(&values);
                            value_changed = true;
                        }
                        values[input_id] = value;
                    }
                }
            }
            if value_changed {
                add_one_interval(&chrom, last_pos, pos, &prev_values, vect, all_length, totaldepth, depthx1);
                last_pos = pos;
                value_changed = false;
            }
        }
        if last_pos != end {
            add_one_interval(&chrom, last_pos, end, &prev_values, vect, all_length, totaldepth, depthx1);
        }
    }
    0
}

fn add_ranges<R: Read + Seek, I: Iterator<Item = String>>(  
    mut reader: R,
    pattern: Regex,
    track: Option<&str>,
    regions: Option<I>,
    vect: &mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64, 
    depthx1: &mut u64
) -> u32 {
    let mut path_buf = vec![];
    if let Some(track_path) = track {
        path_buf.push(track_path.into());
    } else {
        find_tracks(
            &mut reader,
            |path| {
                let stem = path.and_then(|p| p.file_name()).map_or(Cow::Borrowed(""), |x| x.to_string_lossy());
                pattern.is_match(stem.borrow())
            },
            &mut path_buf,
        );
    }
    
    let file_root = Directory::open_root(reader, 8).unwrap();
    let readers: Vec<_> = path_buf.iter()
        .filter_map(|path| {
            match file_root.open(path).unwrap() {
                OpenResult::SubDir(track_root) => Some(D4TrackReader::from_track_root(track_root).unwrap()),
                _ => None,
            }
        })
        .collect();
    
    if readers.is_empty() {
        return 0;
    }
    
    let regions = parse_region_spec(regions, readers[0].chrom_list()).unwrap();
    add_some_interval(&mut readers.into_iter().collect::<Vec<_>>(), &regions, vect, all_length, totaldepth, depthx1)
}

fn read_regions(region_file: &str) -> Vec<String> {
    BufReader::new(File::open(region_file).unwrap())
        .lines()
        .filter_map(|line| {
            let line = line.unwrap();
            if line.starts_with('#') {
                return None;
            }
            let mut splitted = line.trim().split('\t');
            match (splitted.next(), splitted.next(), splitted.next()) {
                (Some(chr), Some(beg), Some(end)) => {
                    beg.parse::<u32>().ok().and_then(|begin| {
                        end.parse::<u32>().ok().map(|end| {
                            format!("{}:{}-{}", chr, begin, end)
                        })
                    })
                },
                _ => None
            }
        })
        .collect()
}

fn quick_select(interval_depth: &[Range], k: u64) -> u16 {
    if interval_depth.len() == 1 {
        return interval_depth[0].value;
    }

    let mut rng = thread_rng();
    let pivot = *interval_depth.choose(&mut rng).unwrap();
    
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut equal = Vec::new();

    for range in interval_depth {
        if range.value < pivot.value {
            left.push(range);
        } else if range.value > pivot.value {
            right.push(range);
        } else {
            equal.push(range);
        }
    }

    let left_sum: u64 = left.iter().map(|r| r.length as u64).sum();
    let equal_sum: u64 = equal.iter().map(|r| r.length as u64).sum();

    if k < left_sum {
        quick_select(&left, k)
    } else if k < left_sum + equal_sum {
        pivot.value
    } else {
        quick_select(&right, k - left_sum - equal_sum)
    }
}


fn cov_stat(interval_depth: &[Range], all_length: &mut u64, mean: &f64, depthx1: &u64) {
    let all_length_f64 = *all_length as f64;

    let stats = interval_depth.par_iter().map(|range| {
        let value = range.value as f64;
        let length = range.length as u64;
        let dis = value - mean;
        (
            if value >= 1.0 { length } else { 0 },
            if value >= 10.0 { length } else { 0 },
            if value >= 20.0 { length } else { 0 },
            if value >= 30.0 { length } else { 0 },
            if value >= 50.0 { length } else { 0 },
            dis * dis * (length as f64),
            if value > mean * 0.2 { length } else { 0 },
        )
    }).reduce(|| (0, 0, 0, 0, 0, 0.0, 0),
               |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3, a.4 + b.4, a.5 + b.5, a.6 + b.6));

    let (x1, x10, x20, x30, x50, variance_sum, large_avg_20) = stats;

    let variance = variance_sum / all_length_f64;
    let std_deviation = variance.sqrt();
    let x1_cov = (x1 as f64) / all_length_f64 * 100.0;
    let x10_cov = (x10 as f64) / all_length_f64 * 100.0;
    let x20_cov = (x20 as f64) / all_length_f64 * 100.0;
    let x30_cov = (x30 as f64) / all_length_f64 * 100.0;
    let x50_cov = (x50 as f64) / all_length_f64 * 100.0;
    let cv = std_deviation / mean;
    let q20_cov = (large_avg_20 as f64) / all_length_f64 * 100.0;

    let idx_q20 = *depthx1 / 5;
    let q20 = quick_select(interval_depth, idx_q20);

    // let q20 = interval_depth.iter()
    //     .scan(0, |acc, range| {
    //         *acc += range.length as u64;
    //         Some((*acc, range.value))
    //     })
    //     .find(|(acc, _)| *acc >= idx_q20)
    //     .map(|(_, value)| value)
    //     .unwrap_or(0);

    let fold80 = (q20 as f64) / mean;

    println!("TotalBases\tCovBases\tCovRatio\tAve_Depth(X)\tDepth>=1X\tDepth>=10X\tDepth>=20X\tDepth>=30X\tDepth>=50X\tFold80\tCV\t>=20%X");
    println!("{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}", 
             all_length, x1, x1_cov, mean, x1_cov, x10_cov, x20_cov, x30_cov, x50_cov, fold80, cv, q20_cov);
}

fn main() {
    let start_time = Instant::now();

    let args = App::new("BamCov")
        .version("0.1")
        .author("liuqingshan")
        .about("bases covered calculation for WGS, exome, or targeted sequencing from mosdepth d4 file")
        .arg(Arg::with_name("d4-file")
            .short("d")
            .long("d4-format")
            .takes_value(true)
            .required(true)
            .help("d4file https://github.com/38/d4-format"))
        .arg(Arg::with_name("region-file")
            .long("region")
            .short("r")
            .takes_value(true)
            .required(true)
            .help("input is in bed file format"))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .takes_value(true)
            .default_value("4")
            .help("Number of threads to use"))
        .get_matches();

    let input_filename = args.value_of("d4-file").unwrap();
    let region_file = args.value_of("region-file").unwrap();
    let thread_count: usize = args.value_of("threads").unwrap().parse().expect("Invalid thread count");

    let regions = Arc::new(read_regions(region_file));
    let chunk_size = (regions.len() + thread_count - 1) / thread_count;

    let vect = Arc::new(Mutex::new(Vec::new()));
    let all_length = Arc::new(Mutex::new(0u64));
    let totaldepth = Arc::new(Mutex::new(0u64));
    let depthx1 = Arc::new(Mutex::new(0u64));

    let handles: Vec<_> = (0..thread_count).map(|i| {
        let start = i * chunk_size;
        let end = std::cmp::min((i + 1) * chunk_size, regions.len());
        
        let regions = Arc::clone(&regions);
        let vect = Arc::clone(&vect);
        let all_length = Arc::clone(&all_length);
        let totaldepth = Arc::clone(&totaldepth);
        let depthx1 = Arc::clone(&depthx1);
        let input_filename = input_filename.to_string();

        thread::spawn(move || {
            let mut local_vect = Vec::new();
            let mut local_all_length = 0u64;
            let mut local_totaldepth = 0u64;
            let mut local_depthx1 = 0u64;

            let mut reader = File::open(&input_filename).unwrap();
            let track_pattern = Regex::new(".*").unwrap();

            for region in &regions[start..end] {
                add_ranges(
                    &mut reader,
                    track_pattern.clone(),
                    None,
                    Some(std::iter::once(region.clone())),
                    &mut local_vect,
                    &mut local_all_length,
                    &mut local_totaldepth,
                    &mut local_depthx1,
                );
                reader.rewind().unwrap();
            }

            let mut vect = vect.lock().unwrap();
            vect.extend(local_vect);
            
            let mut all_length = all_length.lock().unwrap();
            *all_length += local_all_length;
            
            let mut totaldepth = totaldepth.lock().unwrap();
            *totaldepth += local_totaldepth;
            
            let mut depthx1 = depthx1.lock().unwrap();
            *depthx1 += local_depthx1;
        })
    }).collect();

    for handle in handles {
        handle.join().unwrap();
    }

    let vect = Arc::try_unwrap(vect).unwrap().into_inner().unwrap();
    let mut all_length = Arc::try_unwrap(all_length).unwrap().into_inner().unwrap();
    let totaldepth = Arc::try_unwrap(totaldepth).unwrap().into_inner().unwrap();
    let depthx1 = Arc::try_unwrap(depthx1).unwrap().into_inner().unwrap();

    let mean = totaldepth as f64 / all_length as f64;

    cov_stat(&vect, &mut all_length, &mean, &depthx1);

    let elapsed_time = start_time.elapsed();
    println!("Total execution time: {:?}", elapsed_time);
}
