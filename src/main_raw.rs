use std::time::Instant;
use std::{
    borrow::{Borrow, Cow},
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Error, Read, Seek},
    path::Path,
};
use d4::{
    find_tracks,
    ssio::D4TrackReader,
    Chrom,
};
use d4_framefile::{Directory, OpenResult};
use regex::Regex;
use clap::{Arg, App};
use rayon::prelude::*;

struct Range {
    length: u32,
    value: u16
}

fn parse_region_spec<T: Iterator<Item = String>>(
    regions: Option<T>,
    chrom_list: &[Chrom],
) -> std::io::Result<Vec<(usize, u32, u32)>> {
    let region_pattern = Regex::new(r"^(?P<CHR>[^:]+)((:(?P<FROM>\d+)-)?(?P<TO>\d+)?)?$").unwrap();
    let mut ret = Vec::new();

    let chr_map: HashMap<_, _> = chrom_list
        .iter()
        .enumerate()
        .map(|(a, b)| (b.name.to_string(), a))
        .collect();

    if let Some(regions) = regions {
        for region_spec in regions {
            if let Some(captures) = region_pattern.captures(&region_spec) {
                let chr = captures.name("CHR").unwrap().as_str();
                let start: u32 = captures
                    .name("FROM")
                    .map_or(0u32, |x| x.as_str().parse().unwrap_or(0));
                let end: u32 = captures
                    .name("TO")
                    .map_or_else(|| {
                        chr_map.get(chr).map_or(!0, |&id| chrom_list[id].size as u32)
                    }, |x| {
                        x.as_str().parse().unwrap_or(!0)
                    });
                if let Some(&chr) = chr_map.get(chr) {
                    ret.push((chr, start, end));
                } else {
                    eprintln!("Warning: ignore chromosome {} which is not defined in d4 file", chr);
                }
                continue;
            } else {
                return Err(Error::new(std::io::ErrorKind::Other, "Invalid region spec"));
            }
        }
    } else {
        for (id, chrom) in chrom_list.iter().enumerate() {
            ret.push((id, 0, chrom.size as u32));
        }
    }

    ret.sort_unstable();

    Ok(ret)
}


fn add_one_interval(
    chr: & str,
    left: u32,
    right: u32,
    values: &[i32],
    vect:  &mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64,
    depthx1: &mut u64
) {
    // println!("{:?}", values);
    for value in values.iter()
    {
        let vaue_t = *value as u32;
        let length = right - left;
        *all_length += length as u64;
        let interval_depths = (length*vaue_t) as u64;
        *totaldepth += interval_depths;
        if vaue_t >=1 {*depthx1 +=length as u64; };
        vect.push(Range{length:length,  value: value.clone() as u16 });
    }

}


fn add_some_interval <R: Read + Seek>( 
    inputs: &mut [D4TrackReader<R>],
    regions: &[(usize, u32, u32)],
    vect: & mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64,
    depthx1: &mut u64
    ) -> u32
    {
        if inputs.is_empty() {
            return 1;
        }
        for &(cid, begin, end) in regions {
        let chrom_list = inputs[0].chrom_list();
        let chrom = chrom_list[cid].name.as_str().to_string();

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
                let Ok((reported_pos, value)) = input.next().unwrap()  else { todo!() };
                assert_eq!(reported_pos, pos);
                if values[input_id] != value {
                    if !value_changed {
                        prev_values.clone_from(&values);
                        value_changed = true;
                    }
                    values[input_id] = value;
                }
            }
            if value_changed {
                add_one_interval(&chrom,last_pos,pos, prev_values.as_slice(), vect, all_length, totaldepth, depthx1);
                last_pos = pos;
                value_changed = false;
            }
        }
        if last_pos != end {
            add_one_interval(&chrom, last_pos, end, prev_values.as_slice(), vect, all_length, totaldepth, depthx1);
        }
    }
    return  0;
    }

fn add_ranges <R: Read + Seek, I: Iterator<Item = String>> (  
    mut reader: R,
    pattern: Regex,
    track: Option<&str>,
    regions: Option<I>,
    vect: & mut Vec<Range>,
    all_length: &mut u64,
    totaldepth: &mut u64, 
    depthx1: &mut u64)  -> u32 {
        let first =false;
        let mut path_buf = vec![];
        let mut first_found = false;
        if let Some(track_path) = track {
            path_buf.push(track_path.into());
        } else {
            find_tracks(
                &mut reader,
                |path| {
                    let stem = path
                        .map(|what: &Path| {
                            what.file_name()
                                .map(|x| x.to_string_lossy())
                                .unwrap_or_else(|| Cow::<str>::Borrowed(""))
                        })
                        .unwrap_or_default();
                    if pattern.is_match(stem.borrow()) {
                        if first && first_found {
                            false
                        } else {
                            first_found = true;
                            true
                        }
                    } else {
                        false
                    }
                },
                &mut path_buf,
            );
        }
        let file_root = Directory::open_root(reader, 8).unwrap();
        let mut readers = vec![];
        for path in path_buf.iter() {
        let  track_root = match file_root.open(path).unwrap() {
                OpenResult::SubDir(track_root) => track_root,
                _ => {
                    return 0;
                }
            };
            let reader = D4TrackReader::from_track_root(track_root).unwrap();
            readers.push(reader);
        }
        let regions = parse_region_spec(regions, readers[0].chrom_list()).unwrap();
        add_some_interval(& mut readers, &regions, vect, all_length, totaldepth, depthx1);
        return  0;

    }

    
fn all_d4_vec(input_filename: & str,region_file:&str, vect: & mut Vec<Range>, all_length: &mut u64,totaldepth: &mut u64, depthx1: &mut u64 ) -> u32 {
    let regions =  {
        let mut file = BufReader::new(File::open(region_file).unwrap());
        let mut buf = String::new();
        let mut region_list = Vec::new();
        while file.read_line(&mut buf).unwrap() > 0 {
            if &buf[..1] == "#" {
                continue;
            }
            let mut splitted = buf.trim().split('\t');
            let (raw_chr, raw_beg, raw_end) = (splitted.next(), splitted.next(), splitted.next());
            if raw_chr.is_some() && raw_beg.is_some() && raw_end.is_some() {
                if let Ok(begin) = raw_beg.unwrap().parse::<u32>() {
                    if let Ok(end) = raw_end.unwrap().parse::<u32>() {
                        region_list.push(format!("{}:{}-{}", raw_chr.unwrap(), begin, end));
                    }
                }
                buf.clear();
                continue;
            }
            panic!("Invalid bed file");
        }
        Some(region_list.into_iter())
    };
    let track_path = None;
    let track_pattern =  regex::Regex::new(".*").unwrap();
    let reader = File::open(input_filename).unwrap();
    add_ranges( reader, track_pattern, track_path, regions,vect, all_length, totaldepth, depthx1);
    return 0;
}



fn cov_stat(interval_depth: &Vec<Range>, all_length: &mut u64, mean: &f64, depthx1: &mut u64)
{
    let mut idx :u64 = 0;
    let mut  x1:u64 = 0;
    let mut  x10:u64 = 0;
    let mut  x20:u64 = 0;
    let mut  x30:u64 = 0;
    let mut  x50:u64 = 0;

    let idx_q20: u64 = (*depthx1)*2 /10;
    let mut variance:f64= 0.0;
    let mut  q20 :u16 = 0;
    let mut jixu: bool = true;

    // 计算 q20% 
    let mut large_avg_20: u64 = 0;
    let mut less_avg_20: u64 = 0;
    let avg_p20 = (*mean as f64) * 0.2;


    for range_value in interval_depth.iter()
    {
        let value = range_value.value;
        let length = range_value.length as u64;
        if value >=1 {x1+=length; idx += length;}
        if value >=10 {x10+=length; }
        if value >=20 {x20+=length;}
        if value >=30 {x30+=length;}
        if value >=50 {x50+=length;}
        let dis = (value as f64) - mean;
        let pow2 = dis*dis ;
        variance +=  pow2 * (length as f64);
        if jixu {
            if idx >= idx_q20{
                q20 = value; 
                jixu =false;}
        }

        if (value as f64)  > avg_p20 { large_avg_20 +=length;}
        else {less_avg_20 +=length;}
    }


    let all_lenth_f64 = all_length.clone() as f64;
    variance /= all_lenth_f64;
    let  std_deviation = variance.sqrt();
    let x1_cov = (x1 as f64)  / all_lenth_f64 * 100.0;
    let x10_cov = (x10 as f64)  / all_lenth_f64 * 100.0 ;
    let x20_cov = (x20 as f64)  / all_lenth_f64* 100.0;
    let x30_cov = (x30 as f64)  / all_lenth_f64* 100.0;
    let x50_cov = (x50 as f64)  / all_lenth_f64* 100.0;
    let fold80 = (q20 as f64) / mean;
    let cv =  std_deviation / mean;

    let q20_cov = (large_avg_20 as f64 ) / ((large_avg_20 + less_avg_20 ) as f64) *100.0;
    println!("TotalBases\tCovBases\tCovRatio\tAve_Depth(X)\tDepth>=1X\tDepth>=10X\tDepth>=20X\tDepth>=30X\tDepth>=50X\tFold80\tCV\t>=20%X");
    println!("{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}", all_length, x1, x1_cov, mean, x1_cov, x10_cov, x20_cov, x30_cov, x50_cov, fold80, cv, q20_cov);
}


fn main() {
    let args = App::new("BamCov")
        .version("0.1")
        .author("liuqingshan")
        .about("bases covered calculation for WGS, exome, or targeted sequencing from mosdepth  d4 file ")
        .arg(Arg::with_name("d4-file")
            .short("d")
            .long("d4-format")
            .takes_value(true)
            .help("d4file https://github.com/38/d4-format "))
        .arg(Arg::with_name("region-file")
            .long("region")
            .short("r")
            .takes_value(true)
            .help("input is in bed file format"))
        .arg(Arg::with_name("version")
             .long("version")
             .short("V")
             .takes_value(false)
             .help("print version info"))
        .get_matches();

        if args.is_present("version") {
            println!("bamcov v{}", "0.1.0");
            std::process::exit(0)
        };

        let input_filename :&str= args
        .value_of("d4-file")
        .unwrap()
        .trim();
        let region_file :&str = args
        .value_of("region-file")
        .unwrap()
        .trim();

    let start = Instant::now();
    let mut all_length:u64 = 0;
    let mut totaldepth:u64 = 0;
    let mut depthx1: u64 = 0;
    let mut interval_depth: Vec<Range> = Vec::new();
    all_d4_vec(input_filename, region_file, &mut interval_depth, &mut all_length, &mut totaldepth, &mut depthx1);
    let start2 = Instant::now();
    interval_depth.par_sort_unstable_by(|a, b| b.value.cmp(&a.value));
    let duration2 = start2.elapsed();
    //println!("排序用时: {:?}", duration2);
    let mean: f64 = (totaldepth as f64) / (all_length as f64);
    cov_stat(& interval_depth, & mut all_length, & mean, &mut depthx1); 
    let duration = start.elapsed();
    println!("总运行时间: {:?}", duration);
    //println!("finished");
}