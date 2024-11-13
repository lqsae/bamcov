# BamCov

BamCov 是一个基于 mosdepth d4 文件计算 WGS、外显子组或靶向测序的碱基覆盖的 Rust 程序。

## 功能特点

- 支持多线程并行处理，提高计算效率
- 使用 d4 文件格式，兼容 mosdepth 输出
- 可以计算多种覆盖度统计指标
- 支持自定义区域文件进行分析
- 支持自定义深度阈值统计

## 安装

确保您的系统已安装 Rust 环境。克隆此仓库并使用 Cargo 进行编译：

```bash
git clone [repository_url]
cd bamcov
cargo build --release
```

## 使用方法

1. 编译程序：
   ```bash
   cargo build --release
   ```

2. 运行程序：
   ```bash
   cargo run --release -- -h
   ```

3. 示例：
   ```bash
   # 使用默认深度阈值
   ./bamcov -d sample.d4 -r exome.bed

   # 使用自定义深度阈值
   ./bamcov -d sample.d4 -r exome.bed -t 10,20,30,40,50
   ```

## 参数说明

- `-d, --d4-format`: 指定输入的 d4 文件（必需）
- `-r, --region`: 指定输入的 bed 格式区域文件（必需）
- `-t, --threshold`: 指定自定义深度阈值，多个值用逗号分隔（可选）
  - 示例：`-t 10,20,30,40,50`
  - 不指定时使用默认阈值：1X,10X,20X,30X,50X

## 输入文件格式

1. d4 文件：使用 mosdepth 生成的 d4 格式深度文件
2. bed 文件：标准 bed 格式，指定要分析的基因组区域

## 输出结果说明

程序会输出以下统计信息：

- TotalBases: 总碱基数
- CovBases: 覆盖的碱基数（≥1X的碱基数）
- CovRatio: 覆盖率（占总碱基的百分比）
- Ave_Depth(X): 平均测序深度
- Depth>=NX: 深度大于等于 N 的碱基所占比例
  - 默认显示 1X,10X,20X,30X,50X
  - 使用 -t 参数可自定义阈值
- Fold80: 测序均一性参数（80%碱基深度/平均深度）
- CV: 变异系数（深度标准差/平均深度）
- >=20%X: 深度大于等于平均深度 20% 的碱基比例

输出示例：
```
# 默认输出格式
TotalBases  CovBases    CovRatio  Ave_Depth  Depth>=1X  Depth>=10X  Depth>=20X  Depth>=30X  Depth>=50X  Fold80  CV      >=20%X
150000000   149000000  99.333    30.456     99.333     95.234      90.123      85.678      70.456      1.234   0.345   85.678

# 使用自定义阈值 -t 10,20,40
TotalBases  CovBases    CovRatio  Ave_Depth  Depth>=10X  Depth>=20X  Depth>=40X  Fold80  CV      >=20%X
150000000   149000000  99.333    30.456     95.234      90.123      75.678      1.234   0.345   85.678
```

## 性能优化

- 使用多线程并行处理区域
- 采用快速排序算法进行深度排序
- 使用 Rayon 库进行并行计算
- 高效的内存管理和数据结构

## 注意事项

- 确保输入的 d4 文件和区域文件格式正确
- 自定义阈值必须是正整数，多个值用逗号分隔
- 程序会输出处理时间，方便性能分析
- 大文件处理时建议使用 release 模式运行

## 贡献

欢迎提交 issues 和 pull requests 来改进这个项目。

## 许可证

[MIT License](LICENSE)

## 作者

刘青山

## 版本历史

- v0.1.0 (2024-11)
  - 初始版本发布
  - 支持基本的覆盖度统计
  - 添加自定义深度阈值功能