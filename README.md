# BamCov

BamCov 是一个基于 mosdepth d4 文件计算 WGS、外显子组或靶向测序的碱基覆盖的 Rust 程序。

## 功能特点

- 支持多线程并行处理，提高计算效率
- 使用 d4 文件格式，兼容 mosdepth 输出
- 可以计算多种覆盖度统计指标
- 支持自定义区域文件进行分析

## 安装

确保您的系统已安装 Rust 环境。克隆此仓库并使用 Cargo 进行编译：


## 使用方法

1. 编译程序：
   ```bash
   cargo build --release
   ```

2. 运行程序：
   ```bash
   cargo run --release -- -h

3. 示例
   ```bash
   BamCov -d sample.d4 -r exome.bed
   ```
4. 输出示例
   ```
   TotalBases CovBases CovRatio Ave_Depth Depth>=1X Depth>=10X Depth>=20X Depth>=30X Depth>=50X Fold80 CV >=20%X
   150000000 149000000 99.333 30.456 99.333 95.234 90.123 85.678 70.456 1.234 0.345 85.678
   ```

参数说明：
- `-d, --d4-format`: 指定输入的 d4 文件（必需）
- `-r, --region`: 指定输入的 bed 格式区域文件（必需）
- `-t, --threads`: 指定使用的线程数（可选，默认为 4）

## 输入文件格式
1. d4 文件：使用 mosdepth 生成的 d4 格式深度文件
2. bed 文件：标准 bed 格式，指定要分析的基因组区域

## 输出结果

程序会输出以下统计信息：

- TotalBases: 总碱基数
- CovBases: 覆盖的碱基数
- CovRatio: 覆盖率
- Ave_Depth(X): 平均深度
- Depth>=1X: 深度大于等于 1X 的比例
- Depth>=10X: 深度大于等于 10X 的比例
- Depth>=20X: 深度大于等于 20X 的比例
- Depth>=30X: 深度大于等于 30X 的比例
- Depth>=50X: 深度大于等于 50X 的比例
- Fold80: 测序均一性参数
- CV: 变异系数
- 20%X: 深度大于等于平均深度 20% 的比例


## 依赖库

- clap: 命令行参数解析
- d4: d4 文件格式处理
- rayon: 并行计算
- regex: 正则表达式处理
- rand: 随机数生成

## 性能优化

- 使用多线程并行处理区域
- 采用快速选择算法计算分位数
- 使用 Rayon 库进行并行计算

## 注意事项

- 确保输入的 d4 文件和区域文件格式正确
- 程序会输出各个步骤的运行时间，方便性能分析

## 贡献

欢迎提交 issues 和 pull requests 来改进这个项目。

## 许可证

[在此处添加许可证信息]