# Nickelate RKKY

该项目计算双层镍酸盐紧束缚模型的静态磁化率及其对应的 RKKY
相互作用，并保留后续加入基态和激发态计算的模块位置。

## 目录

```text
src/
  NickelateRKKY.jl                 主加载模块
  models/                          哈密顿量与模型参数
  susceptibility/                  磁化率核心算法
  ground_state/                    后续基态模块
  excitations/                     激发态模块与实验代码
scripts/
  susceptibility/                  数值计算入口
  plotting/                        绘图入口
data/
  model/                            模型输入表
  rkky/                             表格数据
  results/                          JLD2 数值结果
results/
  susceptibility/                  磁化率和能带图片
  rkky/                             RKKY 图片
paper/                              TeX、PDF 与演示文稿
references/                         参考论文及提取文本
tests/                              Julia 测试
```

## 环境

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## 常用命令

计算总 RKKY 相互作用：

```powershell
julia --project=. scripts/susceptibility/run_total_rkky.jl
```

绘制高对称路径磁化率：

```powershell
julia --project=. scripts/plotting/plot_chi_highsym.jl
```

运行测试：

```powershell
Get-ChildItem tests/test_*.jl | ForEach-Object {
    julia --project=. $_.FullName
}
```

## 模块边界

`src/models/` 只负责模型和能带哈密顿量；`src/susceptibility/`
负责费米分布、粒子数、化学势、动量空间磁化率和实空间变换。
未来基态和激发态算法分别放入 `src/ground_state/` 与
`src/excitations/`，参数扫描入口放入对应的 `scripts/` 子目录。

Sunny 当前未列入 `Project.toml`，因此
`src/excitations/spin_excitation.jl` 不由主模块自动加载。
# Nickelate RKKY

## Bilayer x2 RKKY scan

The production scan extracts `J0`, `J1`, `J2`, `J3`, and `J1p` in eV from
the layer-resolved x2-y2 susceptibility:

```powershell
$env:JULIA_NUM_THREADS = "20"
julia --project=. scripts/susceptibility/run_x2_bilayer_rkky.jl
```

The default `Nk=Nq=100` calculation contains 61 interaction values and is a
long-running job. It writes an atomically replaceable, resumable result to
`data/results/rkky_x2_bilayer_U0-6_Nk100.jld2`.
