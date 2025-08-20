# GMOS フォルダ

## 概要
GMOSフォルダには、円制限三体問題（CR3BP）における準周期軌道（Quasi-Periodic Orbit）の計算を行うためのJuliaコードが含まれています。特に、GMOS（Generalized Multiple Shooting for Ordinary differential equations）アルゴリズムを用いて、周期軌道からの分岐する準周期軌道の計算を行います。

## ファイル構成

### 主要ファイル
- **main_computeQNRHOs.jl**: メインプログラム。9:2 NRHO（近距離後方ハロー軌道）周辺の準周期軌道を計算
- **gmos.jl**: GMOSアルゴリズムの中核となる機能を提供
- **cr3bp.jl**: 円制限三体問題の方程式とベクトル場を定義
- **po_ms.jl**: 周期軌道のマルチプルシューティング法による計算機能
- **qpo_general.jl**: 一般的な準周期軌道計算機能
- **pseudoarclength.jl**: 擬似弧長継続法の実装
- **genfuncs.jl**: 一般的なユーティリティ関数

### データファイル
- **NRHO92.mat**: 9:2 NRHO軌道の初期条件データ（MATLAB形式）
- **orbits.mat**: 計算結果の軌道データ（MATLAB形式）

## 機能

### 主要機能
1. **周期軌道の計算**: マルチプルシューティング法を用いた周期軌道の精密化
2. **準周期軌道の計算**: GMOSアルゴリズムによる準周期軌道の計算
3. **継続法**: 擬似弧長継続法によるパラメータ変化に伴う軌道族の追跡
4. **多様体計算**: 軌道の安定・不安定多様体の計算

### アルゴリズム
- **GMOS**: Generalized Multiple Shooting for Ordinary differential equations
- **マルチプルシューティング法**: 境界値問題の数値解法
- **Fourier級数展開**: 準周期軌道の表現

## 使用方法
メインプログラムの実行:
```julia
julia main_computeQNRHOs.jl
```

## 依存関係
- DifferentialEquations.jl
- LinearAlgebra.jl
- MAT.jl
- Plots.jl
- その他の数値計算ライブラリ

## 理論的背景
N. Baresi et al. の論文 "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" に基づいたGMOSアルゴリズムの実装です。

## 作成者
Damennick Henry（2021-2022年）