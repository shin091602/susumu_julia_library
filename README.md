# Orbital_Design_Library

## 概要
このコードは円制限３体問題の軌道設計を行うためのJuliaコードである。

## 動作確認
### 基本情報
- OS : Windows 11 Pro
- CPU : AMD Ryzen Threadripper PRO 5975WX 32-Cores        3.60 GHz
- RAM : 256 GB
- GPU : NVIDIA RTX A6000

### 実行環境に関する詳細
- juliaのバージョン：julia 1.11.3+0
- VScodeのバージョン：

### 動作内容に対する負荷情報
- シミュレーション規模

- 計算時間

- CPU負荷

- GPU負荷

- 消費メモリ

- 

## 計算手法はフォルダ内の以下を参照
Dynamical_Systems_Theory_in_CR3BP.pdf
※このpdfは現在、MATLABのコードで説明されたものをアップロードしています。修正し次第、Juliaのものに変更します。
※サイズの影響でpdfファイルしかアップロードできませんでした。

## 環境構築　(Windowsの場合)
1. Juliaのインストール
Microsoft StoreからJuliaをインストールする。または、[Julia公式ウェブサイト](https://julialang.org/downloads/)からJuliaをダウンロードする。
2. ライブラリのインストール
インストールしたJuliaアプリを実行し、REPLを開く。
3. "]"を入力し、パッケージモードに入り、以下を実行する。

- CSVのインストール（CSVファイルの読み書きを行うライブラリ）
  
```add CSV```

- DataFramesのインストール（データを表形式で操作するためのライブラリ）

```add DataFrames```
- DifferentialEquationsのインストール（常微分方程式などを数値的に解くライブラリ）

```add DifferentialEquations```

- StaticArraysのインストール（高速な固定サイズ配列を提供するライブラリ）

```add StaticArrays```

- LinearAlgebraのインストール（線形代数計算用の標準ライブラリ）

```add LinearAlgebra```

- GLMakieのインストール（インタラクティブで高速な可視化を行うためのライブラリ）

```add GLMakie```

- CairoMakieのインストール（高品質なベクター出力を含むMakieバックエンド）

```add CairoMakie```

- Polynomialsのインストール（多項式の定義と演算をサポート）

```add Polynomials```

- Printfのインストール（フォーマット付き文字列出力）

```add Printf```

- Interpolationsのインストール（補間処理用のライブラリ）

```add Interpolations```

- ColorSchemesのインストール（色スキームの管理と可視化用）

```add ColorSchemes```

- Colorsのインストール（色の操作と変換を行うライブラリ）

```add Colors```
