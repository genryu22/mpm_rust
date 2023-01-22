# MPMシミュレータ

## 実行方法

```sh
cargo run --release
```

## 操作方法

| キー | 操作 |
:----:|---- 
| Q | 1ステップ計算を進める。 |
| W | 10ステップ計算を進める。 |
| R | 100ステップ計算を進める。 |
| S | 現在の情報を保存する。(tempフォルダー) |
| Space | 再生 / 一時停止 |


## スナップショット

![image](https://user-images.githubusercontent.com/22574297/213900546-a8bc24cf-bd8a-4cd1-934d-38eae89a8a71.png)

色は圧力を表している。

## 参考リポジトリ、文献

* **[Incremental MPM](https://github.com/nialltl/incremental_mpm)**: Unityによる実装。
* **[High-Performance MLS-MPM Solver with Cutting and Coupling (CPIC) (MIT License)](https://github.com/yuanming-hu/taichi_mpm)**: Taichiライブラリを利用した簡潔な実装。
* Jiang, Chenfanfu, et al. "The affine particle-in-cell method." ACM Transactions on Graphics (TOG) 34.4 (2015): 1-10.
* Hu, Yuanming, et al. "A moving least squares material point method with displacement discontinuity and two-way rigid body coupling." ACM Transactions on Graphics (TOG) 37.4 (2018): 1-14.
* https://github.com/genryu22/test_mpm

## フォルダ構造

| フォルダ | 役割 |
:----:|---- 
| mlsmpm | シミュレータ本体(lib) |
| viewer | 描画用(bin) |
