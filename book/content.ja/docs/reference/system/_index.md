+++
title = "System"
weight = 500
+++

# System

系に含まれる粒子や、系の境界条件、系全体のパラメータなどを設定します。

複数のシミュレーションボックスを持つシミュレーションのためにテーブルの配列として
定義されていますが、通常のシミュレーションでは一つしか定義しません。

## 例

```toml
[[systems]]
attributes.temperature = 300.0
boundary_shape.lower = [  0.0,   0.0,   0.0] # lower limit of the boundary
boundary_shape.upper = [100.0, 100.0, 100.0] # upper limit of the boundary
particles = [
    {m= 1.0, pos=[ 1.000, 2.000, -1.000], vel=[ 0.100,-0.200, 0.300], name="CA", group="A"},
    # ...
]
```

## 入力

- `boundary_shape`: テーブル型
  - 境界の具体的な大きさを設定します。
  - `Unlimited`な境界条件の場合、省略可能です。
- `attributes`: テーブル型
  - 系全体に適用されるパラメータを設定します。
- `particles`: テーブルの配列型
  - 系に含まれる粒子のパラメータと、初期条件を設定します。
  - `m`または`mass`: 浮動小数点数型
  - `pos`または`position`: 浮動小数点数の配列型（要素数は3）
  - `vel`または`velocity`: 浮動小数点数の配列型（要素数は3、省略可能）
  - `name`: 文字列型（省略可能）
  - `group`: 文字列型（省略可能）

### `"Periodic"`境界条件での`boundary_shape`

境界条件が`PeriodicCuboid`の場合、直方体の各座標の最大と最小の座標を取る頂点を設定します。

- `lower`: 浮動小数点の配列型
  - 各座標で最小の座標を取る頂点の位置を指定します。
- `upper`: 浮動小数点の配列型
  - 各座標で最大の座標を取る頂点の位置を指定します。

```toml
boundary_shape = {lower = [0.0, 0.0, 0.0], upper = [10.0, 10.0, 10.0]}
```

または、

```toml
boundary_shape.lower = [  0.0,   0.0,   0.0]
boundary_shape.upper = [100.0, 100.0, 100.0]
```


### `attribute`

系全体に適用されるパラメータを指定します。
力場やシミュレータによって必要なパラメータが変化します。
現在、要求される変数は以下のどれかです。

- `temperature`: 浮動小数点数型
  - 系の参照温度を指定します。温度一定のシミュレーションを行う際は指定しておく必要があります。
- `ionic_strength`: 浮動小数点数型
  - 系のイオン強度を指定します。単位は[mol/L]です。
  - イオン粒子が陽に含まれない系で静電相互作用を使用する際に必要になることがあります。

### `particles`

系の粒子のパラメータと、初期条件を指定します。
一つの粒子に対して一つのテーブルが割り当てられ、質量、位置、速度、名前、グループ名を指定します。

```toml
particles = [
    {mass= 1.0, position=[ 1.000, 2.000, -1.000], velocity=[ 0.100,-0.200, 0.300], name="CA", group="A"},
    # 短縮形
    {m= 1.0, pos=[ 1.000, 2.000, -1.000], vel=[ 0.100,-0.200, 0.300], name="CA", group="A"},
    # ...
]
```

初期速度を指定しなかった場合、Maxwell-Boltzmann分布に従って指定した温度での
初期速度が`Simulator`で設定した`seed`を使ってランダムに生成されます。

そのため、省略した場合は`system.attributes.temperature`の値を設定していることが
要求されます。

```toml
particles = [
    {m= 1.0, pos=[ 1.000, 2.000, -1.000], name="CA", group="A"},
    # ...
]
```

## 他のシミュレーションの最終構造をインポートする

Mjolnirは、ファイル出力時に[MsgPack](https://msgpack.org/)で`System`の全状態をエンコードして出力します。
このファイルを指定することで、前回のシミュレーションの続きからの再開や、途中で計算が止まってしまった時にすぐにやり直すことができます。

```toml
[[systems]]
file_name = "restart.msg"
```

