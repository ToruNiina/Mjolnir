# HarmonicPotential

調和振動子ポテンシャルは以下のような形のポテンシャルです。

{% math %}
U(v) = k(v-v_0)^2
{% endmath %}

シンプルで様々な用途で用いられるポテンシャルです。
結合長、結合角、二面角のどれにも用いられることがあります。

## 例

```toml
[[forcefields.local]]
# ここにはBondLengthなどによって要求されるパラメータが入ります。
# ...
parameters = [
    {indices = [0, 1], v0 = 1.0, k = 100.0},
    # 必要に応じて続きます...
]
```

## 入力

上の例では対応する粒子の番号が2つ書かれていますが、この個数は相互作用の種類によって変動します。

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点数型
  - 最安定点を指定します。
- `indices`: 整数の配列型
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
  - 指定する粒子の数は相互作用によって異なります([`BondLength`](BondLengthInteraction.md)なら2個、[`BondAngle`](BondAngleInteraction.md)なら3個などです)。