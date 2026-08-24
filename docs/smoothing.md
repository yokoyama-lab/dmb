# Smoothing discrete Morse-Bott theory: 連続 Morse-Bott 関数の離散化計画

**問い**: 「全ての連続の Morse 関数を離散 Morse 関数で近似できる」
（Gallais 2010・Benedetti 2016）を Morse-Bott に拡張し，
「全ての連続の Morse-Bott 関数を離散 Morse-Bott 関数で近似できる」を示せるか。

**見込み: 示せる。しかも Morse 版より易しくなる可能性が高い。**
理由は下の「なぜ易しくなるか」。この文書は予想の定式化・証明戦略・
既に機械で確認できている証拠（層 A）・残る難所をまとめた作業計画である。

主張の層は [results.md](results.md) と同じ:
**層 A** = このリポジトリのコードで機械確認済み，**層 B** = 紙の証明が書ける，
**要確認** = 文献の現物での裏取りが済んでいない。

---

## 1. 先行研究（何が示されているか）

**Gallais, *Combinatorial realization of the Thom-Smale complex via discrete
Morse theory*, Ann. Scuola Norm. Sup. Pisa (5) IX (2010) 229–252.**
Theorem 3.1（本文 PDF から逐語確認済み）: $M$ を滑らかな閉有向リーマン多様体，
$f$ を generic な Morse 関数とすると，$C^1$-三角形分割 $T$ と組合せ的 Morse
ベクトル場 $V$ が存在して，(1) 臨界セルと臨界点が全単射，(2) $V$-path と
積分曲線が全単射，(3) 滑らかな Thom-Smale 複体と組合せ的 Thom-Smale 複体が
同型になる。証明はコボルディズム分解による: 核心は Theorem 3.13 =
**臨界点 1 個のコボルディズムに，2 段の simplicial collapse
（$T \searrow T_p^s \cup T_0$ と $T_p^s - \sigma_p \searrow T_p^s \cap T_0$）を
持つ三角形分割が存在する**こと。

**Benedetti, *Smoothing discrete Morse theory*, Ann. SNS Pisa (2016);
[arXiv:1212.0885](https://arxiv.org/abs/1212.0885)**（abstract を逐語確認済み）:
PL ハンドル理論と離散 Morse 理論は同値。**滑らかな Morse ベクトルは全次元で
ある PL 三角形分割上の離散 Morse ベクトルになる**（Gallais の改良）。
逆（離散 → 滑らか）は次元 7 以下。次元制限の源泉は PL/DIFF の差と，
**正則コボルディズムを組合せ的に collapse させる**ために必要な細分の制御にある。

どちらも「近似」の内容は **臨界構造の実現**（臨界点 ↔ 臨界セル）であって，
関数値の $\varepsilon$-近似は付随的（三角形分割を細かくすれば PL 補間が
$C^0$ 近似になる）。本計画でも同じ立場を取る。

## 2. 予想の定式化

**予想 (Smoothing discrete Morse-Bott theory).**
$M$ を滑らかな閉多様体，$f\colon M \to \mathbb{R}$ を Morse-Bott 関数，
臨界部分多様体を $C_1, \dots, C_m$（指数 $\lambda_i$）とする。このとき
$M$ の $C^1$-三角形分割 $T$ と，その上の離散 Morse-Bott 関数
$g$（[arXiv:2511.07864](https://arxiv.org/abs/2511.07864) Definition 18）が存在して:

1. $g$ の（$P_t \ne 0$ の）collection は $C_i$ と全単射に対応する。
2. $C_i$ に対応する collection $K_i$ の制限ホモロジーは
   $P_t(K_i) = t^{\lambda_i} P_t(C_i; F)$（負法束 $\nu^-_i$ が
   $F$-向き付け可能なとき。そうでないときは局所係数 / $F = \mathbb{F}_2$）。
3. したがって $\sum_C P_t(C)$ は滑らかな Morse-Bott 多項式に一致し，
   離散 Morse-Bott 不等式（同 Theorem 4.12）は滑らかな不等式と同じ剰余 $R(t)$ を持つ。
4. （精密化）$f$ が有限群 $G$ の作用で不変で $T$ が $G$-不変に取れるなら，
   $g$ も $G$-不変に取れる。PL 補間 $|g|$ は $f$ の $C^0$-近似にできる。

Morse-Bott の退化極限（$f$ が Morse）では Gallais の定理の帰結
（臨界セルの個数の実現）を回復する。逆に $f$ が定数なら $g$ は定数関数でよい。

## 3. なぜ Morse 版より易しくなる見込みか

Gallais・Benedetti の最難関は **正則コボルディズム
$X \times [0,1]$ を組合せ的に collapse させる**ことである（Whitehead 流の
細分が要り，Benedetti の次元制限 dim ≤ 7 もここから来る）。

離散 Morse-Bott 理論では **collapse は要らない**。必要なのは
「正則な帯の collection が $\sum_C P_t(C)$ に寄与しない」ことだけで，
これは制限ホモロジーの消滅 = **純代数的な条件**である。鍵は次の補題:

**補題 L5（制限ホモロジー = 相対ホモロジー; ⭐ 既に Lean 化済み）.**
$g = $ `max_extension`$(h)$（頂点関数 $h$ の最大値拡張）とし，$K'$ を値 $c$ の
collection の 1 つとする。$\overline{N}$ を $K'$ のセルの閉包，
$\partial_- \overline{N} = \overline{N} \setminus K'$（値 $< c$ の face 全体）と
おくと $\partial_- \overline{N}$ は部分複体であり，$K'$ の制限鎖複体
$(C_*(K'), \partial^{K'})$ は**商鎖複体 $C_*(\overline{N}) / C_*(\partial_-\overline{N})$
そのもの**である。ゆえに
$$P_t(K') = P_t\bigl(|\overline{N}|, |\partial_-\overline{N}|; F\bigr)
\quad(\text{相対単体的ホモロジー}).$$

この同一視は，実は Collapse 定理の形式化プロジェクト
（`tetsuo-jp/discrete-morse-bott`, 研究台帳 P-dmb-01）で
**既に Lean で機械証明されている**（「$(C_*(C), \partial^C)$ が対
$(\overline{C}, \overline{C} \setminus C)$ の相対鎖複体であること」まで Lean 化済み，
2026-08-20 時点）。本計画はそれを滑らかな側と接続するだけでよい。

これが効く理由: 三角形分割を $f$ の準位に適合させておけば，
- **正則な帯**は $(X \times [0,1],\ X \times \lbrace 0 \rbrace)$ 型の対なので
  相対ホモロジーは消える — collapse 不要，細分不要，**次元制限なし**。
- **臨界帯**（$C_i$ を含む帯）は Morse-Bott 補題により
  $(D(\nu^-_i)\text{-束},\ \partial_- )$ 型の対なので，Thom 同型で
  $t^{\lambda_i} P_t(C_i)$ — これが条件 2 になる。

つまり Morse 版で一番重かった組合せ論が，Morse-Bott 版では
ホモロジー計算に置き換わる。**これが本計画の中心的な観察である。**

## 4. 証明戦略

### 戦略 B（本命・直接構成）

1. **L1（適合三角形分割）**: 正則値 $a_0 < a_1 < \dots$ を臨界値を 1 つずつ
   挟むように取り，各 $f^{-1}(a_k)$・各 $C_i$ の管状近傍が部分複体になる
   $C^1$-三角形分割を取る。正則な帯は積として，臨界帯は円板束として
   三角形分割する。**（要確認: Whitehead / Munkres の $C^1$-三角形分割が
   有限個の部分多様体に適合して取れること。Gallais §3 が臨界点の場合に
   同種の構成をしているので，そこを踏襲できる見込み）**
2. **L2（DMBF の自動性; 層 A・一般複体で機械確認済み）**: 頂点に帯の番号
   $h(v) = k$（$v \in f^{-1}[a_{k-1}, a_k]$ 側の割当）を与え，
   $g = $ `max_extension`$(h)$ とする。**どんな複体・どんな $h$ でも $g$ は
   DMBF になる**（`dmb_core.max_extension` の docstring の通り。
   $\sigma \prec \tau \Rightarrow g(\sigma) \le g(\tau)$ なので
   $U^{\mathrm{snc}} = D^{\mathrm{snc}} = 0$）。Gallais が離散勾配ベクトル場を
   手作りしていた工程が丸ごと消える。
3. **L5** で各 collection の $P_t$ を相対ホモロジーに読み替える。
4. **L3（正則帯の消滅）**: L5 の系。$(X \times I, X \times \lbrace 0\rbrace)$ の
   相対ホモロジーは 0。
5. **L4（臨界帯 = Thom 同型）**: L5 + 滑らかな Morse-Bott 補題 + Thom 同型。
   負法束が向き付け不可能なら $\mathbb{F}_2$ か局所係数（予想の但し書き）。
6. **r-連結性の管理**: collection は r-path の同値類なので，帯の連結成分ごとに
   1 つの collection になることを確かめる（max_extension では
   値の等しい codim-1 対がすべて r-辺になるので，帯の双対グラフの連結性に帰着。
   臨界値が帯ごとに 1 つになるよう L1 で分ければよい）。

### 戦略 A（保険・Morsification 経由）

$f$ を Morse 関数 $f_\varepsilon = f + \varepsilon \sum_i \rho_i f_i$（$f_i$ は
$C_i$ 上の Morse 関数）に摂動 → Gallais / Benedetti を適用して DMF を得る →
$C_i$ 由来の臨界セルの群を**同じ値に潰して** collection にまとめ直す
（Lemma 4.8 `morsify` の逆向き）。潰した後の (MB2)(MB4) と r-連結性の検証が
必要で，V-path の構造に踏み込むため戦略 B より重い。B が詰まったときの迂回路。

## 5. 機械で確認済みの証拠（層 A）

すべて `python3 tests/…` あるいは下のスニペットで再現できる。
離散側の計算はこのリポジトリの `max_extension` + `MorseBott.report()` のみ。

| モデル | 滑らかな側 | 離散側（機械確認） | 一致 |
|---|---|---|---|
| トーラスの高さ関数 | 臨界円周 2 本, $R = 0$ | $(1{+}t) + t(1{+}t)$, $R = 0$（[results.md](results.md) §2, $n_i,n_j \le 7$） | ✓ |
| 紡錘形 $S^2$（赤道が最小円周・両極が最大点） | $(1{+}t) + 2t^2$, $R = t$ | 同じ（$n = 3,4,6$ の 3 通り） | ✓ 剰余まで |
| 円筒 $S^1 \times [0,3]$ の高さ | 臨界円周は下端のみ, $R=0$ | $(1{+}t)$ + 正則帯 4 つがすべて $P_t = 0$ | ✓ L3 の証拠 |
| クラインの壺の高さ（上端の円周の法束が非向き付け） | $\mathbb{Q}$: $R = t$ ／ $\mathbb{F}_2$: $R = 0$ | 同じ（同一の $g$ で係数だけ変えて再現） | ✓ 係数体依存まで |

とくに最後の 2 行が重要である:
- 円筒は**境界つき**多様体だが正則帯の消滅（L3）はそのまま成り立つ。
- クラインの壺は**負法束が向き付け不可能な臨界円周**を持ち，滑らかな側の
  「$\mathbb{Q}$ では等号にならず $\mathbb{F}_2$ でなる」がそのまま離散側に現れる。
  予想の但し書き（向き付け条件）が実際に必要で，かつそれで十分らしいことを示す。

再現（検査 `tests/test_smoothing.py` にも同じものを入れてある）:

```python
import dmb_core as D, complexes as X
K = X.annulus(4, 4)                      # 円筒
f = D.max_extension(K, lambda v: v[1])   # 高さ。これだけで DMBF になる（L2）
D.MorseBott(K, f).report()               # R_MB = [] すなわち R = 0
```

## 6. 残る難所（正直に）

1. **L1 の文献確認（最大のリスク）**: 「$C^1$-三角形分割を有限個の部分多様体
   （準位集合・管状近傍）に適合させて取れる」の正確な出典。Whitehead 1940 /
   Munkres *Elementary Differential Topology* / Verona（成層空間）あたり。
   Gallais §3.2 が臨界点の場合の管状近傍でこれをやっているので，
   円板束版に書き直すのが実作業になる。
2. **臨界帯の三角形分割**: $D(\nu^-) \oplus D(\nu^+)$-束を，底空間 $C_i$ の
   三角形分割の上の**セル的な束**として切る必要がある。$C_i$ が三角形分割され
   束が単体ごとに自明化されていれば標準的（fiberwise cone）だが，書き下すと長い。
3. **r-連結性**: 帯が連結でも，値 $c$ のセル同士が r-path で 1 つに繋がることは
   自明ではない（§5 の実験ではすべて成立）。max_extension では帯の
   「上端に触れるセル」の双対的連結性に帰着するはずで，L1 の三角形分割を
   十分細かく取れば通る見込み。
4. **「近似」の意味**: Gallais 型（臨界構造の実現）を主定理にし，$C^0$-近似は
   系として添える。逆向き（離散 → 滑らか, Benedetti の dim ≤ 7 に相当）は
   **本計画ではやらない**（別の問題）。
5. **多重度**: results.md §9 の旧記述「素朴な引き戻しは (MB2) を破る」は，
   準位集合上に頂点を置かない粗い引き戻しの話。L1 のように準位を部分複体に
   すれば max_extension（L2）で回避される — これが §5 の実験がすべて
   一発で DMBF になっている理由である。

## 7. 論文の骨組み（案）

1. Introduction: Gallais–Benedetti の系譜と Morse-Bott 版の主張。
2. 離散 Morse-Bott 理論の復習（arXiv:2511.07864）。
3. L5 と正則帯の消滅（collapse が不要になる，という方法論的な寄与）。
4. 適合三角形分割（L1）と主定理の証明。
5. 対称性つきの精密化（トーラスの $\mathbb{Z}_{n_i}$ 例 = results.md §2–3 が実例）。
6. 計算例と実装（本リポジトリ; 紡錘形・クライン・円筒の表）。

## 8. 次の一手

- [x] L5 は `tetsuo-jp/discrete-morse-bott`（P-dmb-01）で Lean 化済みと判明（2026-08-25）。紙の論文には引用または 5 行の証明を書けばよい
- [ ] L1 の出典確定（Whitehead / Munkres / Gallais §3.2 の再利用）→ 要確認の解消
- [ ] 臨界帯の局所モデル（円板束の fiberwise cone 三角形分割）で L4 を書く
- [ ] $S^2 \times S^1$ の回転など，臨界多様体がトーラスになる 3 次元例を機械で追加
- [ ] 戦略 A の実装実験: `morsify` の逆（クラスタを潰す）が (MB2)(MB4) を
      保つかを乱択で測る（保険の実現可能性の見積り）
