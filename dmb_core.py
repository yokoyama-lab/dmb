#!/usr/bin/env python3
"""トーラスの三角形分割上の離散モースボット理論 (discrete Morse-Bott theory) の計算。

このモジュールは，dmf.py が扱っていた「トーラスの三角形分割上の離散モース関数
(discrete Morse function, DMF)」の計算を，離散モースボット関数
(discrete Morse-Bott function, DMBF) の計算に書き換えたものである。

定義は Y. Nishikawa and T. Yokoyama, *On discrete Morse-Bott theory*,
arXiv:2511.07864 v2 に従う:

  U(σ)      = #{τ^(k+1) : σ ≺ τ, f(σ) ≥ f(τ)}          (noncritical, nc)
  D(σ)      = #{ν^(k-1) : ν ≺ σ, f(ν) ≥ f(σ)}
  U^snc(σ)  = #{τ^(k+1) : σ ≺ τ, f(σ) > f(τ)}          (strictly noncritical, snc)
  D^snc(σ)  = #{ν^(k-1) : ν ≺ σ, f(ν) > f(σ)}

  DMF  (Definition 12): (M1) (M2) U(σ) ≤ 1      (M3) (M4) D(σ) ≤ 1
  DMBF (Definition 18): (M1) (MB2) U^snc(σ) ≤ 1 (M3) (MB4) D^snc(σ) ≤ 1

  (M1)/(M3) は「余次元によらず irregular な face では値が真に増える」という条件で，
  単体的複体では全ての face が regular なので自動的に成り立つ（ここでの複体は
  すべて単体的なので，実装でも vacuous として扱う）。

  collection         : 値が等しく r-path で繋がるセルの同値類（= f の「臨界多様体」の離散版）
  reduced collection : collection のうち weakly critical (U^snc = D^snc = 0) なセルからなる部分
  Theorem 4.12       : Σ_C P_t(C) = P_t(K) + (1 + t) R(t),  R(t) ≥ 0
                       ここで P_t(C) は境界作用素を C に制限した複体 (C_*(C), ∂^C) の
                       ポアンカレ多項式。

DMF は「すべての collection の大きさが 1 か 2」であるような DMBF なので (Theorem 3.2)，
dmf.py が計算していたものはこの計算の特別な場合として復元される。

使い方:
    python3 dmb_core.py          # T(4,4) 上の 7 種類の関数について報告を出す
    python3 dmb_core.py 5 4      # ni=5, nj=4 のトーラスで
    python3 dmb_core.py --table  # 対称性を課したときの DMT と DMBT の差の表
"""

import ast
import os
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd

# ============================================================ 単体的複体


class Complex:
    """有限単体的複体。セルは頂点の tuple（昇順）で表す。"""

    def __init__(self, facets):
        cells = set()
        for fc in facets:
            fc = tuple(sorted(fc))
            for k in range(1, len(fc) + 1):
                cells.update(combinations(fc, k))
        self.cells = sorted(cells, key=lambda s: (len(s), s))
        self.below = {c: [] for c in self.cells}   # 余次元 1 の face
        self.above = {c: [] for c in self.cells}   # 余次元 1 の coface
        for c in self.cells:
            if len(c) == 1:
                continue
            for k in range(len(c)):
                d = c[:k] + c[k + 1:]
                self.below[c].append(d)
                self.above[d].append(c)

    @staticmethod
    def dim(c):
        return len(c) - 1

    def cells_of_dim(self, k):
        return [c for c in self.cells if self.dim(c) == k]

    def counts(self):
        out = {}
        for c in self.cells:
            out[self.dim(c)] = out.get(self.dim(c), 0) + 1
        return dict(sorted(out.items()))

    # -- CWComplex と共通のインタフェース（単体的複体では face はすべて regular）

    @staticmethod
    def lt(a, b):
        """a < b（任意余次元の face 関係）。"""
        return set(a) < set(b)

    @staticmethod
    def is_regular(a, b):
        """単体的複体では，どの face も regular（論文 §2 の意味で）。"""
        return True

    def irregular_above(self, c):
        """c を irregular な face に持つセル（任意余次元）。単体的複体では常に空。"""
        return []

    def irregular_below(self, c):
        return []

    def incidence(self, a, b):
        return simplicial_incidence(a, b)


def simplicial_incidence(a, b):
    """余次元 1 の関係 a ≺ b に対する接続係数 [b:a]（頂点の順序による標準の向き）。"""
    sa, sb = set(a), set(b)
    if len(sb) != len(sa) + 1 or not sa < sb:
        return 0
    missing = (sb - sa).pop()
    return (-1) ** b.index(missing)


incidence = simplicial_incidence   # 後方互換の別名


# ------------------------------------------------------------ 有理数上の階数


def rank_dense(matrix, p=0):
    """密行列の素朴な Gauss 消去による階数。`rank` の参照実装（検査用）。"""
    if p:
        m = [[x % p for x in row] for row in matrix]
    else:
        m = [[Fraction(x) for x in row] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if m else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        if p:
            invp = pow(pv, -1, p)
            m[r] = [(x * invp) % p for x in m[r]]
        else:
            m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                if p:
                    m[i] = [(x - factor * y) % p for x, y in zip(m[i], m[r])]
                else:
                    m[i] = [x - factor * y for x, y in zip(m[i], m[r])]
        r += 1
        if r == rows:
            break
    return r


def rank(matrix, p=0):
    """行列の階数。p = 0 なら有理数体上，p が素数なら F_p 上。

    境界作用素の行列は非常に疎（列ごとに非零が高々 dim+1 個）なので，行を
    「列 -> 値」の辞書で持って消去する。密行列の Gauss 消去と同じ答えを返すが，
    T(16,16) の ∂_2（768x512）で 118 秒 → 0.1 秒になる（`rank_dense` と
    突き合わせて検査してある）。

    係数体を変えられるようにしてあるのは，ねじれのある空間で
    「離散モース不等式の差が係数体に依る」ことを見るため（RP^2 が典型）。"""
    pivots = {}
    r = 0
    for row in matrix:
        if p:
            cur = {j: v % p for j, v in enumerate(row) if v % p}
        else:
            cur = {j: Fraction(v) for j, v in enumerate(row) if v}
        while cur:
            c = min(cur)
            piv = pivots.get(c)
            if piv is None:
                if p:
                    inv = pow(cur[c], -1, p)
                    pivots[c] = {j: (v * inv) % p for j, v in cur.items()}
                else:
                    inv = 1 / cur[c]
                    pivots[c] = {j: v * inv for j, v in cur.items()}
                r += 1
                break
            factor = cur[c]
            for j, v in piv.items():
                nv = cur.get(j, 0) - factor * v
                if p:
                    nv %= p
                if nv:
                    cur[j] = nv
                elif j in cur:
                    del cur[j]
    return r


def smith_diagonal(matrix):
    """整数行列の Smith 標準形の対角成分（0 でないものを約数の順に）。

    行・列の基本変形だけで対角化する素朴な実装。ここで扱う複体は小さいので十分。
    ∂_{k+1} の対角成分のうち 1 より大きいものが H_k のねじれ Z/d を与える。"""
    m = [list(row) for row in matrix]
    rows, cols = len(m), (len(m[0]) if m else 0)
    out = []
    t = 0
    while t < rows and t < cols:
        # 残りの部分行列で絶対値最小の非零成分を (t, t) に持ってくる
        best = None
        for i in range(t, rows):
            for j in range(t, cols):
                if m[i][j] and (best is None or abs(m[i][j]) < abs(m[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        bi, bj = best
        m[t], m[bi] = m[bi], m[t]
        for row in m:
            row[t], row[bj] = row[bj], row[t]
        # (t, t) で残りの行・列を掃き出す。割り切れないうちは繰り返す
        while True:
            done = True
            for i in range(t + 1, rows):
                if m[i][t]:
                    q = m[i][t] // m[t][t]
                    m[i] = [a - q * b for a, b in zip(m[i], m[t])]
                    if m[i][t]:
                        m[t], m[i] = m[i], m[t]
                        done = False
            for j in range(t + 1, cols):
                if m[t][j]:
                    q = m[t][j] // m[t][t]
                    for row in m:
                        row[j] -= q * row[t]
                    if m[t][j]:
                        for row in m:
                            row[t], row[j] = row[j], row[t]
                        done = False
            if done:
                break
        out.append(abs(m[t][t]))
        t += 1
    # 対角化しただけでは d_1 | d_2 | ... にならないので，gcd/lcm で整える
    # （これで初めて invariant factor になる）
    changed = True
    while changed:
        changed = False
        for i in range(len(out) - 1):
            a, b = out[i], out[i + 1]
            if b % a:
                g = gcd(a, b)
                out[i], out[i + 1] = g, a * b // g
                changed = True
    return out


def homology_z(cells, inc=None, dim=None):
    """整数係数ホモロジー。次数ごとに (自由部分の階数, ねじれ係数のリスト)。

    例: RP^2 なら [(1, []), (0, [2]), (0, [])] すなわち H_0 = Z, H_1 = Z/2, H_2 = 0。"""
    inc = simplicial_incidence if inc is None else inc
    dim = (lambda c: len(c) - 1) if dim is None else dim
    by_dim = {}
    for c in cells:
        by_dim.setdefault(dim(c), []).append(c)
    if not by_dim:
        return []
    maxd = max(by_dim)

    def matrix(k):
        """∂_k : C_k → C_{k-1} の行列。"""
        rowsp, colsp = by_dim.get(k - 1, []), by_dim.get(k, [])
        if not rowsp or not colsp:
            return []
        return [[inc(a, b) for b in colsp] for a in rowsp]

    ranks, divisors = {}, {}
    for k in range(maxd + 2):
        mat = matrix(k)
        d = [x for x in smith_diagonal(mat) if x] if mat else []
        divisors[k] = d
        ranks[k] = len(d)
    out = []
    for k in range(maxd + 1):
        nk = len(by_dim.get(k, []))
        free = nk - ranks.get(k, 0) - ranks.get(k + 1, 0)
        torsion = [d for d in divisors.get(k + 1, []) if d > 1]
        out.append((free, torsion))
    return out


def betti(cells, inc=None, dim=None, p=0):
    """セル集合 `cells` に境界作用素を制限した鎖複体の有理ベッチ数。

    `cells` が複体全体なら通常のベッチ数，reduced collection C なら
    (C_*(C), ∂^C) のベッチ数（Theorem 4.12 の P_t(C)）になる。

    inc は接続係数 [b:a] を返す関数，dim はセルの次元を返す関数（どちらも省略時は
    単体的複体のもの）。CW 複体ではその複体の `incidence` と `dim` を渡す。
    p = 0 なら有理数体上，p が素数なら F_p 上のベッチ数。"""
    inc = simplicial_incidence if inc is None else inc
    dim = (lambda c: len(c) - 1) if dim is None else dim
    by_dim = {}
    for c in cells:
        by_dim.setdefault(dim(c), []).append(c)
    if not by_dim:
        return []
    maxd = max(by_dim)
    ranks = {}
    for k in range(maxd + 1):
        rowsp = by_dim.get(k - 1, [])
        colsp = by_dim.get(k, [])
        if not rowsp or not colsp:
            ranks[k] = 0
            continue
        ranks[k] = rank([[inc(a, b) for b in colsp] for a in rowsp], p)
    out = []
    for k in range(maxd + 1):
        ck = len(by_dim.get(k, []))
        out.append(ck - ranks.get(k, 0) - ranks.get(k + 1, 0))
    while out and out[-1] == 0:
        out.pop()
    return out


# ------------------------------------------------------------ 多項式（係数列）


def homology_str(hz):
    """homology_z() の結果を "Z, Z/2, 0" のような文字列にする。"""
    out = []
    for free, tors in hz:
        parts = ([f"Z^{free}" if free > 1 else "Z"] if free else []) + [f"Z/{d}" for d in tors]
        out.append(" ⊕ ".join(parts) if parts else "0")
    return ", ".join(f"H_{k} = {t}" for k, t in enumerate(out))


def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def poly_sub(a, b):
    return poly_add(a, [-c for c in b])


def poly_trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def poly_str(p):
    p = poly_trim(p)
    if not p:
        return "0"
    terms = []
    for k, c in enumerate(p):
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        elif k == 1:
            terms.append("t" if c == 1 else f"{c}t")
        else:
            terms.append(f"t^{k}" if c == 1 else f"{c}t^{k}")
    return " + ".join(terms) if terms else "0"


def poly_div_1pt(p):
    """p = (1 + t) q なる q を返す。割り切れないときは None。"""
    p = poly_trim(p)
    if not p:
        return []
    q = []
    carry = 0
    for c in p:
        carry = c - carry
        q.append(carry)
    if q[-1] != 0:
        return None
    q.pop()
    if poly_trim(poly_add(q, [0] + q)) != p:
        return None
    return q


# ============================================================ CW 複体


class CWComplex:
    """正則とは限らない有限 CW 複体（抽象的な face 関係として扱う）。

    論文 §2 の意味の face 関係 σ < τ（任意余次元），facet σ ≺ τ（余次元 1），
    regular / irregular の区別を持つ。単体的複体はすべての face が regular な
    特別な場合で，(M1)(M3) が vacuous になる。ここでは irregular な face を
    持つ複体（最小 CW 構造など）を扱えるようにする。

    引数:
        dims       {cell: 次元}
        faces      (σ, τ) の列。σ < τ を表す（推移閉包を取るので生成元でよい）
        regular    regular な (σ, τ) の集合（省略すると「すべて irregular」）
        incidence  {(σ, τ): 整数} 余次元 1 の接続係数 [τ:σ]（省略すると 0）

    接続係数は自分で与える。∂∘∂ = 0 になっているかは `check_boundary()` で
    確かめられる（ホモロジーを計算する前に呼ぶこと）。
    """

    def __init__(self, dims, faces, regular=(), incidence=None, name=""):
        self.name = name
        self._dim = dict(dims)
        self.cells = sorted(self._dim, key=lambda c: (self._dim[c], repr(c)))
        # 推移閉包（生成元だけ渡してよいように）
        rel = {(a, b) for a, b in faces}
        changed = True
        while changed:
            changed = False
            for a, b in list(rel):
                for c, d in list(rel):
                    if b == c and (a, d) not in rel:
                        rel.add((a, d))
                        changed = True
        for a, b in rel:
            if self._dim[a] >= self._dim[b]:
                raise ValueError(f"face 関係が次元を上げていない: {a} < {b}")
        self.faces = rel
        self.regular_pairs = {(a, b) for a, b in regular}
        if not self.regular_pairs <= rel:
            raise ValueError("regular な組が face 関係に入っていない")
        self._inc = dict(incidence or {})
        self.below = {c: [] for c in self.cells}
        self.above = {c: [] for c in self.cells}
        for a, b in sorted(rel, key=repr):
            if self._dim[b] == self._dim[a] + 1:
                self.below[b].append(a)
                self.above[a].append(b)

    def dim(self, c):
        return self._dim[c]

    def cells_of_dim(self, k):
        return [c for c in self.cells if self._dim[c] == k]

    def counts(self):
        out = {}
        for c in self.cells:
            out[self._dim[c]] = out.get(self._dim[c], 0) + 1
        return dict(sorted(out.items()))

    def lt(self, a, b):
        return (a, b) in self.faces

    def is_regular(self, a, b):
        return (a, b) in self.regular_pairs

    def irregular_above(self, c):
        """c を irregular な face に持つセル（任意余次元）。(M1) が効く相手。"""
        return [b for (a, b) in sorted(self.faces, key=repr)
                if a == c and not self.is_regular(a, b)]

    def irregular_below(self, c):
        """c の irregular な face（任意余次元）。(M3) が効く相手。"""
        return [a for (a, b) in sorted(self.faces, key=repr)
                if b == c and not self.is_regular(a, b)]

    def incidence(self, a, b):
        return self._inc.get((a, b), 0)

    def check_boundary(self):
        """∂∘∂ = 0 を確かめる。破っている (ν, τ) の一覧を返す（空なら健全）。"""
        bad = []
        for t in self.cells:
            for v in self.cells:
                if self._dim[v] != self._dim[t] - 2:
                    continue
                total = sum(self.incidence(v, s) * self.incidence(s, t)
                            for s in self.below[t])
                if total != 0:
                    bad.append((v, t, total))
        return bad

    def betti(self, p=0):
        return betti(self.cells, self.incidence, self.dim, p)

    def homology_z(self):
        return homology_z(self.cells, self.incidence, self.dim)


# ------------------------------------------------------- 最小 CW 構造の例


def cw_sphere2_minimal():
    """S^2 の最小 CW 構造（0-セル 1 個，2-セル 1 個）。

    v は e の余次元 2 の irregular な face。余次元 1 の組が 1 つも無いので，
    論文 v1 の Definition 12 を字義通り（irregular *facet* にだけ (M1)(M3) を課す）
    読むと**どんな関数も離散モース関数になってしまう**。これが v1 Lemma 3.1
    （DMF ⇒ DMBF）の反例（Lean の `Counterexample.dmf_not_imp_dmb`）。
    v2 の Definition 12 は任意余次元に強めてあるので問題は起きない。"""
    return CWComplex({"v": 0, "e": 2}, [("v", "e")], name="S^2 (最小 CW)")


def cw_circle_minimal():
    """S^1 の最小 CW 構造（0-セル 1 個，1-セル 1 個）。∂e = v - v = 0。"""
    return CWComplex({"v": 0, "e": 1}, [("v", "e")], incidence={("v", "e"): 0},
                     name="S^1 (最小 CW)")


def cw_torus_minimal():
    """トーラスの最小 CW 構造（0-セル 1，1-セル 2，2-セル 1）。

    貼り付け語は a b a^-1 b^-1 なので接続係数はすべて 0。有理数係数のベッチ数は
    (1, 2, 1) で，セル数も (1, 2, 1)。すべての face が irregular。"""
    return CWComplex({"v": 0, "a": 1, "b": 1, "F": 2},
                     [("v", "a"), ("v", "b"), ("a", "F"), ("b", "F")],
                     incidence={("v", "a"): 0, ("v", "b"): 0,
                                ("a", "F"): 0, ("b", "F"): 0},
                     name="T^2 (最小 CW)")


def cw_projective_plane_minimal():
    """RP^2 の最小 CW 構造（各次元 1 個）。貼り付け写像の次数は 2。

    ∂e^2 = 2 e^1，∂e^1 = 0 なので H_1 = Z/2（ねじれ）。有理数係数では P_t = 1。"""
    return CWComplex({"v": 0, "a": 1, "F": 2},
                     [("v", "a"), ("v", "F"), ("a", "F")],
                     incidence={("v", "a"): 0, ("a", "F"): 2},
                     name="RP^2 (最小 CW)")


def cw_from_simplicial(K, name=""):
    """単体的複体を CW 複体として見る（すべての face が regular）。"""
    dims = {c: K.dim(c) for c in K.cells}
    faces = [(a, b) for a in K.cells for b in K.cells if K.lt(a, b)]
    inc = {(a, b): simplicial_incidence(a, b) for a, b in faces}
    return CWComplex(dims, faces, regular=faces, incidence=inc, name=name)


# ============================================================ Morse-Bott 構造


class MorseBott:
    """複体 K 上の関数 f に対する（離散モース／モースボット）データ。"""

    def __init__(self, K, f, p=0):
        """p は係数体（0 = 有理数体，素数 p = F_p）。ホモロジーの計算だけに効く。"""
        self.K = K
        self.f = f
        self.p = p

    # -- 近傍の数え上げ ------------------------------------------------

    def up_nc(self, s):
        """noncritical な coface: f(σ) ≥ f(τ)。"""
        return [t for t in self.K.above[s] if self.f[t] <= self.f[s]]

    def down_nc(self, s):
        """noncritical な face: f(ν) ≥ f(σ)。"""
        return [v for v in self.K.below[s] if self.f[v] >= self.f[s]]

    def up_snc(self, s):
        """strictly noncritical な coface: f(σ) > f(τ)。"""
        return [t for t in self.K.above[s] if self.f[t] < self.f[s]]

    def down_snc(self, s):
        """strictly noncritical な face: f(ν) > f(σ)。"""
        return [v for v in self.K.below[s] if self.f[v] > self.f[s]]

    # -- 定義の検査 ----------------------------------------------------

    # -- (M1)/(M3): irregular な face では値が真に増える ------------------

    def m1_violations(self, strong=True):
        """(M1) σ <^irr τ ⇒ f(σ) < f(τ) を破る組。

        strong=True は論文 v2 の Definition 12/18（**任意余次元**の irregular face）。
        strong=False は v1 の Definition 12 を字義通り読んだ場合（余次元 1 だけ）で，
        `cw_sphere2_minimal()` がその読みの反例になる。単体的複体ではどちらも vacuous。"""
        bad = []
        for c in self.K.cells:
            for t in self.K.irregular_above(c):
                if not strong and self.K.dim(t) != self.K.dim(c) + 1:
                    continue
                if not self.f[c] < self.f[t]:
                    bad.append(("M1", c, [t]))
        return bad

    def m3_violations(self, strong=True):
        """(M3) ν <^irr σ ⇒ f(ν) < f(σ) を破る組。"""
        bad = []
        for c in self.K.cells:
            for v in self.K.irregular_below(c):
                if not strong and self.K.dim(v) != self.K.dim(c) - 1:
                    continue
                if not self.f[v] < self.f[c]:
                    bad.append(("M3", c, [v]))
        return bad

    def dmf_violations(self, strong=True):
        """離散モース関数 (Definition 12) を破る組。

        (M1)(M3) は irregular な face についての条件で，単体的複体では vacuous。
        (M2) U(σ) ≤ 1，(M4) D(σ) ≤ 1。"""
        bad = self.m1_violations(strong) + self.m3_violations(strong)
        for c in self.K.cells:
            u, d = self.up_nc(c), self.down_nc(c)
            if len(u) > 1:
                bad.append(("M2", c, u))
            if len(d) > 1:
                bad.append(("M4", c, d))
        return bad

    def dmb_violations(self):
        """離散モースボット関数 (Definition 18) を破る組。

        (M1)(M3) は任意余次元の irregular face に課される（v1・v2 とも）。
        (MB2) U^snc(σ) ≤ 1，(MB4) D^snc(σ) ≤ 1。"""
        bad = self.m1_violations(True) + self.m3_violations(True)
        for c in self.K.cells:
            u, d = self.up_snc(c), self.down_snc(c)
            if len(u) > 1:
                bad.append(("MB2", c, u))
            if len(d) > 1:
                bad.append(("MB4", c, d))
        return bad

    def is_dmf(self, strong=True):
        return not self.dmf_violations(strong)

    def is_dmb(self):
        return not self.dmb_violations()

    # -- セルの分類 ----------------------------------------------------

    def is_critical(self, s):
        return not self.up_nc(s) and not self.down_nc(s)

    def is_weakly_critical(self, s):
        return not self.up_snc(s) and not self.down_snc(s)

    def critical(self):
        return [c for c in self.K.cells if self.is_critical(c)]

    def weakly_critical(self):
        return [c for c in self.K.cells if self.is_weakly_critical(c)]

    def arrows(self):
        """strictly noncritical pair (σ, τ)，σ ≺ τ かつ f(τ) < f(σ)。

        DMF のときは Forman の組合せベクトル場（V-path の矢印）に一致する。"""
        return [(s, t) for s in self.K.cells for t in self.up_snc(s)]

    # -- collection ----------------------------------------------------

    def collections(self):
        """r-path による同値類。値が等しい余次元 1 のペアを union-find で繋ぐ。"""
        parent = {c: c for c in self.K.cells}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for c in self.K.cells:
            for d in self.K.above[c]:
                if self.f[c] == self.f[d]:
                    rx, ry = find(c), find(d)
                    if rx != ry:
                        parent[rx] = ry
        classes = {}
        for c in self.K.cells:
            classes.setdefault(find(c), []).append(c)
        return [sorted(v, key=lambda s: (len(s), s)) for v in classes.values()]

    def reduced_collections(self):
        """reduced collection の一覧（Definition 21-23 / v2 の reduced collection）。

        定義は「(1) 値が共通 (2) 任意の 2 元が c-path で繋がる (3) すべて weakly
        critical」を満たす極大族。c-path は $f^{-1}(c)$ の中を通ればよく，C の中を
        通る必要はないので，条件 (2) は「同じ collection に属する」ことと同値。
        したがって reduced collection = collection の weakly critical なセル全体
        であり，weakly critical でないセルで分断されていても 1 つと数える
        （Lean 形式化の `Reduced` もセル単位の述語で，分割していない）。

        なお C を facet 関係で連結成分に分けても (C_*(C), ∂^C) は直和に
        分かれるだけなので，Σ_C P_t(C) は分け方に依らない。"""
        wc = set(self.weakly_critical())
        out = []
        for L in self.collections():
            C = [c for c in L if c in wc]
            if C:
                out.append(sorted(C, key=lambda s: (len(s), s)))
        return out

    # -- 多項式 --------------------------------------------------------

    def betti(self, cells, p=None):
        """セル集合の（制限された）ベッチ数。複体自身の接続係数と次元を使う。"""
        return betti(cells, self.K.incidence, self.K.dim,
                     self.p if p is None else p)

    def morse_polynomial(self):
        """M(t) = Σ_{σ critical} t^{dim σ}（離散モース理論の側）。"""
        p = []
        for c in self.critical():
            p = poly_add(p, [0] * self.K.dim(c) + [1])
        return poly_trim(p)

    def morse_bott_polynomial(self):
        """Σ_C P_t(C)（離散モースボット理論の側）。"""
        p = []
        for C in self.reduced_collections():
            p = poly_add(p, self.betti(C))
        return poly_trim(p)

    def report(self):
        """Theorem 4.12（および離散モース不等式）の検査結果をまとめて返す。"""
        pk = poly_trim(self.betti(self.K.cells))
        mb = self.morse_bott_polynomial()
        rb = poly_div_1pt(poly_sub(mb, pk))
        res = {
            "is_dmf": self.is_dmf(),
            "is_dmb": self.is_dmb(),
            "dmb_violations": self.dmb_violations(),
            "dmf_violations": self.dmf_violations(),
            "P_K": pk,
            "collections": len(self.collections()),
            "reduced_collections": self.reduced_collections(),
            "n_weakly_critical": len(self.weakly_critical()),
            "n_critical": len(self.critical()),
            "n_arrows": len(self.arrows()),
            "MB_sum": mb,
            "R_MB": rb,
            "MB_sharp": rb is not None and not poly_trim(rb),
        }
        if self.is_dmf():
            mt = self.morse_polynomial()
            rm = poly_div_1pt(poly_sub(mt, pk))
            res["M"] = mt
            res["R_M"] = rm
        return res


# ============================================================ トーラス


def torus(ni, nj):
    """トーラスの三角形分割。頂点は (i, j) ∈ Z_ni × Z_nj。

    各正方形を dmf.py と同じ向きに 2 つの三角形に分ける:
        上三角 {(i,j), (i,j+1), (i+1,j+1)}
        下三角 {(i,j), (i+1,j), (i+1,j+1)}
    共有される対角線は {(i,j), (i+1,j+1)}。"""
    if ni < 3 or nj < 3:
        raise ValueError("単体的複体になるためには ni, nj ≥ 3 が必要")
    facets = []
    for i in range(ni):
        for j in range(nj):
            a = (i, j)
            b = (i, (j + 1) % nj)
            c = ((i + 1) % ni, (j + 1) % nj)
            d = ((i + 1) % ni, j)
            facets.append((a, b, c))
            facets.append((a, d, c))
    return Complex(facets)


def torus_names(ni, nj):
    """dmf.py と同じセル名 (v_k / e_k / f_k) への辞書 cell -> name を返す。

    ni == nj == grid_size のとき dmf.py の calcDMF が返す添字と一致する。"""
    n = ni * nj
    name = {}
    for i in range(ni):
        for j in range(nj):
            k = i * nj + j
            ip, jp = (i + 1) % ni, (j + 1) % nj
            name[cell((i, j))] = f"v_{k}"
            name[cell((i, j), (i, jp))] = f"e_{k}"                 # 縦（j 方向）
            name[cell((i, j), (ip, j))] = f"e_{n + k}"             # 横（i 方向）
            name[cell((i, j), (ip, jp))] = f"e_{2 * n + k}"        # 斜め
            name[cell((i, j), (i, jp), (ip, jp))] = f"f_{k}"       # 上三角
            name[cell((i, j), (ip, j), (ip, jp))] = f"f_{n + k}"   # 下三角
    return name


def cell(*vertices):
    return tuple(sorted(vertices))


def rotate(ni, nj, c, k=1):
    """i 方向の回転 Z_ni の作用（トーラスの対称性）。"""
    return cell(*[((i + k) % ni, j) for (i, j) in c])


def is_invariant(ni, nj, f):
    """f が i 方向の回転 Z_ni で不変か。"""
    return all(f[c] == f[rotate(ni, nj, c)] for c in f)


# ============================================================ 描画用の座標


def lifted_cells(ni, nj):
    """各セルに，基本領域 [0,ni] x [0,nj] の中での描画位置（持ち上げ）を与える。

    戻り値は cell -> [placement, ...]，placement は [(x, y), ...]（頂点は 1 点，
    辺は 2 点，三角形は 3 点）。最初の placement が標準の持ち上げ。

    トーラスの貼り合わせで図を閉じて見せるため，(ni,0), (0,nj), (ni,nj) だけ
    平行移動しても基本領域に収まるセルは，その複製も持つ。これは dmf.py が
    アポストロフィ付きの名前で複製していた境界のセルに対応する（例: j=0 の
    横の辺は上端 y=nj にも描かれる）。"""
    base = {}
    for i in range(ni):
        for j in range(nj):
            ip, jp = (i + 1) % ni, (j + 1) % nj
            base[cell((i, j))] = [(i, j)]
            base[cell((i, j), (i, jp))] = [(i, j), (i, j + 1)]
            base[cell((i, j), (ip, j))] = [(i, j), (i + 1, j)]
            base[cell((i, j), (ip, jp))] = [(i, j), (i + 1, j + 1)]
            base[cell((i, j), (i, jp), (ip, jp))] = [(i, j), (i, j + 1), (i + 1, j + 1)]
            base[cell((i, j), (ip, j), (ip, jp))] = [(i, j), (i + 1, j), (i + 1, j + 1)]
    out = {}
    for c, pts in base.items():
        places = []
        for a in (0, 1):
            for b in (0, 1):
                q = [(x + a * ni, y + b * nj) for x, y in pts]
                if all(0 <= x <= ni and 0 <= y <= nj for x, y in q):
                    places.append(q)
        out[c] = places
    return out


def centroid(points):
    return (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points))


# ============================================================ 関数の例


def constant_fn(K, value=0):
    """定数関数。自明な DMBF（collection は K 全体ひとつ）。"""
    return dict.fromkeys(K.cells, value)


def max_extension(K, h, refine=False, span_key=None):
    """頂点上の関数 h を「頂点での最大値」でセルに拡張する。

    σ ≺ τ なら常に f(σ) ≤ f(τ) なので U^snc = D^snc = 0，すなわち
    **どんな複体・どんな h でも離散モースボット関数**になり，全セルが weakly
    critical。collection は f の準位集合の連結成分で，滑らかな Morse-Bott 関数の
    準位集合に対応する。

    refine=True では，2 つ以上の準位にまたがるセルの値を 1 つ持ち上げる
    （準位集合が細かく割れる）。「準位」は span_key で測る（既定は h 自身）。
    h が単射でないときは span_key を分けると意味が変わるので注意。"""
    span_key = span_key or h
    f = {}
    for c in K.cells:
        m = max(h(v) for v in c)
        f[c] = 2 * m + (1 if len({span_key(v) for v in c}) > 1 else 0) if refine else m
    return f


def height_fn(K, ni, nj, refine=False):
    """トーラス上の回転対称な離散モースボット関数（φ 方向の高さ）。

    h(j) = min(j, nj - j) を頂点に与えて max_extension したもの。滑らかな
    Morse-Bott 関数（回転対称なトーラスの，回転軸からの距離）の離散化にあたり，
    θ 方向の回転 Z_ni で不変。"""
    return max_extension(K, lambda v: min(v[1] % nj, (-v[1]) % nj), refine,
                         span_key=lambda v: v[1])


def arrowed_dmbf(K, ni, nj):
    """矢印を持つ，回転対称で鋭い離散モースボット関数（帯を 1 つ持ち上げる）。

    最上位の帯（準位 nj-1 と 0 の間）にある **1 次元以上の**セルにだけ値 1 を与え，
    残り（頂点をすべて含む）は 0 とする。高さ関数（`height_fn`）と同じ
    Σ_C P_t(C) = P_t(T^2)，R(t) = 0 を与えるが，**矢印を ni 本持つ**点が違う。

    構成の要点（`ni, nj >= 3` で成り立つ。(M1)(M3) は単体的複体なので vacuous）:

    * (MB2): 値が下がる coface があるのは最上位の帯の**横の辺** ie_{nj-1} だけで，
      その coface は tl_{nj-1}（値 1）と tu_{nj-2}（値 0）なので U^snc = 1。
      縦の辺 je_{nj-1}・斜めの辺 dg_{nj-1} の coface はどちらも同じ帯の三角形
      （値 1）なので U^snc = 0。
    * (MB4): 値が上がる face を持つのは tu_{nj-2} だけで，その辺は
      je_{nj-2}, dg_{nj-2}（値 0）と ie_{nj-1}（値 1）だから D^snc = 1。
    * collection は値 0 の部分と値 1 の部分の 2 つ。reduced collection は
      前者から tu_{nj-2} を，後者から ie_{nj-1} を除いたもので，
      P_t = 1 + t と t + t^2（`ni, nj <= 7` の範囲で計算して確認）。

    「対称性を保ったまま鋭くする」やり方が高さ関数の形しかないわけではないことを
    示す例（docs/results.md §3.6）。"""
    top = nj - 1

    def lifted(c):
        js = {v[1] % nj for v in c}
        return len(c) > 1 and top in js and js <= {top, 0}

    return {c: (1 if lifted(c) else 0) for c in K.cells}


def morsify(K, f):
    """Morsification: f'(σ) = (D + 1) f(σ) + dim σ  (D = dim K)。

    arXiv:2511.07864 の Lemma 4.8 の整数版。DMBF f に対して f' は DMF になり，
    f' の臨界セルはちょうど f の reduced collection の和集合になる。"""
    D = max(K.dim(c) for c in K.cells)
    return {c: (D + 1) * f[c] + K.dim(c) for c in K.cells}


# ------------------------------------------------------- acyclic matching から


def function_from_matching(K, matching, orbit=None):
    """acyclic matching から離散モース関数を作る（修正ハッセ図の位相ソート）。

    matching は σ ≺ τ なる (σ, τ) の集合で，各セルは高々 1 回現れること。
    修正ハッセ図は，matching のペアでは τ → σ，それ以外の σ ≺ τ では σ → τ とした
    有向グラフ。これを位相ソートして順位を値にすると，matching のペアで
    f(τ) < f(σ)，それ以外で f(σ) < f(τ) となり，(M2)(M4) が成り立つ。

    orbit を与えると，同じ軌道のセルに同じ値を与える（対称な DMF）。軌道で
    割った有向グラフに閉路があれば ValueError を投げる。"""
    partner = {}
    for s, t in matching:
        if s in partner or t in partner:
            raise ValueError(f"matching が単射でない: {s}, {t}")
        partner[s] = t
        partner[t] = s

    key = orbit if orbit is not None else (lambda c: c)
    nodes = {key(c) for c in K.cells}
    succ = {n: set() for n in nodes}
    indeg = dict.fromkeys(nodes, 0)

    def link(a, b):
        if b not in succ[a]:
            succ[a].add(b)
            indeg[b] += 1

    for s in K.cells:
        for t in K.above[s]:
            if partner.get(s) == t:
                link(key(t), key(s))
            else:
                link(key(s), key(t))

    order = []
    ready = sorted([n for n in nodes if indeg[n] == 0], key=repr)
    while ready:
        n = ready.pop(0)
        order.append(n)
        for m in sorted(succ[n], key=repr):
            indeg[m] -= 1
            if indeg[m] == 0:
                ready.append(m)
    if len(order) != len(nodes):
        raise ValueError("修正ハッセ図に閉路がある（matching が acyclic でない）")
    level = {n: k for k, n in enumerate(order)}
    return {c: level[key(c)] for c in K.cells}


def canonical_dmf(K):
    """tree-cotree 構成による離散モース関数（連結な 2 次元複体なら何でもよい）。

    1-骨格の全域木で頂点と辺を，双対グラフ（全域木に入っていない辺のみを使う）の
    全域木で三角形と辺を対応させる。臨界セルは
        頂点 1 個，辺 (#E - #V - #F + 2) 本，三角形 1 個
    で，閉曲面なら m_1 = 2 - χ 本になる。トーラスでは (1, 2, 1) でベッチ数に一致し
    離散モース不等式は等号だが，全域木の選び方が対称性を破るので
    Z_ni 不変にはならない。"""
    verts = K.cells_of_dim(0)
    faces = K.cells_of_dim(2)
    matching = []

    # 1-骨格の全域木（BFS）: 頂点 ↔ 親への辺
    root = verts[0]
    seen = {root}
    queue = [root]
    tree_edges = set()
    while queue:
        v = queue.pop(0)
        for e in sorted(K.above[v]):
            w = cell(e[0]) if cell(e[1]) == v else cell(e[1])
            if w not in seen:
                seen.add(w)
                tree_edges.add(e)
                matching.append((w, e))
                queue.append(w)
    if len(seen) != len(verts):
        raise ValueError("1-骨格が連結でない")

    # 双対グラフの全域木（全域木に入っていない辺のみを使う）: 三角形 ↔ 親への辺
    root_f = faces[0]
    seen_f = {root_f}
    queue = [root_f]
    while queue:
        t = queue.pop(0)
        for e in sorted(K.below[t]):
            if e in tree_edges:
                continue
            for u in sorted(K.above[e]):
                if u not in seen_f:
                    seen_f.add(u)
                    tree_edges.add(e)
                    matching.append((e, u))
                    queue.append(u)
    if len(seen_f) != len(faces):
        raise ValueError("双対グラフ（全域木の補集合）が連結でない")

    assert len(matching) == len(verts) - 1 + len(faces) - 1
    return function_from_matching(K, matching)


def invariant_dmf(K, ni, nj):
    """Z_ni（θ 方向の回転）で不変な離散モース関数。臨界セルは 4·ni 個。

    j = 0 から順に帯を下へ潰す，回転で不変な matching:
        頂点 (i,j+1)                    ↔ 縦の辺 {(i,j),(i,j+1)}
        横の辺 {(i,j+1),(i+1,j+1)}      ↔ 上三角 {(i,j),(i,j+1),(i+1,j+1)}
        斜めの辺 {(i,j),(i+1,j+1)}      ↔ 下三角 {(i,j),(i+1,j),(i+1,j+1)}
    ただし最後の帯 (j = nj-1) では上 2 つを行わない。臨界セルは
        j=0 の頂点 ni 個，j=0 の横の辺 ni 個，j=nj-1 の縦の辺 ni 個，
        j=nj-1 の上三角 ni 個
    の 4·ni 個で，M(t) = ni (1 + 2t + t^2) = ni · P_t(T^2)。

    命題: Z_ni が自由に作用するので，不変な DMF の臨界セルの個数は各次元で ni の
    倍数であり，離散モース不等式 m_0 ≥ 1, m_2 ≥ 1 と χ = m_0 - m_1 + m_2 = 0 から
    m_0, m_2 ≥ ni, m_1 = m_0 + m_2 ≥ 2·ni，よって総数は 4·ni 以上。この構成は
    その下限を達成する。"""
    matching = []
    for i in range(ni):
        for j in range(nj):
            ip, jp = (i + 1) % ni, (j + 1) % nj
            diag = cell((i, j), (ip, jp))
            lower = cell((i, j), (ip, j), (ip, jp))
            matching.append((diag, lower))
            if j == nj - 1:
                continue
            matching.append((cell((i, jp)), cell((i, j), (i, jp))))
            upper = cell((i, j), (i, jp), (ip, jp))
            matching.append((cell((i, jp), (ip, jp)), upper))
    return function_from_matching(K, matching,
                                  orbit=lambda c: orbit_key(c, ni, nj))


def orbit_key(c, ni, nj):
    """Z_ni の回転軌道の代表元（i を回して辞書順最小にする）。"""
    return min(rotate(ni, nj, c, k) for k in range(ni))


# ---------------------------------------------------- dmf.py の関数の読み込み


def dmf_from_dmf_py(K, ni, nj, path=None):
    """dmf.py の calcDMF が与える離散モース関数を読み込む（dash 等に依存しない）。

    dmf.py はモジュール読み込み時に Dash サーバを起動するので import せず，
    ast で calcDMF の定義だけを取り出して実行する。"""
    if ni != nj:
        raise ValueError("dmf.py の calcDMF は正方格子 (ni == nj) 専用")
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "dmf.py")
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "calcDMF":
            ns = {}
            exec(compile(ast.Module(body=[node], type_ignores=[]), path, "exec"), ns)
            values, _ = ns["calcDMF"](ni)
            name = torus_names(ni, nj)
            return {c: values[name[c]] for c in K.cells}
    raise ValueError(f"{path} に calcDMF が見つからない")


# ============================================================ ランダムな関数


def random_function(K, rng, spread=3):
    """完全にランダムな整数値（多くは DMBF ではない）。"""
    return {c: rng.randrange(spread) for c in K.cells}


def random_max_extension(K, rng, spread=4):
    """頂点にランダムな値を置いて最大値で拡張（必ず DMBF）。"""
    h = {v[0]: rng.randrange(spread) for v in K.cells_of_dim(0)}
    return max_extension(K, h.__getitem__, refine=rng.random() < 0.5)


def random_matching_dmf(K, rng, tries=30):
    """ランダムな acyclic matching から作った離散モース関数（必ず DMF）。"""
    for _ in range(tries):
        cells = list(K.cells)
        rng.shuffle(cells)
        used, matching = set(), []
        for c in cells:
            if c in used:
                continue
            partners = [t for t in K.above[c] if t not in used]
            if partners and rng.random() < 0.8:
                t = rng.choice(partners)
                used.add(c)
                used.add(t)
                matching.append((c, t))
        try:
            return function_from_matching(K, matching)
        except ValueError:
            continue          # 閉軌道ができたら引き直す
    return constant_fn(K)


def hash_seed(*parts):
    """PYTHONHASHSEED に依らない決定的なシード。"""
    h = 1469598103934665603
    for p in parts:
        for ch in str(p).encode():
            h = ((h ^ ch) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


# ============================================================ 報告


def describe(K, ni, nj, title, f, names=None, p=0):
    X = MorseBott(K, f, p=p)
    r = X.report()
    print(f"=== {title}")
    print(f"    セル数 {K.counts()},  χ = "
          f"{sum((-1) ** k * v for k, v in K.counts().items())},  "
          f"値域 {min(f.values())}..{max(f.values())}"
          + (f",  係数体 F_{p}" if p else ""))
    if ni is not None:
        print(f"    Z_{ni} 回転で不変: "
              f"{'はい' if is_invariant(ni, nj, f) else 'いいえ'}")
    if r["is_dmb"]:
        irr = any(K.irregular_above(c) for c in K.cells)
        print("    (M1)(MB2)(M3)(MB4): OK"
              + ("   [irregular な face があるので (M1)(M3) は実質的]" if irr
                 else "   [(M1)(M3) は単体的複体なので vacuous]"))
    else:
        k, cell_, w = r["dmb_violations"][0]
        print(f"    離散モースボット関数ではない: {len(r['dmb_violations'])} 件の違反, "
              f"最初は {k} at {names[cell_] if names else cell_} (witness {len(w)} 個)")
        return r
    print(f"    離散モース関数か: {'はい' if r['is_dmf'] else 'いいえ'}"
          f"（違反 {len(r['dmf_violations'])} 件）")
    print(f"    collection {r['collections']} 個, reduced collection "
          f"{len(r['reduced_collections'])} 個, weakly critical {r['n_weakly_critical']} セル, "
          f"critical {r['n_critical']} セル, 矢印 (snc pair) {r['n_arrows']} 本")
    kinds = {}
    for C in r["reduced_collections"]:
        b = poly_trim(X.betti(C))
        if not b:
            continue
        key = (len(C), tuple(sorted({K.dim(c) for c in C})), tuple(b))
        kinds[key] = kinds.get(key, 0) + 1
    for (size, dims, b), mult in sorted(kinds.items(), reverse=True):
        print(f"      C: {size} セル (次元 {list(dims)}), P_t(C) = {poly_str(b)}"
              f"{f'   × {mult} 個' if mult > 1 else ''}")
    print(f"    P_t(K)     = {poly_str(r['P_K'])}")
    print(f"    Σ_C P_t(C) = {poly_str(r['MB_sum'])}")
    if r["R_MB"] is None:
        print("    !! Σ_C P_t(C) - P_t(K) が (1+t) で割り切れない")
    else:
        print(f"    R(t)       = {poly_str(r['R_MB'])}"
              f"{'   ← 等号（鋭い）' if r['MB_sharp'] else ''}")
    if r["is_dmf"]:
        print(f"    （離散モース理論）M(t) = {poly_str(r['M'])}, "
              f"R_M(t) = {poly_str(r['R_M']) if r['R_M'] is not None else '割り切れない'}")
    return r


def symmetry_table(sizes):
    """対称性を課したときの離散モース理論と離散モースボット理論の差の表。

    トーラス T(ni, nj) には θ 方向の回転 Z_ni が自由に作用する。

      * 臨界セル 4 個の DMF（tree-cotree）は R(t) = 0 だが Z_ni 不変ではない。
      * Z_ni 不変な DMF の臨界セルは 4·ni 個が最小で（invariant_dmf の docstring の
        命題），そのとき M(t) = ni (1 + 2t + t^2)，R(t) = (ni - 1)(1 + t)。
        つまり対称性を課すと離散モース不等式のずれが ni に比例して悪化する。
      * 一方 Z_ni 不変な離散モースボット関数（高さ関数）は R(t) = 0 のまま。
        collection が 2 つの臨界円周を表し，Σ_C P_t(C) = (1+t) + t(1+t) = P_t(T^2)。
    """
    head = (f"{'T(ni,nj)':>10} | {'DMF 4 cells':>12} | {'inv.DMF #crit':>13} {'R(t)':>10}"
            f" | {'inv.DMBF #C':>11} {'R(t)':>6}")
    print(head)
    print("-" * len(head))
    for ni, nj in sizes:
        K = torus(ni, nj)
        cn = MorseBott(K, canonical_dmf(K)).report()
        iv = MorseBott(K, invariant_dmf(K, ni, nj)).report()
        hb = MorseBott(K, height_fn(K, ni, nj)).report()
        nC = sum(1 for C in hb["reduced_collections"] if poly_trim(betti(C)))
        print(f"{f'T({ni},{nj})':>10} | {poly_str(cn['R_MB']):>12}"
              f" | {iv['n_critical']:>13} {poly_str(iv['R_MB']):>10}"
              f" | {nC:>11} {poly_str(hb['R_MB']):>6}")
    print("\n  列: 「DMF 4 cells」= 非対称な最小 DMF の R(t)，")
    print("      「inv.DMF #crit」= Z_ni 不変な DMF の臨界セル数（最小）と R(t)，")
    print("      「inv.DMBF #C」= Z_ni 不変な DMBF の非自明な reduced collection の個数")
    print("      （= 臨界円周の個数）と R(t)。")


def report_json(K, f, p=0, name=""):
    """機械可読な報告（JSON にできる辞書）。"""
    M = MorseBott(K, f, p=p)
    r = M.report()
    shapes = []
    for C in M.reduced_collections():
        b = poly_trim(M.betti(C))
        if b:
            shapes.append({"cells": len(C),
                           "dims": sorted({K.dim(c) for c in C}),
                           "P_t": b})
    out = {
        "complex": name or getattr(K, "name", ""),
        "cells": K.counts(),
        "field": ("Q" if p == 0 else f"F_{p}"),
        "is_dmf": r["is_dmf"],
        "is_dmb": r["is_dmb"],
        "violations": [{"condition": k, "cell": str(c)} for k, c, _ in r["dmb_violations"]],
        "collections": r["collections"],
        "reduced_collections": len(r["reduced_collections"]),
        "weakly_critical": r["n_weakly_critical"],
        "critical": r["n_critical"],
        "arrows": r["n_arrows"],
        "P_K": r["P_K"],
        "MB_sum": r["MB_sum"],
        "R_MB": r["R_MB"],
        "sharp": r["MB_sharp"],
        "contributing_collections": shapes,
    }
    if r["is_dmf"]:
        out["M"] = r["M"]
        out["R_M"] = r["R_M"]
    return out


def main(argv=None):
    import argparse  # noqa: PLC0415

    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0].endswith(".py"):                # main(sys.argv) でも動くように
        argv = argv[1:]
    ap = argparse.ArgumentParser(
        prog="dmb_core.py",
        description="トーラス（や自分で与えた複体）の上の離散モースボット理論を計算する",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""例:
  python3 dmb_core.py                 T(4,4) 上の 7 種類の関数の報告
  python3 dmb_core.py 5 4             ni=5, nj=4 のトーラスで
  python3 dmb_core.py --table         対称性を課したときの DMT と DMBT の差
  python3 dmb_core.py --json          機械可読な出力
  python3 dmb_core.py --field 2       F_2 係数で（ねじれが見える）
  python3 dmb_core.py --complex my.json --json     自分の複体を読む
""")
    ap.add_argument("ni", nargs="?", type=int, default=4, help="θ 方向の分割数")
    ap.add_argument("nj", nargs="?", type=int, default=4, help="φ 方向の分割数")
    ap.add_argument("--ni", dest="ni_opt", type=int, default=None,
                    help="θ 方向の分割数（位置引数の代わり）")
    ap.add_argument("--nj", dest="nj_opt", type=int, default=None,
                    help="φ 方向の分割数（位置引数の代わり）")
    ap.add_argument("--table", action="store_true",
                    help="対称性を課したときの DMT と DMBT の差の表")
    ap.add_argument("--complex", dest="complex_file", default=None,
                    help="複体（と関数）を JSON から読む。complexes.load_json の形式")
    ap.add_argument("--field", type=int, default=0, metavar="P",
                    help="係数体（0 = 有理数体，素数 p = F_p）")
    ap.add_argument("--json", action="store_true", help="JSON で出力する")
    args = ap.parse_args(argv)

    if args.table:
        symmetry_table([(n, 4) for n in range(3, 9)])
        return 0

    if args.complex_file:
        import json  # noqa: PLC0415

        from complexes import load_json  # noqa: PLC0415
        K, f = load_json(args.complex_file)
        if f is None:
            ap.error("JSON に \"f\"（セルごとの値）が入っていない")
        out = report_json(K, f, args.field, name=args.complex_file)
        if args.json:
            print(json.dumps(out, ensure_ascii=False, indent=1))
        else:
            describe(K, None, None, args.complex_file, f, None, p=args.field)
        return 0 if out["is_dmb"] else 1

    ni = args.ni_opt if args.ni_opt is not None else args.ni
    nj = args.nj_opt if args.nj_opt is not None else args.nj
    if ni < 3 or nj < 3:
        ap.error("ni, nj は 3 以上")
    K = torus(ni, nj)
    names = torus_names(ni, nj)
    catalogue = function_catalogue(K, ni, nj)

    if args.json:
        import json  # noqa: PLC0415

        payload = {"torus": [ni, nj],
                   "functions": {name: report_json(K, f, args.field, f"T({ni},{nj})")
                                 for name, f in catalogue.items()}}
        print(json.dumps(payload, ensure_ascii=False, indent=1))
        return 0

    print(f"トーラス T({ni}, {nj}) の三角形分割"
          + (f"（係数体 F_{args.field}）" if args.field else "") + "\n")
    ok = True
    for title, f in catalogue.items():
        r = describe(K, ni, nj, title, f, names, p=args.field)
        ok &= bool(r["is_dmb"])
        print()
    print("すべての検査を通過" if ok else "検査に失敗したものがある")
    return 0 if ok else 1


def function_catalogue(K, ni, nj):
    """報告に載せる関数の一覧（名前 -> f）。"""
    f_h = height_fn(K, ni, nj)
    out = {
        "定数関数（自明な DMBF）": constant_fn(K),
        "モースボット高さ関数 h(j) = min(j, nj-j)（回転対称）": f_h,
        "同上・細分版 (refine=True)": height_fn(K, ni, nj, True),
        "モースボット高さ関数の Morsification（DMF になる）": morsify(K, f_h),
        "臨界セル 4 個の離散モース関数（tree-cotree, 非対称）": canonical_dmf(K),
        f"回転対称な離散モース関数（臨界セル 4·{ni} 個）": invariant_dmf(K, ni, nj),
        f"矢印 {ni} 本を持つ回転対称で鋭い DMBF（帯を 1 つ持ち上げる）":
            arrowed_dmbf(K, ni, nj),
    }
    if ni == nj:
        try:
            out["dmf.py の calcDMF が与える離散モース関数"] = dmf_from_dmf_py(K, ni, nj)
        except Exception:                                # noqa: BLE001, S110
            pass
    return out


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
