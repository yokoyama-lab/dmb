#!/usr/bin/env python3
"""検査・比較用の有限単体的複体。外部依存なし。

dmb_core.py の理論の計算はトーラスに限らず任意の有限単体的複体で動くので，
ホモロジーが分かっている複体を並べておき，検査に使う。

    python3 complexes.py     # 一覧とベッチ数・オイラー標数を出す
"""

import json
import os
import sys
from itertools import combinations

from dmb_core import Complex, CWComplex, betti, poly_str, torus


def simplex(n):
    """n 単体（可縮）。P_t = 1。"""
    return Complex([tuple(range(n + 1))])


def sphere(n):
    """n 次元球面 = (n+1) 単体の境界。P_t = 1 + t^n。"""
    verts = tuple(range(n + 2))
    return Complex(list(combinations(verts, n + 1)))


def circle(n):
    """n 角形（n ≥ 3）。P_t = 1 + t。"""
    if n < 3:
        raise ValueError("n ≥ 3")
    return Complex([(i, (i + 1) % n) for i in range(n)])


def two_circles(n=3):
    """交わらない 2 つの円周。P_t = 2 + 2t（連結でない複体の検査用）。"""
    a = [(i, (i + 1) % n) for i in range(n)]
    b = [(i + n, (i + 1) % n + n) for i in range(n)]
    return Complex(a + b)


def annulus(ni, nj):
    """円筒 S^1 × [0,1]（ni 周 × nj 段）。P_t = 1 + t。境界を持つ曲面。"""
    facets = []
    for i in range(ni):
        for j in range(nj):
            ip = (i + 1) % ni
            facets.append(((i, j), (i, j + 1), (ip, j + 1)))
            facets.append(((i, j), (ip, j), (ip, j + 1)))
    return Complex(facets)


def moebius():
    """メビウスの帯（5 頂点の最小三角形分割）。P_t = 1 + t，χ = 0，向き付け不可能。"""
    return Complex([(i, (i + 1) % 5, (i + 2) % 5) for i in range(5)])


def klein_bottle(ni, nj):
    """クラインの壺（正方形の i 方向の貼り合わせで j を反転）。P_t = 1 + t（有理数係数）。

    頂点は (i, j) ∈ Z_ni × Z_nj。i = ni での貼り合わせだけ (0, -j) に飛ばす。
    単体的複体になるには ni, nj ≥ 5 程度が必要（`is_closed_surface` で確認できる）。"""
    def v(i, j):
        return (0, (-j) % nj) if i == ni else (i % ni, j % nj)

    facets = []
    for i in range(ni):
        for j in range(nj):
            facets.append((v(i, j), v(i, j + 1), v(i + 1, j + 1)))
            facets.append((v(i, j), v(i + 1, j), v(i + 1, j + 1)))
    return Complex(facets)


def klein_bottle_sym(ni, nj):
    """回転で不変な三角形分割のクラインの壺（各正方形に中心を入れて 4 分割）。

    `klein_bottle` の対角線による分割は，貼り合わせの反転 j → −j と両立しないので
    i 方向の平行移動が自己同型にならない。中心を入れて 4 分割すると反転で不変になり，
    平行移動 `klein_shift` が自己同型になる。

    頂点は ("v", i, j)（格子点）と ("c", i, j)（正方形 (i,j) の中心）。
    セル数は (2·ni·nj, 6·ni·nj, 4·ni·nj)，χ = 0，H_* = (Z, Z ⊕ Z/2, 0)。"""
    def v(i, j):
        return ("v", 0, (-j) % nj) if i == ni else ("v", i % ni, j % nj)

    facets = []
    for i in range(ni):
        for j in range(nj):
            a, b = v(i, j), v(i + 1, j)
            c, d = v(i + 1, j + 1), v(i, j + 1)
            m = ("c", i % ni, j % nj)
            facets += [(p, q, m) for p, q in ((a, b), (b, c), (c, d), (d, a))]
    return Complex(facets)


def klein_shift(cell, ni, nj, times=1):
    """`klein_bottle_sym` の i 方向の平行移動（貼り合わせのところで j を反転）。

    頂点は (ni, j) ~ (0, −j)，正方形は (ni, j) ~ (0, −j−1) に貼り合わさることに注意。
    この写像 g は位数 2·ni（ni 回まわると反転 j → −j が残る）で，**自由ではない**。
    ni が奇数のときに限り，g^2 が生成する Z_ni が自由に作用する
    （反転は位数 2 で，奇数位数の部分群には入らないから）。"""
    for _ in range(times):
        out = []
        for kind, i, j in cell:
            if i == ni - 1:
                out.append((kind, 0, (-j) % nj) if kind == "v"
                           else (kind, 0, (-j - 1) % nj))
            else:
                out.append((kind, i + 1, j))
        cell = tuple(sorted(out))
    return cell


def klein_free_action(ni, nj):
    """ni が奇数のときの自由な Z_ni 作用（g^2 が生成する）。作用を返す。"""
    if ni % 2 == 0:
        raise ValueError("自由になるのは ni が奇数のときだけ（g^ni = 反転 が入るため）")
    return lambda cell, k=1: klein_shift(cell, ni, nj, 2 * k)


def projective_plane():
    """実射影平面 RP^2（6 頂点の最小三角形分割）。有理数係数で P_t = 1，χ = 1。"""
    return Complex([
        (1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 6), (1, 2, 6),
        (2, 3, 5), (3, 4, 6), (4, 5, 2), (5, 6, 3), (6, 2, 4),
    ])


# ------------------------------------------------------------------ 検査補助


def euler(K):
    return sum((-1) ** k * v for k, v in K.counts().items())


def boundary_edges(K):
    """ちょうど 1 つの三角形に含まれる辺（2 次元複体の境界）。"""
    return [e for e in K.cells_of_dim(1) if len(K.above[e]) == 1]


def is_closed_surface(K):
    """閉曲面の組合せ的条件: どの辺もちょうど 2 つの三角形に含まれ，
    どの頂点のリンクも 1 本の閉路（＝頂点まわりの三角形が輪になる）。"""
    if any(K.dim(c) != 2 for c in K.cells if not K.above[c]):
        return False
    for e in K.cells_of_dim(1):
        if len(K.above[e]) != 2:
            return False
    for v in K.cells_of_dim(0):
        star = [t for t in K.cells_of_dim(2) if set(v) < set(t)]
        link = [tuple(sorted(set(t) - set(v))) for t in star]
        deg = {}
        for a, b in link:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        if not deg or any(d != 2 for d in deg.values()) or len(link) != len(deg):
            return False
    return True


# ------------------------------------------------------------ JSON 入出力


def _to_cell(x):
    """JSON の配列を tuple に戻す（頂点が (i, j) のような組でもよいように）。"""
    if isinstance(x, list):
        return tuple(_to_cell(v) for v in x)
    return x


def load_json(path_or_obj):
    """JSON から複体（と関数）を読む。自分の複体を持ち込むための入口。

    セルの名前は数・文字列のほか，(i, j) のような配列でもよい（JSON の配列は
    tuple に戻す）。単体的複体:
        {"type": "simplicial",
         "facets": [[0,1,2], [0,1,3]],
         "f": [[[0,1,2], 3], [[0], 0], ...]}     # [セル, 値] の対の並び。省略可

    CW 複体（正則とは限らない。セルは文字列の名前）:
        {"type": "cw",
         "dims": {"v": 0, "a": 1, "F": 2},
         "faces": [["v","a"], ["a","F"], ["v","F"]],
         "regular": [["v","a"]],                 # 省略すると全部 irregular
         "incidence": [["v","a",0], ["a","F",2]],
         "f": {"v": 0, "a": 1, "F": 2}}          # 辞書でも対の並びでもよい

    戻り値は (複体, f または None)。値が一部のセルにしか無ければエラーにする
    （黙って 0 を埋めると間違った結論が出るので）。"""
    if isinstance(path_or_obj, (str, bytes, os.PathLike)):
        with open(path_or_obj, encoding="utf-8") as fh:
            obj = json.load(fh)
    else:
        obj = path_or_obj
    kind = obj.get("type", "simplicial")
    if kind == "simplicial":
        K = Complex([_to_cell(fc) for fc in obj["facets"]])
    elif kind == "cw":
        raw_dims = obj["dims"]
        dims = {_to_cell(c): d for c, d in
                (raw_dims.items() if isinstance(raw_dims, dict) else raw_dims)}
        K = CWComplex(dims,
                      [(_to_cell(a), _to_cell(b)) for a, b in obj.get("faces", [])],
                      regular=[(_to_cell(a), _to_cell(b))
                               for a, b in obj.get("regular", [])],
                      incidence={(_to_cell(a), _to_cell(b)): v
                                 for a, b, v in obj.get("incidence", [])},
                      name=obj.get("name", ""))
    else:
        raise ValueError(f"未知の type: {kind}")

    raw = obj.get("f")
    if raw is None:
        return K, None
    pairs = raw.items() if isinstance(raw, dict) else raw
    known = set(K.cells)
    f = {}
    for cell_, val in pairs:
        c = _to_cell(cell_)
        if c not in known:
            raise ValueError(f"複体に無いセル: {cell_!r}")
        f[c] = val
    missing = [c for c in K.cells if c not in f]
    if missing:
        raise ValueError(f"値が与えられていないセルがある（{len(missing)} 個。"
                         f"例: {missing[0]!r}）")
    return K, f


def _to_json(x):
    """tuple を JSON の配列にする（_to_cell の逆）。"""
    if isinstance(x, tuple):
        return [_to_json(v) for v in x]
    return x


def dump_json(K, f=None, path=None):
    """複体（と関数）を JSON にする。load_json で読み戻せる。"""
    if isinstance(K, CWComplex):
        # セル名が全部文字列なら読みやすい辞書に，そうでなければ対の並びにする
        dims = ({c: K.dim(c) for c in K.cells} if all(isinstance(c, str) for c in K.cells)
                else [[_to_json(c), K.dim(c)] for c in K.cells])
        obj = {"type": "cw", "name": K.name, "dims": dims,
               "faces": [[_to_json(a), _to_json(b)]
                         for a, b in sorted(K.faces, key=repr)],
               "regular": [[_to_json(a), _to_json(b)]
                           for a, b in sorted(K.regular_pairs, key=repr)],
               "incidence": [[_to_json(a), _to_json(b), K.incidence(a, b)]
                             for a in K.cells for b in K.above[a]]}
    else:
        obj = {"type": "simplicial",
               "facets": [_to_json(c) for c in K.cells if not K.above[c]]}
    if f is not None:
        obj["f"] = [[_to_json(c), f[c]] for c in K.cells]
    text = json.dumps(obj, ensure_ascii=False, indent=1)
    if path:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
    return text


CATALOGUE = {
    "point (Δ^0)": (lambda: simplex(0), [1]),
    "segment (Δ^1)": (lambda: simplex(1), [1]),
    "triangle (Δ^2, 円板)": (lambda: simplex(2), [1]),
    "S^1 (3 頂点)": (lambda: circle(3), [1, 1]),
    "S^1 (7 頂点)": (lambda: circle(7), [1, 1]),
    "S^1 ⊔ S^1": (two_circles, [2, 2]),
    "S^2 = ∂Δ^3": (lambda: sphere(2), [1, 0, 1]),
    "S^3 = ∂Δ^4": (lambda: sphere(3), [1, 0, 0, 1]),
    "円筒 (5×2)": (lambda: annulus(5, 2), [1, 1]),
    "メビウスの帯": (moebius, [1, 1]),
    "RP^2 (6 頂点)": (projective_plane, [1]),
    "クラインの壺 (5×5)": (lambda: klein_bottle(5, 5), [1, 1]),
    "クラインの壺・対称 (3×4)": (lambda: klein_bottle_sym(3, 4), [1, 1]),
    "トーラス T(4,4)": (lambda: torus(4, 4), [1, 2, 1]),
}


def main(argv=None):
    del argv                      # 引数は取らない（署名を他のスクリプトと揃えるため）
    ok = True
    print(f"{'複体':<22} {'セル数':<22} {'χ':>3}  {'P_t':<14} {'期待と一致'}")
    for name, (build, expected) in CATALOGUE.items():
        K = build()
        b = betti(K.cells)
        match = b == expected
        ok &= match
        print(f"{name:<22} {str(K.counts()):<22} {euler(K):>3}  "
              f"{poly_str(b):<14} {'OK' if match else 'MISMATCH ' + poly_str(expected)}")
    print()
    print("すべて期待どおり" if ok else "一致しないものがある")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
