#!/usr/bin/env python3
"""トーラス上の Z_ni 不変な離散モースボット関数の全数探索。外部依存なし。

Z_ni（θ 方向の回転）はトーラス T(ni, nj) のセルに自由に作用するので，不変な関数は
軌道（6·nj 個）の上の関数と同じである。値の範囲を区切れば全数探索できる。

何が分かるか:
  * 不変な DMBF がどれだけあるか，そのうち R(t) = 0（鋭い）ものはどれか
  * 鋭いものの reduced collection の形（＝「臨界円周」の並び方）の分類
  * 「不変な DMF は 4·ni 個の臨界セルが要る」のに対し，不変な DMBF は
    ずっと少ない情報で鋭くなる，という強み 2 の全数版

    python3 search.py                    # T(4,3), 値 3 通りで全数探索
    python3 search.py --ni 5 --nj 3 --values 3
    python3 search.py --nj 4 --values 2 --sharp-only
"""

import argparse
import sys
from collections import Counter

import dmb_core as core

TYPE_ORDER = {1: 0, 2: 1, 3: 2}       # 頂点 → 辺 → 三角形


class SearchTruncated(Exception):
    """--max-nodes に当たって全探索を終えられなかった。

    「全部見た」と「途中で止めた」を取り違えないよう，例外にして区別する。"""

    def __init__(self, nodes, n_dmb, n_found):
        super().__init__(f"節点 {nodes} で打ち切り（DMBF {n_dmb} 個まで）")
        self.nodes, self.n_dmb, self.n_found = nodes, n_dmb, n_found


def invariant_orbits(K, ni, nj):
    """Z_ni の軌道を，探索に都合のよい順（下の準位から，低い次元から）に並べる。"""
    groups = {}
    for c in K.cells:
        groups.setdefault(core.orbit_key(c, ni, nj), []).append(c)

    def sort_key(item):
        key, _cells = item
        return (base_level(key, nj), TYPE_ORDER[len(key)], orbit_label(key, nj))

    return [cells for _, cells in sorted(groups.items(), key=sort_key)]


def prune_predicate(K, orbits, value):
    """割り当て済みのセルだけで見て (MB2)/(MB4) が既に破れているか。

    値を足しても U^snc / D^snc は減らないので，これが真なら以後どう埋めても
    離散モースボット関数にならない ＝ 枝を捨てて safe（健全性は
    `survives_pruning` を使って検査してある）。"""

    def bad(t):
        touched = set(orbits[t])
        for c in orbits[t]:
            touched.update(K.above[c])
            touched.update(K.below[c])
        for c in touched:
            if c not in value:
                continue
            if sum(1 for d in K.above[c] if d in value and value[d] < value[c]) > 1:
                return True
            if sum(1 for d in K.below[c] if d in value and value[d] > value[c]) > 1:
                return True
        return False

    return bad


def survives_pruning(ni, nj, vals):
    """軌道ごとの値 `vals` が枝刈りを通り抜けるか（枝刈りの健全性検査用）。

    離散モースボット関数に対応する `vals` は必ず True になるはずである。"""
    K = core.torus(ni, nj)
    orbits = invariant_orbits(K, ni, nj)
    value = {}
    bad = prune_predicate(K, orbits, value)
    for t, cells in enumerate(orbits):
        for c in cells:
            value[c] = vals[t]
        if bad(t):
            return False
    return True


def search(ni, nj, nvalues=3, sharp_only=False, limit=None, progress=None,
           sink=None, max_nodes=None):
    """不変な離散モースボット関数を全数探索する。

    値は 0..nvalues-1。値の平行移動は自明なので探索空間から外していないが，
    順序しか効かないので nvalues を変えると別の族が出ることに注意。

    戻り値は (見つかった関数の一覧, 調べた節点数, DMBF の個数)。一覧の各要素は
    (軌道ごとの値のタプル, report の辞書)。

    sink を渡すと，一覧に貯める代わりに sink(vals, rep) を呼ぶ（一覧は空のまま
    返る）。値の個数を増やすと見つかる関数が数百万になり，全部を持つと
    メモリが尽きるので，集計だけしたいときに使う。"""
    K = core.torus(ni, nj)
    orbits = invariant_orbits(K, ni, nj)
    n = len(orbits)
    value = {}                                   # cell -> 値（割り当て済みのみ）
    found, nodes, n_dmb, n_found = [], 0, 0, 0

    locally_bad = prune_predicate(K, orbits, value)

    pk = core.poly_trim(core.betti(K.cells))       # P_t(K) は毎回同じなので 1 度だけ
    reps = [cells[0] for cells in orbits]          # 各軌道の代表元

    def is_dmb_invariant(M):
        """不変な f の (MB2)/(MB4) 検査。軌道の代表元だけ見れば十分。

        回転は f を保つ複体の自己同型なので U^snc(gσ) = U^snc(σ)。
        全セルを見る is_dmb() と同値であることは検査で確かめてある。"""
        for c in reps:
            if len(M.up_snc(c)) > 1 or len(M.down_snc(c)) > 1:
                return False
        return True

    def record():
        nonlocal n_dmb, n_found
        M = core.MorseBott(K, dict(value))
        if not is_dmb_invariant(M):                # 高い report() の前に安く弾く
            return
        n_dmb += 1
        mb = M.morse_bott_polynomial()
        rt = core.poly_div_1pt(core.poly_sub(mb, pk))
        sharp = rt is not None and not core.poly_trim(rt)
        if sharp_only and not sharp:
            return
        # 矢印は軌道ごとに同じ本数なので，代表元だけ数えて ni 倍すればよい
        n_arrows = ni * sum(len(M.up_snc(c)) for c in reps)
        vals = tuple(value[orbits[t][0]] for t in range(n))
        rep = {"MB_sum": mb, "R_MB": rt, "MB_sharp": sharp, "P_K": pk,
               "shape": shape_of(M), "n_arrows": n_arrows}
        n_found += 1
        if sink is None:
            found.append((vals, rep))
        else:
            sink(vals, rep)

    def rec(t):
        nonlocal nodes
        if limit is not None and n_found >= limit:
            return
        if max_nodes is not None and nodes > max_nodes:
            raise SearchTruncated(nodes, n_dmb, n_found)
        if t == n:
            record()
            return
        for v in range(nvalues):
            nodes += 1
            for c in orbits[t]:
                value[c] = v
            if not locally_bad(t):
                if progress and nodes % progress == 0:
                    print(f"  ... 節点 {nodes}，見つかった {n_found}", file=sys.stderr)
                rec(t + 1)
            for c in orbits[t]:
                del value[c]

    rec(0)
    return found, nodes, n_dmb


def shape_of(M):
    """非自明な reduced collection の（P_t, セル数）の多重集合。臨界円周の並び。"""
    out = []
    for C in M.reduced_collections():
        b = tuple(core.poly_trim(M.betti(C)))
        if b:
            out.append(b)
    return tuple(sorted(out))


def shape_str(shape):
    if not shape:
        return "（寄与する collection なし）"
    counts = Counter(shape)
    return "  ".join(f"({core.poly_str(list(p))})"
                     + (f" ×{k}" if k > 1 else "")
                     for p, k in sorted(counts.items()))


def classify(found):
    """見つかった関数を「臨界円周の並び」と R(t) で分類する。"""
    shapes = Counter()
    for _, r in found:
        shapes[(r["shape"], tuple(core.poly_trim(r["R_MB"] or [])))] += 1
    return shapes


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Z_ni 不変な離散モースボット関数の全数探索",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--ni", type=int, default=4)
    ap.add_argument("--nj", type=int, default=3)
    ap.add_argument("--values", type=int, default=3, help="値の個数（既定 3）")
    ap.add_argument("--sharp-only", action="store_true", help="R(t) = 0 のものだけ")
    ap.add_argument("--limit", type=int, default=None, help="見つける個数の上限")
    ap.add_argument("--progress", type=int, default=None, help="節点数ごとに進捗を出す")
    ap.add_argument("--show", type=int, default=3, help="代表例をいくつ表示するか")
    ap.add_argument("--max-nodes", type=int, default=None,
                    help="この節点数を超えたら打ち切る（打ち切りはそう明記する）")
    args = ap.parse_args(argv)

    if args.ni < 3 or args.nj < 3:
        ap.error("ni, nj は 3 以上")
    K = core.torus(args.ni, args.nj)
    orbits = invariant_orbits(K, args.ni, args.nj)
    print(f"トーラス T({args.ni}, {args.nj}): セル {len(K.cells)} 個，"
          f"Z_{args.ni} の軌道 {len(orbits)} 個（各 {args.ni} セル）")
    print(f"値は 0..{args.values - 1} の {args.values} 通り "
          f"→ 素朴には {args.values}^{len(orbits)} 通り\n")

    # 見つかった関数を貯めずに集計する（値を増やすと数百万になるため）
    tally = Counter()          # (shape, R(t), 矢印あり) -> 個数
    sharp_tally = Counter()    # 鋭いものだけ: (shape, 矢印あり) -> 個数
    n_sharp = n_arrowed = n_sharp_arrowed = 0
    examples = []              # 鋭くて矢印のあるもの，鋭いもの，の順に少しだけ

    def sink(vals, rep):
        nonlocal n_sharp, n_arrowed, n_sharp_arrowed
        arrowed = rep["n_arrows"] > 0
        rt = None if rep["R_MB"] is None else tuple(core.poly_trim(rep["R_MB"]))
        tally[(rep["shape"], rt, arrowed)] += 1
        if rep["MB_sharp"]:
            sharp_tally[(rep["shape"], arrowed)] += 1
        n_sharp += rep["MB_sharp"]
        n_arrowed += arrowed
        if rep["MB_sharp"] and arrowed:
            n_sharp_arrowed += 1
            if len(examples) < args.show:
                examples.append((vals, rep))

    truncated = None
    try:
        _, nodes, n_dmb = search(args.ni, args.nj, args.values,
                                 sharp_only=args.sharp_only, limit=args.limit,
                                 progress=args.progress, sink=sink,
                                 max_nodes=args.max_nodes)
    except SearchTruncated as exc:
        truncated, nodes, n_dmb = exc, exc.nodes, exc.n_dmb

    naive = args.values ** len(orbits)
    print(f"調べた節点 {nodes}（枝刈り。素朴な全数は {naive} 通り）")
    if truncated:
        print("!! --max-nodes で打ち切った。以下は全数ではない")
    else:
        print("全探索を終えた（打ち切りなし）")
    print(f"不変な離散モースボット関数: {n_dmb} 個")
    print(f"  そのうち R(t) = 0（鋭い）        : {n_sharp} 個")
    print(f"  そのうち矢印を持つ               : {n_arrowed} 個")
    print(f"  鋭くて，かつ矢印を持つ           : {n_sharp_arrowed} 個")

    print("\n「臨界円周の並び」ごとの個数"
          "（非自明な reduced collection の P_t を並べたもの）:")
    for (shape, r, arrowed), cnt in tally.most_common(12):
        rt = "割り切れない" if r is None else core.poly_str(list(r))
        print(f"  {cnt:>8} 個   {shape_str(shape):<44} R(t) = {rt:<10} "
              f"矢印 {'あり' if arrowed else 'なし'}")

    if sharp_tally:
        print("\n鋭いもの（R(t) = 0）だけの内訳:")
        for (shape, arrowed), cnt in sharp_tally.most_common():
            print(f"  {cnt:>8} 個   {shape_str(shape):<44} "
                  f"矢印 {'あり' if arrowed else 'なし'}")
        print(f"  鋭いものの「臨界円周の並び」は {len({s for s, _ in sharp_tally})} 種類")

    if examples:
        print(f"\n鋭くて矢印を持つ例（先頭 {len(examples)} 件）の軌道ごとの値:")
        names = [orbit_label(orbits[t][0], args.nj) for t in range(len(orbits))]
        print("    " + "  ".join(f"{s:>6}" for s in names))
        for vals, rep in examples:
            print("    " + "  ".join(f"{v:>6}" for v in vals)
                  + f"   矢印 {rep['n_arrows']} 本")
    return 2 if truncated else 0


def base_level(c, nj):
    """セルが乗っている「下の準位」。j の集合が {j} か {j, j+1 mod nj} になることを使う。"""
    js = {v[1] for v in c}
    if len(js) == 1:
        return next(iter(js))
    a, b = sorted(js)
    return a if (a + 1) % nj == b else b


def orbit_label(c, nj):
    """軌道の見出し（型と準位）。v=頂点, ie=横の辺, je=縦の辺, dg=斜めの辺,
    tu=上三角, tl=下三角。"""
    j = base_level(c, nj)
    if len(c) == 1:
        return f"v{j}"
    if len(c) == 2:
        if len({v[1] for v in c}) == 1:
            return f"ie{j}"
        return f"je{j}" if len({v[0] for v in c}) == 1 else f"dg{j}"
    mult = Counter(v[1] for v in c)[j]
    return f"{'tl' if mult == 2 else 'tu'}{j}"


if __name__ == "__main__":
    raise SystemExit(main())
