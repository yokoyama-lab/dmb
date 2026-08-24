#!/usr/bin/env python3
"""離散モースボット理論の例集（読んで理解するための小さい例と，DMBT ならではの強み）。

    python3 examples.py            # 全部
    python3 examples.py tutorial   # 手で確かめられる小さい例だけ
    python3 examples.py strength   # 離散モース理論にはできないことだけ

各例は最後に自分で検算する（食い違えば終了コード 1）。
"""

import itertools
import sys

import complexes as X
import dmb_core as D

RULE = "─" * 78


def cell_name(c):
    """セルの表示名。単体は頂点を並べ，CW 複体のセルはその ID をそのまま。"""
    if isinstance(c, tuple):
        return "".join(str(v) for v in c)
    return str(c)


def table(K, f, names=None, title=""):
    """セルごとに U, D, U^snc, D^snc と分類を並べる。定義を手で追うための表。"""
    M = D.MorseBott(K, f)
    index = {}
    for k, L in enumerate(sorted(M.collections(), key=lambda L: (-len(L), L[0]))):
        for c in L:
            index[c] = k
    if title:
        print(f"  {title}")
    print(f"    {'cell':<14}{'f':>3} {'U':>3}{'D':>3} {'Usnc':>5}{'Dsnc':>5}  "
          f"{'臨界':<6}{'弱臨界':<8}{'collection'}")
    for c in K.cells:
        nm = names[c] if names else cell_name(c)
        print(f"    {nm:<14}{f[c]:>3} {len(M.up_nc(c)):>3}{len(M.down_nc(c)):>3} "
              f"{len(M.up_snc(c)):>5}{len(M.down_snc(c)):>5}  "
              f"{'○' if M.is_critical(c) else '·':<6}"
              f"{'○' if M.is_weakly_critical(c) else '·':<8}{index[c]}")
    r = M.report()
    print(f"    ⇒ DMF: {'はい' if r['is_dmf'] else 'いいえ'} / "
          f"DMBF: {'はい' if r['is_dmb'] else 'いいえ'} / "
          f"collection {[len(L) for L in sorted(M.collections(), key=len, reverse=True)]} / "
          f"Σ_C P_t(C) = {D.poly_str(r['MB_sum'])} / P_t(K) = {D.poly_str(r['P_K'])} / "
          f"R(t) = {D.poly_str(r['R_MB']) if r['R_MB'] is not None else '割り切れない'}")
    return r


# ============================================================ 理解のための例


def ex1_one_edge():
    print(RULE)
    print("例 1: 1 本の辺（頂点 0, 1 と辺 01）。同じ複体・3 通りの関数で違いを見る。")
    print("  U(σ)      = 値が σ 以下の余次元 1 の coface の数    （離散モース理論）")
    print("  U^snc(σ)  = 値が σ 未満の余次元 1 の coface の数    （離散モースボット理論）")
    print("  違いは「値が等しいものを数えるかどうか」だけ。ここが collection を生む。\n")
    K = D.Complex([(0, 1)])
    v0, v1, e = (0,), (1,), (0, 1)
    ok = True

    r = table(K, {v0: 0, v1: 1, e: 2}, title="(a) 単射: f = 0, 1, 2 → 全部が臨界")
    ok &= r["is_dmf"] and r["M"] == [2, 1]

    print()
    r = table(K, {v0: 0, v1: 1, e: 1},
              title="(a') 単射でない DMF: f(v1) = f(e) = 1 → {v1, e} が対になる")
    ok &= r["is_dmf"] and r["M"] == [1] and D.poly_trim(r["R_MB"]) == []
    print("      v1 は U = 1（e の値が等しい）だが U^snc = 0 なので weakly critical。")
    print("      これが Theorem 3.2 の #L = 2 で reduced な collection。")

    print()
    r = table(K, {v0: 0, v1: 0, e: 0}, title="(b) 定数: collection は複体全体 1 つ")
    ok &= r["is_dmb"] and not r["is_dmf"] and r["MB_sum"] == [1]
    print("      定数関数はいつでも DMBF で，Σ_C P_t(C) = P_t(K)（R = 0）。")
    print("      鋭くはなるが collection が大きすぎて何も分解していない —")
    print("      「R = 0 かどうか」だけでなく「collection がどれだけ細かいか」が大事。")
    return ok


def ex2_theorem32():
    print(RULE)
    print("例 2: DMF ⊊ DMBF（Theorem 3.2）。パス 0 - 1 - 2 の上で。")
    print("  Theorem 3.2: f が DMF ⇔ f が DMBF かつ各 collection L が")
    print("               (i) #L = 1，または (ii) #L = 2 かつ L が reduced。\n")
    K = D.Complex([(0, 1), (1, 2)])
    v0, v1, v2, e01, e12 = (0,), (1,), (2,), (0, 1), (1, 2)
    ok = True

    r = table(K, {v0: 0, v1: 1, v2: 0, e01: 1, e12: 0},
              title="(a) #L = 2 だが reduced でない: L = {v1, e01}，v1 は U^snc = 1")
    ok &= r["is_dmb"] and not r["is_dmf"]
    print("      DMBF だが DMF ではない（U(v1) = 2: e01 は値が等しく，e12 は値が小さい）。")
    print("      ⇒ (ii) の「reduced」は落とせない。")

    print()
    r = table(K, {v0: 0, v1: 0, v2: 0, e01: 0, e12: 0},
              title="(b) #L = 5 の collection（定数関数）: DMBF だが DMF ではない")
    ok &= r["is_dmb"] and not r["is_dmf"]
    return ok


def ex3_mb2_violation():
    print(RULE)
    print("例 3: (MB2) を破る例。頂点 v から値の小さい辺が 2 本出ていると U^snc = 2。")
    print("  論文の Example（値 2 のセルへ 2 本の辺が出る図）と同じ形。\n")
    K = D.Complex([(0, 1), (0, 2)])
    f = {(0,): 1, (1,): 0, (2,): 0, (0, 1): 0, (0, 2): 0}
    M = D.MorseBott(K, f)
    table(K, f, title="f(v0) = 1，他は 0")
    bad = M.dmb_violations()
    print(f"      違反: {[(k, ''.join(map(str, c))) for k, c, _ in bad]}")
    print("      U^snc(v0) = 2 なので (MB2) を破る。片方の辺の値を v0 と等しくすれば DMBF になる。")
    f2 = dict(f)
    f2[(0, 2)] = 1
    ok = M.dmb_violations() and D.MorseBott(K, f2).is_dmb()
    print(f"      f(e02) を 1 に上げると DMBF: {D.MorseBott(K, f2).is_dmb()}")
    return bool(ok)


def ex4_critical_circle():
    print(RULE)
    print("例 4: 一番小さい「臨界円周」。S^1（3 頂点）に定数関数。\n")
    K = X.circle(3)
    r = table(K, D.constant_fn(K), title="collection は S^1 全体，P_t(C) = 1 + t")
    print("      滑らかな Morse-Bott 理論で臨界部分多様体が S^1 のとき，")
    print("      その寄与は t^λ P_t(S^1) = t^λ (1 + t)。その離散版がこれ。")
    print()
    K2 = X.circle(6)
    f = D.max_extension(K2, lambda v: min(v, 6 - v))
    r2 = table(K2, f, names=None,
               title="S^1（6 頂点）の高さ関数: 臨界点 2 つ（極小・極大）に分かれる")
    print("      collection は 4 つだが，寄与するのは極小の点と極大の辺の 2 つ:")
    print("      Σ_C P_t(C) = 1 + t = P_t(S^1)。")
    return r["MB_sum"] == [1, 1] and r2["MB_sum"] == [1, 1]


def ex5_morsification():
    print(RULE)
    print("例 5: Morsification f'(σ) = (dim K + 1) f(σ) + dim σ。")
    print("  DMBF を DMF に摂動する。臨界セルはちょうど ⋃ reduced collection。\n")
    K = X.circle(3)
    f = D.constant_fn(K)
    M, Y = D.MorseBott(K, f), D.MorseBott(K, D.morsify(K, f))
    table(K, D.morsify(K, f), title="S^1 の定数関数を Morsify すると全セルが臨界")
    red = {c for C in M.reduced_collections() for c in C}
    print(f"      ⋃ reduced collection = {len(red)} セル，"
          f"f' の臨界セル = {len(Y.critical())} セル → 一致: {red == set(Y.critical())}")
    print("      DMT では M(t) = 3 + 3t で R = 2（鈍い）。DMBT では collection ごとに")
    print("      束ね直して Σ_C P_t(C) = 1 + t で R = 0（鋭い）。")
    print("      ＝ Morse-Bott 理論の御利益は「束ね直し」そのものである。")
    return red == set(Y.critical()) and M.morse_bott_polynomial() == [1, 1]


def ex6_cw_complexes():
    print(RULE)
    print("例 6: 正則とは限らない CW 複体と，条件 (M1)(M3)。")
    print("  単体的複体では face はすべて regular なので (M1)(M3) は自動的に成り立つ。")
    print("  一般の CW 複体には irregular な face があり，そこでは値が真に増えることを")
    print("  要求する。これが (M1)(M3) で，**任意の余次元**に課される点が要点。\n")
    ok = True

    print("  最小 CW 構造では，セル数がベッチ数と一致するほど小さくできる:")
    for build in (D.cw_circle_minimal, D.cw_sphere2_minimal,
                  D.cw_torus_minimal, D.cw_projective_plane_minimal):
        K = build()
        print(f"    {K.name:<16} セル数 {K.counts()}  ∂∘∂ = 0: "
              f"{'はい' if not K.check_boundary() else 'いいえ'}  "
              f"P_t(K) = {D.poly_str(K.betti())}")
        ok &= not K.check_boundary()

    print("\n  最小 CW 構造では face が全部 irregular なので，(M1)(M3) が f を次元に沿って")
    print("  真に増やす。すると U^snc = D^snc = 0 で全セルが weakly critical になり，")
    print("  Σ_C P_t(C) = Σ_σ t^{dim σ}（セルの数え上げ）になる:")
    for build, sharp in ((D.cw_torus_minimal, True), (D.cw_sphere2_minimal, True),
                         (D.cw_projective_plane_minimal, False)):
        K = build()
        f = {c: K.dim(c) for c in K.cells}          # 次元に沿って増える関数
        r = D.MorseBott(K, f).report()
        print(f"    {K.name:<16} DMF: {'○' if r['is_dmf'] else '×'}  "
              f"DMBF: {'○' if r['is_dmb'] else '×'}  "
              f"Σ_C P_t(C) = {D.poly_str(r['MB_sum'])}  "
              f"P_t(K) = {D.poly_str(r['P_K'])}  R(t) = {D.poly_str(r['R_MB'])}")
        ok &= r["is_dmb"] and (r["MB_sharp"] == sharp)
    print("    トーラスと S^2 はセル数がベッチ数と等しいので R = 0（完全）。")
    print("    RP^2 だけ R = t ≠ 0 — 有理数係数で b_2 = 0 なのに m_2 = 1 だから（強み 1）。")

    print("\n  (M1)(M3) を「余次元 1 の irregular facet だけ」と読むと壊れる:")
    K = D.cw_sphere2_minimal()
    f = {"v": 1, "e": 0}
    M = D.MorseBott(K, f)
    print(f"    S^2 の最小 CW 構造には余次元 1 の組が 1 つも無い（"
          f"{[(a, b) for a in K.cells for b in K.above[a]]}）。")
    print(f"    f(v) = 1 > 0 = f(e) は，v1 の字義通りの読みでは DMF: "
          f"{M.is_dmf(strong=False)}")
    print(f"    しかし v は e の余次元 2 の irregular な face なので DMBF ではない: "
          f"{M.is_dmb()}")
    print(f"    違反 = {[(k, cell_name(c)) for k, c, _ in M.dmb_violations()]}")
    print("    ⇒ これが arXiv:2511.07864 v1 の Lemma 3.1（DMF ⇒ DMBF）の反例。")
    print("      v2 の Definition 12 は任意余次元に強めてあるので，この f は DMF でもなくなり")
    print(f"      （v2 の読みで DMF: {M.is_dmf(strong=True)}），齟齬は消える。")
    ok &= M.is_dmf(strong=False) and not M.is_dmb() and not M.is_dmf(strong=True)
    return ok


# ============================================================ DMBT の強み


def st1_nonorientable():
    print(RULE)
    print("強み 1: 閉非向き付け曲面では，離散モース不等式は**決して**等号にならないが，")
    print("        離散モースボット不等式は等号になる。\n")
    print("  理由（証明）: 閉曲面は H_2(K; Z/2) ≠ 0 なので，どんな離散モース関数でも")
    print("  m_2 ≥ 1。一方 RP^2・クラインの壺は有理数係数で b_2 = 0。よって")
    print("  M(t) ≠ P_t(K)，すなわち R(t) ≠ 0 が常に成り立つ。")
    print("  （最小 CW 構造でも同じ: RP^2 は 3 セルだが m_2 = 1 が要る。例 6 参照。）\n")
    ok = True
    for name, K in (("RP^2", X.projective_plane()), ("クラインの壺", X.klein_bottle(5, 5))):
        r = D.MorseBott(K, D.canonical_dmf(K)).report()
        print(f"  {name}: P_t(K) = {D.poly_str(r['P_K'])},  "
              f"最小級の DMF は M(t) = {D.poly_str(r['M'])} → R(t) = {D.poly_str(r['R_M'])} ≠ 0")
        ok &= bool(D.poly_trim(r["R_M"]))
        c = D.MorseBott(K, D.constant_fn(K)).report()
        print(f"        定数関数（DMBF）: Σ_C P_t(C) = {D.poly_str(c['MB_sum'])} "
              f"→ R(t) = {D.poly_str(c['R_MB'])}")
        ok &= c["MB_sharp"]

    print("\n  自明な定数関数だけではない。RP^2 の頂点に 0..2 の高さを与えて最大値で")
    print("  拡張した 3^6 = 729 個の関数のうち，鋭い (R = 0) ものを数えると:")
    K = X.projective_plane()
    verts = sorted(v[0] for v in K.cells_of_dim(0))
    sharp = []
    for vals in itertools.product(range(3), repeat=len(verts)):
        h = dict(zip(verts, vals))
        M = D.MorseBott(K, D.max_extension(K, h.__getitem__))
        rr = M.report()
        if rr["is_dmb"] and rr["MB_sharp"]:
            sharp.append((len(M.collections()), h))
    best = max(sharp, key=lambda t: t[0])[1]
    print(f"        鋭い DMBF: {len(sharp)} 個 / 729（collection 数の最大は "
          f"{max(n for n, _ in sharp)}）")
    print(f"        例: h = {best}")
    M = D.MorseBott(K, D.max_extension(K, best.__getitem__))
    for C in M.reduced_collections():
        print(f"          C: {len(C):>2} セル（次元 {sorted({K.dim(c) for c in C})}）"
              f"  P_t(C) = {D.poly_str(D.betti(C))}")
    ok &= len(sharp) > 0
    print("\n  この差は**係数体に依る**。ねじれ Z/2 は F_2 係数ではベッチ数として見えるので，")
    print("  F_2 上では同じ離散モース関数がちょうど鋭くなる:")
    print(f"    {'複体':<14}{'係数':<7}{'P_t(K)':<16}{'M(t)':<16}{'R_DMF':<8}{'R_DMBF'}")
    for name, K in (("RP^2", X.projective_plane()),
                    ("クラインの壺", X.klein_bottle(5, 5)),
                    ("T(4,4)", D.torus(4, 4))):
        for q, lab in ((0, "Q"), (2, "F_2"), (3, "F_3")):
            a = D.MorseBott(K, D.canonical_dmf(K), p=q).report()
            b = D.MorseBott(K, D.constant_fn(K), p=q).report()
            print(f"    {name:<14}{lab:<7}{D.poly_str(a['P_K']):<16}"
                  f"{D.poly_str(a['M']):<16}{D.poly_str(a['R_M']):<8}"
                  f"{D.poly_str(b['R_MB'])}")
            ok &= b["MB_sharp"]
    print("    整数係数のホモロジー（Smith 標準形で計算）:")
    for name, K in (("RP^2", X.projective_plane()),
                    ("クラインの壺", X.klein_bottle(5, 5))):
        print(f"      {name:<12} {D.homology_str(D.homology_z(K.cells))}")

    print("\n  ⇒ ねじれのある空間では，ねじれの位数を割らない標数の係数体で")
    print("     DMT の不等式は原理的に埋まらない差を持つ（m_2 ≥ 1 は係数に依らないのに")
    print("     b_2 = 0 になるため）。DMBT は collection の（制限された）ホモロジーで")
    print("     その差を吸収するので，**どの係数体でも**鋭くできる。")
    return ok


def st2_symmetry():
    print(RULE)
    print("強み 2: 対称性を課すと DMT は悪化するが，DMBT は悪化しない。")
    print("        トーラス T(n, 4) に θ 方向の回転 Z_n（セルに自由に作用）。\n")
    print("  命題: Z_n 不変な DMF の臨界セルは 4n 個以上（臨界セル集合が自由軌道の和 +")
    print("  離散モース不等式 + χ = 0）。この下限は達成でき，そのとき R = (n-1)(1+t)。\n")
    print(f"  {'n':>3} {'非対称な最小 DMF':>18} "
          f"{'Z_n 不変な DMF':>22} {'Z_n 不変な DMBF':>24}")
    ok = True
    for n in (3, 4, 5, 6, 7):
        K = D.torus(n, 4)
        cn = D.MorseBott(K, D.canonical_dmf(K)).report()
        iv = D.MorseBott(K, D.invariant_dmf(K, n, 4)).report()
        hb = D.MorseBott(K, D.height_fn(K, n, 4)).report()
        ncirc = sum(1 for C in hb["reduced_collections"] if D.poly_trim(D.betti(C)))
        col_a = "臨界 4, R = 0"
        col_b = "臨界 {}, R = {}".format(iv["n_critical"], D.poly_str(iv["R_MB"]))
        col_c = "臨界円周 {} 本, R = {}".format(ncirc, D.poly_str(hb["R_MB"]))
        print(f"  {n:>3} {col_a:>18} {col_b:>22} {col_c:>24}")
        ok &= (cn["MB_sharp"] and iv["n_critical"] == 4 * n
               and iv["R_MB"] == [n - 1, n - 1] and hb["MB_sharp"])
    print("\n  非対称な最小 DMF は鋭いが対称でない。対称にすると R が n に比例して悪化する。")
    print("  DMBF なら対称なまま R = 0 のまま。＝ 対称性のある対象での御利益。")
    return ok


def st3_torus_height():
    print(RULE)
    print("強み 3: 滑らかな Morse-Bott 関数の再現。")
    print("        回転対称なトーラスの「軸からの距離」は臨界円周 2 本（指数 0, 1）を持ち")
    print("        Σ_i t^λi P_t(S^1) = (1+t) + t(1+t) = 1 + 2t + t^2 = P_t(T^2)。\n")
    ok = True
    for ni, nj in ((4, 4), (5, 6)):
        K = D.torus(ni, nj)
        M = D.MorseBott(K, D.height_fn(K, ni, nj))
        r = M.report()
        print(f"  T({ni},{nj}): collection {r['collections']} 個のうち，寄与するのは 2 つ")
        for C in sorted(M.reduced_collections(), key=len):
            b = D.poly_trim(D.betti(C))
            if b:
                print(f"      C: {len(C):>3} セル（次元 {sorted({K.dim(c) for c in C})}）"
                      f"  P_t(C) = {D.poly_str(b)}")
        print(f"      Σ_C P_t(C) = {D.poly_str(r['MB_sum'])} = P_t(T^2)，R = 0")
        ok &= r["MB_sharp"] and r["MB_sum"] == [1, 2, 1]
    print("\n  離散モース関数では臨界「点」しか持てないので，この円周は表現できない。")
    return ok


def st4_klein_bottle():
    print(RULE)
    print("強み 4: ねじれと対称性が同時に効く例（クラインの壺）。")
    print("        強み 1（ねじれ）と強み 2（対称性）は独立な障害で，足し算になる。\n")
    print("  クラインの壺は χ = 0 の閉非向き付け曲面で，回転で不変な三角形分割を取ると")
    print("  i 方向の平行移動 g が自己同型になる。g は位数 2·ni（ni 回まわると貼り合わせの")
    print("  反転 j → −j が残る）で自由ではないが，**ni が奇数なら g^2 が生成する Z_ni は")
    print("  自由**（反転は位数 2 で，奇数位数の部分群に入らないから）。\n")
    ok = True
    for ni, nj in ((3, 4), (5, 4)):
        K = X.klein_bottle_sym(ni, nj)
        act = X.klein_free_action(ni, nj)
        free = all(act(c, k) != c for c in K.cells for k in range(1, ni))
        print(f"  Kl({ni},{nj}): セル {K.counts()}，"
              f"H_* = {D.homology_str(D.homology_z(K.cells))}")
        print(f"    Z_{ni} は自由に作用するか: {'はい' if free else 'いいえ'}")
        ok &= free
        print(f"      {'係数':<6}{'P_t(K)':<16}{'非対称な最小 DMF':<22}"
              f"{'不変な DMF（下限）':<24}{'不変 DMBF'}")
        for q, lab in ((0, "Q"), (2, "F_2")):
            cn = D.MorseBott(K, D.canonical_dmf(K), p=q).report()
            cs = D.MorseBott(K, D.constant_fn(K), p=q).report()
            bound = [ni, 2 * ni, ni]          # m_0 = m_2 = ni が最小（下限の等号）
            r_inv = D.poly_div_1pt(D.poly_sub(bound, cn["P_K"]))
            print(f"      {lab:<6}{D.poly_str(cn['P_K']):<16}"
                  f"{'R = ' + D.poly_str(cn['R_M']):<22}"
                  f"{'R ≥ ' + D.poly_str(r_inv):<24}"
                  f"{'R = ' + D.poly_str(cs['R_MB'])}")
            ok &= cs["MB_sharp"]
        print()
    print("  読み方（ni = 3 の場合）:")
    print("    Q 上  : 対称性を捨てても R = t が残る（ねじれのぶん）。")
    print("            対称にすると R ≥ 2 + 3t = (ni−1)(1+t) + t —— 対称性のぶん + ねじれのぶん。")
    print("    F_2 上: ねじれがベッチ数として見えるので R = 0 にできる。")
    print("            対称にすると R ≥ 2 + 2t = (ni−1)(1+t) —— 対称性のぶんだけが残る。")
    print("    どちらの係数でも，離散モースボット関数なら R = 0。")
    print("  ＝ 2 つの障害は独立で，係数体を変えるとねじれのぶんだけを消せる。")
    print("\n  注意: 「不変な DMF」の行は**下限**である（臨界セルが 4·ni 個以上，")
    print("  §3.2 と同じ議論）。トーラスと違い，下限を達成する不変 DMF は構成していない。")
    return ok


def st5_arrows_and_sharpness():
    print(RULE)
    print("強み 5: 矢印（勾配）と臨界円周は同居できる。しかも鋭いまま。")
    print("        最上位の帯の 1 次元以上のセルだけ値 1 にする（dmb_core.arrowed_dmbf）。\n")
    print("  高さ関数（強み 3）は鋭いが矢印を 1 本も持たない。不変な DMF（強み 2）は")
    print("  矢印を持つが鋭くない。両方を兼ねられるか？ → 兼ねられる。\n")
    print(f"  {'T(ni,nj)':>10} {'矢印':>6} {'collection':>11} "
          f"{'Σ_C P_t(C)':>14} {'R(t)':>6}  DMF か")
    ok = True
    for ni, nj in ((3, 3), (5, 4), (4, 6)):
        K = D.torus(ni, nj)
        r = D.MorseBott(K, D.arrowed_dmbf(K, ni, nj)).report()
        print(f"  T({ni},{nj})".rjust(12)
              + f" {r['n_arrows']:>5} {r['collections']:>11} "
                f"{D.poly_str(r['MB_sum']):>14} {D.poly_str(r['R_MB']):>6}"
                f"  {'はい' if r['is_dmf'] else 'いいえ'}")
        ok &= (r["is_dmb"] and not r["is_dmf"] and r["MB_sharp"]
               and r["n_arrows"] == ni and D.is_invariant(ni, nj,
                                                          D.arrowed_dmbf(K, ni, nj)))
    print("\n  reduced collection は高さ関数と同じ 1 + t と t + t^2 の 2 本。")
    print("  違いは矢印が ni 本あること。全数探索でも，鋭い不変 DMBF の大多数は")
    print("  矢印を持つ（T(3,3) の値 2 通りで 689 個中 558 個）。docs/results.md §3.6。")
    return ok


TUTORIAL = [ex1_one_edge, ex2_theorem32, ex3_mb2_violation, ex4_critical_circle,
            ex5_morsification, ex6_cw_complexes]
STRENGTH = [st1_nonorientable, st2_symmetry, st3_torus_height, st4_klein_bottle,
            st5_arrows_and_sharpness]


def main(argv=None):
    argv = sys.argv if argv is None else argv
    which = argv[1] if len(argv) > 1 else "all"
    groups = {"tutorial": TUTORIAL, "strength": STRENGTH,
              "all": TUTORIAL + STRENGTH}
    if which not in groups:
        print("使い方: python3 examples.py [tutorial|strength|all]")
        return 2
    ok = True
    for fn in groups[which]:
        ok &= bool(fn())
        print()
    print(RULE)
    print("すべての例が期待どおり" if ok else "期待と違う例がある")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
