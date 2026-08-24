#!/usr/bin/env python3
"""dmb_core.py の検査。外部依存なし。

    python3 -m unittest test_dmb_core -v
    python3 test_dmb_core.py

検出力を測るために，正しい関数を壊す変異注入（負の対照）も入れてある。
"""

import unittest

import dmb_core as D

SIZES = [(3, 3), (3, 5), (4, 4), (5, 4), (4, 6)]


class TestComplex(unittest.TestCase):
    def test_torus_counts_and_euler(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            self.assertEqual(K.counts(), {0: ni * nj, 1: 3 * ni * nj, 2: 2 * ni * nj})
            chi = sum((-1) ** k * v for k, v in K.counts().items())
            self.assertEqual(chi, 0)

    def test_boundary_squared_is_zero(self):
        """∂∘∂ = 0（接続係数の向きの整合性）。"""
        K = D.torus(3, 4)
        for t in K.cells_of_dim(2):
            for v in K.cells_of_dim(0):
                self.assertEqual(
                    sum(D.incidence(v, e) * D.incidence(e, t) for e in K.below[t]), 0)

    def test_betti_of_torus(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            self.assertEqual(D.betti(K.cells), [1, 2, 1])

    def test_names_are_a_bijection(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            names = D.torus_names(ni, nj)
            self.assertEqual(set(names), set(K.cells))
            self.assertEqual(len(set(names.values())), len(K.cells))


class TestPolynomials(unittest.TestCase):
    def test_div_by_one_plus_t(self):
        self.assertEqual(D.poly_div_1pt([]), [])
        self.assertEqual(D.poly_div_1pt([1, 1]), [1])
        self.assertEqual(D.poly_div_1pt([3, 6, 3]), [3, 3])
        self.assertIsNone(D.poly_div_1pt([1]))
        self.assertIsNone(D.poly_div_1pt([1, 2, 1, 5]))

    def test_str(self):
        self.assertEqual(D.poly_str([1, 2, 1]), "1 + 2t + t^2")
        self.assertEqual(D.poly_str([0, 0]), "0")


class TestMorseBottHeight(unittest.TestCase):
    """回転対称な離散モースボット関数（滑らかな Morse-Bott 高さ関数の離散化）。"""

    def test_is_dmb_and_sharp(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            f = D.height_fn(K, ni, nj)
            X = D.MorseBott(K, f)
            r = X.report()
            self.assertTrue(r["is_dmb"], (ni, nj))
            self.assertTrue(D.is_invariant(ni, nj, f))
            self.assertEqual(r["P_K"], [1, 2, 1])
            self.assertEqual(r["MB_sum"], [1, 2, 1])
            self.assertEqual(D.poly_trim(r["R_MB"]), [])          # R(t) = 0

    def test_two_critical_circles(self):
        """非自明な reduced collection は 2 つで，P_t = (1+t) と t(1+t)。

        滑らかな側の「指数 0 と 1 の臨界円周」に対応する:
        Σ_i t^{λ_i} P_t(S^1) = (1+t) + t(1+t) = 1 + 2t + t^2 = P_t(T^2)。"""
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            X = D.MorseBott(K, D.height_fn(K, ni, nj))
            polys = sorted(tuple(D.poly_trim(D.betti(C)))
                           for C in X.reduced_collections() if D.poly_trim(D.betti(C)))
            self.assertEqual(polys, [(0, 1, 1), (1, 1)], (ni, nj))

    def test_min_circle_is_the_level_zero_circle(self):
        """P_t = 1 + t の reduced collection は j = 0 の円周（頂点 ni 個 + 横の辺 ni 個）。"""
        ni, nj = 5, 4
        K = D.torus(ni, nj)
        X = D.MorseBott(K, D.height_fn(K, ni, nj))
        C = next(C for C in X.reduced_collections() if D.poly_trim(D.betti(C)) == [1, 1])
        self.assertEqual(len(C), 2 * ni)
        self.assertTrue(all(v[1] == 0 for c in C for v in c))


class TestDiscreteMorseFunctions(unittest.TestCase):
    def test_canonical_dmf_has_four_critical_cells(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            f = D.canonical_dmf(K)
            r = D.MorseBott(K, f).report()
            self.assertTrue(r["is_dmf"], (ni, nj))
            self.assertTrue(r["is_dmb"])                # Theorem 3.2: DMF ⊂ DMBF
            self.assertEqual(r["n_critical"], 4)
            self.assertEqual(r["M"], [1, 2, 1])
            self.assertEqual(D.poly_trim(r["R_M"]), [])
            self.assertFalse(D.is_invariant(ni, nj, f))  # 対称性は破れている

    def test_invariant_dmf_attains_4n(self):
        """Z_ni 不変な DMF の臨界セルは 4·ni 個（下限を達成）。"""
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            f = D.invariant_dmf(K, ni, nj)
            r = D.MorseBott(K, f).report()
            self.assertTrue(r["is_dmf"], (ni, nj))
            self.assertTrue(D.is_invariant(ni, nj, f))
            self.assertEqual(r["n_critical"], 4 * ni)
            self.assertEqual(r["M"], [ni, 2 * ni, ni])
            self.assertEqual(r["R_M"], [ni - 1, ni - 1])   # (ni - 1)(1 + t)

    def test_arrowed_dmbf_is_sharp_and_has_arrows(self):
        """矢印を持つ鋭い不変 DMBF（docs/results.md §3.6）。

        高さ関数と同じ Σ_C P_t(C) = P_t(T^2)・R(t) = 0 を与えながら，
        矢印をちょうど ni 本持つ（＝「鋭い不変 DMBF は矢印を持たない」は偽）。"""
        for ni, nj in SIZES:
            with self.subTest(ni=ni, nj=nj):
                K = D.torus(ni, nj)
                f = D.arrowed_dmbf(K, ni, nj)
                r = D.MorseBott(K, f).report()
                self.assertTrue(r["is_dmb"])
                self.assertFalse(r["is_dmf"])            # 矢印と collection が同居する
                self.assertTrue(D.is_invariant(ni, nj, f))
                self.assertEqual(r["n_arrows"], ni)      # 軌道 1 本ぶん
                self.assertEqual(r["MB_sum"], [1, 2, 1])
                self.assertEqual(D.poly_trim(r["R_MB"]), [])
                shapes = sorted(tuple(D.poly_trim(D.betti(C)))
                                for C in r["reduced_collections"]
                                if D.poly_trim(D.betti(C)))
                self.assertEqual(shapes, [(0, 1, 1), (1, 1)])   # t + t^2 と 1 + t

    def test_arrowed_dmbf_arrows_are_one_free_orbit(self):
        """矢印は最上位の帯の横の辺 1 軌道ぶんだけ（帯の外には出ない）。"""
        ni, nj = 5, 4
        K = D.torus(ni, nj)
        X = D.MorseBott(K, D.arrowed_dmbf(K, ni, nj))
        arrows = X.arrows()
        self.assertEqual(len(arrows), ni)
        for s, t in arrows:
            self.assertEqual(K.dim(s), 1)
            self.assertEqual(K.dim(t), 2)
            self.assertEqual({v[1] for v in s}, {nj - 1})   # 最上位の横の辺 ie
            # 相手は 1 つ下の帯の上三角 tu（準位 nj-2 の頂点をちょうど 1 つ持つ）
            self.assertEqual({v[1] for v in t}, {nj - 2, nj - 1})
            self.assertEqual(sum(1 for v in t if v[1] == nj - 2), 1)
        # 矢印の集合も Z_ni 不変
        pairs = {(s, t) for s, t in arrows}
        self.assertTrue(all((D.rotate(ni, nj, s), D.rotate(ni, nj, t)) in pairs
                            for s, t in pairs))

    def test_arrowed_dmbf_over_the_documented_range(self):
        """docs/results.md §3.6 の「ni, nj ∈ {3,…,7} の 25 通りすべて」の裏付け。"""
        n = 0
        for ni in range(3, 8):
            for nj in range(3, 8):
                K = D.torus(ni, nj)
                r = D.MorseBott(K, D.arrowed_dmbf(K, ni, nj)).report()
                self.assertTrue(r["is_dmb"], (ni, nj))
                self.assertEqual(r["n_arrows"], ni, (ni, nj))
                self.assertEqual(r["MB_sum"], [1, 2, 1], (ni, nj))
                self.assertEqual(D.poly_trim(r["R_MB"]), [], (ni, nj))
                n += 1
        self.assertEqual(n, 25)

    def test_section_36_comparison_table(self):
        """docs/results.md §3.6 の比較表の数値（腐らせないための検査）。"""
        for ni, nj in SIZES:
            with self.subTest(ni=ni, nj=nj):
                K = D.torus(ni, nj)
                rows = {name: (len(set(f.values())), D.MorseBott(K, f).report())
                        for name, f in (
                            ("height", D.height_fn(K, ni, nj)),
                            ("arrowed", D.arrowed_dmbf(K, ni, nj)),
                            ("dmf", D.invariant_dmf(K, ni, nj)))}
                nvals, r = rows["height"]
                self.assertEqual((nvals, r["collections"], r["n_arrows"]),
                                 (nj // 2 + 1, 2 * (nj // 2), 0))
                nvals, r = rows["arrowed"]
                self.assertEqual((nvals, r["collections"], r["n_arrows"]), (2, 2, ni))
                nvals, r = rows["dmf"]
                self.assertEqual((nvals, r["collections"], r["n_arrows"]),
                                 (6 * nj, 6 * ni * nj, ni * (3 * nj - 2)))

    def test_any_band_may_be_lifted(self):
        """持ち上げる帯はどれでもよい（j 方向の平行移動で移り合う）。

        帯を 2 つ持ち上げると DMBF のままだが鋭さは失われる（§3.6 の「近くの反例」）。"""
        ni, nj = 5, 5
        K = D.torus(ni, nj)

        def lift(bands):
            def sel(c):
                js = {v[1] % nj for v in c}
                return len(c) > 1 and any(b in js and js <= {b, (b + 1) % nj}
                                          for b in bands)
            return {c: (1 if sel(c) else 0) for c in K.cells}

        self.assertEqual(lift([nj - 1]), D.arrowed_dmbf(K, ni, nj))
        for b in range(nj):
            r = D.MorseBott(K, lift([b])).report()
            self.assertTrue(r["is_dmb"], b)
            self.assertEqual(r["n_arrows"], ni)
            self.assertEqual(D.poly_trim(r["R_MB"]), [])
        r2 = D.MorseBott(K, lift([0, 2])).report()          # 隣り合わない 2 本
        self.assertTrue(r2["is_dmb"])
        self.assertEqual(r2["MB_sum"], [2, 4, 2])
        self.assertEqual(D.poly_trim(r2["R_MB"]), [1, 1])   # 1 + t，鋭くない

    def test_critical_cells_of_invariant_dmf_form_free_orbits(self):
        """不変な DMF の臨界セルは Z_ni の軌道の和なので，各次元で ni の倍数。"""
        ni, nj = 5, 4
        K = D.torus(ni, nj)
        X = D.MorseBott(K, D.invariant_dmf(K, ni, nj))
        crit = set(X.critical())
        self.assertTrue(all(D.rotate(ni, nj, c) in crit for c in crit))
        for k in (0, 1, 2):
            self.assertEqual(len([c for c in crit if K.dim(c) == k]) % ni, 0)

    def test_theorem_32_collection_sizes(self):
        """Theorem 3.2: DMF なら各 collection の大きさは 1 か 2 で reduced。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        for f in (D.canonical_dmf(K), D.invariant_dmf(K, ni, nj)):
            X = D.MorseBott(K, f)
            self.assertTrue(X.is_dmf())
            self.assertTrue(all(len(L) in (1, 2) for L in X.collections()))

    def test_refined_height_loses_sharpness(self):
        """細分版の高さ関数は DMBF だが鋭くない: reduced collection が 2·nj 個になり
        R(t) = (nj - 1)(1 + t)。collection の取り方が鋭さを決めることの反例。"""
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            r = D.MorseBott(K, D.height_fn(K, ni, nj, refine=True)).report()
            self.assertTrue(r["is_dmb"], (ni, nj))
            self.assertEqual(len(r["reduced_collections"]), 2 * nj, (ni, nj))
            self.assertEqual(r["MB_sum"], [nj, 2 * nj, nj], (ni, nj))
            self.assertEqual(r["R_MB"], [nj - 1, nj - 1], (ni, nj))

    def test_constant_function_is_the_trivial_dmbf(self):
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        X = D.MorseBott(K, D.constant_fn(K))
        self.assertTrue(X.is_dmb())
        self.assertFalse(X.is_dmf())
        self.assertEqual(len(X.collections()), 1)
        self.assertEqual(X.morse_bott_polynomial(), [1, 2, 1])


class TestMorsification(unittest.TestCase):
    """Morsification 定理: f'(σ) = (D+1) f(σ) + dim σ は DMF で，
    その臨界セルはちょうど f の reduced collection の和集合。"""

    def test_morsify(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            for f in (D.height_fn(K, ni, nj), D.height_fn(K, ni, nj, True),
                      D.constant_fn(K), D.canonical_dmf(K),
                      D.invariant_dmf(K, ni, nj)):
                X = D.MorseBott(K, f)
                self.assertTrue(X.is_dmb())
                Y = D.MorseBott(K, D.morsify(K, f))
                self.assertTrue(Y.is_dmf(), (ni, nj))
                union = {c for C in X.reduced_collections() for c in C}
                self.assertEqual(set(Y.critical()), union, (ni, nj))


class TestAgainstDmfPy(unittest.TestCase):
    """dmf.py（元の離散モース理論の計算）が離散モースボット理論の特別な場合として
    再現されること。"""

    def test_calcdmf_is_a_dmf_and_a_dmbf(self):
        for n in (3, 4, 5, 6, 7):
            K = D.torus(n, n)
            try:
                f = D.dmf_from_dmf_py(K, n, n)
            except Exception as exc:                     # noqa: BLE001
                self.skipTest(f"dmf.py の calcDMF を読めない: {exc}")
            r = D.MorseBott(K, f).report()
            self.assertTrue(r["is_dmf"], f"grid_size={n}")
            self.assertTrue(r["is_dmb"], f"grid_size={n}")
            self.assertEqual(r["n_critical"], 4, f"grid_size={n}")
            self.assertEqual(D.poly_trim(r["R_MB"]), [], f"grid_size={n}")


class TestNegativeControls(unittest.TestCase):
    """変異注入: 正しい関数を壊したら検出できるか（検出力の測定）。"""

    def test_lowering_one_edge_breaks_mb4_everywhere(self):
        """辺 1 本の値を全体より下げると，その辺で D^snc = 2 となり (MB4) が破れる。
        どの辺で壊しても検出される（検出率 100%）。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = D.height_fn(K, ni, nj)
        edges = K.cells_of_dim(1)
        detected = 0
        for c in edges:
            g = dict(f)
            g[c] = min(f.values()) - 1
            X = D.MorseBott(K, g)
            if not X.is_dmb() and {k for k, _, _ in X.dmb_violations()} == {"MB4"}:
                detected += 1
        self.assertEqual(detected, len(edges))

    def test_mb2_violation_is_detected(self):
        """1 つの頂点の 2 本の余次元 1 の coface を下げると (MB2) が破れる。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = D.height_fn(K, ni, nj)
        v = K.cells_of_dim(0)[0]
        g = dict(f)
        for e in K.above[v][:2]:
            g[e] = f[v] - 1
        X = D.MorseBott(K, g)
        self.assertFalse(X.is_dmb())
        self.assertIn("MB2", {k for k, _, _ in X.dmb_violations()})
        self.assertEqual(len(X.up_snc(v)), 2)

    def test_broken_dmf_is_detected_but_may_stay_dmb(self):
        """DMF の値を 1 つ潰すと (M2)/(M4) は破れるが (MB2)/(MB4) は残りうる
        （まさにこれが離散モースボット関数が離散モース関数より広い理由）。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = D.canonical_dmf(K)
        v = K.cells_of_dim(0)[0]
        g = dict(f)
        for e in K.above[v]:
            g[e] = f[v]                       # v のまわりの辺を全部同じ値にする
        X = D.MorseBott(K, g)
        self.assertFalse(X.is_dmf())          # U(v) = 6 > 1
        self.assertTrue(X.is_dmb())           # U^snc(v) = 0
        self.assertGreater(len(X.collections()[0]), 1)

    def test_mb4_violation_is_detected(self):
        """1 つの辺の 2 つの端点を辺より高くすると (MB4) が破れる。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = D.constant_fn(K)
        e = K.cells_of_dim(1)[0]
        g = dict(f)
        g[e] = -1
        X = D.MorseBott(K, g)
        self.assertFalse(X.is_dmb())
        self.assertEqual({k for k, _, _ in X.dmb_violations()}, {"MB4"})

    def test_wrong_incidence_breaks_betti(self):
        """接続係数の符号を落とすと（∂∘∂ ≠ 0）ベッチ数が変わること
        ＝ betti() が向きを本当に使っていることの確認。"""
        K = D.torus(3, 4)

        def unsigned(a, b):
            return 1 if set(a) < set(b) and len(b) == len(a) + 1 else 0

        self.assertNotEqual(D.betti(K.cells, unsigned), [1, 2, 1])
        self.assertEqual(D.betti(K.cells), [1, 2, 1])
        self.assertEqual(D.betti(K.cells, K.incidence), [1, 2, 1])


class TestLayout(unittest.TestCase):
    """描画用の持ち上げ座標の整合性（図が正しい三角形分割を描いているか）。"""

    def test_lift_covers_every_cell(self):
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            pos = D.lifted_cells(ni, nj)
            self.assertEqual(set(pos), set(K.cells))
            for c, places in pos.items():
                self.assertTrue(places)
                for p in places:
                    self.assertEqual(len(p), len(c))
                    for x, y in p:
                        self.assertTrue(0 <= x <= ni and 0 <= y <= nj)

    def test_every_placement_reduces_to_the_cell(self):
        """どの持ち上げも mod (ni, nj) で元のセルの頂点集合に戻る。"""
        for ni, nj in SIZES:
            for c, places in D.lifted_cells(ni, nj).items():
                for p in places:
                    self.assertEqual({(x % ni, y % nj) for x, y in p}, set(c))

    def test_triangles_are_unit_halves_and_tile_the_square(self):
        """三角形はすべて面積 1/2 で，合計面積は ni·nj（基本領域をちょうど覆う）。"""
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            pos = D.lifted_cells(ni, nj)
            total = 0.0
            for t in K.cells_of_dim(2):
                self.assertEqual(len(pos[t]), 1)      # 三角形は複製を持たない
                (x0, y0), (x1, y1), (x2, y2) = pos[t][0]
                area = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) / 2
                self.assertAlmostEqual(area, 0.5)
                total += area
            self.assertAlmostEqual(total, ni * nj)

    def test_every_triangle_edge_is_drawn_on_its_corners(self):
        """三角形の 3 辺は，どれかの持ち上げでその三角形の 3 頂点を結ぶ
        （＝図の中で辺が三角形の縁として現れる）。"""
        for ni, nj in SIZES:
            K = D.torus(ni, nj)
            pos = D.lifted_cells(ni, nj)
            for t in K.cells_of_dim(2):
                corners = set(pos[t][0])
                for e in K.below[t]:
                    self.assertTrue(any(set(p) <= corners for p in pos[e]), (t, e))
                    for v in K.below[e]:
                        self.assertTrue(any(set(p) <= corners for p in pos[v]), (t, v))

    def test_boundary_copies_close_the_picture(self):
        """基本領域の 4 辺すべてに，辺のセルが描かれている。"""
        ni, nj = 4, 5
        K = D.torus(ni, nj)
        pos = D.lifted_cells(ni, nj)
        segs = [p for c in K.cells_of_dim(1) for p in pos[c]]
        for coord, val in ((0, 0), (0, ni), (1, 0), (1, nj)):
            on = [p for p in segs if all(q[coord] == val for q in p)]
            self.assertEqual(len(on), ni if coord == 1 else nj, (coord, val))

    def test_centroid(self):
        self.assertEqual(D.centroid([(0, 0), (2, 0)]), (1.0, 0.0))
        self.assertEqual(D.centroid([(0, 0), (3, 0), (0, 3)]), (1.0, 1.0))


class TestMatchingGuards(unittest.TestCase):
    def test_cyclic_matching_is_rejected(self):
        """閉軌道を含む matching は位相ソートで弾かれる（acyclic の検査）。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        matching = []
        for i in range(ni):
            for j in range(nj):
                ip, jp = (i + 1) % ni, (j + 1) % nj
                matching.append((D.cell((i, j), (ip, jp)),
                                 D.cell((i, j), (ip, j), (ip, jp))))
                matching.append((D.cell((i, jp)), D.cell((i, j), (i, jp))))
                matching.append((D.cell((i, jp), (ip, jp)),
                                 D.cell((i, j), (i, jp), (ip, jp))))
        with self.assertRaises(ValueError):
            D.function_from_matching(K, matching)


if __name__ == "__main__":
    unittest.main(verbosity=2)
