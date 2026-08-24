#!/usr/bin/env python3
"""complexes.py の検査と，トーラス以外の複体の上での離散モースボット理論の検査。

理論の計算はトーラス専用ではないので，ホモロジーが分かっている複体で確かめる。

    python3 -m unittest test_complexes -v
"""

import random
import unittest

import complexes as X
import dmb_core as D


class TestCatalogue(unittest.TestCase):
    def test_betti_numbers(self):
        for name, (build, expected) in X.CATALOGUE.items():
            with self.subTest(name):
                self.assertEqual(D.betti(build().cells), expected)

    def test_euler_characteristics(self):
        expected = {
            "point (Δ^0)": 1, "segment (Δ^1)": 1, "triangle (Δ^2, 円板)": 1,
            "S^1 (3 頂点)": 0, "S^1 (7 頂点)": 0, "S^1 ⊔ S^1": 0,
            "S^2 = ∂Δ^3": 2, "S^3 = ∂Δ^4": 0, "円筒 (5×2)": 0,
            "メビウスの帯": 0, "RP^2 (6 頂点)": 1, "クラインの壺 (5×5)": 0,
            "クラインの壺・対称 (3×4)": 0,
            "トーラス T(4,4)": 0,
        }
        for name, (build, betti) in X.CATALOGUE.items():
            with self.subTest(name):
                K = build()
                self.assertEqual(X.euler(K), expected[name])
                # χ = Σ (-1)^k b_k
                self.assertEqual(X.euler(K),
                                 sum((-1) ** k * b for k, b in enumerate(betti)))

    def test_boundary_squared_is_zero(self):
        for name, (build, _) in X.CATALOGUE.items():
            with self.subTest(name):
                K = build()
                for c in K.cells:
                    if K.dim(c) < 2:
                        continue
                    for g in K.cells:
                        if K.dim(g) != K.dim(c) - 2:
                            continue
                        self.assertEqual(
                            sum(D.incidence(g, h) * D.incidence(h, c) for h in K.below[c]),
                            0, (name, c, g))

    def test_closed_surfaces(self):
        for K in (X.sphere(2), X.torus(4, 4), X.klein_bottle(5, 5), X.projective_plane()):
            self.assertTrue(X.is_closed_surface(K))
            self.assertEqual(X.boundary_edges(K), [])

    def test_surfaces_with_boundary(self):
        for K, nboundary in ((X.moebius(), 5), (X.annulus(5, 2), 10), (X.simplex(2), 3)):
            self.assertFalse(X.is_closed_surface(K))
            self.assertEqual(len(X.boundary_edges(K)), nboundary)

    def test_klein_bottle_is_nonorientable(self):
        """クラインの壺は χ = 0 かつ b_2 = 0（向き付け不可能）でトーラスと区別される。"""
        Kl = X.klein_bottle(5, 5)
        T = X.torus(5, 5)
        self.assertEqual(X.euler(Kl), X.euler(T))
        self.assertEqual(D.betti(Kl.cells), [1, 1])
        self.assertEqual(D.betti(T.cells), [1, 2, 1])


class TestComplexStructure(unittest.TestCase):
    def test_face_relations_are_symmetric(self):
        for name, (build, _) in X.CATALOGUE.items():
            with self.subTest(name):
                K = build()
                for c in K.cells:
                    for d in K.below[c]:
                        self.assertIn(c, K.above[d])
                        self.assertEqual(K.dim(d), K.dim(c) - 1)
                        self.assertTrue(set(d) < set(c))
                    for d in K.above[c]:
                        self.assertIn(c, K.below[d])

    def test_incidence_is_nonzero_exactly_on_facets(self):
        K = X.sphere(2)
        for a in K.cells:
            for b in K.cells:
                facet = set(a) < set(b) and K.dim(b) == K.dim(a) + 1
                self.assertEqual(D.incidence(a, b) != 0, facet, (a, b))

    def test_closure_is_taken(self):
        """facet だけ渡してもすべての面が入る。"""
        K = D.Complex([(0, 1, 2)])
        self.assertEqual(K.counts(), {0: 3, 1: 3, 2: 1})

    def test_betti_is_additive_over_disjoint_cell_sets(self):
        """接続のないセル集合に分けるとベッチ数は和になる。

        reduced collection を facet 関係の連結成分に分けても Σ_C P_t(C) が
        変わらないことの根拠。"""
        K = X.two_circles(4)
        A = [c for c in K.cells if max(c) < 4]
        B = [c for c in K.cells if min(c) >= 4]
        self.assertEqual(len(A) + len(B), len(K.cells))
        self.assertFalse(any(set(a) < set(b) or set(b) < set(a) for a in A for b in B))
        self.assertEqual(D.betti(A), [1, 1])
        self.assertEqual(D.betti(B), [1, 1])
        self.assertEqual(D.poly_trim(D.betti(A + B)),
                         D.poly_trim(D.poly_add(D.betti(A), D.betti(B))))


class TestEquivariantKleinBottle(unittest.TestCase):
    """回転で不変な三角形分割のクラインの壺と，その自由な Z_ni 作用。

    ねじれ（強み 1）と対称性（強み 2）が同時に効く例。"""

    def test_is_a_klein_bottle(self):
        for ni, nj in ((3, 4), (4, 4), (5, 4), (3, 6)):
            K = X.klein_bottle_sym(ni, nj)
            with self.subTest(ni=ni, nj=nj):
                self.assertEqual(K.counts(),
                                 {0: 2 * ni * nj, 1: 6 * ni * nj, 2: 4 * ni * nj})
                self.assertEqual(X.euler(K), 0)
                self.assertTrue(X.is_closed_surface(K))
                self.assertEqual(D.betti(K.cells), [1, 1])
                self.assertEqual(D.homology_z(K.cells), [(1, []), (1, [2]), (0, [])])

    def test_shift_is_an_automorphism(self):
        """対角線で分けた `klein_bottle` では自己同型にならないが，
        中心を入れた `klein_bottle_sym` ではなる。"""
        for ni, nj in ((3, 4), (4, 4), (5, 6)):
            K = X.klein_bottle_sym(ni, nj)
            cells = set(K.cells)
            with self.subTest(ni=ni, nj=nj):
                for c in K.cells:
                    self.assertIn(X.klein_shift(c, ni, nj), cells)

    def test_shift_has_order_2ni(self):
        """ni 回まわると貼り合わせの反転が残るので位数は 2·ni。"""
        ni, nj = 3, 4
        K = X.klein_bottle_sym(ni, nj)
        for c in K.cells:
            self.assertEqual(X.klein_shift(c, ni, nj, 2 * ni), c)
        self.assertTrue(any(X.klein_shift(c, ni, nj, ni) != c for c in K.cells))

    def test_free_only_for_odd_ni(self):
        """g^2 が生成する Z_ni が自由なのは ni が奇数のときだけ。"""
        for ni in (3, 5, 7):
            nj = 4
            K = X.klein_bottle_sym(ni, nj)
            act = X.klein_free_action(ni, nj)
            with self.subTest(ni=ni):
                for c in K.cells:
                    for k in range(1, ni):
                        self.assertNotEqual(act(c, k), c)
                    self.assertEqual(act(c, ni), c)
        for ni in (4, 6):
            with self.subTest(ni=ni):
                with self.assertRaises(ValueError):
                    X.klein_free_action(ni, 4)

    def test_both_obstructions_add_up(self):
        """R(不変 DMF の下限) = (ni−1)(1+t) + t（対称性のぶん + ねじれのぶん）。

        F_2 上ではねじれのぶんが消えて (ni−1)(1+t) だけになる。"""
        for ni in (3, 5):
            K = X.klein_bottle_sym(ni, 4)
            bound = [ni, 2 * ni, ni]                 # m_0 = m_2 = ni（下限の等号）
            with self.subTest(ni=ni):
                pq = D.betti(K.cells)                # Q: 1 + t
                p2 = D.betti(K.cells, p=2)           # F_2: 1 + 2t + t^2
                self.assertEqual(pq, [1, 1])
                self.assertEqual(p2, [1, 2, 1])
                r_q = D.poly_div_1pt(D.poly_sub(bound, pq))
                r_2 = D.poly_div_1pt(D.poly_sub(bound, p2))
                self.assertEqual(r_2, [ni - 1, ni - 1])          # (ni−1)(1+t)
                self.assertEqual(r_q, [ni - 1, ni])              # 上に t を足したもの
                self.assertEqual(D.poly_sub(r_q, r_2), [0, 1])   # 差はちょうど t

    def test_dmb_is_sharp_over_both_fields(self):
        for ni in (3, 5):
            K = X.klein_bottle_sym(ni, 4)
            for p in (0, 2, 3):
                with self.subTest(ni=ni, p=p):
                    self.assertTrue(
                        D.MorseBott(K, D.constant_fn(K), p=p).report()["MB_sharp"])

    def test_minimal_dmf_shows_the_torsion_gap_over_q_only(self):
        K = X.klein_bottle_sym(3, 4)
        f = D.canonical_dmf(K)
        self.assertEqual(D.MorseBott(K, f, p=0).report()["R_M"], [0, 1])   # R = t
        self.assertEqual(D.poly_trim(D.MorseBott(K, f, p=2).report()["R_M"]), [])


class TestTheoryOnOtherComplexes(unittest.TestCase):
    """dmb_core の理論の計算がトーラス以外でも正しいこと。"""

    def test_constant_function_gives_connected_components(self):
        """定数関数の collection は連結成分そのもの，Σ_C P_t(C) = P_t(K)，R = 0。"""
        for name, (build, expected) in X.CATALOGUE.items():
            with self.subTest(name):
                K = build()
                M = D.MorseBott(K, D.constant_fn(K))
                r = M.report()
                self.assertTrue(r["is_dmb"])
                self.assertEqual(r["MB_sum"], expected)
                self.assertEqual(D.poly_trim(r["R_MB"]), [])
                self.assertEqual(len(M.collections()), expected[0])  # b_0 = 連結成分の数

    def test_tree_cotree_dmf_on_closed_surfaces(self):
        """tree-cotree の離散モース関数は，どの閉曲面でも m = (1, 2-χ, 1)。

        トーラス (1,2,1)・S^2 (1,0,1)・クラインの壺 (1,2,1)・RP^2 (1,1,1)。"""
        cases = [(X.sphere(2), [1, 0, 1], [1, 0, 1]),
                 (X.torus(4, 4), [1, 2, 1], [1, 2, 1]),
                 (X.klein_bottle(5, 5), [1, 2, 1], [1, 1]),
                 (X.projective_plane(), [1, 1, 1], [1])]
        for K, morse, betti in cases:
            with self.subTest(K.counts()):
                r = D.MorseBott(K, D.canonical_dmf(K)).report()
                self.assertTrue(r["is_dmf"])
                self.assertEqual(r["M"], morse)
                self.assertEqual(r["P_K"], betti)
                self.assertEqual(r["M"][1], 2 - X.euler(K) + morse[0] + morse[2] - 2)
                self.assertIsNotNone(r["R_M"])
                self.assertTrue(all(c >= 0 for c in r["R_M"]))

    def test_rp2_morse_inequality_is_strict_over_the_rationals(self):
        """RP^2 は有理数係数で P_t = 1 だが，どんな離散モース関数でも
        m_1 ≥ 1, m_2 ≥ 1（ねじれのため）。tree-cotree は R(t) = t を与える。"""
        K = X.projective_plane()
        r = D.MorseBott(K, D.canonical_dmf(K)).report()
        self.assertEqual(r["M"], [1, 1, 1])
        self.assertEqual(r["P_K"], [1])
        self.assertEqual(r["R_M"], [0, 1])          # R(t) = t

    def test_max_extension_is_dmb_on_every_complex(self):
        """max_extension は複体・h によらず必ず離散モースボット関数（全セル weakly
        critical）で，Theorem 4.12 が成り立つ（R(t) ≥ 0）。

        ただし**鋭い（R = 0）とは限らない**: h が滑らかな Morse-Bott 関数の
        離散化になっているときだけ等号になる（トーラスの高さ関数がその例）。"""
        for name, (build, _) in X.CATALOGUE.items():
            for seed in (0, 1, 2):
                with self.subTest(name=name, seed=seed):
                    K = build()
                    # hash() は PYTHONHASHSEED で変わるので使わない（再現性のため）
                    rng = random.Random(seed)
                    h = {v[0]: rng.randrange(5) for v in K.cells_of_dim(0)}
                    f = D.max_extension(K, h.__getitem__)
                    M = D.MorseBott(K, f)
                    r = M.report()
                    self.assertTrue(r["is_dmb"])
                    self.assertEqual(len(M.weakly_critical()), len(K.cells))
                    self.assertIsNotNone(r["R_MB"], (name, seed))
                    self.assertTrue(all(c >= 0 for c in r["R_MB"]), (name, seed))

    def test_torus_height_is_the_sharp_case(self):
        """同じ max_extension でも，トーラスの高さ関数は R(t) = 0 で鋭い。"""
        for ni, nj in ((3, 4), (4, 4), (5, 6)):
            K = X.torus(ni, nj)
            r = D.MorseBott(K, D.height_fn(K, ni, nj)).report()
            self.assertEqual(r["MB_sum"], [1, 2, 1])
            self.assertEqual(D.poly_trim(r["R_MB"]), [])

    def test_sphere_height_has_two_critical_circles(self):
        """S^2 上の「高さ」（頂点を 2 段に分ける）は，最小・最大の 2 つの
        非自明な reduced collection を持つ。"""
        K = X.sphere(2)
        f = D.max_extension(K, lambda v: 0 if v in (0, 1) else 1)
        M = D.MorseBott(K, f)
        r = M.report()
        self.assertTrue(r["is_dmb"])
        self.assertEqual(r["MB_sum"], [1, 0, 1])
        self.assertEqual(D.poly_trim(r["R_MB"]), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
