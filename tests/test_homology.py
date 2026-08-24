#!/usr/bin/env python3
"""整数係数ホモロジー（Smith 標準形）と，係数体を変えたベッチ数の検査。

ねじれのある空間では「離散モース不等式が鋭いかどうか」が係数体に依る。
その現象を検査で固定する。
"""

import random
import unittest

import complexes as X
import dmb_core as D


class TestSmithNormalForm(unittest.TestCase):
    def test_diagonal_of_small_matrices(self):
        self.assertEqual(D.smith_diagonal([[2]]), [2])
        self.assertEqual(D.smith_diagonal([[0]]), [])
        self.assertEqual(D.smith_diagonal([[1, 0], [0, 1]]), [1, 1])
        self.assertEqual(D.smith_diagonal([[2, 0], [0, 3]]), [1, 6])   # 約数列 1 | 6
        self.assertEqual(D.smith_diagonal([[6, 4], [4, 6]]), [2, 10])
        self.assertEqual(D.smith_diagonal([]), [])

    def test_divisibility_chain_and_determinant(self):
        """対角成分は約数の列をなし，正方行列なら積が行列式の絶対値。"""
        rng = random.Random(D.hash_seed("snf"))
        for _ in range(60):
            n = rng.randrange(1, 5)
            m = [[rng.randrange(-4, 5) for _ in range(n)] for _ in range(n)]
            d = D.smith_diagonal(m)
            for a, b in zip(d, d[1:]):
                self.assertEqual(b % a, 0, (m, d))
            self.assertEqual(len(d), D.rank(m))          # 非零の個数 = 階数
            if len(d) == n:
                prod = 1
                for x in d:
                    prod *= x
                self.assertEqual(prod, abs(round(determinant(m))), (m, d))


def determinant(m):
    n = len(m)
    if n == 1:
        return m[0][0]
    total = 0
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in m[1:]]
        total += (-1) ** j * m[0][j] * determinant(minor)
    return total


class TestRank(unittest.TestCase):
    """疎な消去 (rank) と素朴な密の消去 (rank_dense) が一致すること。

    rank は T(16,16) の ∂_2 で 118 秒 → 0.1 秒にするための実装なので，
    参照実装との突き合わせで正しさを担保する。"""

    def test_matches_dense_on_random_matrices(self):
        rng = random.Random(D.hash_seed("rank"))
        for _ in range(120):
            rows, cols = rng.randrange(1, 8), rng.randrange(1, 8)
            m = [[rng.randrange(-3, 4) for _ in range(cols)] for _ in range(rows)]
            for p in (0, 2, 3, 5):
                with self.subTest(p=p, m=m):
                    self.assertEqual(D.rank(m, p), D.rank_dense(m, p))

    def test_matches_dense_on_boundary_matrices(self):
        for name, build in (("T(5,4)", lambda: D.torus(5, 4)),
                            ("RP^2", X.projective_plane),
                            ("Klein", lambda: X.klein_bottle(5, 5))):
            K = build()
            by_dim = {}
            for c in K.cells:
                by_dim.setdefault(K.dim(c), []).append(c)
            for k in (1, 2):
                mat = [[D.simplicial_incidence(a, b) for b in by_dim[k]]
                       for a in by_dim[k - 1]]
                for p in (0, 2, 3):
                    with self.subTest(name=name, k=k, p=p):
                        self.assertEqual(D.rank(mat, p), D.rank_dense(mat, p))

    def test_edge_cases(self):
        self.assertEqual(D.rank([]), 0)
        self.assertEqual(D.rank([[0, 0], [0, 0]]), 0)
        self.assertEqual(D.rank([[2, 0], [0, 2]], p=2), 0)   # F_2 では退化する
        self.assertEqual(D.rank([[2, 0], [0, 2]]), 2)


class TestIntegralHomology(unittest.TestCase):
    """ねじれを含む整数係数ホモロジー。"""

    CASES = [
        ("S^2", lambda: X.sphere(2), [(1, []), (0, []), (1, [])]),
        ("S^1", lambda: X.circle(5), [(1, []), (1, [])]),
        ("円板", lambda: X.simplex(2), [(1, []), (0, []), (0, [])]),
        ("トーラス", lambda: X.torus(4, 4), [(1, []), (2, []), (1, [])]),
        ("RP^2", X.projective_plane, [(1, []), (0, [2]), (0, [])]),
        ("クラインの壺", lambda: X.klein_bottle(5, 5), [(1, []), (1, [2]), (0, [])]),
        ("メビウスの帯", X.moebius, [(1, []), (1, []), (0, [])]),
    ]

    def test_homology_with_torsion(self):
        for name, build, want in self.CASES:
            with self.subTest(name):
                self.assertEqual(D.homology_z(build().cells), want)

    def test_free_part_matches_rational_betti(self):
        for name, build, _ in self.CASES:
            with self.subTest(name):
                K = build()
                self.assertEqual(D.poly_trim([free for free, _ in D.homology_z(K.cells)]),
                                 D.poly_trim(D.betti(K.cells)))

    def test_minimal_cw_structures(self):
        for build, want in ((D.cw_circle_minimal, [(1, []), (1, [])]),
                            (D.cw_sphere2_minimal, [(1, []), (0, []), (1, [])]),
                            (D.cw_torus_minimal, [(1, []), (2, []), (1, [])]),
                            (D.cw_projective_plane_minimal, [(1, []), (0, [2]), (0, [])])):
            K = build()
            with self.subTest(K.name):
                self.assertEqual(K.homology_z(), want)

    def test_homology_str(self):
        self.assertEqual(D.homology_str([(1, []), (0, [2]), (0, [])]),
                         "H_0 = Z, H_1 = Z/2, H_2 = 0")
        self.assertEqual(D.homology_str([(1, []), (1, [2])]),
                         "H_0 = Z, H_1 = Z ⊕ Z/2")


class TestBettiOverFiniteFields(unittest.TestCase):
    def test_universal_coefficients_on_examples(self):
        """F_p 係数のベッチ数は，自由部分 + p で割れるねじれ 2 つ分。"""
        for name, build, _ in TestIntegralHomology.CASES:
            K = build()
            hz = D.homology_z(K.cells)
            for p in (2, 3, 5):
                with self.subTest(name=name, p=p):
                    want = []
                    for k, (free, tors) in enumerate(hz):
                        prev = hz[k - 1][1] if k > 0 else []
                        want.append(free + sum(1 for d in tors if d % p == 0)
                                    + sum(1 for d in prev if d % p == 0))
                    self.assertEqual(D.poly_trim(D.betti(K.cells, p=p)),
                                     D.poly_trim(want), (name, p))

    def test_rp2_is_a_homology_sphere_over_f2_only(self):
        K = X.projective_plane()
        self.assertEqual(D.betti(K.cells), [1])              # Q
        self.assertEqual(D.betti(K.cells, p=2), [1, 1, 1])   # F_2
        self.assertEqual(D.betti(K.cells, p=3), [1])         # F_3


class TestSharpnessDependsOnTheField(unittest.TestCase):
    """本題: 「離散モース不等式が鋭いか」は係数体に依るが，
    離散モースボット不等式はどの係数体でも鋭くできる。"""

    NONORIENTABLE = [("RP^2", X.projective_plane), ("Klein", lambda: X.klein_bottle(5, 5))]

    def test_dmf_is_not_sharp_over_q_but_is_over_f2(self):
        for name, build in self.NONORIENTABLE:
            K = build()
            f = D.canonical_dmf(K)
            with self.subTest(name):
                over_q = D.MorseBott(K, f, p=0).report()
                over_f2 = D.MorseBott(K, f, p=2).report()
                self.assertTrue(D.poly_trim(over_q["R_M"]), "Q 上では鋭くないはず")
                self.assertEqual(D.poly_trim(over_f2["R_M"]), [], "F_2 上では鋭いはず")
                self.assertEqual(over_q["M"], over_f2["M"])   # M(t) は係数に依らない

    def test_dmb_is_sharp_over_every_field(self):
        for name, build in self.NONORIENTABLE:
            K = build()
            for p in (0, 2, 3, 5):
                with self.subTest(name=name, p=p):
                    r = D.MorseBott(K, D.constant_fn(K), p=p).report()
                    self.assertTrue(r["MB_sharp"])

    def test_theorem_412_holds_over_finite_fields(self):
        """Theorem 4.12 は体上の鎖複体の議論なので，F_p 上でも成り立つはず。"""
        for name, build in [("RP^2", X.projective_plane), ("T(3,3)", lambda: D.torus(3, 3)),
                            ("メビウス", X.moebius), ("S^2", lambda: X.sphere(2))]:
            K = build()
            for t in range(4):
                rng = random.Random(D.hash_seed(name, "field", t))
                for gen in (D.random_max_extension, D.random_matching_dmf):
                    f = gen(K, rng)
                    for p in (0, 2, 3):
                        M = D.MorseBott(K, f, p=p)
                        if not M.is_dmb():
                            continue
                        with self.subTest(name=name, gen=gen.__name__, t=t, p=p):
                            r = M.report()
                            self.assertIsNotNone(r["R_MB"])
                            self.assertTrue(all(c >= 0 for c in r["R_MB"]))

    def test_minimal_cw_projective_plane(self):
        """最小 CW 構造（3 セル）でも同じ: Q 上で R = t，F_2 上で R = 0。"""
        K = D.cw_projective_plane_minimal()
        f = {c: K.dim(c) for c in K.cells}
        self.assertEqual(D.MorseBott(K, f, p=0).report()["R_M"], [0, 1])
        self.assertEqual(D.poly_trim(D.MorseBott(K, f, p=2).report()["R_M"]), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
