#!/usr/bin/env python3
"""正則とは限らない CW 複体 (dmb_core.CWComplex) の検査。

単体的複体では (M1)(M3) が vacuous なので，そこだけを見ていると
「任意余次元の irregular face に課す」という条件の意味が検査されない。
最小 CW 構造でそこを突く。
"""

import random
import unittest

import complexes as X
import dmb_core as D

MINIMAL = [D.cw_circle_minimal, D.cw_sphere2_minimal,
           D.cw_torus_minimal, D.cw_projective_plane_minimal]


class TestCWComplex(unittest.TestCase):
    def test_minimal_cw_structures(self):
        want = {
            "S^1 (最小 CW)": ({0: 1, 1: 1}, [1, 1]),
            "S^2 (最小 CW)": ({0: 1, 2: 1}, [1, 0, 1]),
            "T^2 (最小 CW)": ({0: 1, 1: 2, 2: 1}, [1, 2, 1]),
            "RP^2 (最小 CW)": ({0: 1, 1: 1, 2: 1}, [1]),
        }
        for build in MINIMAL:
            K = build()
            with self.subTest(K.name):
                counts, betti = want[K.name]
                self.assertEqual(K.counts(), counts)
                self.assertEqual(K.check_boundary(), [], "∂∘∂ = 0 でない")
                self.assertEqual(K.betti(), betti)

    def test_transitive_closure(self):
        """生成元だけ渡しても face 関係の推移閉包が取られる。"""
        K = D.CWComplex({"a": 0, "b": 1, "c": 2}, [("a", "b"), ("b", "c")])
        self.assertTrue(K.lt("a", "c"))
        self.assertEqual(sorted(K.irregular_above("a")), ["b", "c"])
        self.assertEqual(K.above["a"], ["b"])          # 余次元 1 は b だけ

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):            # 次元を上げていない
            D.CWComplex({"a": 1, "b": 1}, [("a", "b")])
        with self.assertRaises(ValueError):            # regular が face でない
            D.CWComplex({"a": 0, "b": 1}, [("a", "b")], regular=[("b", "a")])

    def test_all_faces_irregular_by_default(self):
        K = D.cw_torus_minimal()
        for a in K.cells:
            for b in K.cells:
                if K.lt(a, b):
                    self.assertFalse(K.is_regular(a, b))
                    self.assertIn(b, K.irregular_above(a))
                    self.assertIn(a, K.irregular_below(b))


class TestSimplicialAsCW(unittest.TestCase):
    """単体的複体を CW 複体として見ても結論が変わらないこと（往復検査）。"""

    CASES = [("S^2", X.sphere(2)), ("メビウス", X.moebius()),
             ("T(3,3)", X.torus(3, 3)), ("RP^2", X.projective_plane())]

    def test_same_betti(self):
        for name, K in self.CASES:
            with self.subTest(name):
                self.assertEqual(D.cw_from_simplicial(K).betti(), D.betti(K.cells))

    def test_same_verdicts_and_polynomials(self):
        for name, K in self.CASES:
            C = D.cw_from_simplicial(K)
            for t in range(4):
                rng = random.Random(D.hash_seed(name, t))
                for gen in (D.random_function, D.random_max_extension,
                            D.random_matching_dmf):
                    f = gen(K, rng)
                    with self.subTest(name=name, gen=gen.__name__, t=t):
                        a = D.MorseBott(K, f).report()
                        b = D.MorseBott(C, f).report()
                        for key in ("is_dmf", "is_dmb", "P_K", "MB_sum", "R_MB",
                                    "collections", "n_critical", "n_weakly_critical"):
                            self.assertEqual(a[key], b[key], (name, gen.__name__, key))

    def test_m1_m3_are_vacuous_on_simplicial_complexes(self):
        for name, K in self.CASES:
            C = D.cw_from_simplicial(K)
            f = D.random_function(K, random.Random(D.hash_seed(name, "m1")))
            with self.subTest(name):
                M = D.MorseBott(C, f)
                self.assertEqual(M.m1_violations(), [])
                self.assertEqual(M.m3_violations(), [])
                self.assertEqual(K.irregular_above(K.cells[0]), [])


class TestCounterexampleToLemma31(unittest.TestCase):
    """arXiv:2511.07864 v1 の Lemma 3.1（DMF ⇒ DMBF）への反例。

    Lean の `Counterexample.dmf_not_imp_dmb` と同じ例。S^2 の最小 CW 構造は
    余次元 1 の組を 1 つも持たないので，(M1)(M3) を「irregular *facet* だけ」と
    読むとどんな関数も DMF になってしまう。"""

    def setUp(self):
        self.K = D.cw_sphere2_minimal()
        self.f = {"v": 1, "e": 0}
        self.M = D.MorseBott(self.K, self.f)

    def test_no_codimension_one_pairs(self):
        self.assertEqual([(a, b) for a in self.K.cells for b in self.K.above[a]], [])
        self.assertTrue(self.K.lt("v", "e"))
        self.assertFalse(self.K.is_regular("v", "e"))
        self.assertEqual(self.K.dim("e") - self.K.dim("v"), 2)

    def test_weak_reading_makes_everything_a_dmf(self):
        self.assertTrue(self.M.is_dmf(strong=False))
        for f in ({"v": 1, "e": 0}, {"v": 0, "e": 1}, {"v": 5, "e": 5}):
            self.assertTrue(D.MorseBott(self.K, f).is_dmf(strong=False))

    def test_not_dmb(self):
        self.assertFalse(self.M.is_dmb())
        self.assertEqual({k for k, _, _ in self.M.dmb_violations()}, {"M1", "M3"})

    def test_strong_reading_removes_the_mismatch(self):
        """v2 の Definition 12（任意余次元）なら DMF でもなくなる。"""
        self.assertFalse(self.M.is_dmf(strong=True))

    def test_increasing_function_is_fine(self):
        M = D.MorseBott(self.K, {"v": 0, "e": 1})
        self.assertTrue(M.is_dmf(strong=True))
        self.assertTrue(M.is_dmb())


class TestPerfectFunctionsOnMinimalCW(unittest.TestCase):
    """最小 CW 構造では，次元に沿って増える関数がそのまま「完全な」関数になる。"""

    def test_torus_and_sphere_are_perfect(self):
        for build, poly in ((D.cw_torus_minimal, [1, 2, 1]),
                            (D.cw_sphere2_minimal, [1, 0, 1])):
            K = build()
            r = D.MorseBott(K, {c: K.dim(c) for c in K.cells}).report()
            with self.subTest(K.name):
                self.assertTrue(r["is_dmf"])
                self.assertTrue(r["is_dmb"])
                self.assertEqual(r["M"], poly)
                self.assertEqual(r["P_K"], poly)
                self.assertEqual(D.poly_trim(r["R_MB"]), [])     # R = 0（完全）

    def test_projective_plane_cannot_be_perfect(self):
        """RP^2 は最小 CW 構造（3 セル）でも R(t) = t ≠ 0。

        有理数係数で b_2 = 0 なのに m_2 = 1 が要るため。強み 1 の最小の証拠。"""
        K = D.cw_projective_plane_minimal()
        r = D.MorseBott(K, {c: K.dim(c) for c in K.cells}).report()
        self.assertTrue(r["is_dmf"])
        self.assertEqual(r["M"], [1, 1, 1])
        self.assertEqual(r["P_K"], [1])
        self.assertEqual(r["R_MB"], [0, 1])                      # R(t) = t


if __name__ == "__main__":
    unittest.main(verbosity=2)
