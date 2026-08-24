#!/usr/bin/env python3
"""Z_ni 不変な離散モースボット関数の全数探索 (search.py) の検査。

全数探索そのものは重いので，既定では軽い検査だけを走らせる。
`DMB_SLOW=1` を付けると，総当たりとの突き合わせと分類まで確かめる。
"""

import itertools
import os
import unittest

import dmb_core as D
import search as S

SLOW = os.environ.get("DMB_SLOW")


class TestOrbits(unittest.TestCase):
    def test_orbit_structure(self):
        """軌道は 6·nj 個，どれも ni セルで，全体をちょうど覆う。"""
        for ni, nj in ((3, 3), (4, 3), (5, 4), (4, 6)):
            with self.subTest(ni=ni, nj=nj):
                K = D.torus(ni, nj)
                orbits = S.invariant_orbits(K, ni, nj)
                self.assertEqual(len(orbits), 6 * nj)
                self.assertEqual({len(o) for o in orbits}, {ni})
                flat = [c for o in orbits for c in o]
                self.assertEqual(sorted(flat), sorted(K.cells))

    def test_orbits_are_rotation_orbits(self):
        ni, nj = 5, 4
        K = D.torus(ni, nj)
        for orbit in S.invariant_orbits(K, ni, nj):
            members = set(orbit)
            for c in orbit:
                self.assertEqual({D.rotate(ni, nj, c, k) for k in range(ni)}, members)

    def test_labels_are_unique_and_typed(self):
        ni, nj = 4, 5
        K = D.torus(ni, nj)
        labels = [S.orbit_label(o[0], nj) for o in S.invariant_orbits(K, ni, nj)]
        self.assertEqual(len(set(labels)), len(labels))
        for j in range(nj):
            for kind in ("v", "ie", "je", "dg", "tu", "tl"):
                self.assertIn(f"{kind}{j}", labels)

    def test_base_level(self):
        nj = 5
        self.assertEqual(S.base_level(D.cell((0, 2)), nj), 2)
        self.assertEqual(S.base_level(D.cell((0, 2), (0, 3)), nj), 2)
        self.assertEqual(S.base_level(D.cell((0, 4), (0, 0)), nj), 4)   # 貼り合わせ
        self.assertEqual(S.base_level(D.cell((0, 1), (1, 1)), nj), 1)


class TestPruningIsSound(unittest.TestCase):
    """枝刈りが離散モースボット関数を取りこぼさないこと。"""

    def test_known_invariant_dmbf_survive(self):
        for ni, nj in ((3, 3), (4, 3), (5, 4)):
            K = D.torus(ni, nj)
            orbits = S.invariant_orbits(K, ni, nj)
            functions = {
                "constant": D.constant_fn(K),
                "height": D.height_fn(K, ni, nj),
                "height refined": D.height_fn(K, ni, nj, refine=True),
                "invariant DMF": D.invariant_dmf(K, ni, nj),
                "morsified": D.morsify(K, D.height_fn(K, ni, nj)),
            }
            for name, f in functions.items():
                with self.subTest(ni=ni, nj=nj, f=name):
                    self.assertTrue(D.is_invariant(ni, nj, f))
                    self.assertTrue(D.MorseBott(K, f).is_dmb())
                    vals = tuple(f[o[0]] for o in orbits)
                    self.assertTrue(S.survives_pruning(ni, nj, vals))

    def test_obvious_violation_is_pruned(self):
        """頂点のまわりの縦の辺と斜めの辺を両方下げると U^snc = 2 なので枝刈りされる。"""
        ni, nj = 4, 3
        K = D.torus(ni, nj)
        orbits = S.invariant_orbits(K, ni, nj)
        vals = []
        for o in orbits:
            label = S.orbit_label(o[0], nj)
            vals.append(0 if label in ("je0", "dg0") else 1)
        f = {c: vals[k] for k, o in enumerate(orbits) for c in o}
        self.assertFalse(D.MorseBott(K, f).is_dmb())
        self.assertFalse(S.survives_pruning(ni, nj, tuple(vals)))


class TestOrbitRepresentativeCheck(unittest.TestCase):
    """不変な f では，軌道の代表元だけ見た (MB2)/(MB4) の検査が全セル版と同値。

    探索はこれで 5 倍速くなっているので，同値性を固定しておく。"""

    def test_agrees_with_full_check(self):
        import random
        for ni, nj in ((3, 3), (4, 3), (5, 4)):
            K = D.torus(ni, nj)
            orbits = S.invariant_orbits(K, ni, nj)
            reps = [o[0] for o in orbits]
            for t in range(120):
                rng = random.Random(D.hash_seed("reps", ni, nj, t))
                vals = [rng.randrange(3) for _ in orbits]
                f = {c: vals[k] for k, o in enumerate(orbits) for c in o}
                M = D.MorseBott(K, f)
                by_reps = all(len(M.up_snc(c)) <= 1 and len(M.down_snc(c)) <= 1
                              for c in reps)
                with self.subTest(ni=ni, nj=nj, t=t):
                    self.assertEqual(by_reps, M.is_dmb())

    def test_counts_per_cell_are_constant_along_orbits(self):
        """そもそも U^snc, D^snc が軌道の上で一定であること（同値性の理由）。"""
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = D.height_fn(K, ni, nj)
        M = D.MorseBott(K, f)
        for orbit in S.invariant_orbits(K, ni, nj):
            self.assertEqual(len({len(M.up_snc(c)) for c in orbit}), 1)
            self.assertEqual(len({len(M.down_snc(c)) for c in orbit}), 1)


class TestSearch(unittest.TestCase):
    def test_results_are_invariant_dmbf(self):
        ni, nj = 4, 3
        found, nodes, n_dmb = S.search(ni, nj, 2, limit=40)
        self.assertTrue(found)
        self.assertGreater(nodes, 0)
        K = D.torus(ni, nj)
        orbits = S.invariant_orbits(K, ni, nj)
        for vals, r in found:
            f = {c: vals[k] for k, o in enumerate(orbits) for c in o}
            self.assertTrue(D.is_invariant(ni, nj, f))
            M = D.MorseBott(K, f)
            self.assertTrue(M.is_dmb())
            self.assertEqual(r["MB_sum"], M.morse_bott_polynomial())
            self.assertEqual(r["MB_sharp"], not D.poly_trim(r["R_MB"] or []))

    def test_sharp_only_filters(self):
        found, _, _ = S.search(4, 3, 2, sharp_only=True, limit=20)
        self.assertTrue(found)
        for _, r in found:
            self.assertTrue(r["MB_sharp"])
            self.assertEqual(r["MB_sum"], [1, 2, 1])

    def test_shape_str(self):
        self.assertEqual(S.shape_str(()), "（寄与する collection なし）")
        self.assertEqual(S.shape_str(((1, 1), (0, 1, 1))), "(t + t^2)  (1 + t)")
        self.assertEqual(S.shape_str(((1, 1), (1, 1))), "(1 + t) ×2")

    def test_cli(self):
        self.assertEqual(S.main(["--ni", "4", "--nj", "3", "--values", "2",
                                 "--limit", "5", "--show", "1"]), 0)
        with self.assertRaises(SystemExit):
            S.main(["--ni", "2"])


@unittest.skipUnless(SLOW, "DMB_SLOW=1 のときだけ（総当たりとの突き合わせで ~30 秒）")
class TestExhaustive(unittest.TestCase):
    def test_matches_brute_force(self):
        """枝刈りありの探索と，枝刈りなしの総当たりで個数が一致すること。

        T(3,3)・値 2 通りなら 2^18 = 262144 通りを全部試せる。"""
        ni, nj, nv = 3, 3, 2
        K = D.torus(ni, nj)
        orbits = S.invariant_orbits(K, ni, nj)
        brute = 0
        for vals in itertools.product(range(nv), repeat=len(orbits)):
            f = {c: vals[k] for k, o in enumerate(orbits) for c in o}
            if D.MorseBott(K, f).is_dmb():
                brute += 1
        _, _, n_dmb = S.search(ni, nj, nv, sharp_only=True)
        self.assertEqual(n_dmb, brute)
        self.assertEqual(brute, 6522)

    def test_classification_of_sharp_functions(self):
        """**鋭い不変 DMBF は 2 種類しかない**（T(3,3) と T(4,3)，値 2 通り）:

          * (1 + t) と t(1 + t) の 2 本の臨界円周 — 滑らかな Morse-Bott 関数と同じ形
          * collection 1 つが P_t(T^2) を丸ごと担う自明な形（定数関数など）
        """
        for ni, nj in ((3, 3), (4, 3)):
            with self.subTest(ni=ni, nj=nj):
                found, _, _ = S.search(ni, nj, 2, sharp_only=True)
                shapes = {shape for (shape, _) in S.classify(found)}
                self.assertEqual(shapes, {((0, 1, 1), (1, 1)), ((1, 2, 1),)})
                counts = S.classify(found)
                trivial = counts[(((1, 2, 1),), ())]
                self.assertEqual(trivial, 2, "自明な形は定数関数の類だけ")
                self.assertGreater(counts[(((0, 1, 1), (1, 1)), ())], 600)


if __name__ == "__main__":
    unittest.main(verbosity=2)
