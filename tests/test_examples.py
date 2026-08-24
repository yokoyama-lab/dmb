#!/usr/bin/env python3
"""examples.py / complexes.py / dmb_core.py の実行そのものの検査（自己検算つき）。

これらは「走らせて読む」ための台本なので，落ちないことと，中の検算が
すべて通ることを確かめる。
"""

import contextlib
import io
import unittest

import complexes
import dmb_core
import examples


def run(main, argv=None):
    """標準出力を捨てて main を走らせ，終了コードを返す。"""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = main(argv)
    return code, buf.getvalue()


class TestRunnableScripts(unittest.TestCase):
    def test_examples_tutorial(self):
        code, out = run(examples.main, ["examples.py", "tutorial"])
        self.assertEqual(code, 0, out[-2000:])
        self.assertIn("すべての例が期待どおり", out)

    def test_examples_strength(self):
        code, out = run(examples.main, ["examples.py", "strength"])
        self.assertEqual(code, 0, out[-2000:])
        for must in ("RP^2", "クラインの壺", "Z_n 不変な DMF", "臨界円周"):
            self.assertIn(must, out)

    def test_examples_bad_argument(self):
        code, _ = run(examples.main, ["examples.py", "no-such-group"])
        self.assertEqual(code, 2)

    def test_complexes_catalogue(self):
        code, out = run(complexes.main)
        self.assertEqual(code, 0, out[-2000:])
        self.assertIn("すべて期待どおり", out)

    def test_dmb_core_report(self):
        code, out = run(dmb_core.main, ["dmb_core.py", "4", "4"])
        self.assertEqual(code, 0, out[-3000:])
        self.assertIn("すべての検査を通過", out)

    def test_dmb_core_table(self):
        code, out = run(dmb_core.main, ["dmb_core.py", "--table"])
        self.assertEqual(code, 0)
        self.assertIn("inv.DMF", out)


class TestStrengthClaims(unittest.TestCase):
    """examples.py が主張している「DMBT ならでは」の内容を，出力とは独立に確かめる。"""

    def test_no_dmf_is_sharp_on_nonorientable_closed_surfaces(self):
        """閉非向き付け曲面では m_2 ≥ 1 > b_2 = 0 なので M(t) ≠ P_t(K)。

        ランダムな離散モース関数 60 個で R(t) ≠ 0 を確かめる。"""
        import random
        for name, K in (("RP^2", complexes.projective_plane()),
                        ("Klein", complexes.klein_bottle(5, 5))):
            b = dmb_core.betti(K.cells)
            self.assertLess(len(b), 3, name)           # b_2 = 0（有理数係数）
            for t in range(60):
                f = dmb_core.random_matching_dmf(K, random.Random(dmb_core.hash_seed(name, 'sharp', t)))
                r = dmb_core.MorseBott(K, f).report()
                if not r["is_dmf"]:
                    continue
                self.assertGreaterEqual(len(r["M"]), 3, (name, t))   # m_2 ≥ 1
                self.assertTrue(dmb_core.poly_trim(r["R_M"]), (name, t))

    def test_dmb_is_sharp_on_nonorientable_closed_surfaces(self):
        """一方 DMBF は鋭くできる（定数関数だけでなく非自明なものも）。"""
        K = complexes.projective_plane()
        self.assertTrue(dmb_core.MorseBott(K, dmb_core.constant_fn(K)).report()["MB_sharp"])
        h = {1: 2, 2: 2, 3: 2, 4: 2, 5: 1, 6: 0}
        M = dmb_core.MorseBott(K, dmb_core.max_extension(K, h.__getitem__))
        r = M.report()
        self.assertTrue(r["is_dmb"])
        self.assertTrue(r["MB_sharp"])
        self.assertGreater(len(M.collections()), 1, "非自明（collection が 2 つ以上）")


if __name__ == "__main__":
    unittest.main(verbosity=2)
