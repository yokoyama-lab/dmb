#!/usr/bin/env python3
"""滑らかな Morse-Bott 関数の離散化計画 (docs/smoothing.md) の層 A 証拠。外部依存なし。

準位に適合した三角形分割の上で max_extension が，滑らかな側の Morse-Bott
多項式（剰余 R(t)・係数体依存を含めて）を再現することを，モデルケースで検査する。

    python3 -m unittest tests.test_smoothing -v
"""

import unittest

import complexes as X
import dmb_core as D


def spindle(n):
    """紡錘形 S²: 赤道 n 頂点（h=0）と両極 N, S（h=1）。

    滑らかな側の対応物は「赤道が指数 0 の臨界円周・両極が指数 2 の非退化な
    最大点」の Morse-Bott 関数。MB_t = (1+t) + 2t², R(t) = t（等号でない）。"""
    eq = [f"e{i}" for i in range(n)]
    tris = [("N", eq[i], eq[(i + 1) % n]) for i in range(n)] + \
           [("S", eq[i], eq[(i + 1) % n]) for i in range(n)]
    K = D.Complex(tris)
    f = D.max_extension(K, lambda v: 0 if str(v).startswith("e") else 1)
    return K, f


class TestSmoothingEvidence(unittest.TestCase):
    def test_spindle_reproduces_smooth_morse_bott_polynomial(self):
        """紡錘形 S²: Σ_C P_t(C) = (1+t) + 2t²，R(t) = t（滑らかな剰余まで一致）。"""
        for n in (3, 4, 6):
            with self.subTest(n=n):
                K, f = spindle(n)
                r = D.MorseBott(K, f).report()
                self.assertTrue(r["is_dmb"])
                self.assertEqual(r["MB_sum"], [1, 1, 2])       # 1 + t + 2t²
                self.assertEqual(r["P_K"], [1, 0, 1])          # S²
                self.assertEqual(D.poly_trim(r["R_MB"]), [0, 1])   # R = t
                shapes = sorted(tuple(D.poly_trim(D.betti(C)))
                                for C in r["reduced_collections"]
                                if D.poly_trim(D.betti(C)))
                # 赤道の円周 (1+t) と，極ごとの t² が 2 つ
                self.assertEqual(shapes, [(0, 0, 1), (0, 0, 1), (1, 1)])

    def test_cylinder_regular_bands_vanish(self):
        """円筒の高さ: 正則な帯（臨界点を含まない準位帯）の寄与はすべて 0（補題 L3）。

        境界つき多様体でもそのまま成り立つ。collapse は使っていない。"""
        K = X.annulus(4, 4)
        f = D.max_extension(K, lambda v: v[1])
        M = D.MorseBott(K, f)
        r = M.report()
        self.assertTrue(r["is_dmb"])
        self.assertEqual(D.poly_trim(r["R_MB"]), [])           # R = 0（鋭い）
        shapes = sorted(tuple(D.poly_trim(M.betti(C))) for C in M.collections())
        self.assertEqual(shapes.count(()), 4)                  # 正則帯 4 つが 0
        self.assertIn((1, 1), shapes)                          # 下端の円周 1+t

    def test_klein_height_matches_smooth_coefficient_dependence(self):
        """クラインの壺の高さ: 上端の臨界円周の負法束が向き付け不可能なので，
        滑らかな側では Q で R = t・F₂ で R = 0。離散側も同じになる。

        予想（docs/smoothing.md §2）の向き付けの但し書きが必要かつ十分である
        ことの証拠。"""
        K = X.klein_bottle_sym(3, 4)
        f = D.max_extension(K, lambda v: min(v[1] % 4, (-v[1]) % 4))
        r_q = D.MorseBott(K, f).report()
        r_f2 = D.MorseBott(K, f, p=2).report()
        self.assertTrue(r_q["is_dmb"])
        self.assertEqual(r_q["MB_sum"], [1, 2, 1])             # Σ は係数体に依らない
        self.assertEqual(r_f2["MB_sum"], [1, 2, 1])
        self.assertEqual(D.poly_trim(r_q["R_MB"]), [0, 1])     # Q: R = t
        self.assertEqual(D.poly_trim(r_f2["R_MB"]), [])        # F₂: R = 0


if __name__ == "__main__":
    unittest.main()
