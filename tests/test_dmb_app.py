#!/usr/bin/env python3
"""可視化 (dmb.py) の検査。dash / plotly / numpy が無ければ skip する。

    python3 -m unittest test_dmb_app -v
    DMB_RENDER=1 python3 -m unittest test_dmb_app   # 画像に書き出して形も見る
                                                    # （kaleido と Chrome が要る）
"""

import os
import unittest

import dmb_core as D

try:
    import dmb
    HAVE_DASH = True
except ImportError as exc:                                   # pragma: no cover
    HAVE_DASH = False
    IMPORT_ERROR = exc


@unittest.skipUnless(HAVE_DASH, "dash / plotly / numpy が入っていない")
class TestCallback(unittest.TestCase):
    """コールバックが全ての操作の組合せで落ちないこと。"""

    def test_all_option_combinations(self):
        checks = [[], ["showColor"], ["showCollection"], ["showArrow"], ["showWeak"],
                  ["showColor", "showCollection", "showArrow", "showWeak"]]
        n = 0
        for fkey in dmb.FUNCTIONS:
            for radio in ("cell", "value", "collection", "none"):
                for opts in checks:
                    with self.subTest(fkey=fkey, radio=radio, opts=opts):
                        out = dmb.update(fkey, "YlGn", radio, 4, 4, 2, opts)
                        self.assertEqual(len(out), 5)
                        fig2, fig3, scale, style, report = out
                        self.assertTrue(fig2.data)
                        self.assertTrue(fig3.data)
                        self.assertIn("display", style)
                        self.assertIn("P_t(K)", report)
                        n += 1
        self.assertGreater(n, 100)

    def test_grid_sizes(self):
        for ni in (3, 4, 7):
            for nj in (3, 5, 6):
                for smooth in (1, 3):
                    with self.subTest(ni=ni, nj=nj, smooth=smooth):
                        dmb.update("height", "YlGn", "value", ni, nj, smooth,
                                   ["showCollection", "showArrow"])

    def test_input_guards(self):
        """入力欄が空 (None) や範囲外でも落ちず，範囲に丸められる。"""
        for ni, nj, smooth in ((None, None, None), (0, 0, 0), (999, 999, 999),
                               (2, 2, -5), ("4", "4", "2")):
            with self.subTest(ni=ni, nj=nj, smooth=smooth):
                *_, report = dmb.update("height", "YlGn", "value", ni, nj, smooth, [])
                self.assertIn("T(", report)

    def test_unavailable_function_falls_back(self):
        """正方格子でない場合の dmf.py の関数など，使えない選択でも落ちない。"""
        *_, report = dmb.update("dmf_py", "YlGn", "value", 5, 4, 2, [])
        self.assertIn("使えない", report)

    def test_unknown_function_key(self):
        *_, report = dmb.update("no-such-function", "YlGn", "value", 4, 4, 2, [])
        self.assertIn("使えない", report)

    def test_report_matches_the_core(self):
        """画面に出る報告が dmb_core の計算と一致すること。"""
        K = D.torus(5, 4)
        r = D.MorseBott(K, D.height_fn(K, 5, 4)).report()
        *_, report = dmb.update("height", "YlGn", "value", 5, 4, 2, [])
        self.assertIn(f"Σ_C P_t(C) = {D.poly_str(r['MB_sum'])}", report)
        self.assertIn(f"P_t(K)     = {D.poly_str(r['P_K'])}", report)
        self.assertIn("R(t)       = 0", report)


@unittest.skipUnless(HAVE_DASH, "dash / plotly / numpy が入っていない")
class TestDrawingGeometry(unittest.TestCase):
    """図形そのものの整合性（目視の代わりに数値で見る）。"""

    def test_arrows_are_between_adjacent_cells(self):
        """矢印は隣り合うセルの重心を結ぶので短い（貼り合わせをまたいでも）。"""
        for ni, nj in ((3, 3), (5, 4), (4, 6)):
            K = D.torus(ni, nj)
            pos = D.lifted_cells(ni, nj)
            M = D.MorseBott(K, D.canonical_dmf(K))
            for s, t in M.arrows():
                for p1 in [D.centroid(q) for q in pos[t]]:
                    p0 = dmb.nearest_centroid(pos[s], p1)
                    self.assertLessEqual(abs(p1[0] - p0[0]) + abs(p1[1] - p0[1]), 1.0)

    def test_all_drawn_points_are_inside_the_padded_frame(self):
        """描画点が軸の範囲からはみ出さない（ラベルの見切れ対策）。"""
        ni, nj = 5, 4
        fig2, *_ = dmb.update("height", "YlGn", "value", ni, nj, 2,
                              ["showCollection", "showWeak", "showArrow"])
        xr = fig2.layout.xaxis.range
        yr = fig2.layout.yaxis.range
        self.assertLess(xr[0], 0)
        self.assertGreater(xr[1], ni)
        self.assertLess(yr[0], 0)
        self.assertGreater(yr[1], nj)
        for tr in fig2.data:
            for v, lo, hi in ((tr.x, xr[0], xr[1]), (tr.y, yr[0], yr[1])):
                for q in (v or ()):
                    if q is not None:
                        self.assertTrue(lo <= q <= hi, (tr.mode, q, lo, hi))

    def test_3d_points_lie_on_or_just_above_the_torus(self):
        """3 次元の点が本当にトーラス上（辺・頂点は法線方向にわずかに浮く）にある。"""
        ni, nj = 5, 4
        _, fig3, *_ = dmb.update("height", "YlGn", "none", ni, nj, 3, ["showCollection"])
        for tr in fig3.data:
            xs, ys, zs = tr.x, tr.y, tr.z
            for x, y, z in zip(xs, ys, zs):
                if x is None:
                    continue
                # (√(x²+y²) − R)² + z² = r²  （浮かせた分だけ半径が大きくなる）
                rad = ((x ** 2 + y ** 2) ** 0.5 - dmb.R) ** 2 + z ** 2
                self.assertGreaterEqual(rad ** 0.5, dmb.r - 1e-6)
                self.assertLessEqual(rad ** 0.5, dmb.r + 4 * dmb.LIFT)

    def test_mesh_is_visible_over_the_fill(self):
        """着色していても三角形分割が見えること: 塗りは半透明で，
        全ての辺に不透明な縁取りのトレースがある。"""
        for opts in (["showColor"], ["showCollection"]):
            fig2, *_ = dmb.update("height", "YlGn", "none", 4, 4, 2, opts)
            fills = [tr.fillcolor for tr in fig2.data if tr.fill == 'toself']
            self.assertTrue(fills)
            for col in fills:
                self.assertTrue(col.startswith("rgba"), col)
                self.assertLess(float(col.rsplit(",", 1)[1].rstrip(")")), 1.0)
            casings = [tr for tr in fig2.data
                       if tr.mode == 'lines' and tr.line.color == dmb.CASING]
            self.assertEqual(len(casings), 1)   # 着色時は必ず 1 本引かれる
            # 縁取りは全ての辺（の全ての複製）を含む
            n_seg = sum(len(p) for p in D.lifted_cells(4, 4).values()
                        if len(p[0]) == 2)
            self.assertEqual(len([q for q in casings[0].x if q is None]), n_seg)

    def test_no_casing_when_not_coloured(self):
        """着色していないときは縁取りを引かない（元の dmf.py と同じ軽い線）。"""
        fig2, *_ = dmb.update("dmf_min", "YlGn", "none", 4, 4, 2, ["showArrow"])
        self.assertFalse([tr for tr in fig2.data
                          if tr.mode == 'lines' and tr.line.color == dmb.CASING])

    def test_label_positions_are_cell_centroids(self):
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        pos = D.lifted_cells(ni, nj)
        fig2, *_ = dmb.update("cellname-check", "YlGn", "cell", ni, nj, 2, [])
        texts = [tr for tr in fig2.data if tr.mode == 'text']
        self.assertEqual(len(texts), 1)
        want = {D.centroid(pos[c][0]) for c in K.cells}
        got = set(zip(texts[0].x, texts[0].y))
        self.assertEqual(got, want)


@unittest.skipUnless(HAVE_DASH and os.environ.get("DMB_RENDER"),
                     "DMB_RENDER=1 のときだけ（kaleido と Chrome が要る）")
class TestRendering(unittest.TestCase):
    """実際に画像に書き出せること（形の目視確認に使ったのと同じ経路）。"""

    def test_write_png(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            fig2, fig3, *_ = dmb.update("height", "YlGn", "value", 5, 4, 3,
                                        ["showCollection", "showWeak"])
            for i, fig in enumerate((fig2, fig3)):
                path = os.path.join(d, f"fig{i}.png")
                fig.write_image(path, width=600, height=600)
                self.assertGreater(os.path.getsize(path), 5000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
