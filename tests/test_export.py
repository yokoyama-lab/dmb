#!/usr/bin/env python3
"""図の書き出し (export.py) の検査。

TikZ 出力は依存ゼロなので常に検査する。LaTeX が入っていれば実際に組んでみる
（`DMB_LATEX=1`）。plotly を使う形式は無ければ skip。
"""

import os
import shutil
import subprocess
import tempfile
import unittest

import dmb_core as D
import export


def tikz_for(fkey="height", ni=4, nj=4, **kw):
    K = D.torus(ni, nj)
    f = export.FUNCTIONS[fkey](K, ni, nj)
    return export.tikz(K, f, ni, nj, **kw)


class TestTikz(unittest.TestCase):
    def test_well_formed(self):
        for fkey in export.FUNCTIONS:
            with self.subTest(fkey):
                src = tikz_for(fkey, arrows=True, weak=True)
                self.assertEqual(src.count("\\begin{tikzpicture}"), 1)
                self.assertEqual(src.count("\\end{tikzpicture}"), 1)
                self.assertNotIn("nan", src.lower())
                self.assertNotIn("None", src)
                for line in src.splitlines():
                    if line.strip().startswith(("\\fill", "\\draw", "\\node")):
                        self.assertTrue(line.rstrip().endswith(";"), line)

    def test_counts_match_the_complex(self):
        """三角形・辺・頂点の描画数が，持ち上げの数と一致する。"""
        ni, nj = 4, 5
        K = D.torus(ni, nj)
        pos = D.lifted_cells(ni, nj)
        src = tikz_for("height", ni, nj)
        want_fill = sum(len(pos[c]) for c in K.cells_of_dim(2))
        want_draw = sum(len(pos[c]) for c in K.cells_of_dim(1))
        want_point = sum(len(pos[c]) for c in K.cells_of_dim(0))
        self.assertEqual(src.count("\\fill["), want_fill)
        self.assertEqual(len([ln for ln in src.splitlines()
                              if ln.strip().startswith("\\draw[")
                              and "dmbarrow" not in ln]), want_draw)
        self.assertEqual(src.count("dmbpoint"), want_point + 1)   # +1 は style 定義

    def test_arrows_match_the_snc_pairs(self):
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = export.FUNCTIONS["dmf_min"](K, ni, nj)
        M = D.MorseBott(K, f)
        pos = D.lifted_cells(ni, nj)
        src = export.tikz(K, f, ni, nj, arrows=True)
        want = sum(len(pos[t]) for _, t in M.arrows())
        self.assertEqual(src.count("dmbarrow]"), want)

    def test_weak_circles_match(self):
        ni, nj = 4, 4
        K = D.torus(ni, nj)
        f = export.FUNCTIONS["height"](K, ni, nj)
        pos = D.lifted_cells(ni, nj)
        src = export.tikz(K, f, ni, nj, weak=True)
        want = sum(len(pos[c]) for c in D.MorseBott(K, f).weakly_critical())
        self.assertEqual(src.count("dmbweak]"), want)   # style 定義は "dmbweak/"

    def test_options(self):
        self.assertNotIn("\\fill[", tikz_for(fill="none"))
        self.assertIn("dmbval", tikz_for(fill="value"))
        self.assertNotIn("$", tikz_for(label="none"))
        self.assertIn("v_{0}", tikz_for(label="name"))
        self.assertIn("これは説明", tikz_for(caption="これは説明"))

    def test_deterministic(self):
        self.assertEqual(tikz_for(arrows=True, weak=True),
                         tikz_for(arrows=True, weak=True))

    def test_standalone_wrapper(self):
        src = export.tikz_standalone(tikz_for())
        self.assertTrue(src.startswith("\\documentclass"))
        self.assertIn("\\begin{document}", src)
        self.assertIn("\\end{document}", src)


class TestCLI(unittest.TestCase):
    def test_writes_a_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "fig.tex")
            self.assertEqual(export.main(["--format", "tikz", "--ni", "3", "--nj", "3",
                                          "-o", path]), 0)
            with open(path, encoding="utf-8") as fh:
                self.assertIn("\\begin{tikzpicture}", fh.read())

    def test_rejects_small_grids(self):
        with self.assertRaises(SystemExit):
            export.main(["--ni", "2"])

    def test_rejects_unknown_function(self):
        with self.assertRaises(SystemExit):
            export.main(["--function", "no-such"])


@unittest.skipUnless(os.environ.get("DMB_LATEX") and shutil.which("pdflatex"),
                     "DMB_LATEX=1 かつ pdflatex があるときだけ")
class TestLatexCompiles(unittest.TestCase):
    """書き出した TikZ が本当に組めること。"""

    def test_compiles(self):
        cases = [{"arrows": True, "weak": True},
                 {"fill": "none", "label": "none", "arrows": True},
                 {"fill": "value", "label": "name"}]
        with tempfile.TemporaryDirectory() as d:
            for i, kw in enumerate(cases):
                src = export.tikz_standalone(tikz_for("dmf_min", 4, 4, **kw))
                path = os.path.join(d, f"f{i}.tex")
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(src)
                r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                                    "-halt-on-error", f"f{i}.tex"],
                                   cwd=d, capture_output=True, text=True, timeout=180)
                pdf = os.path.join(d, f"f{i}.pdf")
                self.assertTrue(os.path.exists(pdf), r.stdout[-3000:])
                self.assertGreater(os.path.getsize(pdf), 1000)


class TestPlotlyFormats(unittest.TestCase):
    def test_html(self):
        try:
            import dmb  # noqa: F401
        except ImportError:
            self.skipTest("dash / plotly が無い")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "fig.html")
            self.assertEqual(export.main(["--format", "html", "--ni", "3", "--nj", "3",
                                          "-o", path]), 0)
            self.assertGreater(os.path.getsize(path), 10000)

    def test_reports_missing_dependency(self):
        real = export.plotly_figures

        def boom(*a, **k):
            raise ImportError("dash が無い（模擬）")

        export.plotly_figures = boom
        try:
            self.assertEqual(export.main(["--format", "html", "-o", "-"]), 2)
        finally:
            export.plotly_figures = real


if __name__ == "__main__":
    unittest.main(verbosity=2)
