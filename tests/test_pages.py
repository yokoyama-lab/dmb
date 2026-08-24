#!/usr/bin/env python3
"""GitHub Pages 用の書き出し (pages.py) の検査。

数式を伏せる処理は依存ゼロなので常に検査する。サイト全体の書き出しは
dash / plotly / numpy が無ければ skip し，markdown の有無で results.html の
有無が変わることも見る。

    python3 -m unittest tests.test_pages -v
"""

import re
import tempfile
import unittest
from pathlib import Path

import pages

try:
    import dmb  # noqa: F401
    HAVE_DASH = True
except ImportError:                                          # pragma: no cover
    HAVE_DASH = False

try:
    import markdown  # noqa: F401
    HAVE_MARKDOWN = True
except ImportError:                                          # pragma: no cover
    HAVE_MARKDOWN = False


@unittest.skipUnless(HAVE_MARKDOWN, "markdown が入っていない")
class TestMarkdown(unittest.TestCase):
    """数式が Markdown に壊されないこと（添字の `_` 対・`*`・`\\{` の escape）。"""

    def test_math_survives(self):
        """実際に docs/results.md で壊れた形（添字の `_` 対・`*`・`\\{`）。"""
        body = pages.render_markdown(
            "$\\mathbb{Z}_{n_i}$ が $\\mathbb{Z}_{n_j}$ に，"
            "$H_*(\\mathbb{RP}^2)$，$n_i \\in \\{3,4\\}$。\n")
        self.assertIn("$\\mathbb{Z}_{n_i}$", body)
        self.assertIn("$\\mathbb{Z}_{n_j}$", body)
        self.assertIn("$H_*(\\mathbb{RP}^2)$", body)
        self.assertIn("\\{3,4\\}", body)
        self.assertNotIn("<em>", body)
        self.assertNotIn("@@MATH", body)

    def test_angle_bracket_in_math_is_escaped(self):
        """数式の中の不等号がタグの始まりと読まれないこと。"""
        body = pages.render_markdown("$m_2 = 1 > 0 = b_2$ かつ $a < b$ かつ $x \\& y$。\n")
        self.assertIn("&lt;", body)
        self.assertNotIn("$a < b$", body)
        self.assertIn("&amp;", body)

    def test_display_math_and_emphasis(self):
        body = pages.render_markdown("$$\\sum_C P_t(C) = P_t(K) + (1 + t) R(t)$$\n\n"
                                     "これは *強調* である。\n")
        self.assertIn("$$\\sum_C P_t(C) = P_t(K) + (1 + t) R(t)$$", body)
        self.assertIn("<em>強調</em>", body)

    def test_real_results_document(self):
        page = pages.results_page()
        self.assertIsNotNone(page)
        self.assertIn("MathJax", page)
        self.assertNotIn("@@MATH", page)


class TestMarkdownSource(unittest.TestCase):
    """原稿の側の検査（GitHub 上での見え方。依存ゼロなので常に走る）。

    GitHub の Markdown は数式より先に backslash escape を解くので，数式の中の
    `\\{` `\\#` は `{` `#` になってから KaTeX に渡る。`\\#` は
    「You can't use 'macro parameter character #' in math mode」で落ち，
    `\\{` は**黙って波括弧が消える**（集合の記法が崩れる）。Pages 側は
    pages.py が数式を退避するので壊れず，GitHub 上だけが壊れるので気づきにくい。
    英字だけの綴り（`\\lbrace` `\\lvert`）を使えばどちらでも読める。
    """

    #: CommonMark が escape として解く文字
    PUNCT = set("!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~")
    DOCS = sorted(Path(__file__).resolve().parent.parent.glob("*.md")) + \
        sorted((Path(__file__).resolve().parent.parent / "docs").glob("*.md"))

    def test_no_backslash_escaped_punctuation_in_math(self):
        bad = []
        for path in self.DOCS:
            text = path.read_text(encoding="utf-8")
            for m in pages.MATH_RE.finditer(text):
                for esc in re.finditer(r"\\(.)", m.group(0)):
                    if esc.group(1) in self.PUNCT:
                        line = text[:m.start()].count("\n") + 1
                        bad.append(f"{path.name}:{line} {m.group(0)[:40]}")
        self.assertEqual(bad, [], "GitHub 上で数式が壊れる書き方: " + "; ".join(bad))


class TestFigureSpecs(unittest.TestCase):
    """並べる図の指定が dmb.py の語彙と合っていること（skip でも常に見る）。"""

    def test_slugs_unique(self):
        slugs = [s["slug"] for s in pages.FIGURES]
        self.assertEqual(len(slugs), len(set(slugs)))

    def test_options_are_known(self):
        known = {"showColor", "showCollection", "showArrow", "showWeak"}
        for spec in pages.FIGURES:
            with self.subTest(spec["slug"]):
                self.assertLessEqual(set(spec["options"]), known)
                self.assertGreaterEqual(spec["ni"], 3)
                self.assertGreaterEqual(spec["nj"], 3)

    @unittest.skipUnless(HAVE_DASH, "dash / plotly / numpy が入っていない")
    def test_fkeys_exist(self):
        for spec in pages.FIGURES:
            self.assertIn(spec["fkey"], dmb.FUNCTIONS)


@unittest.skipUnless(HAVE_DASH, "dash / plotly / numpy が入っていない")
class TestBuild(unittest.TestCase):
    """書き出したサイトが，開ける形になっていること。"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.out = Path(cls.tmp.name) / "_site"
        pages.build(cls.out)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_index_links_every_figure(self):
        index = (self.out / "index.html").read_text(encoding="utf-8")
        for spec in pages.FIGURES:
            self.assertIn(f'href="{spec["slug"]}.html"', index)

    def test_each_figure_page_has_both_plots(self):
        for spec in pages.FIGURES:
            with self.subTest(spec["slug"]):
                page = (self.out / f"{spec['slug']}.html").read_text(encoding="utf-8")
                # 2 次元と 3 次元で 2 枚
                self.assertEqual(page.count("Plotly.newPlot"), 2)
                self.assertEqual(page.count("plotly-graph-div"), 2)
                # plotly.js は CDN から 1 回だけ読む
                self.assertEqual(page.count(pages.PLOTLY_CDN), 1)
                # Theorem 4.12 の検算が載っている
                self.assertIn("P_t(K)", page)
                self.assertIn("Σ_C P_t(C)", page)

    def test_no_absolute_paths(self):
        """相対リンクだけ（サブディレクトリ配信でも壊れない）。"""
        for path in self.out.glob("*.html"):
            with self.subTest(path.name):
                text = path.read_text(encoding="utf-8")
                self.assertNotIn('href="/', text)
                self.assertNotIn('src="/', text)

    def test_results_page_follows_markdown(self):
        self.assertEqual((self.out / "results.html").exists(), HAVE_MARKDOWN)

    def test_nojekyll_and_images(self):
        self.assertTrue((self.out / ".nojekyll").exists())
        self.assertTrue((self.out / "img" / "collections-2d.png").exists())


if __name__ == "__main__":
    unittest.main()
