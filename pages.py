#!/usr/bin/env python3
"""GitHub Pages 用の静的サイトを書き出す。

dmb.py は Dash（Flask サーバ）なので Pages には置けない。ここでは同じ図を
plotly の自己完結 HTML として**あらかじめ書き出し**，索引を付けて並べる。
回転・ズーム・ホバーはブラウザ側で効くが，関数や分割数の切り替えはできない
（切り替えたい場合は手元で `python3 dmb.py` を動かす）。

    python3 pages.py -o _site        # _site/index.html ができる
    python3 -m http.server -d _site  # 手元で確認

要 dash / plotly / numpy（`pip install -e ".[app]"`）。
docs/results.md の HTML 化には markdown（`pip install -e ".[pages]"`）も使うが，
無ければその項目だけ落として GitHub 上の版へリンクする。
"""

import argparse
import html
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
REPO = "https://github.com/yokoyama-lab/dmb"
PAPER = "https://arxiv.org/abs/2511.07864"

# 並べる図。dmb.py のコントロールのうち「見せたい組み合わせ」を固定したもの。
# options は dmb.update の check-options と同じ語彙。
FIGURES = [
    {"slug": "height", "fkey": "height", "ni": 5, "nj": 4,
     "options": ["showCollection", "showWeak"],
     "blurb": "回転対称な高さ関数。値の等しいセルが collection をなす典型例。"},
    {"slug": "height-refined", "fkey": "height_refined", "ni": 5, "nj": 4,
     "options": ["showCollection", "showWeak"],
     "blurb": "高さ関数の細分版。collection が細かく割れる。"},
    {"slug": "constant", "fkey": "constant", "ni": 4, "nj": 4,
     "options": ["showCollection", "showWeak"],
     "blurb": "定数関数。複体全体がただ 1 つの collection になる自明な DMBF。"},
    {"slug": "morsified", "fkey": "morsified", "ni": 5, "nj": 4,
     "options": ["showCollection", "showArrow", "showWeak"],
     "blurb": "高さ関数の Morsification。DMBF から DMF へ落とすと矢印が現れる。"},
    {"slug": "dmf-min", "fkey": "dmf_min", "ni": 5, "nj": 4,
     "options": ["showCollection", "showArrow", "showWeak"],
     "blurb": "tree-cotree で作る臨界セル 4 個の DMF（回転対称ではない）。"},
    {"slug": "dmf-invariant", "fkey": "dmf_invariant", "ni": 5, "nj": 4,
     "options": ["showCollection", "showArrow", "showWeak"],
     "blurb": "回転対称な DMF。対称性を保つと臨界セルは 4·n_i 個まで増える。"},
    {"slug": "arrowed", "fkey": "arrowed", "ni": 5, "nj": 4,
     "options": ["showCollection", "showArrow", "showWeak"],
     "blurb": "矢印を持つ鋭い DMBF。高さ関数と同じ R(t) = 0 だが，"
              "最上位の帯を持ち上げることで矢印が n_i 本現れる（§3.6）。"},
    {"slug": "dmf-py", "fkey": "dmf_py", "ni": 4, "nj": 4,
     "options": ["showCollection", "showArrow", "showWeak"],
     "blurb": "dmf.py の calcDMF が作る DMF（正方格子のみ）。"},
]

SMOOTH = 3          # 3 次元表示の細分数
PLOTLY_CDN = "https://cdn.plot.ly/plotly-3.7.0.min.js"

CSS = """
:root { --fg:#1a1a1a; --muted:#666; --line:#e2e2e2; --bg:#fff; --code:#f6f6f6; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--fg); line-height:1.7;
       font-family: -apple-system, "Hiragino Sans", "Noto Sans JP", sans-serif; }
.wrap { max-width: 1180px; margin:0 auto; padding: 0 20px 64px; }
header.site { border-bottom:1px solid var(--line); margin-bottom:32px; }
header.site .wrap { padding-top:20px; padding-bottom:16px; }
header.site a { color:inherit; text-decoration:none; }
nav a { margin-right:18px; font-size:14px; color:var(--muted); }
nav a:hover { color:var(--fg); }
h1 { font-size:26px; margin:8px 0; }
h2 { font-size:20px; margin-top:40px; border-bottom:1px solid var(--line);
     padding-bottom:6px; }
p, li { font-size:15px; }
.muted { color:var(--muted); font-size:14px; }
.figs { display:flex; flex-wrap:wrap; gap:24px; }
.figs > div { flex:1 1 480px; min-width:0; }
.figs figcaption, .cap { text-align:center; font-size:14px; color:var(--muted);
                         margin-bottom:4px; }
pre.report { background:var(--code); padding:14px; border-radius:6px;
             font-size:13px; line-height:1.55; overflow-x:auto; white-space:pre-wrap; }
table.gallery { border-collapse:collapse; width:100%; }
table.gallery td, table.gallery th { border-bottom:1px solid var(--line);
                                     padding:10px 8px; text-align:left; vertical-align:top; }
table.gallery td:first-child { white-space:nowrap; }
img { max-width:100%; height:auto; }
code { background:var(--code); padding:1px 5px; border-radius:4px; font-size:13px; }
"""


def document(title, body, *, plotly=False, mathjax=False):
    """共通の枠。生成物は単体で開けるようにしておく（相対リンクのみ）。"""
    head = [f"<title>{html.escape(title)}</title>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            f"<style>{CSS}</style>"]
    if plotly:
        head.append(f'<script src="{PLOTLY_CDN}" charset="utf-8"></script>')
    if mathjax:
        head.append("<script>window.MathJax={tex:{inlineMath:[['$','$'],"
                    "['\\\\(','\\\\)']],displayMath:[['$$','$$'],"
                    "['\\\\[','\\\\]']]}};</script>")
        head.append('<script id="MathJax-script" async '
                    'src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">'
                    "</script>")
    return ("<!doctype html>\n<html lang=\"ja\">\n<head>\n" + "\n".join(head) +
            "\n</head>\n<body>\n"
            '<header class="site"><div class="wrap">'
            '<h1><a href="index.html">離散モースボット理論 — トーラス上の例</a></h1>'
            '<nav><a href="index.html">図の一覧</a>'
            '<a href="results.html">計算結果</a>'
            f'<a href="{REPO}">GitHub</a><a href="{PAPER}">arXiv:2511.07864</a></nav>'
            "</div></header>\n"
            f'<div class="wrap">\n{body}\n</div>\n</body>\n</html>\n')


def figure_page(spec, titles):
    """図 1 枚ぶんのページ（2 次元・3 次元・検算の報告）を作る。"""
    import dmb  # noqa: PLC0415  plotly / dash が要るのでここで読む

    fig2, fig3, _, _, report = dmb.update(
        spec["fkey"], "YlGn", "value", spec["ni"], spec["nj"], SMOOTH, spec["options"])
    divs = []
    for fig, caption in ((fig2, "基本領域（三角形分割）"), (fig3, "トーラス")):
        fig.update_layout(width=None, autosize=True)
        divs.append(
            f'<div><div class="cap">{caption}</div>' +
            fig.to_html(full_html=False, include_plotlyjs=False,
                        default_width="100%", default_height="560px",
                        config={"responsive": True}) +
            "</div>")
    title = titles.get(spec["fkey"], spec["fkey"])
    body = (f"<h2>{html.escape(title)}</h2>\n"
            f'<p>{html.escape(spec["blurb"])}</p>\n'
            f'<p class="muted">T({spec["ni"]}, {spec["nj"]})　'
            f'表示: {", ".join(spec["options"])}　'
            f'再現: <code>python3 dmb.py</code> で同じ設定を選ぶ</p>\n'
            f'<div class="figs">{"".join(divs)}</div>\n'
            f'<pre class="report">{html.escape(report)}</pre>\n')
    return document(title, body, plotly=True)


MATH_RE = re.compile(r"(\$\$.+?\$\$|\$[^$\n]+?\$)", re.S)


def render_markdown(text):
    """数式を伏せてから Markdown に通す。

    Markdown は数式を知らないので，そのまま通すと壊れる。実際に docs/results.md で
    起きたもの: `\\mathbb{Z}_{n_i}` と `\\mathbb{Z}_{n_j}` の `_` が対になって
    `<em>` に化ける，`H_*` の `*` も同様，`\\{` の backslash が食われる。
    そこで `$...$` を先に取り出して placeholder に置き換え，変換後に戻す。
    戻すときに `<` `&` を実体参照にするのは，数式の中の不等号を
    タグの始まりと読まれないため（MathJax は textContent を見るので影響しない）。"""
    import markdown  # noqa: PLC0415

    spans = []

    def stash(m):
        spans.append(m.group(0))
        return f"@@MATH{len(spans) - 1}@@"

    body = markdown.markdown(MATH_RE.sub(stash, text),
                             extensions=["fenced_code", "tables", "toc"])
    for k, span in enumerate(spans):
        body = body.replace(f"@@MATH{k}@@", html.escape(span, quote=False))
    return body


def results_page():
    """docs/results.md を HTML 化する。markdown が無ければ None。"""
    src = DOCS / "results.md"
    if not src.exists():
        return None
    try:
        body = render_markdown(src.read_text(encoding="utf-8"))
    except ImportError:
        print("markdown が無いので results.html は作らない "
              '（pip install -e ".[pages]"）', file=sys.stderr)
        return None
    return document("計算結果", body, mathjax=True)


def index_page(specs, titles, has_results):
    rows = []
    for spec, title in ((s, titles.get(s["fkey"], s["fkey"])) for s in specs):
        rows.append(f'<tr><td><a href="{spec["slug"]}.html">{html.escape(title)}</a></td>'
                    f'<td>T({spec["ni"]}, {spec["nj"]})</td>'
                    f'<td>{html.escape(spec["blurb"])}</td></tr>')
    results_link = ('<a href="results.html">計算結果（層 A / 層 B）</a>' if has_results
                    else f'<a href="{REPO}/blob/main/docs/results.md">'
                         "計算結果（層 A / 層 B）</a>")
    imgs = ""
    if (DOCS / "img").exists():
        imgs = ('<div class="figs">'
                '<div><img src="img/collections-2d.png" alt="collections 2D"></div>'
                '<div><img src="img/collections-3d.png" alt="collections 3D"></div>'
                "</div>")
    body = (
        "<p>Nishikawa–Yokoyama, <em>On discrete Morse-Bott theory</em> "
        f'(<a href="{PAPER}">arXiv:2511.07864</a>) の定義に沿って，'
        "三角形分割したトーラス $T(n_i, n_j)$ 上の離散モースボット関数を計算し，"
        "図示したもの。</p>"
        f"{imgs}"
        "<h2>図の一覧</h2>"
        "<p>各ページは plotly の静的な書き出しで，回転・ズーム・ホバーが効く。"
        "関数や分割数を切り替えたい場合は手元で "
        f'<code>python3 dmb.py</code> を動かす（<a href="{REPO}">導入手順</a>）。</p>'
        '<table class="gallery"><tr><th>関数</th><th>分割</th><th>説明</th></tr>'
        + "".join(rows) + "</table>"
        "<h2>ほかの資料</h2><ul>"
        f"<li>{results_link}</li>"
        f'<li><a href="{REPO}/blob/main/README.md">README（日本語）</a> / '
        f'<a href="{REPO}/blob/main/README.en.md">README (English)</a></li>'
        f'<li><a href="{REPO}">ソース一式</a>'
        "（理論の計算 <code>dmb_core.py</code> は標準ライブラリだけで動く）</li>"
        "</ul>")
    return document("離散モースボット理論 — トーラス上の例", body, mathjax=True)


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    import dmb  # noqa: PLC0415

    titles = dmb.FUNCTIONS
    built = []
    for spec in FIGURES:
        page = figure_page(spec, titles)
        (out / f"{spec['slug']}.html").write_text(page, encoding="utf-8")
        built.append(spec)
        print(f"書き出した: {out / (spec['slug'] + '.html')}", file=sys.stderr)

    results = results_page()
    if results is not None:
        (out / "results.html").write_text(results, encoding="utf-8")
        print(f"書き出した: {out / 'results.html'}", file=sys.stderr)

    if (DOCS / "img").exists():
        shutil.copytree(DOCS / "img", out / "img", dirs_exist_ok=True)

    (out / "index.html").write_text(
        index_page(built, titles, results is not None), encoding="utf-8")
    (out / ".nojekyll").write_text("", encoding="utf-8")   # _ 始まりの名前を消さない
    print(f"書き出した: {out / 'index.html'}", file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="GitHub Pages 用の静的サイトを書き出す")
    ap.add_argument("-o", "--out", default="_site", help="出力先（既定 _site）")
    args = ap.parse_args(argv)
    try:
        return build(Path(args.out))
    except ImportError as exc:
        print(f"図の書き出しには dash / plotly / numpy が要る: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
