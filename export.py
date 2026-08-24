#!/usr/bin/env python3
"""図の書き出し。TikZ（依存ゼロ）と，plotly があれば HTML / PNG / SVG。

論文の図は tikzpicture なので，**LaTeX にそのまま貼れる TikZ 出力を依存ゼロで**
作れるようにしてある（plotly も Chrome も要らない）。

    python3 export.py --format tikz --function height --ni 5 --nj 4 -o fig.tex
    python3 export.py --format tikz --function dmf_min --arrows --no-fill -o vf.tex
    python3 export.py --format html --function height -o fig.html      # 要 plotly
    python3 export.py --format png  --function height -o fig.png       # 要 kaleido

出力した TikZ は tikz パッケージだけで組める（xcolor は tikz が読み込む）。
"""

import argparse
import sys

import dmb_core as core

FUNCTIONS = {
    "height": lambda K, ni, nj: core.height_fn(K, ni, nj),
    "height_refined": lambda K, ni, nj: core.height_fn(K, ni, nj, refine=True),
    "constant": lambda K, ni, nj: core.constant_fn(K),
    "morsified": lambda K, ni, nj: core.morsify(K, core.height_fn(K, ni, nj)),
    "dmf_min": lambda K, ni, nj: core.canonical_dmf(K),
    "dmf_invariant": lambda K, ni, nj: core.invariant_dmf(K, ni, nj),
    "arrowed": lambda K, ni, nj: core.arrowed_dmbf(K, ni, nj),
}

# collection の塗り分けに使う色（RGB 0-255）。TikZ では \definecolor で定義する
PALETTE = [
    (228, 26, 28), (55, 126, 184), (77, 175, 74), (152, 78, 163),
    (255, 127, 0), (166, 86, 40), (247, 129, 191), (153, 153, 153),
]


# ------------------------------------------------------------------ TikZ


def tikz(K, f, ni, nj, scale=1.2, fill="collection", label="value",
         arrows=False, weak=False, caption=None):
    """基本領域 [0,ni] x [0,nj] の三角形分割を tikzpicture として書き出す。

    fill  : "collection" | "value" | "none"
    label : "value" | "name" | "collection" | "none"
    arrows: strictly noncritical pair を矢印で描く
    weak  : weakly critical なセルを丸で囲む
    """
    M = core.MorseBott(K, f)
    pos = core.lifted_cells(ni, nj)
    names = core.torus_names(ni, nj)
    collections = sorted(M.collections(), key=lambda L: (-len(L), L[0]))
    index = {c: k for k, L in enumerate(collections) for c in L}
    lo, hi = min(f.values()), max(f.values())
    weakset = set(M.weakly_critical())

    def colour(c):
        if fill == "collection":
            return f"dmbc{index[c] % len(PALETTE)}"
        if fill == "value":
            return "dmbval"
        return None

    def opacity(c):
        if fill == "value":
            span = (hi - lo) or 1
            return 0.15 + 0.55 * (f[c] - lo) / span
        return 0.35

    def text(c):
        if label == "value":
            return str(f[c])
        if label == "name":
            return names[c].replace("_", "_{") + "}" if "_" in names[c] else names[c]
        if label == "collection":
            return str(index[c])
        return None

    out = []
    out.append("% dmb: python3 export.py --format tikz "
               f"--ni {ni} --nj {nj} で生成")
    out.append("% 必要なパッケージ: \\usepackage{tikz}")
    for k, (r, g, b) in enumerate(PALETTE):
        out.append(f"\\definecolor{{dmbc{k}}}{{RGB}}{{{r},{g},{b}}}")
    out.append("\\definecolor{dmbval}{RGB}{35,110,60}")
    out.append(f"\\begin{{tikzpicture}}[scale={scale},")
    out.append("    dmbpoint/.style={fill,shape=circle,inner sep=1.1pt,outer sep=0pt},")
    out.append("    dmbarrow/.style={->,>=stealth,thick,gray!70!black},")
    out.append("    dmbweak/.style={draw=red,circle,inner sep=1.6pt,thick}]")

    out.append("  % 2-セル")
    for c in K.cells_of_dim(2):
        col = colour(c)
        for pts in pos[c]:
            path = " -- ".join(f"({x},{y})" for x, y in pts)
            if col:
                out.append(f"  \\fill[{col},opacity={opacity(c):.2f}] {path} -- cycle;")

    out.append("  % 1-セル")
    for c in K.cells_of_dim(1):
        col = colour(c) if fill == "collection" else None
        style = f"[{col}!70!black,thick]" if col else "[gray!60]"
        for (x0, y0), (x1, y1) in pos[c]:
            out.append(f"  \\draw{style} ({x0},{y0}) -- ({x1},{y1});")

    out.append("  % 0-セル")
    for c in K.cells_of_dim(0):
        col = colour(c) if fill == "collection" else None
        style = f"dmbpoint,{col}!70!black" if col else "dmbpoint"
        for (x, y), in pos[c]:
            out.append(f"  \\node[{style}] at ({x},{y}) {{}};")

    if arrows:
        out.append("  % 矢印（strictly noncritical pair）")
        for s, t in M.arrows():
            for pts in pos[t]:
                p1 = core.centroid(pts)
                p0 = min((core.centroid(q) for q in pos[s]),
                         key=lambda q: (q[0] - p1[0]) ** 2 + (q[1] - p1[1]) ** 2)
                mid = (p0[0] + 0.88 * (p1[0] - p0[0]), p0[1] + 0.88 * (p1[1] - p0[1]))
                out.append(f"  \\draw[dmbarrow] ({p0[0]:.3f},{p0[1]:.3f}) -- "
                           f"({mid[0]:.3f},{mid[1]:.3f});")

    if weak:
        out.append("  % weakly critical なセル")
        for c in K.cells:
            if c not in weakset:
                continue
            for pts in pos[c]:
                x, y = core.centroid(pts)
                out.append(f"  \\node[dmbweak] at ({x:.3f},{y:.3f}) {{}};")

    if label != "none":
        out.append("  % ラベル")
        for c in K.cells:
            s = text(c)
            if s is None:
                continue
            x, y = core.centroid(pos[c][0])
            out.append(f"  \\node[font=\\scriptsize] at ({x:.3f},{y:.3f}) {{${s}$}};")

    if caption:
        out.append(f"  \\node[below] at ({ni / 2},-0.4) {{{caption}}};")
    out.append("\\end{tikzpicture}")
    return "\n".join(out) + "\n"


def tikz_standalone(body):
    """そのまま latexmk できる最小の文書に包む。"""
    return ("\\documentclass[tikz,border=5pt]{standalone}\n"
            "\\usepackage{tikz}\n"
            "\\begin{document}\n" + body + "\\end{document}\n")


# ------------------------------------------------------ plotly を使う書き出し


def plotly_figures(fkey, ni, nj, smooth, options):
    """dmb.py のコールバックを呼んで figure を作る（plotly と dash が要る）。"""
    import dmb  # noqa: PLC0415
    fig2, fig3, _, _, report = dmb.update(fkey, "YlGn", "value", ni, nj, smooth, options)
    return fig2, fig3, report


# ------------------------------------------------------------------ CLI


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="離散モースボット関数の図を書き出す",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("--format", default="tikz",
                    choices=["tikz", "tikz-standalone", "html", "png", "svg"])
    ap.add_argument("--function", default="height", choices=sorted(FUNCTIONS))
    ap.add_argument("--ni", type=int, default=5, help="θ 方向の分割数（既定 5）")
    ap.add_argument("--nj", type=int, default=4, help="φ 方向の分割数（既定 4）")
    ap.add_argument("--scale", type=float, default=1.2, help="TikZ の scale")
    ap.add_argument("--smooth", type=int, default=3, help="3D の細分数")
    ap.add_argument("--fill", default="collection",
                    choices=["collection", "value", "none"])
    ap.add_argument("--label", default="value",
                    choices=["value", "name", "collection", "none"])
    ap.add_argument("--arrows", action="store_true", help="snc pair を矢印で描く")
    ap.add_argument("--weak", action="store_true", help="weakly critical を丸で囲む")
    ap.add_argument("--caption", default=None)
    ap.add_argument("--three-d", action="store_true", help="3 次元の図を書き出す")
    ap.add_argument("-o", "--out", default="-", help="出力先（既定は標準出力）")
    args = ap.parse_args(argv)

    if args.ni < 3 or args.nj < 3:
        ap.error("ni, nj は 3 以上")
    K = core.torus(args.ni, args.nj)
    f = FUNCTIONS[args.function](K, args.ni, args.nj)

    if args.format in ("tikz", "tikz-standalone"):
        body = tikz(K, f, args.ni, args.nj, scale=args.scale, fill=args.fill,
                    label=args.label, arrows=args.arrows, weak=args.weak,
                    caption=args.caption)
        text = tikz_standalone(body) if args.format == "tikz-standalone" else body
        if args.out == "-":
            sys.stdout.write(text)
        else:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(text)
            print(f"書き出した: {args.out}", file=sys.stderr)
        return 0

    options = (["showCollection"] if args.fill == "collection" else
               ["showColor"] if args.fill == "value" else [])
    if args.arrows:
        options.append("showArrow")
    if args.weak:
        options.append("showWeak")
    try:
        fig2, fig3, _ = plotly_figures(args.function, args.ni, args.nj,
                                       args.smooth, options)
    except ImportError as exc:
        print(f"{args.format} の書き出しには plotly / dash が要る: {exc}", file=sys.stderr)
        return 2
    fig = fig3 if args.three_d else fig2
    out = args.out if args.out != "-" else f"dmb.{args.format}"
    if args.format == "html":
        fig.write_html(out, include_plotlyjs="cdn")
    else:
        fig.write_image(out)
    print(f"書き出した: {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
