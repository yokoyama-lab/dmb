#!/usr/bin/env python3
"""トーラス上の離散モースボット関数 (discrete Morse-Bott function) の可視化。

dmf.py（離散モース関数の図示）を離散モースボット理論に書き換えたもの。
理論の計算は dmb_core.py（外部依存なし）にあり，ここは Dash による表示だけを行う。

図示するもの:
    * セルの値 f
    * collection（値が等しく r-path で繋がるセルの同値類）の色分け
    * weakly critical なセル（U^snc = D^snc = 0）＝ reduced collection に入るセル
    * 矢印 = strictly noncritical pair σ ≺ τ, f(τ) < f(σ)
      （離散モース関数のときは Forman の V-path の矢印に一致する）
    * Theorem 4.12 の検算  Σ_C P_t(C) = P_t(K) + (1 + t) R(t)

    python3 dmb.py       # http://127.0.0.1:8050/
"""

import functools
import logging

import numpy as np
import plotly.colors as pc
import plotly.figure_factory as ff
import plotly.graph_objects as go
from dash import Dash, Input, Output, dcc, html

import dmb_core as core
from dmb_core import centroid

marker_size = 10        # line と marker の size は共通，初期値 10
R = 3
r = 1
m = 0.3                 # 3D の表示範囲の余白（大きいとトーラスが小さく映る）
MAX_GRID = 16           # 分割数を増やしたいとき変更
CASING = 'rgba(60,60,60,0.85)'   # 三角形分割を塗りの上でも見えるようにする縁取り
FILL_ALPHA = 0.5                 # 2 次元の塗りの不透明度（濃すぎるとラベルが読めない）
LIFT = 0.03                      # 3 次元で辺・頂点を面から浮かせる量（z-fighting 対策）


def torus_point(theta, phi):
    x = (R + r * np.cos(phi)) * np.cos(theta)
    y = (R + r * np.cos(phi)) * np.sin(theta)
    z = r * np.sin(phi)
    return (x, y, z)


def to_torus(ni, nj, x, y, lift=0.0):
    """基本領域の座標 (x, y) をトーラス上の点に写す。

    lift > 0 のときは曲面の法線方向に持ち上げる。辺や頂点を面と同一平面に置くと
    z-fighting で途切れて見えるため，少しだけ外側に出して描くのに使う。"""
    theta = 2 * np.pi * x / ni
    phi = 2 * np.pi * y / nj
    px, py, pz = torus_point(theta, phi)
    if not lift:
        return (px, py, pz)
    nx = np.cos(phi) * np.cos(theta)
    ny = np.cos(phi) * np.sin(theta)
    nz = np.sin(phi)
    return (px + lift * nx, py + lift * ny, pz + lift * nz)


def get_colorscales(category):
    return [
        name for name in dir(getattr(pc, category, None))
        if isinstance(getattr(getattr(pc, category, None), name), list) and name != '__all__'
    ]


colorscales = (get_colorscales("sequential") + get_colorscales("diverging")
               + get_colorscales("cyclical"))

# collection の色分けに使う離散パレット
QUALITATIVE = (pc.qualitative.Plotly + pc.qualitative.Set2 + pc.qualitative.Dark24)

FUNCTIONS = {
    "height": "Morse-Bott 高さ関数（回転対称）",
    "height_refined": "Morse-Bott 高さ関数・細分版",
    "constant": "定数関数（自明な DMBF）",
    "morsified": "高さ関数の Morsification（DMF）",
    "dmf_min": "臨界セル 4 個の DMF（tree-cotree・非対称）",
    "dmf_invariant": "回転対称な DMF（臨界セル 4·ni 個）",
    "dmf_py": "dmf.py の calcDMF（正方格子のみ）",
}


@functools.lru_cache(maxsize=64)
def analyse(fkey, ni, nj):
    """(関数, 分割数) ごとの計算結果をまとめて返す（表示の切り替えでは再計算しない）。

    ホモロジーの計算は分割数を上げると効くので，チェックボックスやラベルを
    変えるたびに計算し直さないようにキャッシュする。"""
    K = core.torus(ni, nj)
    try:
        f = build_function(fkey, K, ni, nj)
        note = ""
    except Exception as exc:                                   # noqa: BLE001
        f = core.constant_fn(K)
        note = f"（{FUNCTIONS.get(fkey, fkey)} は使えない: {exc}）"
    X = core.MorseBott(K, f)
    return K, f, X, X.report(), core.torus_names(ni, nj), core.lifted_cells(ni, nj), note


def build_function(key, K, ni, nj):
    """関数名から，セル -> 値 の辞書を作る。"""
    if key == "height":
        return core.height_fn(K, ni, nj)
    if key == "height_refined":
        return core.height_fn(K, ni, nj, refine=True)
    if key == "constant":
        return core.constant_fn(K)
    if key == "morsified":
        return core.morsify(K, core.height_fn(K, ni, nj))
    if key == "dmf_min":
        return core.canonical_dmf(K)
    if key == "dmf_invariant":
        return core.invariant_dmf(K, ni, nj)
    if key == "dmf_py":
        return core.dmf_from_dmf_py(K, ni, nj)
    raise ValueError(key)


# ------------------------------------------------------------------ 色の割当


def value_colors(scale, f):
    lo, hi = min(f.values()), max(f.values())
    span = (hi - lo) or 1
    cache = {}

    def color(cellkey):
        v = f[cellkey]
        if v not in cache:
            cache[v] = pc.sample_colorscale(scale, (v - lo) / span)[0]
        return cache[v]

    return color, lo, hi


def collection_colors(X):
    """collection ごとに色を割り当てる（1 セルだけの collection は灰色）。"""
    out = {}
    k = 0
    for L in sorted(X.collections(), key=len, reverse=True):
        if len(L) == 1:
            out[L[0]] = 'rgba(0,0,0,0.25)'
            continue
        c = QUALITATIVE[k % len(QUALITATIVE)]
        k += 1
        for cellkey in L:
            out[cellkey] = c
    return out


def rgb_of(col):
    """色文字列を (r, g, b) にする（塗りの明るさでラベル色を決めるため）。"""
    if col.startswith('#'):
        return pc.hex_to_rgb(col)
    if col.startswith('rgb'):
        parts = col[col.index('(') + 1:col.index(')')].split(',')
        return tuple(float(x) for x in parts[:3])
    return (255, 255, 255)


def with_alpha(col, a=FILL_ALPHA):
    """塗り用に半透明にする。濃い色で黒いラベルや縁取りが潰れるのを防ぐ。"""
    r_, g_, b_ = rgb_of(col)
    return f'rgba({r_:.0f},{g_:.0f},{b_:.0f},{a})'


def nearest_centroid(places, target):
    """target に最も近い placement の重心（貼り合わせをまたぐ矢印のため）。"""
    return min((centroid(p) for p in places),
               key=lambda q: (q[0] - target[0]) ** 2 + (q[1] - target[1]) ** 2)


# --------------------------------------------------------------- Dash アプリ

app = Dash(__name__)
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)


@app.callback(
    [Output("graph_triangulation", "figure"), Output("graph_torus", "figure"),
     Output("graph_colorscale", "figure"), Output("graph_colorscale", "style"),
     Output("report", "children")],
    [Input("dropdown-function", "value"), Input("dropdown", "value"),
     Input("radio-options", "value"), Input("text-input-ni", "value"),
     Input("text-input-nj", "value"), Input("text-input-smooth", "value"),
     Input("check-options", "value")]
)
def update(fkey, scale, radio_value, ni, nj, smooth, check_options):
    ni = max(3, min(MAX_GRID, int(ni or 4)))
    nj = max(3, min(MAX_GRID, int(nj or 4)))
    smooth = max(1, min(MAX_GRID, int(smooth or 2)))
    K, f, X, rep, names, pos, note = analyse(fkey, ni, nj)
    color_by_value, lo, hi = value_colors(scale, f)
    coll_color = collection_colors(X)
    show_color = 'showColor' in check_options
    show_coll = 'showCollection' in check_options
    weak = set(X.weakly_critical())
    reduced_cells = {c for C in X.reduced_collections() for c in C}
    font_size = 10 if radio_value == "cell" else 16

    def cell_color(cellkey, default):
        if show_coll:
            return coll_color[cellkey]
        if show_color:
            return color_by_value(cellkey)
        return default

    collection_index = {}
    for k, L in enumerate(sorted(X.collections(), key=lambda L: (-len(L), L[0]))):
        for cellkey in L:
            collection_index[cellkey] = k

    def label_of(cellkey):
        if radio_value == "cell":
            return names[cellkey]
        if radio_value == "value":
            return str(f[cellkey])
        if radio_value == "collection":
            return str(collection_index[cellkey])
        return None

    # 色ごとにセルをまとめる（1 セル 1 トレースにすると分割数を上げたとき重い）
    groups = {}
    defaults = {0: 'rgba(0,0,0,0.35)', 1: 'rgba(0,0,0,0.25)', 2: 'white'}
    for c in K.cells:
        d = K.dim(c)
        groups.setdefault(d, {}).setdefault(cell_color(c, defaults[d]), []).append(c)

    # ------------------------------------------------ 2 次元（三角形分割）

    fig2 = go.Figure()
    for col, cells in groups[2].items():
        xs, ys = [], []
        for c in cells:
            for p in pos[c]:
                xs += [q[0] for q in p] + [p[0][0], None]
                ys += [q[1] for q in p] + [p[0][1], None]
        fig2.add_trace(go.Scatter(x=xs, y=ys, mode='lines', hoverinfo='skip',
                                  line=dict(width=0), fill='toself',
                                  fillcolor=with_alpha(col)))

    # 辺: 着色していると塗りと同じ色になって三角形分割が見えなくなるので，
    # まず暗い縁取り（casing）を全辺に引き，その上に各辺の色を細く重ねる。
    # 着色していないときは元から見えるので縁取りは引かない（線が重くなるだけ）。
    if show_color or show_coll:
        edge_segments = []
        for cells in groups[1].values():
            for c in cells:
                edge_segments += pos[c]
        xs, ys = [], []
        for p in edge_segments:
            xs += [p[0][0], p[1][0], None]
            ys += [p[0][1], p[1][1], None]
        fig2.add_trace(go.Scatter(x=xs, y=ys, mode='lines', hoverinfo='skip',
                                  line=dict(color=CASING, width=4)))
    for col, cells in groups[1].items():
        xs, ys = [], []
        for c in cells:
            for p in pos[c]:
                xs += [p[0][0], p[1][0], None]
                ys += [p[0][1], p[1][1], None]
        fig2.add_trace(go.Scatter(x=xs, y=ys, mode='lines', hoverinfo='skip',
                                  line=dict(color=col, width=2)))

    for col, cells in groups[0].items():
        pts = [p[0] for c in cells for p in pos[c]]
        fig2.add_trace(go.Scatter(x=[q[0] for q in pts], y=[q[1] for q in pts],
                                  mode='markers', hoverinfo='skip',
                                  marker=dict(color=col, size=marker_size,
                                              line=dict(color=CASING, width=1.5))))

    # weakly critical なセル（reduced collection に入るセル）を丸で囲む
    if 'showWeak' in check_options:
        pts = [centroid(p) for c in K.cells if c in weak for p in pos[c]]
        fig2.add_trace(go.Scatter(
            x=[q[0] for q in pts], y=[q[1] for q in pts], mode='markers', hoverinfo='skip',
            marker=dict(color='rgba(0,0,0,0)', size=marker_size + 8,
                        line=dict(color='crimson', width=2))))

    # 矢印 = strictly noncritical pair σ ≺ τ, f(τ) < f(σ)
    if 'showArrow' in check_options:
        xs, ys, us, vs = [], [], [], []
        for s, t in X.arrows():
            for p1 in [centroid(q) for q in pos[t]]:
                p0 = nearest_centroid(pos[s], p1)
                xs.append(p0[0])
                ys.append(p0[1])
                us.append(p1[0] - p0[0])
                vs.append(p1[1] - p0[1])
        if xs:
            quiver = ff.create_quiver(x=xs, y=ys, u=us, v=vs, scale=0.9,
                                      arrow_scale=0.3, line_color='rgba(0,0,0,0.65)')
            for trace in quiver.data:
                trace.hoverinfo = 'skip'
                fig2.add_trace(trace)

    if radio_value != "none":
        labels, lx, ly = [], [], []
        for c in K.cells:
            text = label_of(c)
            if text is None:
                continue
            p = centroid(pos[c][0])          # 標準の持ち上げにだけラベルを置く
            labels.append(text)
            lx.append(p[0])
            ly.append(p[1])
        fig2.add_trace(go.Scatter(x=lx, y=ly, mode='text', hoverinfo='skip', text=labels,
                                  textposition='middle center',
                                  textfont=dict(color='black', size=font_size)))

    size = max(ni, nj) * 120
    pad = 0.55
    fig2.update_layout(
        xaxis=dict(showticklabels=False, range=[-pad, ni + pad]),
        yaxis=dict(showticklabels=False, scaleanchor="x", range=[-pad, nj + pad]),
        showlegend=False,
        plot_bgcolor='white',
        height=size,
        width=size,
        margin=dict(l=0, r=0, t=0, b=0),
    )

    # ------------------------------------------------------- 3 次元（トーラス）
    # 平行移動した複製は 3 次元では同じ点に写るので，標準の持ち上げだけを使う。

    fig3 = go.Figure()
    for col, cells in groups[2].items():
        xs, ys, zs, ii, jj, kk = [], [], [], [], [], []
        for c in cells:
            a, b, cc = pos[c][0]
            index = {}
            base = len(xs)
            for p in range(smooth + 1):
                for q in range(smooth + 1 - p):
                    u, v = p / smooth, q / smooth
                    x = a[0] + u * (b[0] - a[0]) + v * (cc[0] - a[0])
                    y = a[1] + u * (b[1] - a[1]) + v * (cc[1] - a[1])
                    px, py, pz = to_torus(ni, nj, x, y)
                    index[(p, q)] = base + len(index)
                    xs.append(px)
                    ys.append(py)
                    zs.append(pz)
            for p in range(smooth):
                for q in range(smooth - p):
                    ii.append(index[(p, q)])
                    jj.append(index[(p + 1, q)])
                    kk.append(index[(p, q + 1)])
                    if p + q <= smooth - 2:
                        ii.append(index[(p + 1, q)])
                        jj.append(index[(p + 1, q + 1)])
                        kk.append(index[(p, q + 1)])
        fig3.add_trace(go.Mesh3d(x=xs, y=ys, z=zs, i=ii, j=jj, k=kk,
                                 color=col, opacity=1, hoverinfo='skip'))

    edge3d = dict(groups[1])
    casing3d = ([(CASING, [c for cs in edge3d.values() for c in cs])]
                if (show_color or show_coll) else [])
    for col, cells in casing3d + list(edge3d.items()):
        width = 7 if col is CASING else 4
        xs, ys, zs = [], [], []
        for c in cells:
            p0, p1 = pos[c][0]
            for s in range(smooth + 1):
                u = s / smooth
                px, py, pz = to_torus(ni, nj, p0[0] + u * (p1[0] - p0[0]),
                                      p0[1] + u * (p1[1] - p0[1]), LIFT)
                xs.append(px)
                ys.append(py)
                zs.append(pz)
            xs.append(None)
            ys.append(None)
            zs.append(None)
        fig3.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                                    line=dict(color=col, width=width),
                                    showlegend=False, hoverinfo='skip'))

    for col, cells in groups[0].items():
        coords = [to_torus(ni, nj, *pos[c][0][0], lift=LIFT * 1.5) for c in cells]
        fig3.add_trace(go.Scatter3d(
            x=[p[0] for p in coords], y=[p[1] for p in coords], z=[p[2] for p in coords],
            mode='markers', marker=dict(color=col, size=3), hoverinfo='skip'))

    if radio_value != "none":
        coords, texts = [], []
        for c in K.cells_of_dim(0):
            text = label_of(c)
            if text is None:
                continue
            coords.append(to_torus(ni, nj, *pos[c][0][0], lift=LIFT * 3))
            texts.append(text)
        fig3.add_trace(go.Scatter3d(
            x=[p[0] for p in coords], y=[p[1] for p in coords], z=[p[2] for p in coords],
            mode='text', text=texts, textposition='top center',
            textfont=dict(color="black", size=font_size), hoverinfo='skip'))

    fig3.update_layout(
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=(r + m) / (R + r + m)),
            xaxis=dict(showticklabels=False, visible=False, showgrid=False,
                       range=[-R - r - m, R + r + m]),
            yaxis=dict(showticklabels=False, visible=False, showgrid=False,
                       range=[-R - r - m, R + r + m]),
            zaxis=dict(showticklabels=False, visible=False, showgrid=False,
                       range=[-r - m, r + m]),
        ),
        showlegend=False,
        plot_bgcolor='white',
        height=size,
        width=size,
        margin=dict(l=0, r=0, t=0, b=0),
    )

    # ------------------------------------------------------------ カラースケール

    steps = 128
    fig_scale = go.Figure(go.Heatmap(
        z=[list(range(steps))],
        x=[lo + (hi - lo) * s / (steps - 1) for s in range(steps)],
        colorscale=scale, showscale=False, hoverinfo='skip'))
    fig_scale.update_layout(
        xaxis=dict(linewidth=1, linecolor='white'),
        yaxis=dict(showticklabels=False),
        showlegend=False, plot_bgcolor='white',
        height=40, width=size, margin=dict(l=0, r=0, t=0, b=0))
    style = {'display': 'block'} if (show_color and not show_coll) else {'display': 'none'}

    # ------------------------------------------------------------------ 報告

    lines = [f"{FUNCTIONS.get(fkey, fkey)}   T({ni}, {nj})   {note}",
             f"セル数 {K.counts()}   値域 {lo}..{hi}   "
             f"Z_{ni} 回転で不変: {'はい' if core.is_invariant(ni, nj, f) else 'いいえ'}"]
    if not rep["is_dmb"]:
        kind, s, w = rep["dmb_violations"][0]
        lines.append(f"離散モースボット関数ではない: 違反 {len(rep['dmb_violations'])} 件, "
                     f"最初は {kind} at {names[s]}（witness {len(w)} 個）")
    else:
        lines.append("(M1)(MB2)(M3)(MB4): OK  [(M1)(M3) は単体的複体なので vacuous]")
        lines.append(f"離散モース関数か: {'はい' if rep['is_dmf'] else 'いいえ'}"
                     f"（(M2)/(M4) の違反 {len(rep['dmf_violations'])} 件）")
        lines.append(f"collection {rep['collections']} 個 / reduced collection "
                     f"{len(rep['reduced_collections'])} 個 / weakly critical "
                     f"{rep['n_weakly_critical']} セル（うち reduced collection に入るもの "
                     f"{len(reduced_cells)}）/ critical {rep['n_critical']} セル / "
                     f"矢印 {rep['n_arrows']} 本")
        kinds = {}
        for C in rep["reduced_collections"]:
            b = core.poly_trim(core.betti(C))
            if not b:
                continue
            key = (len(C), tuple(sorted({K.dim(c) for c in C})), tuple(b))
            kinds[key] = kinds.get(key, 0) + 1
        for (sz, dims, b), mult in sorted(kinds.items(), reverse=True):
            lines.append(f"    C: {sz} セル（次元 {list(dims)}）  P_t(C) = {core.poly_str(b)}"
                         + (f"   × {mult} 個" if mult > 1 else ""))
        lines.append(f"P_t(K)     = {core.poly_str(rep['P_K'])}")
        lines.append(f"Σ_C P_t(C) = {core.poly_str(rep['MB_sum'])}")
        if rep["R_MB"] is None:
            lines.append("!! Σ_C P_t(C) - P_t(K) が (1+t) で割り切れない")
        else:
            lines.append(f"R(t)       = {core.poly_str(rep['R_MB'])}"
                         + ("   ← 等号（鋭い）" if rep["MB_sharp"] else ""))
        if rep["is_dmf"]:
            lines.append(
                f"（離散モース理論）M(t) = {core.poly_str(rep['M'])},  R_M(t) = "
                + (core.poly_str(rep['R_M']) if rep['R_M'] is not None else "割り切れない"))
    return fig2, fig3, fig_scale, style, "\n".join(lines)


def number_input(id_, value, minimum=3, maximum=MAX_GRID):
    return dcc.Input(id=id_, type='number', min=minimum, max=maximum, step=1, value=value,
                     style={'height': '30px', 'width': '45px', 'display': 'inline-block',
                            'vertical-align': 'middle'})


inline = {'margin-left': '10px', 'display': 'inline-block', 'vertical-align': 'middle'}

app.layout = html.Div([
    html.Div([
        html.Div([
            html.P("Function:", style=inline),
            dcc.Dropdown(id='dropdown-function',
                         options=[{'label': v, 'value': k} for k, v in FUNCTIONS.items()],
                         value='height',
                         style={'height': '30px', 'width': '360px',
                                'display': 'inline-block', 'vertical-align': 'middle'}),
            html.P("Grid Size (i, θ 方向):", style=inline),
            number_input('text-input-ni', 4),
            html.P("Grid Size (j, φ 方向):", style=inline),
            number_input('text-input-nj', 4),
            html.P("Smooth:", style=inline),
            number_input('text-input-smooth', 2, minimum=1),
            html.P("Color Scale:", style=inline),
            dcc.Dropdown(id='dropdown', options=colorscales, value='YlGn',
                         style={'height': '30px', 'width': '200px',
                                'display': 'inline-block', 'vertical-align': 'middle'}),
            html.Div([
                dcc.RadioItems(
                    id='radio-options',
                    options=[
                        {'label': 'Cell Name', 'value': 'cell'},
                        {'label': 'Value f', 'value': 'value'},
                        {'label': 'Collection', 'value': 'collection'},
                        {'label': 'No Label', 'value': 'none'},
                    ],
                    value='value', inline=True, style={'display': 'inline-block'}),
            ]),
            html.Div([
                dcc.Checklist(
                    id='check-options',
                    options=[
                        {'label': 'Show Color (値)', 'value': 'showColor'},
                        {'label': 'Show Collection (色分け)', 'value': 'showCollection'},
                        {'label': 'Show Arrow (snc pair)', 'value': 'showArrow'},
                        {'label': 'Show Weakly Critical', 'value': 'showWeak'},
                    ],
                    value=['showCollection'], inline=True, style={'display': 'inline-block'}),
            ]),
        ])
    ], style={'display': 'flex', 'align-items': 'center'}),
    html.Div([
        html.Div([
            html.P('figure : Triangulation of Torus',
                   style={'margin': '0px', 'textAlign': 'center'}),
            dcc.Graph(id='graph_triangulation'),
        ]),
        html.Div([
            html.P('figure : Torus', style={'margin': '0px', 'textAlign': 'center'}),
            dcc.Graph(id='graph_torus'),
        ]),
    ], style={'display': 'flex', 'flex-direction': 'row'}),
    dcc.Graph(id='graph_colorscale', style={'display': 'none'}),
    html.Pre(id='report', style={'fontFamily': 'monospace', 'fontSize': '13px',
                                 'background': '#f6f6f6', 'padding': '10px',
                                 'whiteSpace': 'pre-wrap'}),
])


if __name__ == "__main__":
    run = getattr(app, "run", None) or app.run_server
    run(port=8050)   # port 8050 を使用できない場合は変更
