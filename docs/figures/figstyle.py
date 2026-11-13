# Matplotlib style of ThetaL's documentation figures: serif text and math, boxed axes with inward ticks,
# no grid, no titles, axis labels "description, symbol", italic panel labels (a), (b), lines
# with filled circle markers, black dashed lines for references, legends without a frame. Colours by the data's job
# (QUAL qualitative, SEQ sequential, DIV diverging).
#   import sys; sys.path.insert(0, <this folder>); import figstyle as fs; fs.apply()
# Docs mode (FIGSTYLE_DOCS=light or dark): a transparent background, black ink for the light theme and light grey
# for the dark one; `save` then writes docs/src/assets/<name>-light.png or <name>-dark.png (the site's
# CSS shows the one matching the theme).
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DOCS = os.environ.get("FIGSTYLE_DOCS", "")       # "", "light" or "dark"
INK = "#d4d4d4" if DOCS == "dark" else "black"
GREY = "0.5" if DOCS == "dark" else "0.55"
# Colour by the data's job (one rule each); every role has its own steps for the dark docs theme:
#   qualitative (identity: cases, methods, mesh classes): QUAL, in this fixed order, never cycled — jewel tones
#     indigo, rose, teal, amber, violet, lime, sky, pink (checked: light on
#     white, worst adjacent CVD ΔE 10.2, normal-vision 24.3, all ≥ 3:1 contrast; dark on #1f2424, CVD 8.3,
#     normal-vision 21.8, all ≥ 3:1; the first three slots separate in every pairing in both modes);
#   sequential (an ordered parameter, a magnitude): SEQ, one hue (indigo) with monotone lightness, from the
#     strongest contrast with the background to the weakest; seq(n) samples it;
#   diverging (a signed quantity about a reference): DIV, indigo → neutral grey → rose;
#   fields of a magnitude with white arrows: FIELD, dark → light, monotone lightness.
DARK = DOCS == "dark"
QUAL = (["#6875f5", "#f05470", "#14a394", "#c97a0a", "#9b6cf6", "#4d9e1a", "#2493d1", "#e7559c"] if DARK else
        ["#4f46e5", "#e11d48", "#0d9488", "#d97706", "#7c3aed", "#65a30d", "#0284c7", "#db2777"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#a5b4fc", "#818cf8", "#6366f1", "#4f46e5", "#4338ca"] if DARK else
                                        ["#26215f", "#3730a3", "#4f46e5", "#6366f1", "#8b97f9"])
DIV = LinearSegmentedColormap.from_list("div", ["#312e81", "#4f46e5", "#a5b4fc", "#4b5050" if DARK else "#e5e7eb",
                                                "#fda4af", "#e11d48", "#881337"])
FIELD = LinearSegmentedColormap.from_list("field", ["#0b1026", "#1e1b4b", "#3730a3", "#4f74e0", "#5fb8e8", "#e0f2fe"])
ACCENT = QUAL[0]     # a single series
MESH = {"interior": QUAL[0], "graded": QUAL[2], "weighted": QUAL[1]}

# Application themes (fs.apply("microbial") etc.): the docs page's hue for every role. QUAL becomes four tones of
# that hue, stepped in lightness with small hue shifts (worst CVD ΔE 6.7–11.6, normal-vision 9.5–16, so at
# most four cases and always with a legend); SEQ, FIELD and ACCENT follow the same hue. Light/dark as above.
THEMES = {
    "membranes": dict(qual=["#0361bd", "#22affc", "#012681", "#7870e4"], qual_dark=["#62a7fd", "#b8e1fe", "#3c6ed7", "#b8b9fe"],
                      seq=["#012854", "#054180", "#0961bb", "#3280dd", "#57a0fb"],
                      seq_dark=["#d4e6fe", "#9dc7fe", "#60a7ff", "#3986e4", "#1d6ec9"],
                      field=["#000d22", "#012854", "#034b95", "#3280dd", "#85baff", "#e9f3ff"]),
    "rigid": dict(qual=["#b7162d", "#ff6b91", "#640501", "#da4273"], qual_dark=["#fc7374", "#fecbd4", "#cc3526", "#fea2b7"],
                  seq=["#52030e", "#80041a", "#af2934", "#d24c50", "#f66c6d"],
                  seq_dark=["#ffdbd8", "#ffaaa7", "#fd7273", "#d95256", "#bd373f"],
                  field=["#220002", "#52020d", "#910a20", "#d24c50", "#fd9693", "#ffedec"]),
    "chiral": dict(qual=["#a44c05", "#e99606", "#502603", "#c5684b"], qual_dark=["#fa7c20", "#fed6a9", "#af5b0a", "#fead91"],
                   seq=["#441b00", "#6a2f03", "#9c4804", "#c65f0f", "#eb7c33"],
                   seq_dark=["#ffddca", "#feaf82", "#f2823b", "#cf630d", "#ad5003"],
                   field=["#1b0700", "#441c02", "#7b3700", "#c65f0f", "#f99e6a", "#ffeee5"]),
    "microbial": dict(qual=["#01791e", "#80b934", "#013e1e", "#15997d"], qual_dark=["#54bf5c", "#b5f171", "#138b4b", "#23e0b8"],
                      seq=["#033309", "#035112", "#107823", "#3a9742", "#5bb661"],
                      seq_dark=["#a5fea8", "#81dd85", "#61bd67", "#409d48", "#248430"],
                      field=["#001201", "#023309", "#085e18", "#3a9742", "#85cc87", "#e3f9e3"]),
    "langmuir": dict(qual=["#7541b8", "#9c95fc", "#49046f", "#ba54b9"], qual_dark=["#b589fc", "#d6d6ff", "#924fc4", "#fe93fb"],
                     seq=["#370164", "#522288", "#7347af", "#9065d0", "#ae84f2"],
                     seq_dark=["#e9deff", "#ceb5fd", "#b58bf9", "#966cd7", "#7e53bc"],
                     field=["#140129", "#331556", "#5c2f95", "#9065d0", "#c4a4fd", "#f4efff"]),
}


_BASE = (QUAL, SEQ, FIELD)


def theme(name):
    """Switch QUAL, SEQ, FIELD, ACCENT and MESH to an application theme (see THEMES); None restores the default."""
    global QUAL, SEQ, FIELD, ACCENT, MESH
    if name is None:
        QUAL, SEQ, FIELD = _BASE
        ACCENT = QUAL[0]
        MESH = {"interior": QUAL[0], "graded": QUAL[2], "weighted": QUAL[1]}
        return
    t = THEMES[name]
    QUAL = t["qual_dark"] if DARK else t["qual"]
    SEQ = LinearSegmentedColormap.from_list("seq_" + name, t["seq_dark"] if DARK else t["seq"])
    FIELD = LinearSegmentedColormap.from_list("field_" + name, t["field"])
    ACCENT = QUAL[0]
    MESH = {"interior": QUAL[0], "graded": QUAL[2], "weighted": QUAL[1]}


LW = 1.4             # data lines
MS = 3.5             # markers
WIDTH = 6.5          # text width of doc/numerics.tex (in): draw at print size, include at \textwidth


def apply(name=None):
    """Set the rcParams; with a name, also the application theme (see THEMES)."""
    if name:
        theme(name)
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.labelsize": 10, "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.linewidth": 1.0, "axes.edgecolor": INK, "axes.labelcolor": INK, "axes.grid": False,
        "axes.spines.top": True, "axes.spines.right": True,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "lines.linewidth": LW, "lines.markersize": MS,
        "legend.frameon": False, "legend.handlelength": 2.0, "legend.labelspacing": 0.3,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "figure.dpi": 100,
    })


def seq(n):
    """n colours along the sequential map, dark (first) to light (last)."""
    return [SEQ(x) for x in (np.linspace(0, 1, n) if n > 1 else [0.5])]


def panel(ax, label):
    """Italic panel label (a), (b), … above the top-left corner of the axes, clear of the top tick label."""
    ax.annotate(f"({label})", xy=(0, 1), xycoords="axes fraction", xytext=(-3, 7), textcoords="offset points",
                ha="right", va="bottom", fontstyle="italic", fontsize=10)


def line(ax, x, y, color, marker=True, **kw):
    """Data line with filled circle markers of the same colour."""
    return ax.plot(x, y, "-o" if marker else "-", color=color, mfc=color, mec=color, **kw)


def save(fig, path):
    """Save as path.pdf and path.png (200 dpi); in docs mode as a transparent docs/src/assets/<name>-<mode>.png."""
    if DOCS:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "assets",
                           f"{os.path.basename(path)}-{DOCS}.png")
        fig.savefig(out, dpi=200, transparent=True)
        print(os.path.normpath(out))
        return
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=200)
    print(path + ".pdf")


def sparse(ax, x, y, n=18):
    """Thin a digitized curve to about n points evenly spaced along its arclength in axes units: the points of each
    arclength bin are averaged (which also smooths pixel noise). Call after the axis limits are set."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    ok = (x >= min(x0, x1)) & (x <= max(x0, x1)) & (y >= min(y0, y1)) & (y <= max(y0, y1))
    x, y = x[ok], y[ok]
    ds = np.hypot(np.diff(x) / (x1 - x0), np.diff(y) / (y1 - y0))
    cut = np.nonzero(ds > 0.03)[0] + 1            # gaps (columns dropped by the digitizer): no bin spans one
    ds[cut - 1] = 0
    s = np.r_[0, np.cumsum(ds)]
    k = np.minimum((s / s[-1] * n).astype(int), n - 1)
    seg = np.searchsorted(cut, np.arange(len(x)), side="right")
    big = np.bincount(seg)[seg] >= 3                 # isolated columns are overlap artefacts
    sel = []
    for i in range(n):                               # one point per bin, from its most populated segment
        m = (k == i) & big
        if m.any():
            sel.append(m & (seg == np.bincount(seg[m]).argmax()))
    return np.array([x[m].mean() for m in sel]), np.array([y[m].mean() for m in sel])


def ref_markers(ax, x, y, color, n=18, **kw):
    """Results of another reference (digitized or tabulated, not analytical): open circles, thinned by sparse()."""
    xs, ys = sparse(ax, x, y, n)
    kw = dict(dict(ms=3.5, mew=0.9, ls="none", marker="o", mfc="none", zorder=5), **kw)
    return ax.plot(xs, ys, color=color, mec=color, **kw)
