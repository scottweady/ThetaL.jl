# The figures of ThetaL's documentation, from the CSV files of make_data.jl, in the light and dark renderings of the
# site (docs/src/assets/<name>-light.png, -dark.png):
#   for m in light dark; do FIGSTYLE_DOCS=$m python docs/figures/plot_figures.py; done
# (from the package root). Curves at successive times are colored along the sequential map (time is ordered); schemes and
# variants by the qualitative palette.
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
fs.apply()


def read(name):
    """A CSV of make_data.jl as a dict of columns (numbers where they parse)."""
    with open(os.path.join(DATA, name)) as f:
        head = f.readline().strip().split(",")
        rows = [line.strip().split(",") for line in f if line.strip()]
    cols = {}
    for j, h in enumerate(head):
        col = [r[j] for r in rows]
        try:
            cols[h] = np.array([float(x) for x in col])
        except ValueError:
            cols[h] = np.array(col)
    return cols


def curves_at_times(ax, name, scale_area=False, lw=1.0):
    """Every curve of a gallery CSV, colored by its time (sequential, light = later)."""
    d = read(name)
    times = np.unique(d["t"])
    cols = fs.seq(len(times))
    for t, c in zip(times, cols):
        sel = d["t"] == t
        for cid in np.unique(d["curve"][sel]):
            m = sel & (d["curve"] == cid)
            x, y = d["x"][m], d["y"][m]
            if scale_area:                           # scaled to the area π, about the centroid
                A = 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))
                k = np.sqrt(np.pi / A)
                x, y = (x - x.mean()) * k, (y - y.mean()) * k
            ax.plot(x, y, color=c, lw=lw)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    return times


def gallery():
    fig, axs = plt.subplots(1, 5, figsize=(fs.WIDTH * 1.25, 1.75))
    specs = [("gallery_mcf.csv", "curve shortening", False), ("gallery_sd.csv", "surface diffusion", False),
             ("gallery_aniso.csv", "anisotropic, scaled", True), ("gallery_vortex.csv", "wind-up in a vortex", False),
             ("gallery_coarsening.csv", "coarsening", False)]
    for ax, (name, label, scale) in zip(axs, specs):
        curves_at_times(ax, name, scale_area=scale, lw=0.9)
        ax.set_xlabel(label)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.tight_layout(w_pad=0.5)
    fs.save(fig, "gallery")


def checks():
    fig, axs = plt.subplots(1, 3, figsize=(fs.WIDTH, 2.0))
    d = read("mcf_area.csv")
    fs.line(axs[0], d["t"], d["area"], fs.QUAL[0], label="computed")
    axs[0].plot(d["t"], d["exact"], "k--", lw=1.0, label=r"$A_0 - 2\pi t$")
    axs[0].set_xlabel(r"time, $t$"); axs[0].set_ylabel(r"area, $A$"); axs[0].legend()
    fs.panel(axs[0], "a")
    d = read("aniso_area.csv")
    fs.line(axs[1], d["t"], d["area"], fs.QUAL[0], label="computed")
    axs[1].plot(d["t"], d["exact"], "k--", lw=1.0, label=r"$\pi - 2\pi t$")
    axs[1].set_xlabel(r"time, $t$"); axs[1].set_ylabel(r"area, $A$"); axs[1].legend()
    fs.panel(axs[1], "b")
    d = read("coarsening_radii.csv")
    ids = np.unique(d["curve"])
    order = np.argsort([d["radius"][d["curve"] == i][0] for i in ids])     # colored by initial radius
    cols = fs.seq(len(ids))[::-1]
    for c, i in zip(cols, ids[order]):
        m = d["curve"] == i
        axs[2].plot(d["t"][m], d["radius"][m], color=c, lw=1.0)
    axs[2].set_xlabel(r"time, $t$"); axs[2].set_ylabel(r"radius, $L/2\pi$")
    fs.panel(axs[2], "c")
    fig.tight_layout()
    fs.save(fig, "checks")


def slope(ax, x0, y0, p, label, decades=1.0):
    """A dashed reference line of slope p from (x0, y0) down `decades` in x, labeled at its lower end."""
    x = np.array([x0, x0 * 10 ** (-decades)])
    ax.plot(x, y0 * (x / x0) ** p, "k--", lw=0.8)
    ax.annotate(label, (x[1], y0 * (x[1] / x0) ** p), textcoords="offset points", xytext=(-3, 0), fontsize=8,
                ha="right", va="center")


def below(ax):
    """The legend under the axes, clear of the data."""
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), fontsize=7.5, ncol=1)


def loglog(ax):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())


def convergence():
    fig, axs = plt.subplots(1, 4, figsize=(fs.WIDTH * 1.35, 3.0))
    d = read("conv_rigid.csv")
    for k, name in enumerate(["ETDRK2", "ExpEuler", "ETDRK2 without the rigid split"]):
        m = (d["scheme"] == name) & np.isfinite(d["error"])
        fs.line(axs[0], d["dt"][m], d["error"][m], fs.QUAL[k], label=name.replace("ETDRK2 without the rigid split", "no rigid split"))
    loglog(axs[0])
    axs[0].set_xlabel(r"step, $\Delta t$"); axs[0].set_ylabel("error after a turn"); below(axs[0])
    fs.panel(axs[0], "a")
    d = read("conv_linear.csv")
    for k, name in enumerate(["ETDRK2", "ExpEuler"]):
        m = d["scheme"] == name
        fs.line(axs[1], d["dt"][m], d["error"][m], fs.QUAL[k], label=name)
    slope(axs[1], 0.1, 3e-4, 2, r"$\Delta t^2$", 0.6); slope(axs[1], 0.1, 2.5e-2, 1, r"$\Delta t$", 0.6)
    loglog(axs[1])
    axs[1].set_xlabel(r"step, $\Delta t$"); axs[1].set_ylabel("error, linear flow"); below(axs[1])
    fs.panel(axs[1], "b")
    d = read("conv_sd.csv")
    for k, (name, start) in enumerate([("ETDRK2", "relaxed"), ("ETDRK2", "from the shape"),
                                       ("ExpEuler", "relaxed"), ("ExpEuler", "from the shape")]):
        m = (d["scheme"] == name) & (d["start"] == start)
        style = "-o" if start == "relaxed" else ":o"
        c = fs.QUAL[0 if name == "ETDRK2" else 1]
        axs[2].plot(d["dt"][m], d["error"][m], style, color=c, mfc=c, mec=c, label=f"{name}, {start}")
    loglog(axs[2])
    axs[2].set_xlabel(r"step, $\Delta t$"); axs[2].set_ylabel(r"error in $L$, surface diffusion")
    below(axs[2])
    fs.panel(axs[2], "c")
    d = read("conv_manufactured.csv")
    scales = np.unique(d["symbol_scale"])
    for c, s in zip(fs.seq(len(scales)), scales):
        m = d["symbol_scale"] == s
        fs.line(axs[3], d["dt"][m], d["error"][m], c, label=f"symbol ×{s:g}")
    slope(axs[3], 0.01, 5e-6, 2, r"$\Delta t^2$", 0.6)
    loglog(axs[3])
    axs[3].set_xlabel(r"step, $\Delta t$"); axs[3].set_ylabel("error, manufactured"); below(axs[3])
    fs.panel(axs[3], "d")
    fig.tight_layout()
    fs.save(fig, "convergence")


def step_control():
    d = read("step_control.csv")
    fig, axs = plt.subplots(1, 2, figsize=(fs.WIDTH, 2.0))
    fs.line(axs[0], d["t"], d["dt"], fs.QUAL[0])
    loglog(axs[0])
    axs[0].set_xlabel(r"time, $t$"); axs[0].set_ylabel(r"accepted step, $\Delta t$")
    fs.panel(axs[0], "a")
    fs.line(axs[1], d["t"], d["error_estimate"], fs.QUAL[0])
    axs[1].axhline(1e-4, color="k", ls="--", lw=0.8)
    loglog(axs[1])
    axs[1].set_xlabel(r"time, $t$"); axs[1].set_ylabel("embedded estimate")
    fs.panel(axs[1], "b")
    fig.tight_layout()
    fs.save(fig, "step_control")


gallery()
checks()
convergence()
step_control()
