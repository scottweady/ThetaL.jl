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


def check(size=(3.4, 2.1)):
    """A small figure for one example's section."""
    fig, ax = plt.subplots(figsize=size)
    return fig, ax


def check_mcf():
    d = read("mcf_area.csv")
    fig, ax = check()
    fs.line(ax, d["t"], d["area"], fs.QUAL[0], label="computed")
    ax.plot(d["t"], d["exact"], "k--", lw=1.0, label=r"$A_0 - 2\pi t$")
    ax.set_xlabel(r"time, $t$"); ax.set_ylabel(r"area, $A$"); ax.legend()
    fig.tight_layout()
    fs.save(fig, "check_mcf")


def check_vortex():
    d = read("vortex_resolution.csv")
    fig, axs = plt.subplots(1, 2, figsize=(fs.WIDTH * 0.85, 2.1))
    fs.line(axs[0], d["t"], d["length"], fs.QUAL[0], marker=False)
    axs[0].set_xlabel(r"time, $t$"); axs[0].set_ylabel(r"length, $L$")
    fs.panel(axs[0], "a")
    axs[1].step(d["t"], d["N"], where="post", color=fs.QUAL[0], lw=fs.LW)
    axs[1].set_yscale("log", base=2)
    axs[1].set_xlabel(r"time, $t$"); axs[1].set_ylabel(r"points, $N$")
    fs.panel(axs[1], "b")
    fig.tight_layout()
    fs.save(fig, "check_vortex")


def relative_area_error(ax, d, k=None, **kw):
    m = np.ones(len(d["t"]), bool) if k is None else d["curve"] == k
    ax.plot(d["t"][m], (d["area"][m] - d["exact"][m]) / d["exact"][m], lw=1.2, **kw)


def check_hele_shaw():
    d = read("hele_shaw_area.csv")
    fig, ax = check()
    relative_area_error(ax, d, color=fs.QUAL[0])
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel(r"time, $t$"); ax.set_ylabel(r"$(A - A_0 - Jt)\,/\,(A_0 + Jt)$")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()
    fs.save(fig, "check_hele_shaw")


def check_bubbles():
    d = read("bubbles_area.csv")
    fig, ax = check()
    for k in np.unique(d["curve"]):
        relative_area_error(ax, d, k, color=fs.QUAL[int(k) - 1], label=f"bubble {int(k)}")
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel(r"time, $t$"); ax.set_ylabel(r"$(A - A_0 - Jt)\,/\,(A_0 + Jt)$"); ax.legend()
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()
    fs.save(fig, "check_bubbles")


def check_vesicles():
    d = read("vesicles_energy.csv")
    fig, ax = check()
    for k in np.unique(d["curve"]):
        m = (d["curve"] == k) & (d["t"] > 0)
        fs.line(ax, d["t"][m], d["energy"][m], fs.QUAL[int(k) - 1], marker=False, label=rf"$\nu = {d['nu'][m][0]:.2f}$")
    ax.set_xscale("log"); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"time, $t$"); ax.set_ylabel(r"bending energy, $E$"); ax.legend(ncol=2)
    fig.tight_layout()
    fs.save(fig, "check_vesicles")


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


def movie(data, name, ghost_every=10, fmt=".1f", width=3.6):
    """A run as an animated GIF, docs/src/assets/<name>-<mode>.gif: every curve at each saved time, filled, and the
    curves of every ghost_every-th frame left behind, faint; the time in the format fmt; width inches wide (at 110
    dpi). A GIF has no partial
    transparency: the page's background."""
    from matplotlib.animation import FuncAnimation, PillowWriter
    d = read(data)
    times = np.unique(d["t"])
    bg = "#1f2424" if fs.DARK else "#ffffff"
    (x0, x1), (y0, y1) = (d["x"].min(), d["x"].max()), (d["y"].min(), d["y"].max())
    pad = 0.04 * max(x1 - x0, y1 - y0)
    w = width
    fig, ax = plt.subplots(figsize=(w, w * (y1 - y0 + 4 * pad) / (x1 - x0 + 2 * pad)), facecolor=bg)
    ax.set_facecolor(bg); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + 3 * pad)    # room above for the time
    fig.subplots_adjust(0, 0, 1, 1)
    label = ax.text(0.02, 0.98, "", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    current = []

    def frame(i):
        t = times[i]
        sel = d["t"] == t
        while current:
            current.pop().remove()
        for cid in np.unique(d["curve"][sel]):
            m = sel & (d["curve"] == cid)
            x, y = d["x"][m], d["y"][m]
            if i % ghost_every == 0:                  # leave the curve behind (the start too)
                ax.plot(x, y, color=fs.GRAY, lw=0.6, alpha=0.6)
            current.extend(ax.fill(x, y, color=fs.QUAL[0], alpha=0.15, lw=0))
            current.extend(ax.plot(x, y, color=fs.QUAL[0], lw=1.2))
        label.set_text(f"$t = {t:{fmt}}$")

    hold = 15                                         # frames to pause on the last time
    anim = FuncAnimation(fig, frame, frames=list(range(len(times))) + [len(times) - 1] * hold, blit=False)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "assets",
                       f"{name}-{fs.DOCS or 'light'}.gif")
    anim.save(out, writer=PillowWriter(fps=12), dpi=110, savefig_kwargs={"facecolor": bg})
    plt.close(fig)
    print(os.path.normpath(out))


def hero():
    """The last frame of the Hele-Shaw movie as a still, for the home page: the bubble at t = 8, filled, and the
    interface at each whole time before it, faint."""
    d = read("hele_shaw_movie.csv")
    times = np.unique(d["t"])
    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    for t in times[(times > 0) & np.isclose(times % 1.0, 0.0)][:-1]:
        m = d["t"] == t
        ax.plot(d["x"][m], d["y"][m], color=fs.GRAY, lw=0.5, alpha=0.6)
    m = d["t"] == times[-1]
    ax.fill(d["x"][m], d["y"][m], color=fs.QUAL[0], alpha=0.15, lw=0)
    ax.plot(d["x"][m], d["y"][m], color=fs.QUAL[0], lw=1.2)
    ax.set_aspect("equal"); ax.set_axis_off()
    fs.save(fig, "hero")


check_mcf()
check_vortex()
check_hele_shaw()
check_bubbles()
check_vesicles()
convergence()
step_control()
movie("mcf_movie.csv", "mean_curvature_flow", ghost_every=8, fmt=".3f")
movie("vortex_movie.csv", "vortex_windup", ghost_every=20)
movie("hele_shaw_movie.csv", "hele_shaw")
movie("bubbles_movie.csv", "hele_shaw_bubbles")
movie("vesicles_movie.csv", "vesicles", ghost_every=1000, fmt=".2g", width=5.4)
hero()
