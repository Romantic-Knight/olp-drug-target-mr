"""build_restructured_sc_figures.py — Plan 3 restructure of the single-cell figures (Python side).

Companion to build_restructured_figures.R, which handles the R/ggplot plates. The single-cell
plates are matplotlib and the source scripts are procedural (no reusable panel functions), so the
panels are recovered by running each source script with figure writing neutralised and then
**re-parenting the real Axes objects** into a new composite figure. The drawn content (points,
colours, ticks and any `ax.legend(...)` anchored to the axes) is thus the original panel — nothing
is cropped and no panel is re-drawn from scratch.

Writes ONLY to figures/main_figures_v2/ and figures/supplementary_v2/; no source figure, script or
data file is modified.

STATUS (2026-09-23): the compose approach below was **tested and rejected**. Running the source
scripts with writes redirected into a temp directory works (that is what makes their own
text-collision QA pass), but matplotlib cannot re-parent an Axes across figures: `add_axes(ax)`
refuses it, and bypassing the check by repointing `ax.figure`/`_axstack` renders every panel at its
ORIGINAL position, producing an unreadable pile of overlapping panels. See
results/reports/figure_restructure_plan3_mapping.md ("The single-cell plates ... need a decision").
The reusable part is `run()` + the savefig redirect; the `adopt()`/`write()` composition is dead.

Run in the project's scRNA environment (analysis/scrna/pixi.toml):
    cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py
"""
import importlib
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.figure import Figure       # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, "analysis", "scripts"))

MAINDIR = os.path.join(ROOT, "figures", "main_figures_v2")
SUPPDIR = os.path.join(ROOT, "figures", "supplementary_v2")
for d in (MAINDIR, SUPPDIR):
    os.makedirs(d, exist_ok=True)

# --- redirect figure writing (used ONLY while composing, never when verifying) ---------------
# The redirect must NOT be active when a refactored module's own main() is run for verification:
# otherwise its savefig calls are silently diverted to the scratch directory and the byte-compare
# passes without anything having been written (a false positive — this actually happened).
import tempfile                        # noqa: E402

_SCRATCH = tempfile.mkdtemp(prefix="scfig_scratch_")
_real_savefig = Figure.savefig
_PATCHED = {"on": False}


def enable_redirect():
    """Divert source-script figure writes to a scratch dir. Compose-only."""
    if _PATCHED["on"]:
        return

    def _redirected(self, fname, *a, **k):
        return _real_savefig(self, os.path.join(_SCRATCH, os.path.basename(str(fname))), *a, **k)

    Figure.savefig = _redirected
    plt.close = lambda *a, **k: None
    plt.show = lambda *a, **k: None
    _PATCHED["on"] = True


def disable_redirect():
    Figure.savefig = _real_savefig
    _PATCHED["on"] = False

_RUN = {}


def run(script):
    """Execute a source script; its figures are built but never written."""
    if script in _RUN:
        return _RUN[script]
    path = os.path.join(ROOT, "analysis", "scripts", script)
    ns = {"__name__": "__main__", "__file__": path}
    print(f"  running {script} ...", flush=True)
    exec(compile(open(path, encoding="utf-8").read(), path, "exec"), ns)
    _RUN[script] = ns
    return ns


def adopt(*a, **k):
    """REMOVED — matplotlib cannot re-parent an Axes across figures (see module docstring)."""
    raise NotImplementedError(
        "Axes cannot be moved between figures; refactor the source script into panel "
        "functions instead (route C in results/reports/figure_restructure_plan3_mapping.md)")


def new_fig(w_mm, h_mm):
    """Create a figure already at its final size.

    This matters: any geometry expressed in figure fractions but derived from text extents (panel
    letters are placed relative to their title's rendered top, see letters_grid) depends on the
    canvas size, because a fixed point size occupies a different fraction of a 4.8 in default canvas
    than of the final 9.8 in plate. Building first and resizing in write() therefore produced letter
    positions computed for the wrong canvas.
    """
    return plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))


def write(fig, stem, w_mm, h_mm):
    fig.set_size_inches(w_mm / 25.4, h_mm / 25.4)
    _real_savefig(fig, f"{stem}.png", dpi=300)
    _real_savefig(fig, f"{stem}.pdf")
    print(f"  saved {os.path.basename(stem)} ({w_mm:.0f} x {h_mm:.0f} mm)")


# def build_xxx() functions are added below as each plate is built.


def build_main4():
    """Main 4 - lesion atlas, cluster annotation and composition (6 panels).

    Panel A's and C's legends are anchored OUTSIDE their axes (bbox x > 1), so the grid keeps a
    dedicated middle column to hold them; placing B or D there would bury them.
    """
    print("== Main 4 (single-cell atlas, annotation and composition) ==")
    ero = importlib.import_module("scrna_erosion_gse211630")            # Fig6 A/B
    um = importlib.import_module("scrna_suppfig5_umap_gse211630")       # S6 A/B
    comp = importlib.import_module("scrna_plot_composition_tnk_gse211630")   # S9A
    m7 = importlib.import_module("scrna_suppfig7_marker_dotplot_gse211630")  # S7

    fig = new_fig(180, 250)
    # hspace raised from 0.26 and the plate made taller: the A/B to C/D gap and the E to F gap were
    # reported as too tight (the letters and C's legend need that band).
    outer = fig.add_gridspec(4, 3, width_ratios=[1.0, 1.25, 1.0],
                             height_ratios=[1.50, 1.50, 0.78, 2.22],
                             left=0.032, right=0.975, top=0.972, bottom=0.018,
                             wspace=0.10, hspace=0.40)
    axA = fig.add_subplot(outer[0, 0])
    axB = fig.add_subplot(outer[0, 2])
    axC = fig.add_subplot(outer[1, 0])
    axD = fig.add_subplot(outer[1, 2])
    axE = fig.add_subplot(outer[2, :])
    gsF = outer[3, :].subgridspec(2, 1, height_ratios=[2.10, 0.42], hspace=0.34)
    axF = fig.add_subplot(gsF[0])
    axFc = fig.add_subplot(gsF[1])

    # Letters are placed by this builder (internal tags suppressed) so every plate shares one
    # hierarchy: letter > title > uniform body text. restyle() runs before letter_at().
    ero.panel_A(axA, tag=None)
    # Both big keys are single-column, small-font and tightly spaced so they fit the middle
    # column: a single 28-row column needs ~57 mm, and the original spacing would need ~99 mm.
    # They are legends, so they take the dense body size, not the panel's own body size.
    ero.panel_A_legend(axA, ncol=1, fontsize=BODY_FS_DENSE, labelspacing=0.05,
                       ms=3.4, ms_dbl=3.0)
    # B's key sits BELOW its axes as a two-column strip, matching D: the map's corners all carry
    # points at this panel size, and the band under row 1 is empty in the right column.
    ero.panel_B(axB, tag=None, ms=3.4, legend_loc="upper center",
                legend_bbox=(0.5, -0.03), legend_ncol=2)
    um.panel_a(axC, tag=False, ncol=1, labelspacing=0.05,
               title="Raw Leiden clusters (c0\u2013c20)",
               ms=3.4, ms_dbl=3.0)
    # D's key moves from the source panel's lower-left into the map's empty upper-right corner: at
    # this panel size the lower-left position sits over the embedding's structure.
    # D's embedding fills its panel - every corner has points - so no in-axes position can hold the
    # key without covering data. It is parked in the empty band UNDER the axes instead, as a
    # two-column strip, which also keeps it off the map entirely.
    um.panel_b(axD, tag=False, title="Per-sample coverage",
               legend_loc="upper center", legend_bbox=(0.5, -0.035),
               legend_ncol=2, legend_fontsize=4.8,
               ms=2.8, labelspacing=0.34, borderpad=0.25)
    comp.panel_composition_bars(axE, letter=None, cond_y=1.30)
    m7.panel_dotplot(axF)
    m7.panel_coverage(axFc)

    # A/B and C/D are rows: the two panels in each row must sit on the same horizontal line.
    align_row(fig, (axA, axB))
    align_row(fig, (axC, axD))

    # A and C carry the long legends; their own axes take the dense body size.
    restyle(axA, body=BODY_FS_DENSE)
    restyle(axB)
    restyle(axC, body=BODY_FS_DENSE)
    restyle(axD)
    restyle(axE)
    restyle(axF, body=BODY_FS_DENSE)     # 51 marker labels
    restyle(axFc, body=BODY_FS_DENSE)

    # E's own n= labels sit just above its axes, where the panel title also lands; give the axes
    # headroom so the two do not touch. Its y-axis label is parked outside the axes (at -0.18 axes
    # widths in the source) and needs both a left margin and a nearer anchor to stay on the canvas.
    axE.set_ylim(0, 1.42)
    p = axE.get_position()
    axE.set_position([p.x0 + 0.055, p.y0, p.width - 0.055, p.height])
    axE.yaxis.set_label_coords(-0.085, 0.5)

    # Letters sit level with their title's own line, one column per plate column. Several of these
    # titles are wide enough, when centred, to reach the letter's position; each such title is
    # nudged right so the letter gets a clear gutter. `title_after` moves only a title that would
    # actually collide. NOTE: this must run AFTER restyle(), which re-emits the title.
    xcol_left = 0.014
    xcol_right = axB.get_position().x0 - 0.032
    # Four of the six titles were wide enough, when centred, to reach the letter's position; they
    # are shortened here so the letter has a clear gutter. Re-emitting a shifted title was tried and
    # rendered the text twice, so the titles themselves are trimmed instead. Titles only; no data.
    def retitle(ax, txt):
        for t in (ax._left_title, ax._right_title):
            t.set_text("")
        ax.set_title(txt, fontsize=TITLE_FS, loc="center", fontweight="normal")

    retitle(axA, "OLP atlas (n = 46,478)")
    retitle(axB, "NEOLP vs erosive EOLP")
    retitle(axC, "Raw Leiden clusters (c0\u2013c20)")
    retitle(axE, "Cell-type composition (families)")
    retitle(axF, "Lineage markers used for cluster annotation")
    letters_grid(fig, None, [
        ((axA, "A", xcol_left), (axB, "B", xcol_right)),
        ((axC, "C", xcol_left), (axD, "D", xcol_right)),
    ])
    letter_at(fig, axE, "E", x_fig=xcol_left)
    letter_at(fig, axF, "F", x_fig=xcol_left)


    write(fig, os.path.join(MAINDIR, "fig05_atlas_composition"), 180, 250)


def title_left(fig, ax):
    """Figure-fraction x where the axes' title text begins (0 when there is no title)."""
    t = ax.title
    if not t.get_text().strip():
        return ax.get_position().x0
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    tb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    return inv.transform((tb.x0, tb.y0))[0]


def letters_grid(fig, entries, row_pairs):
    """Place panel letters with a per-row shared height and per-column shared x.

    A letter sits just above its panel's title, and titles differ in height (a two-line title pushes
    its letter higher than a one-line one), so panels side by side end up with letters at different
    heights unless the row's height is computed once and shared. Same for the column: the letters of
    B and D were placed left of their own titles, and different title lengths put them on different
    x. `row_pairs` groups the entries that must line up (each a row of the plate).
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()

    def line_centre(ax):
        """Vertical centre of the title's LAST line, in figure fraction.

        Axes titles are bottom-anchored, so the last line of a two-line title sits at the same
        height as the only line of a one-line title in the panel beside it. Aligning the letter with
        that line therefore puts A and B (or C and D) at the same height even when their titles have
        different line counts, and it places the letter level with the title rather than floating
        above it.
        """
        t = ax.title
        if not t.get_text().strip():
            return ax.get_position().y1, ax.get_position().x1
        tb = t.get_window_extent(renderer=rend)
        x0, y0 = inv.transform((tb.x0, tb.y0))
        x1, y1 = inv.transform((tb.x1, tb.y1))
        n = max(1, t.get_text().count("\n") + 1)
        return y0 + (y1 - y0) / (2.0 * n), x1

    for pair in row_pairs:
        for ax, letter, x in pair:
            y, tx1 = line_centre(ax)
            fig.text(x, y, letter.upper(), fontsize=LETTER_FS, fontweight="bold",
                     ha="left", va="center", zorder=20)


def letter_at(fig, ax, letter, x_fig=None):
    """Place a panel letter at the panel's top-left corner, just above its title.

    x_fig overrides the horizontal position (figure fraction) so several panels can share one
    column of letters.
    """
    fig.canvas.draw()
    bb = ax.get_position()
    # The title may be on any of the three artists; use whichever has text.
    t = ax.title if ax.title.get_text().strip() else (
        ax._left_title if ax._left_title is not None and ax._left_title.get_text().strip()
        else ax._right_title)
    inv = fig.transFigure.inverted()
    if t is not None and t.get_text().strip():
        tb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        _, ty0 = inv.transform((tb.x0, tb.y0))
        _, ty1 = inv.transform((tb.x1, tb.y1))
        n = max(1, t.get_text().count("\n") + 1)
        # Level with the title's LAST line. Centred titles have different line counts side by side,
        # and axes titles are bottom-anchored, so the last lines line up across a row.
        y = ty0 + (ty1 - ty0) / (2.0 * n)
    else:
        y = bb.y1 + 0.006
    # A centred title can reach the letter's column. Where it does - typically a long title whose
    # first line starts further left than the axes - shift the whole block right so the letter gets
    # a clear gutter. Shifting before reading the extent above would have measured the old position,
    # so the measured y stays valid only because the shift is horizontal.
    if t.get_text().strip():
        tx = t.get_window_extent(renderer=fig.canvas.get_renderer())
        x0f = inv.transform((tx.x0, tx.y0))[0]
        need = (bb.x0 if x_fig is None else x_fig) + 0.028
        if x0f < need:
            t.set_x(t.get_position()[0] + (need - x0f) / max(bb.width, 1e-6))
    # x_fig pins the letter to a shared column; panels whose axes are indented (e.g. for a y-axis
    # label) would otherwise put their letter further right than the rest of the column.
    x = bb.x0 if x_fig is None else x_fig
    fig.text(x, y, letter.upper(), fontsize=LETTER_FS,
             fontweight="bold", ha="left", va="center")


# =====================================================================
# STYLE — one hierarchy for every recomposed plate
#   panel letter  >  panel title  >  in-panel body text
# Body text is uniform *within* each panel; a dense panel may take a smaller body size than a
# roomy one, which is why restyle() takes it per panel rather than fixing one global value.
# =====================================================================
LETTER_FS = 12.0     # panel letters (largest level)
TITLE_FS = 8.8       # panel titles
BODY_FS = 6.4        # default in-panel text (labels, ticks, legends, annotations)
BODY_FS_DENSE = 5.4  # dense panels: long label lists / many-row legends


def letters_row(fig, pairs, **kw):
    """Letter every panel of a row, each in its own top-left corner.

    With the letters inside the axes, a row of panels that share a top edge (see align_row) gets
    letters at one height automatically, and panels in the same column share one vertical line.
    The helper is kept so callers read the same way as before.
    """
    for ax, letter in pairs:
        letter_at(fig, ax, letter, **kw)


def restyle(ax, body=BODY_FS, title=TITLE_FS, centre_title=True):
    """Enforce the hierarchy on one panel: centred title + one uniform body size.

    Everything the panel draws that is not its title or its letter is set to `body` — axis
    labels, tick labels, legends, and free text. Call this AFTER the panel function and BEFORE
    adding the panel letter, so the letter keeps LETTER_FS.
    """
    # matplotlib keeps SEPARATE title artists for loc="left"/"center"/"right", so switching a panel
    # from a left title to a centred one leaves the left one drawn as well - which rendered C's and
    # E's titles twice, overlapping. Clear the other two before setting the wanted one.
    keep = "center" if centre_title else "left"
    # Read the text BEFORE blanking the other artists: a panel whose title was set with loc="left"
    # keeps it on ax._left_title, so blanking first would destroy the only copy.
    txt = ax.title.get_text() or (ax._left_title.get_text() if ax._left_title else "")
    for side, artist in (("left", ax._left_title), ("center", ax.title),
                         ("right", ax._right_title)):
        if artist is not None and side != keep:
            artist.set_text("")
    if txt:
        # weight is normalised too: some source panels set bold titles and others do not, which made
        # C and D stand out as heavier than the rest of the plate.
        ax.set_title(txt, fontsize=title, loc=keep, pad=4.0, fontweight="normal")
    for t in ax.texts:
        t.set_fontsize(body)
    for t in (ax.xaxis.label, ax.yaxis.label):
        t.set_fontsize(body)
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        t.set_fontsize(body)
    lg = ax.get_legend()
    if lg is not None:
        for t in lg.get_texts():
            t.set_fontsize(body)
        if lg.get_title():
            lg.get_title().set_fontsize(body)
    for cb in getattr(ax, "_colorbars", []):
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(body)
        if cb.ax.get_ylabel():
            cb.ax.set_ylabel(cb.ax.get_ylabel(), fontsize=body)
    return ax


def title_after(ax, x_min, fs=TITLE_FS):
    """Shift a panel's title right if it would start before figure fraction x_min.

    Panels whose centred title is wide enough to reach the plate's left edge leave no room for the
    panel letter beside it, unless the title moves. Only used where a title actually collides.
    """
    fig = ax.figure
    if title_left(fig, ax) >= x_min:
        return ax
    return shift_title_right(ax, shift=0.0, fs=fs, x_min=x_min)


def centre_title(ax, fs=TITLE_FS):
    """Centre a panel's title without touching any other text.

    For panels whose text deliberately uses more than one size (a dense label list plus a smaller
    annotation), restyle() would flatten them onto one value, so only the title is adjusted here.
    """
    # A panel whose title was set with loc="left" keeps it on ax._left_title, where get_title()
    # cannot see it - so move it onto the centred artist and blank the others.
    txt = ax.get_title() or (ax._left_title.get_text() if ax._left_title is not None else "")
    if txt:
        for art in (ax._left_title, ax._right_title):
            art.set_text("")
        ax.set_title(txt, fontsize=fs, loc="center", pad=4.0, fontweight="normal")
    return ax


def shift_title_right(ax, shift=0.022, fs=TITLE_FS, x_min=None):
    """Keep a panel's centred title but start it a little right of the axes' left edge.

    Some full-width panels (C and F here) have a centred title that begins exactly at the axes' left
    edge, leaving no room for the panel letter: above the axes it lands on the title, left of the
    axes it overlaps the first word. The title is re-emitted as figure text shifted right by roughly
    the letter's width, which opens that corner. Title text is unchanged; only its x moves.
    """
    txt = ax.get_title()
    if not txt:
        return ax
    ax.set_title("")
    bb = ax.get_position()
    fig = ax.figure
    fs_pt = fs
    # centre the shifted title on the axes so it still reads as centred
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    tb = ax.text(0, 0, txt, fontsize=fs_pt)
    w = tb.get_window_extent(renderer=fig.canvas.get_renderer()).width
    tb.remove()
    w_fig = inv.transform((w, 0))[0]
    x = bb.x0 + (bb.width - w_fig) / 2.0 + shift
    if x_min is not None:
        x = max(x, x_min)
    fig.text(x, bb.y1 + 0.012, txt, fontsize=fs_pt, ha="left", va="bottom")
    return ax


def align_row(fig, axes):
    """Put a row of axes back on one horizontal line by restoring their gridspec cells.

    A panel that lets matplotlib reserve space inside itself (a colourbar attached with `ax=`)
    ends up shorter than its neighbours, so the row's tops and bottoms no longer match. The cell
    bounds are authoritative: `SubplotSpec.get_position` returns the cell, not the stolen rect.
    Axes placed by explicit rect (no subplotspec) are left alone.
    """
    fixed = []
    for a in axes:
        if a is None:
            continue
        ss = a.get_subplotspec()
        if ss is None:
            fixed.append(a)
            continue
        cell = ss.get_position(fig)
        a.set_position([cell.x0, cell.y0, cell.width, cell.height])
        fixed.append(a)
    return fixed


def build_main5():
    """Main 5 - MR candidate-gene expression across lesional cell types (10 panels).

    A-E the five candidate genes on the OLP UMAP (GDF11, SUOX, EGFL8, TAPBP, HLA-DPA1), sharing one
    colour scale; F-J the same five genes as per-cell-type dotplots.

    Arrangement: two blocks over a shared 3-column grid. The UMAP block runs 3 across over two rows
    (A B C / D E + the note), so each map is ~55 mm instead of the ~33 mm it got when five ran across
    a single row. The dotplot block stays 2 across over three rows (F G / H I / J + keys): the
    longest cell-type label ("c16 CD8+ T (GZMB+/CXCL13+ cytotoxic-exhausted-like)") is ~57 mm at the
    label size, so a 3-across dotplot would push the label column over its neighbours.

    The bulleted caveat paragraph is kept on the plate (its content is not stated elsewhere); it
    sits in the UMAP block's sixth cell, immediately above the dotplots it describes, which frees
    the dotplot keys cell to hold the two keys alone.
    """
    print("== Main 5 (candidate-gene expression across cell types) ==")
    cand = run("plot_scrna_candidates_gse211630.py")   # candidate UMAP panels
    ero = importlib.import_module("scrna_erosion_gse211630")   # erosion dotplot grid

    fig = new_fig(180, 250)
    outer = fig.add_gridspec(2, 1, height_ratios=[1.86, 3.30],
                             left=0.045, right=0.975, top=0.968, bottom=0.016,
                             hspace=0.15)
    gsTop = outer[0].subgridspec(2, 3, wspace=0.14, hspace=0.46)
    gsBot = outer[1].subgridspec(3, 2, wspace=0.50, hspace=0.55)

    # ---- A-E: the five candidate-gene maps, 3 across over two rows
    up_axes = []
    for i, g in enumerate(cand["CAND"]):
        ax = fig.add_subplot(gsTop[i // 3, i % 3])
        if cand["panel_candidate_umap"](ax, g) is None:
            ax.set_title(g, loc="left", fontsize=8)
        up_axes.append(ax)
    for ax in up_axes:
        restyle(ax)                       # centred gene title, one body size
    letters_row(fig, ((up_axes[0], "A"), (up_axes[1], "B"), (up_axes[2], "C")))
    letters_row(fig, ((up_axes[3], "D"), (up_axes[4], "E")))

    # The UMAP block's sixth cell stays empty: five maps in a 3x2 block leave one free.
    # The caveat paragraph lives in the dotplot keys cell below, which is ~55 mm wide (2-across)
    # and the only cell wide enough to hold its 45-character lines without wrapping.

    # ---- F-J: the five dotplots, then the keys in the grid's free sixth cell
    dot_axes = [fig.add_subplot(gsBot[0, 0]), fig.add_subplot(gsBot[0, 1]),
                fig.add_subplot(gsBot[1, 0]), fig.add_subplot(gsBot[1, 1]),
                fig.add_subplot(gsBot[2, 0])]
    axKeys = fig.add_subplot(gsBot[2, 1])
    # A flat 1-D array of six axes: erosion_grid indexes it with .flat. A ragged np.array here would
    # silently collapse into a list.
    axes6 = np.empty(6, dtype=object)
    axes6[:5] = dot_axes
    axes6[5] = axKeys
    # notes at their source size: this cell is the widest one on the plate.
    ero.erosion_grid(axes6, label_fontsize=5.4)
    for ax in axes6[:5]:
        restyle(ax, body=5.4, title=9.4)
    restyle(axes6[5], title=0)
    letters_row(fig, ((dot_axes[0], "F"), (dot_axes[1], "G")))
    letters_row(fig, ((dot_axes[2], "H"), (dot_axes[3], "I")))
    letters_row(fig, ((dot_axes[4], "J"),))

    write(fig, os.path.join(MAINDIR, "fig06_candidate_expression"), 180, 250)


def build_main6():
    """Main 6 - tissue-level target mapping, T/NK state and communication (6 panels).

    A reverse-mapping dotplot (target message per subtype)
    B Oaxaca tissue decomposition
    C T/NK module scores: the seven signature panels, spanning as a 2x4 band
    D presence-flux heatmap | E strength shift | F ligand extremes
    """
    print("== Main 6 (target mapping, T/NK state, communication) ==")
    f8 = importlib.import_module("build_final_sc_fig8")
    comp = importlib.import_module("scrna_plot_composition_tnk_gse211630")

    fig = new_fig(240, 270)
    # Wider canvas than the other plates: at 180 mm this plate had to stack four bands and ran
    # ~320 mm tall, and every panel was narrow. At 240 mm the same bands give each panel ~33% more
    # width and the plate reads landscape instead of a column strip. NOTE: 240 mm exceeds the
    # usual 180 mm single-page limit, so this plate would be scaled to 75% for submission.
    # D | E sit on one row and F spans below; row 0 keeps A | B, and the C band is 4 columns of
    # 2 rows, the widest arrangement whose per-panel labels still fit.
    # hspace 0.75 left ~34 mm between every band, which read as dead space. 0.34 -> ~17 mm, which
    # still holds A's two-row key below its axes and the C band's row-to-row title gap.
    outer = fig.add_gridspec(4, 1, height_ratios=[2.30, 1.42, 0.90, 0.72],
                             left=0.170, right=0.950, top=0.968, bottom=0.018,
                             hspace=0.34)
    # B starts at the same x as E (see gsDE below): the two right-column panels then share one
    # letter column, which is what makes B and E line up vertically.
    gsAB = outer[0].subgridspec(1, 2, width_ratios=[1.676, 1.0], wspace=0.37)
    # 0.80 of a panel height between the two signature rows: the lower row's title has to clear
    # the upper row's "module score" x-label.
    gsC = outer[1].subgridspec(2, 4, wspace=0.34, hspace=0.80)
    gsDE = outer[2].subgridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.95)

    axA = fig.add_subplot(gsAB[0, 0])
    axB = fig.add_subplot(gsAB[0, 1])
    # A is one column of a 180 mm plate, not half of a 12.2 in one: shrink its headers and row
    # labels so the five gene badges and the 21 row names stay inside the panel.
    f8.panel_A(axA, tag=None, badge_fontsize=4.0, gene_fontsize=6.6, ylabel_fontsize=5.2,
               size_key=False, legend_bbox=(-0.02, -0.055),
               legend_fontsize=5.0, legend_labelspacing=0.30, legend_ncol=2,
               title="Target-message direction inside each\ncell subtype (EOLP vs NEOLP)")
    # Two-line titles: at ~43 mm a descriptive one-line title cannot fit at 8.2 pt.
    # B's key goes inside the axes (empty lower-right); below the axes it reached into the C band.
    f8.panel_B(axB, tag=None,
               title="Tissue target message:\ncomposition vs per-cell shift",
               legend_loc="lower right", legend_bbox=(0.985, 0.02), legend_fontsize=5.4)

    cells = [(r, c) for r in range(2) for c in range(4)]
    c_axes = []
    for i, sig in enumerate(comp.SIG_ORDER):
        ax = fig.add_subplot(gsC[cells[i]])
        comp.panel_tnk_signature(ax, sig, ylabels=(i % 4 == 0))
        c_axes.append(ax)
    # The band's last cell is free (7 signatures in a 4x2 grid). The signature panels colour their
    # points by condition and nothing else on the plate keys that colour, so it goes there rather
    # than leaving the corner empty.
    axCKey = fig.add_subplot(gsC[cells[7]])
    restyle(axCKey, title=0)
    axCKey.axis("off")            # plain legend host: no frame, no ticks, no labels
    axCKey.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=4.6,
                                      color=comp.COL["NEOLP"], label="NEOLP (n = 3)"),
                           plt.Line2D([], [], marker="o", ls="", ms=4.6,
                                      color=comp.COL["EOLP"], label="EOLP (n = 2)")],
                  loc="upper left", frameon=True, framealpha=0.95, edgecolor="#DDDDDD",
                  fontsize=BODY_FS_DENSE, handletextpad=0.5, labelspacing=0.55,
                  borderpad=0.6, title="condition", title_fontsize=BODY_FS_DENSE)

    axD = fig.add_subplot(gsDE[0, 0])
    axE = fig.add_subplot(gsDE[0, 1])
    # The colour bar gets its own axes: attached with `ax=` it would take its space out of panel D,
    # leaving D shorter than E on the same row.
    caxD = fig.add_axes([0.0, 0.0, 0.012, 0.08])       # provisional; positioned after D is aligned
    axF = fig.add_subplot(outer[3, 0])
    f8.panel_C(axD, tag=None, colorbar=True, cbar_ax=caxD,
               title="Incoming focus signalling\nper receiver")
    f8.panel_D(axE, tag=None, title="Sustained T/NK signalling\ndeclines in EOLP")
    f8.panel_E(axF, tag=None)

    # D | E and A | B are rows: same top and bottom for each pair.
    align_row(fig, (axA, axB))
    align_row(fig, (axD, axE))
    # Put the colour bar just right of D, then make room for it deterministically: E's wide
    # cell-type labels are measured, and E is pushed right (narrowing it) until it clears the bar.
    bD, bE = axD.get_position(), axE.get_position()
    caxD.set_position([bD.x1 + 0.014, bD.y0, 0.013, bD.height])
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    lab_w = max([inv.transform((bb.x1, bb.y1))[0] - inv.transform((bb.x0, bb.y0))[0]
                 for bb in (t.get_window_extent(renderer=rend) for t in axE.get_yticklabels())
                 if bb.width > 0] or [0.0])
    need = caxD.get_position().x1 + 0.008 + lab_w
    if bE.x0 < need:
        axE.set_position([need, bE.y0, max(0.10, bE.x1 - need), bE.height])

    # Hierarchy: A keeps its own two label sizes (see centre_title); everything else gets one
    # body size per panel.
    centre_title(axA)
    restyle(axB)
    for ax in c_axes:
        restyle(ax, body=BODY_FS_DENSE)
    restyle(axD, body=BODY_FS_DENSE)
    restyle(axE)
    restyle(axF)
    restyle(caxD, title=0)          # colour bar label/ticks follow the dense body size

    # Letters: A | B on one line, D | E on one line, C above the signature band, F above its band.
    letters_row(fig, ((axA, "A"), (axB, "B")))
    letters_row(fig, ((axD, "D"), (axE, "E")))
    # A, C, D and F all start at the plate's left margin, so their letters share one column; B and
    # E share the right column (their axes' left edges are aligned above).
    letter_at(fig, c_axes[0], "C")
    letter_at(fig, axF, "F")

    write(fig, os.path.join(MAINDIR, "fig07_target_mapping"), 240, 270)


def build_supp1():
    """Supp 1 - single-cell dataset provenance and per-cell quality control (3 panels).

    Mirrors the source plate's five-band structure (panel, colour key, panel, colour key, panel):
    the keys are separate functions in the module, so they need their own gridspec cells.
    """
    print("== Supp 1 (dataset provenance and per-cell QC) ==")
    qc = importlib.import_module("scrna_suppfig5_dataset_qc_gse211630")
    fig = new_fig(180, 240)
    gs = fig.add_gridspec(5, 1, height_ratios=[1.00, 0.13, 2.05, 0.13, 2.45],
                          left=0.075, right=0.975, top=0.972, bottom=0.012,
                          hspace=0.32)
    a_ax = qc.panel_a(fig, gs[0], with_letter=False)
    qc.legend_a(fig, gs[1])
    b_ax = qc.panel_b(fig, gs[2], with_letter=False)
    qc.legend_b(fig, gs[3])
    c_ax = qc.panel_c(fig, gs[4], with_letter=False)
    # Dense QC panels: every panel's own text takes the dense body size, titles are centred, and
    # the letters are placed here so they are uppercase and sized by this plate's hierarchy.
    for ax in list(a_ax) + list(b_ax) + list(c_ax):
        restyle(ax, body=BODY_FS_DENSE)
    # Letters go in the plate's left margin, above their panel: this module's panels are built
    # from multi-axes gridspecs, so anchoring to one axes inside a panel lands mid-figure.
    def plate_letter(group, letter):
        top = max(a.get_position().y1 for a in group)
        fig.text(0.014, top + 0.008, letter, fontsize=LETTER_FS, fontweight="bold",
                 ha="left", va="bottom")

    plate_letter(a_ax, "A")
    plate_letter(b_ax, "B")
    plate_letter(c_ax, "C")
    write(fig, os.path.join(SUPPDIR, "suppfig01_dataset_qc"), 180, 240)


def build_supp2():
    """Supp 2 - candidate-gene expression landscape in the OLP atlas (A dot plot, B UMAP facets)."""
    print("== Supp 2 (candidate expression landscape) ==")
    ca = importlib.import_module("scrna_suppfig8_candidate_expression_gse211630")
    fig = new_fig(180, 215)
    # title=None suppresses the module's left-flush subtitle line; a centred one is set below so
    # this plate follows the same "title centred" rule as the rest.
    A = ca.panel_A(fig, [0.300, 0.480, 0.420, 0.465],
                   [0.745, 0.620, 0.020, 0.180], [0.720, 0.430, 0.255, 0.170],
                   title=None, with_letter=False)
    B = ca.panel_B(fig, [0.045, 0.180, 0.865, 0.190], [0.930, 0.195, 0.016, 0.150],
                   with_letter=False)
    restyle(A["axA"], body=BODY_FS_DENSE)
    A["axA"].set_title("Candidate expression across 19 OLP cell types\n"
                       "(c12/c18 doublet clusters excluded)", fontsize=TITLE_FS,
                       loc="center", pad=5.0)
    for ax in B["b_ax"]:
        restyle(ax)
    letter_at(fig, A["axA"], "A")
    letter_at(fig, B["b_ax"][0], "B")
    write(fig, os.path.join(SUPPDIR, "suppfig02_candidate_landscape"), 180, 215)


def build_supp3():
    """Supp 3 - per-gene detection summary (A) and erosive composition shift (B)."""
    print("== Supp 3 (detection summary and composition shift) ==")
    ca = importlib.import_module("scrna_suppfig8_candidate_expression_gse211630")
    comp = importlib.import_module("scrna_plot_composition_tnk_gse211630")
    fig = new_fig(180, 200)
    # Panel A keeps the source's own letter here (its axes-anchored "C" is left in place), so the
    # plate reads as one continuous pair of panels rather than two renumbered ones.
    C = ca.panel_C(fig, [0.050, 0.680, 0.925, 0.250], with_letter=False)
    axB = fig.add_subplot(fig.add_gridspec(
        1, 1, left=0.300, right=0.975, top=0.560, bottom=0.085)[0, 0])
    comp.panel_subtype_dots(axB, letter=None)              # B: composition shift
    restyle(C["axC"], body=BODY_FS_DENSE)                  # A: summary table, no title
    restyle(axB)
    # A's letter goes at the top-left of its own rect (letter_at is axes-relative, and the
    # table axes has no useful data transform to anchor against).
    fig.text(0.048, 0.938, "A", fontsize=LETTER_FS, fontweight="bold", ha="left", va="bottom")
    letter_at(fig, axB, "B")
    write(fig, os.path.join(SUPPDIR, "suppfig03_detection_composition"), 180, 200)


BUILDERS = {"main4": build_main4, "main5": build_main5, "main6": build_main6,
            "supp1": build_supp1, "supp2": build_supp2, "supp3": build_supp3}

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which.startswith("__run_main:"):
        # Verification hook: run a refactored figure script's own main() UNPATCHED so its
        # canonical output is really written and can be byte-compared. The savefig redirect must
        # stay off here, otherwise the compare passes against an untouched file.
        import importlib
        mod = importlib.import_module(which.split(":", 1)[1])
        mod.main()
    else:
        enable_redirect()
        todo = BUILDERS.values() if which == "all" else [BUILDERS[which]]
        for fn in todo:
            fn()
        disable_redirect()
    print("done")
