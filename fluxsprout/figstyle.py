"""Publication figure style for FluxSprout.

Metabolic Engineering: full width 190 mm (we use 180 mm), single column 90 mm. Sans-serif type throughout
(Arial or Helvetica when installed, otherwise the metrically identical Liberation Sans), smallest text 6.5 pt at
print size, panel letters 10 pt bold. Layout is handled by matplotlib's constrained layout, and panel letters are
placed after layout so they never collide with titles, tick labels or neighbouring panels.
"""
import matplotlib as mpl, matplotlib.pyplot as plt

MM = 1 / 25.4
FULL, HALF = 180 * MM, 88 * MM
INK, MUTED, GRID = '#1f2328', '#5b6470', '#e5e7eb'
C = dict(tip='#d97706', cyc='#059669', art='#dc2626', ven='#2563eb', other='#9ca3af', nonec='#7c3aed',
         ok='#059669', no='#9ca3af', bad='#b91c1c', accent='#2563eb', orange='#ea580c')
STATE_COL = {'aDA': '#dc2626', 'aCran': '#f87171', 'aTrunk_mixed': '#fca5a5', 'ISV_tip': '#d97706', 'ISV_cyc': '#059669',
             'cranV': '#93c5fd', 'cranV_cyc': '#60a5fa', 'vPCV': '#2563eb', 'vPCV_cyc': '#1d4ed8', 'vCV': '#1e3a8a',
             'angioblast': '#7c3aed', 'angioblast_cyc': '#a78bfa', 'endocard': '#7f1d1d', 'HE_EMP': '#6b7280', 'unresolved': '#d1d5db'}
FONT_STACK = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']
FS = dict(base=7, small=6.5, title=7.5, letter=10)


def setup():
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': FONT_STACK, 'mathtext.default': 'regular', 'mathtext.fontset': 'custom', 'mathtext.rm': 'sans', 'mathtext.it': 'sans:italic', 'mathtext.bf': 'sans:bold', 'mathtext.cal': 'sans', 'mathtext.sf': 'sans',
        'font.size': FS['base'], 'axes.titlesize': FS['title'], 'axes.labelsize': FS['base'],
        'xtick.labelsize': FS['small'], 'ytick.labelsize': FS['small'], 'legend.fontsize': FS['small'],
        'legend.title_fontsize': FS['small'], 'figure.titlesize': FS['title'], 'figure.labelsize': FS['base'],
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4,
        'ytick.minor.width': 0.4, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.major.pad': 2,
        'ytick.major.pad': 2, 'axes.labelpad': 3, 'axes.titlepad': 5, 'lines.linewidth': 0.9,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': INK, 'axes.labelcolor': INK,
        'xtick.color': INK, 'ytick.color': INK, 'text.color': INK, 'axes.titleweight': 'bold',
        'axes.titlelocation': 'left', 'legend.frameon': False, 'legend.handlelength': 1.2,
        'legend.handletextpad': 0.4, 'legend.columnspacing': 1.0, 'legend.borderaxespad': 0.2,
        'figure.dpi': 150, 'savefig.dpi': 600, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03,
        'figure.constrained_layout.h_pad': 0.04, 'figure.constrained_layout.w_pad': 0.04,
        'figure.constrained_layout.hspace': 0.04, 'figure.constrained_layout.wspace': 0.04,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})


def new_figure(height_mm, width=FULL):
    return plt.figure(figsize=(width, height_mm * MM), layout='constrained')


def panel(ax, letter, **_):
    """Register a panel letter; it is drawn by label_panels() once the layout is final."""
    f = ax.figure
    while getattr(f, '_parent', None) is not None and f._parent is not f: f = f._parent
    f._fs_panels = getattr(f, '_fs_panels', []) + [(ax, letter)]


def label_panels(fig, gap_pt=3):
    """Place each registered letter on the baseline of the panel title, left of everything in that panel."""
    pairs = getattr(fig, '_fs_panels', [])
    if not pairs: return
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    fig.set_layout_engine('none')                   # freeze the constrained layout before adding letters
    px = fig.dpi / 72.0
    for ax, letter in pairs:
        tb = ax.get_tightbbox(r); title = ax._left_title if ax.get_title(loc='left') else (ax.title if ax.get_title() else None)
        t = fig.text(0, 0, letter, fontsize=FS['letter'], fontweight='bold', ha='left', va='baseline')
        w = t.get_window_extent(r).width
        if title is not None:
            tbb = title.get_window_extent(r)
            x = min(tb.x0, tbb.x0 - w - gap_pt * px); y = tbb.y0 + title.get_fontsize() * 0.22 * px
        else:
            x = tb.x0 - w - gap_pt * px if tb.x0 - w - gap_pt * px > 0 else tb.x0; y = tb.y1 - FS['letter'] * px
        t.set_position(inv.transform((x, y)))


def save(fig, name, outdir='figures'):
    import os; os.makedirs(outdir, exist_ok=True)
    label_panels(fig)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/{name}.{ext}')
