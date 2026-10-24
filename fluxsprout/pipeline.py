"""FluxSprout pipeline: context models, E-Flux and knockout screens (HiGHS engine).
The notebook calls these functions and caches their outputs in results/."""
import json, pickle, time, numpy as np, pandas as pd
from fslib import *
from ctxlib import eval_gpr, localT2
from hlp import Net, fastcore, gimme, gapfill_support

THRESH = {'default': (25, 75), 'loose': (15, 65), 'strict': (40, 85)}
EC_STATES = ['aDA', 'aCran', 'aTrunk_mixed', 'ISV_tip', 'ISV_cyc', 'cranV', 'cranV_cyc', 'vPCV', 'vPCV_cyc', 'vCV',
             'angioblast', 'angioblast_cyc', 'endocard', 'HE_EMP']
SPROUT = ['ISV_tip', 'ISV_cyc']                                 # intersegmental sprouting population
CYCLING = ['ISV_cyc', 'vPCV_cyc', 'cranV_cyc', 'angioblast_cyc']
QUIESCENT = ['aDA', 'aCran', 'vPCV', 'vCV', 'cranV']             # established, non-sprouting vessels
NONEC = ['Erythroid', 'Macrophage', 'Myeloid', 'Muscle', 'Neural prog.', 'Neuron', 'Fibroblast', 'Mesenchyme',
         'Post. mes.', 'Epithelium', 'Periderm', 'Xantho.', 'LPM/cardiac']

def load_ref():
    ref = pickle.load(open('cache/ref_medium_consistent.pkl', 'rb'))
    return ref, Net(ref)

def expression():
    M = pd.read_csv('results/mean_cp10k_by_group.csv', index_col=0)
    return M, np.log1p(M)

def protect_ids(ref):
    p = {BIOMASS, NGAM, 'MAR10062', 'MAR10063', 'MAR10064', 'MAR10065'} | {r.id for r in ref.exchanges}
    return p & {r.id for r in ref.reactions}

MU_REF = {}
CORE_PROTECT = {BIOMASS, 'MAR10062', 'MAR10063', 'MAR10064', 'MAR10065'}   # forced into every FASTCORE core

def rxn_scores(ref, gene_vals):
    d = gene_vals.to_dict()
    return np.array([eval_gpr(r.gene_reaction_rule, d) for r in ref.reactions], float)

def build_context(ref, net, setting='default', algo='fastcore', groups=None, L=None, verbose=True):
    """Return {group: sorted reaction ids}. FASTCORE on local-T2 cores, or GIMME on log expression."""
    if L is None: _, L = expression()
    act, thr = localT2(L, *THRESH[setting])
    net.reset_bounds(); MU_REF['v'] = net.fba()
    prot = protect_ids(ref); pidx = {net.idx[r] for r in prot}
    out = {}
    for g in (groups or list(L.columns)):
        t = time.time()
        if algo == 'fastcore':
            sc = rxn_scores(ref, act[g].astype(float))
            core = (set(np.where(sc == 1.0)[0]) - {i for i in range(net.n) if net.exch[i]}) | {net.idx[r] for r in CORE_PROTECT if r in net.idx}
            ngam = net.idx[NGAM]; lo = net.lb[ngam]; net.lb[ngam] = 0.0   # no forced flux during extraction
            A, dropped = fastcore(net, sorted(core))
            net.lb[ngam] = lo
            # closure: support biomass at 10% of the reference optimum with maintenance on (micronutrient routes)
            A |= gapfill_support(net, A | pidx, net.obj, 0.1 * MU_REF.get('v', net.fba()))
            keep = A | pidx
        else:
            sc = rxn_scores(ref, L[g]); sc[net.exch] = np.nan
            keep = gimme(net, sc, thr[1]) | set(np.where(sc >= thr[1])[0]) | pidx
        out[g] = sorted(net.ids[i] for i in keep)
        net.reset_bounds()
        if verbose: print(f'{algo}/{setting} {g:16s} {len(keep):5d} rxns ({time.time()-t:.0f}s)', flush=True)
    return out

def mask_net(net, keep_ids, medium=None):
    """Restrict the shared reference network to a context model by bounds (and optionally change medium)."""
    net.reset_bounds()
    keep = np.zeros(net.n, bool); keep[[net.idx[r] for r in keep_ids if r in net.idx]] = True
    net.lb = np.where(keep, net.lb0, 0.0); net.ub = np.where(keep, net.ub0, 0.0)
    if medium is not None:
        for k, v in medium.items():
            if EXCH[k] in net.idx: net.lb[net.idx[EXCH[k]]] = -float(v) if keep[net.idx[EXCH[k]]] else 0.0
    net._h = None
    return net

def ko_screen(net, ctx, medium=None, genes=None):
    """Gene knockout growth ratios for every context model -> DataFrame genes x groups (+ WT row)."""
    res, wt = {}, {}
    for g, keep in ctx.items():
        mask_net(net, keep, medium)
        mu, r = net.gene_screen(genes)
        wt[g] = mu; res[g] = pd.Series(r)
    R = pd.DataFrame(res).clip(0, 1)
    R.loc['__WT__'] = pd.Series(wt)
    return R

MEDIA = {'baseline': {}, 'hypoxia': {'glucose': 30, 'O2': 1, 'glutamine': 20},
         'ischemia': {'glucose': 2, 'O2': 0.8, 'glutamine': 2}, 'hyperglycemia': {'glucose': 40},
         '+carnitine': {'carnitine': 1}, 'low amino acids': {a: 0.1 for a in AA}}

def full_medium(name):
    return baseline_medium(extra=MEDIA[name])

# ------------------------------------------------------------------ E-Flux (expression-scaled bounds)
EXEMPT_EFLUX = {BIOMASS, NGAM, 'MAR10062', 'MAR10063', 'MAR10064', 'MAR10065'}

def eflux_global(net, ref, gene_cp10k, U, eps=0.05, norm='p95'):
    """E-Flux with a global capacity U: bound_j = s_hat_j * U for GPR reactions (exchanges/biomass/NGAM/pools exempt).
    norm='p95': s/95th percentile, capped at 1 (default); norm='max': s/max (conventional E-Flux)."""
    net.reset_bounds()
    s = rxn_scores(ref, gene_cp10k); has = ~np.isnan(s)
    den = np.nanpercentile(s[has], 95) if norm == 'p95' else np.nanmax(s)
    shat = np.full(net.n, np.nan); shat[has] = np.clip(s[has] / den, eps, 1.0)
    ex = np.array([net.exch[i] or net.ids[i] in EXEMPT_EFLUX for i in range(net.n)])
    sel = has & ~ex
    b = shat * U
    net.ub = np.where(sel & (net.ub0 > 0), np.minimum(net.ub0, b), net.ub0)
    net.lb = np.where(sel & (net.lb0 < 0), np.maximum(net.lb0, -b), net.lb0)
    net._h = None
    return shat

RXI = dict(PGK='MAR04368', PK='MAR04358', ATPS='MAR06916', CI='MAR06921', CIII='MAR06918', CIV='MAR06914',
           PHGDH='MAR03839', PPAT='MAR04406', CPSII='MAR04034', CACT='MAR02629', CARNm='MAR02590', IDH2='MAR03958',
           GLC='MAR09034', O2='MAR09048', LAC='MAR09135', GLN='MAR09063', SDH='MAR04652')

def readouts_vec(net, v):
    g = lambda k: float(v[net.idx[RXI[k]]]) if RXI[k] in net.idx else 0.0
    mu = float(v[net.obj]); gly = max(g('PGK'), 0) + max(g('PK'), 0); ox = max(g('ATPS'), 0)
    return dict(growth=mu, glc_uptake=-g('GLC'), o2_uptake=-g('O2'), lac_secr=g('LAC'), gln_uptake=-g('GLN'),
                glyc_ATP=gly, oxphos_ATP=ox, glyc_ATP_frac=gly / (gly + ox) if gly + ox > 0 else np.nan,
                lac_per_glc=g('LAC') / -g('GLC') if -g('GLC') > 1e-9 else np.nan,
                PHGDH=g('PHGDH'), PPAT_per_growth=g('PPAT') / mu if mu > 1e-9 else np.nan,
                CPSII_per_growth=g('CPSII') / mu if mu > 1e-9 else np.nan, CACT=g('CACT'), SDH=-g('SDH'), IDH2=g('IDH2'))

def solve_unit(net, loopcheck=True):
    """pFBA + CycleFreeFlux loop check; returns (readouts, loop_gap, n_at_bound)."""
    v = net.pfba()
    if v is None: return None
    gap = np.nan
    if loopcheck:
        ll = net.cycle_free(v); gap = float(np.abs(v - ll).sum()) if ll is not None else np.nan
    atb = int(np.sum((np.abs(v) > 1e-6) & ((np.abs(v - net.ub) < 1e-6) | (np.abs(v - net.lb) < 1e-6)) & ~net.exch))
    return readouts_vec(net, v), gap, atb, v
