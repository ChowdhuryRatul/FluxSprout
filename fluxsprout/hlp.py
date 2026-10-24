"""hlp: a small, fast constraint-based engine on top of HiGHS (highspy), used by FluxSprout.

Why: COBRApy+GLPK needs ~2 s per LP on this 11k-reaction network; warm-started HiGHS needs ~0.03 s,
which makes full knockout screens, FVA, loop audits and FASTCORE over many models practical.
All routines were checked against COBRApy on the E. coli core model (see notebook, section 0).
"""
import numpy as np, scipy.sparse as sp, highspy
from cobra.util.array import create_stoichiometric_matrix

INF = highspy.kHighsInf

def _highs(quiet=True):
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.setOptionValue('primal_feasibility_tolerance', 1e-8)
    h.setOptionValue('dual_feasibility_tolerance', 1e-8)
    return h

def _pass(h, A, c, lb, ub, rlo, rhi, sense_max=False):
    lp = highspy.HighsLp()
    A = sp.csc_matrix(A)
    lp.num_col_, lp.num_row_ = A.shape[1], A.shape[0]
    lp.col_cost_ = np.asarray(c, float); lp.col_lower_ = np.asarray(lb, float); lp.col_upper_ = np.asarray(ub, float)
    lp.row_lower_ = np.asarray(rlo, float); lp.row_upper_ = np.asarray(rhi, float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    lp.sense_ = highspy.ObjSense.kMaximize if sense_max else highspy.ObjSense.kMinimize
    h.passModel(lp)

class Net:
    """Stoichiometric network snapshot of a COBRA model (bounds, S, GPR map)."""
    def __init__(self, model, objective=None):
        self.ids = [r.id for r in model.reactions]
        self.idx = {r: i for i, r in enumerate(self.ids)}
        self.S = create_stoichiometric_matrix(model, array_type='dok').tocsc()
        self.m, self.n = self.S.shape
        self.lb0 = np.array([r.lower_bound for r in model.reactions], float)
        self.ub0 = np.array([r.upper_bound for r in model.reactions], float)
        self.lb, self.ub = self.lb0.copy(), self.ub0.copy()
        self.exch = np.array([r.boundary for r in model.reactions])
        obj = objective or [r.id for r in model.reactions if r.objective_coefficient != 0][0]
        self.obj = self.idx[obj]
        # gene -> reactions disabled by its knockout
        self.gene_ko = {}
        for g in model.genes:
            dead = [self.idx[r.id] for r in g.reactions if not r.gpr.eval(knockouts={g.id})]
            if dead: self.gene_ko[g.id] = dead
        self._h = None

    # ------------------------------------------------------------ basic LP
    def reset_bounds(self):
        self.lb, self.ub = self.lb0.copy(), self.ub0.copy(); self._h = None

    def set_bounds(self, rid, lo, hi):
        i = self.idx[rid]; self.lb[i], self.ub[i] = lo, hi
        if self._h is not None: self._h.changeColBounds(i, float(lo), float(hi))

    def _base(self):
        if self._h is None:
            self._h = _highs()
            c = np.zeros(self.n); c[self.obj] = 1.0
            _pass(self._h, self.S, c, self.lb, self.ub, np.zeros(self.m), np.zeros(self.m), sense_max=True)
        return self._h

    def fba(self, return_x=False):
        h = self._base(); h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            return (0.0, None) if return_x else 0.0
        v = h.getInfo().objective_function_value
        return (v, np.array(h.getSolution().col_value)) if return_x else v

    def ko_growth(self, rxn_idx):
        """Growth with reactions `rxn_idx` set to zero (warm start, bounds restored afterwards)."""
        h = self._base()
        for i in rxn_idx: h.changeColBounds(int(i), 0.0, 0.0)
        h.run()
        ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
        v = h.getInfo().objective_function_value if ok else 0.0
        for i in rxn_idx: h.changeColBounds(int(i), float(self.lb[i]), float(self.ub[i]))
        return max(v, 0.0)

    def gene_screen(self, genes=None):
        wt = self.fba(); out = {}
        for g in (genes or self.gene_ko):
            out[g] = self.ko_growth(self.gene_ko.get(g, [])) / wt if wt > 1e-9 else np.nan
        return wt, out

    def fva(self, rxn_idx=None, frac=0.0, lb_obj=None):
        """Flux variability (min,max) for reactions; optimum constraint obj >= frac*max."""
        h = _highs(); _pass(h, self.S, np.zeros(self.n), self.lb, self.ub, np.zeros(self.m), np.zeros(self.m))
        if frac > 0 or lb_obj is not None:
            mu = self.fba() if lb_obj is None else lb_obj / max(frac, 1e-12)
            h.changeColBounds(self.obj, float(frac * mu), float(self.ub[self.obj]))
        rxn_idx = range(self.n) if rxn_idx is None else rxn_idx
        lo, hi = {}, {}
        for i in rxn_idx:
            for sense, store in ((highspy.ObjSense.kMinimize, lo), (highspy.ObjSense.kMaximize, hi)):
                h.changeColCost(int(i), 1.0); h.changeObjectiveSense(sense); h.run()
                ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
                store[i] = h.getInfo().objective_function_value if ok else np.nan
                h.changeColCost(int(i), 0.0)
        return lo, hi

    # ------------------------------------------------------------ pFBA
    def pfba(self, frac=0.9999):
        mu, _ = self.fba(return_x=True)
        n, m = self.n, self.m
        # variables [v (n), a (n)]; rows: S v = 0 ; a - v >= 0 ; a + v >= 0
        I = sp.identity(n, format='csc')
        A = sp.vstack([sp.hstack([self.S, sp.csc_matrix((m, n))]), sp.hstack([-I, I]), sp.hstack([I, I])]).tocsc()
        lb = np.r_[self.lb, np.zeros(n)]; ub = np.r_[self.ub, np.full(n, INF)]
        lb[self.obj] = frac * mu
        rlo = np.r_[np.zeros(m), np.zeros(2 * n)]; rhi = np.r_[np.zeros(m), np.full(2 * n, INF)]
        c = np.r_[np.zeros(n), np.ones(n)]
        h = _highs(); _pass(h, A, c, lb, ub, rlo, rhi); h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return None
        return np.array(h.getSolution().col_value[:n])

    # ------------------------------------------------------------ loop removal (CycleFreeFlux, Desouki 2015)
    def cycle_free(self, v, fix=None):
        """Nearest loop-free flux: exchanges fixed at v, internal reactions confined to [min(0,v),max(0,v)],
        minimize sum |v|. `fix`: optional extra reaction index held at its value. Returns flux vector."""
        lb = np.where(self.exch, v, np.minimum(0, v)); ub = np.where(self.exch, v, np.maximum(0, v))
        if fix is not None: lb[fix] = ub[fix] = v[fix]
        c = np.sign(v) * (~self.exch)
        h = _highs(); _pass(h, self.S, c, lb, ub, np.zeros(self.m), np.zeros(self.m)); h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return None
        return np.array(h.getSolution().col_value)

    def loopless_fva(self, rxn_idx, frac=0.0):
        """Loop-corrected FVA (CycleFreeFlux-style): optimize, then test whether the optimum survives loop removal;
        if not, re-optimize inside the sign pattern of the loop-free solution. Returns {i:(min,max)} and a flag
        dict telling whether the plain FVA value was loop-dependent."""
        h = _highs(); _pass(h, self.S, np.zeros(self.n), self.lb, self.ub, np.zeros(self.m), np.zeros(self.m))
        if frac > 0:
            mu = self.fba(); h.changeColBounds(self.obj, float(frac * mu), float(self.ub[self.obj]))
        res, loopdep = {}, {}
        for i in rxn_idx:
            vals = []; dep = False
            for sense in (highspy.ObjSense.kMinimize, highspy.ObjSense.kMaximize):
                h.changeColCost(int(i), 1.0); h.changeObjectiveSense(sense); h.run()
                if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
                    vals.append(np.nan); h.changeColCost(int(i), 0.0); continue
                v = np.array(h.getSolution().col_value); opt = v[i]
                ll0 = self.cycle_free(v)
                if ll0 is not None and abs(ll0[i] - opt) <= 1e-6 * max(1.0, abs(opt)):
                    vals.append(opt)
                else:
                    dep = True
                    if ll0 is None: vals.append(np.nan)
                    else:
                        # re-optimize reaction i inside the loop-free sign pattern
                        lb = np.where(ll0 < -1e-9, self.lb, 0.0); ub = np.where(ll0 > 1e-9, self.ub, 0.0)
                        lb[i], ub[i] = self.lb[i], self.ub[i]
                        lb = np.where(self.exch, self.lb, lb); ub = np.where(self.exch, self.ub, ub)
                        h2 = _highs(); c = np.zeros(self.n); c[i] = 1
                        _pass(h2, self.S, c, lb, ub, np.zeros(self.m), np.zeros(self.m), sense_max=(sense == highspy.ObjSense.kMaximize))
                        h2.run()
                        vals.append(h2.getInfo().objective_function_value if h2.getModelStatus() == highspy.HighsModelStatus.kOptimal else ll0[i])
                h.changeColCost(int(i), 0.0)
            res[i] = tuple(vals); loopdep[i] = dep
        return res, loopdep

# ------------------------------------------------------------ FASTCORE (Vlassis 2014) on HiGHS
def _lp7(S, lb, ub, J, sign, eps, relaxed=False):
    m, n = S.shape; k = len(J)
    # vars [v (n), z (k)] ; rows: S v = 0 ; z_j - s_j v_j <= 0
    rows = sp.csc_matrix((np.r_[-np.array([sign[j] for j in J]), np.ones(k)],
                          (np.r_[np.arange(k), np.arange(k)], np.r_[J, n + np.arange(k)])), shape=(k, n + k))
    A = sp.vstack([sp.hstack([S, sp.csc_matrix((m, k))]), rows]).tocsc()
    c = np.r_[np.zeros(n), np.ones(k)]
    zlo = np.full(k, -INF) if relaxed else np.zeros(k)
    h = _highs(); _pass(h, A, c, np.r_[lb, zlo], np.r_[ub, np.full(k, eps)],
                        np.r_[np.zeros(m), np.full(k, -INF)], np.r_[np.zeros(m), np.zeros(k)], sense_max=True)
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return np.zeros(n)
    return np.array(h.getSolution().col_value[:n])

def _lp10(S, lb, ub, K, P, sign, eps, scale=1e3):
    m, n = S.shape; p = len(P)
    lb = lb * scale; ub = ub * scale
    lb = lb.copy(); ub = ub.copy()
    for k in K:
        if sign[k] > 0: lb[k] = max(lb[k], eps * scale)
        else: ub[k] = min(ub[k], -eps * scale)
    if p == 0:
        h = _highs(); _pass(h, S, np.zeros(n), lb, ub, np.zeros(m), np.zeros(m)); h.run()
        return np.array(h.getSolution().col_value) / scale if h.getModelStatus() == highspy.HighsModelStatus.kOptimal else np.zeros(n)
    P = np.asarray(P)
    E = sp.csc_matrix((np.ones(p), (np.arange(p), P)), shape=(p, n)); Ip = sp.identity(p, format='csc')
    A = sp.vstack([sp.hstack([S, sp.csc_matrix((m, p))]), sp.hstack([-E, Ip]), sp.hstack([E, Ip])]).tocsc()
    c = np.r_[np.zeros(n), np.ones(p)]
    h = _highs(); _pass(h, A, c, np.r_[lb, np.zeros(p)], np.r_[ub, np.full(p, INF)],
                        np.r_[np.zeros(m), np.zeros(2 * p)], np.r_[np.zeros(m), np.full(2 * p, INF)])
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return np.zeros(n)
    return np.array(h.getSolution().col_value[:n]) / scale

def fastcore(net, core_idx, eps=1e-4):
    S, lb, ub = net.S, net.lb, net.ub
    n = net.n; allr = set(range(n))
    irrev = {i for i in range(n) if lb[i] >= 0}
    sign = {i: 1.0 for i in range(n)}
    C = set(core_idx); J = C & irrev; P = allr - C; A = set()
    if J:
        v = _lp7(S, lb, ub, sorted(J), sign, eps); A = set(np.where(np.abs(v) > 0.99 * eps)[0])
    J = C - A; flipped = False; singleton = False; dropped = []
    while J:
        P = P - A
        Ji = [min(J)] if singleton else sorted(J)
        v = _lp7(S, lb, ub, Ji, sign, eps)
        K = [j for j in Ji if sign[j] * v[j] >= 0.99 * eps]
        supp = set()
        if K:
            v2 = _lp10(S, lb, ub, K, sorted(P), sign, eps)
            supp = set(np.where(np.abs(v2) > 0.99 * eps)[0])
        A |= supp
        if J & A:
            J = J - A; flipped = False
        else:
            JiRev = [j for j in Ji if j not in irrev]
            if flipped or not JiRev:
                if singleton:
                    dropped += Ji; J = J - set(Ji); flipped = False
                else:
                    flipped = False; singleton = True
            else:
                for j in JiRev: sign[j] = -sign[j]
                flipped = True
    return A, dropped

def gimme(net, score, threshold, frac=0.9):
    """GIMME (Becker & Palsson 2008): min sum (thr - s_j)|v_j| over low-score reactions, objective >= frac*max."""
    mu = net.fba(); n, m = net.n, net.m
    w = np.where(np.isnan(score) | (score >= threshold), 0.0, threshold - np.nan_to_num(score))
    I = sp.identity(n, format='csc')
    A = sp.vstack([sp.hstack([net.S, sp.csc_matrix((m, n))]), sp.hstack([-I, I]), sp.hstack([I, I])]).tocsc()
    lb = np.r_[net.lb, np.zeros(n)]; ub = np.r_[net.ub, np.full(n, INF)]; lb[net.obj] = frac * mu
    h = _highs(); _pass(h, A, np.r_[np.zeros(n), w], lb, ub, np.r_[np.zeros(m), np.zeros(2 * n)], np.r_[np.zeros(m), np.full(2 * n, INF)])
    h.run()
    v = np.array(h.getSolution().col_value[:n])
    return set(np.where(np.abs(v) > 1e-7)[0])

def can_carry(net, eps=1e-5, verbose=False, probe_patience=5):
    """Indices of reactions that can carry non-zero flux under the current bounds (FASTCC-style batched LP7):
    repeatedly maximize the number of still-unknown reactions carrying >= eps, forward then backward,
    then test the stubborn remainder one by one."""
    S, lb, ub = net.S, net.lb, net.ub; n = net.n
    unknown = set(np.where((ub > 0) | (lb < 0))[0]); carry = set()
    for direction in (1.0, -1.0):
        while True:
            J = sorted(i for i in unknown if (ub[i] > 0 if direction > 0 else lb[i] < 0))
            if not J: break
            sign = {i: direction for i in range(n)}
            v = _lp7(S, lb, ub, J, sign, eps, relaxed=True)
            new = set(np.where(np.abs(v) >= 0.99 * eps)[0]) & unknown
            if verbose: print('dir', direction, 'batch', len(J), 'new', len(new), flush=True)
            if not new: break
            carry |= new; unknown -= new
    # random-objective probing: vertices of the flux polytope expose remaining flux-capable reactions
    rng = np.random.default_rng(0); misses = 0
    h = _highs(); _pass(h, S, np.zeros(n), lb, ub, np.zeros(net.m), np.zeros(net.m), sense_max=True)
    while unknown and misses < probe_patience:
        c = np.zeros(n); U = np.array(sorted(unknown)); c[U] = rng.choice([-1.0, 1.0], len(U)) * rng.random(len(U))
        h.changeColsCost(n, np.arange(n, dtype=np.int32), c); h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: misses += 1; continue
        v = np.array(h.getSolution().col_value)
        new = set(np.where(np.abs(v) >= 0.99 * eps)[0]) & unknown
        if verbose: print('probe new', len(new), flush=True)
        if new: carry |= new; unknown -= new; misses = 0
        else: misses += 1
    return carry

def gapfill_support(net, have, target, min_flux, extra_lb=None):
    """Smallest-flux set of reactions outside `have` needed for reaction `target` >= min_flux
    (min sum |v_j| over j not in have, reactions in have are free). Support threshold 1e-9."""
    n, m = net.n, net.m
    lb = net.lb.copy(); ub = net.ub.copy(); lb[target] = max(lb[target], min_flux)
    if extra_lb:
        for i, v in extra_lb.items(): lb[i] = max(lb[i], v)
    P = np.array(sorted(set(range(n)) - set(have)))
    p = len(P)
    E = sp.csc_matrix((np.ones(p), (np.arange(p), P)), shape=(p, n)); Ip = sp.identity(p, format='csc')
    A = sp.vstack([sp.hstack([net.S, sp.csc_matrix((m, p))]), sp.hstack([-E, Ip]), sp.hstack([E, Ip])]).tocsc()
    h = _highs(); _pass(h, A, np.r_[np.zeros(n), np.ones(p)], np.r_[lb, np.zeros(p)], np.r_[ub, np.full(p, INF)],
                        np.r_[np.zeros(m), np.zeros(2 * p)], np.r_[np.zeros(m), np.full(2 * p, INF)])
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return set()
    v = np.array(h.getSolution().col_value[:n])
    return set(np.where(np.abs(v) > 1e-9)[0])
