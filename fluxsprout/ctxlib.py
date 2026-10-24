"""GPR evaluation and local-T2 activity calls for FluxSprout.

GPR rules are evaluated with AND = min and OR = max; local-T2 thresholds follow Richelle et al. (2019).
FASTCORE, GIMME and E-Flux are implemented on the HiGHS engine in hlp.py and pipeline.py.
"""
import re, numpy as np, pandas as pd

# ---------------------------------------------------------------- GPR evaluation
_tok = re.compile(r'\(|\)|\band\b|\bor\b|[^\s()]+')

def eval_gpr(rule, val, missing=np.nan):
    """Evaluate a GPR string with AND->min, OR->max. Genes absent from `val` take `missing`
    and are ignored inside AND/OR (so an undetected gene neither switches a reaction off nor on)."""
    if not rule or not rule.strip():
        return np.nan
    toks = _tok.findall(rule)
    pos = 0
    def expr():
        nonlocal pos
        left = term(); vals = [left]
        while pos < len(toks) and toks[pos] == 'or':
            pos += 1; vals.append(term())
        v = [x for x in vals if not np.isnan(x)]
        return max(v) if v else np.nan
    def term():
        nonlocal pos
        left = factor(); vals = [left]
        while pos < len(toks) and toks[pos] == 'and':
            pos += 1; vals.append(factor())
        v = [x for x in vals if not np.isnan(x)]
        return min(v) if v else np.nan
    def factor():
        nonlocal pos
        t = toks[pos]; pos += 1
        if t == '(':
            v = expr(); pos += 1; return v
        return val.get(t, missing)
    return expr()

def reaction_scores(model, gene_val):
    d = gene_val.to_dict() if hasattr(gene_val, 'to_dict') else gene_val
    return pd.Series({r.id: eval_gpr(r.gene_reaction_rule, d) for r in model.reactions})

# ---------------------------------------------------------------- local T2 activity
def localT2(expr_df, lower_q=25, upper_q=75):
    """expr_df: genes x samples (log1p mean CP10k). Returns boolean genes x samples activity.
    Global thresholds from all values; between them a gene is active if above its own mean across samples."""
    v = expr_df.values.ravel()
    lo, hi = np.percentile(v, lower_q), np.percentile(v, upper_q)
    gmean = expr_df.mean(axis=1)
    act = (expr_df >= hi) | ((expr_df >= lo) & expr_df.ge(gmean, axis=0))
    act &= expr_df > 0
    return act, (lo, hi)
