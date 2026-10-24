"""FluxSprout helper library: medium, objective and reaction sets."""
import pickle, numpy as np, pandas as pd, cobra
from cobra.flux_analysis import pfba
cobra.Configuration().solver='glpk'
EXCH={'glucose':'MAR09034','O2':'MAR09048','alanine':'MAR09061','arginine':'MAR09066','asparagine':'MAR09062','aspartate':'MAR09070',
'cysteine':'MAR09065','glutamine':'MAR09063','glutamate':'MAR09071','glycine':'MAR09067','histidine':'MAR09038','isoleucine':'MAR09039',
'leucine':'MAR09040','lysine':'MAR09041','methionine':'MAR09042','phenylalanine':'MAR09043','proline':'MAR09068','serine':'MAR09069',
'threonine':'MAR09044','tryptophan':'MAR09045','tyrosine':'MAR09064','valine':'MAR09046','palmitate':'MAR00611','oleate':'MAR00650',
'linoleate':'MAR09035','linolenate':'MAR09036','DHA':'MAR00573','arachidonate':'MAR00568','cholesterol':'MAR09285','choline':'MAR09083',
'ethanolamine':'MAR09084','folate':'MAR09146','thiamin':'MAR09159','biotin':'MAR09109','riboflavin':'MAR09143','nicotinamide':'MAR09378',
'pyridoxine':'MAR09144','pantothenate':'MAR09145','inositol':'MAR09361','retinol':'MAR09147','lipoic acid':'MAR09167','ascorbate':'MAR09158',
'alpha-tocopherol':'MAR09151','Pi':'MAR09072','Fe2+':'MAR09076','sulfate':'MAR09074','H2O':'MAR09047','H+':'MAR09079','zinc':'MAR09080','gamma-tocopherol':'MAR09153','B12':'MAR09269','carnitine':'MAR09292'}
AA=['alanine','arginine','asparagine','aspartate','cysteine','glutamate','glycine','histidine','isoleucine','leucine','lysine','methionine',
    'phenylalanine','proline','serine','threonine','tryptophan','tyrosine','valine']
FA=['palmitate','oleate','linoleate','linolenate','DHA','arachidonate']
VIT=['folate','thiamin','biotin','riboflavin','nicotinamide','pyridoxine','pantothenate','inositol','retinol','lipoic acid','ascorbate','alpha-tocopherol','gamma-tocopherol','B12']
IONS=['Pi','Fe2+','sulfate','H2O','H+','zinc']
BIOMASS='MAR13082'   # generic human cell biomass (45 ATP GAM)
NGAM='MAR03964'
PPI_LEGACY=['MAR03975','MAR03977','MAR03979','MAR06912']
OXPHOS6=['MAR06921','MAR06918','MAR06914','MAR06916','MAR06911','MAR13081']
OXPHOS10=OXPHOS6+PPI_LEGACY

def baseline_medium(glc=10,o2=20,gln=10,aa=1,fa=1,vit=0.1,chol=1,extra=None):
    med={'glucose':glc,'O2':o2,'glutamine':gln,'cholesterol':chol,'choline':1,'ethanolamine':1}
    med.update({a:aa for a in AA}); med.update({f:fa for f in FA}); med.update({v:vit for v in VIT}); med.update({i:1000 for i in IONS})
    if extra: med.update(extra)
    return med

def load_model(path='cache/zgem.pkl'):
    m=pickle.load(open(path,'rb'))
    for r in m.reactions:   # open biomass bounds
        if r.id in ('MAR13082','MAR04413','MAR00021'): r.bounds=(0,1000)
    return m

def apply_medium(m,med,unbounded=False):
    """Close all uptakes, open all secretion; set uptake bounds from `med` (names->value)."""
    for r in m.exchanges:
        r.lower_bound=-1000 if unbounded else 0; r.upper_bound=1000
    if not unbounded:
        for k,v in med.items():
            if EXCH[k] in m.reactions: m.reactions.get_by_id(EXCH[k]).lower_bound=-float(v)
    return m

def set_objective(m,rid=BIOMASS):
    m.objective=rid; return m
