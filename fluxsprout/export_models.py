"""Write the FluxSprout models as standalone COBRApy-readable files.

Run from the repository root after the notebook (or with results/ and cache/ present):
    python fluxsprout/export_models.py

Output (models/):
    FluxSprout_reference.xml.gz / .json      reference network, embryonic medium, NGAM = 1, objective MAR13082
    context/<member>/<group>.xml.gz          108 context models (4 ensemble members x 27 groups), SBML L3 FBC2
    model_index.csv                          one row per file with size and the growth COBRApy finds on reading it back
"""
import sys, os, json, pickle, gzip, shutil
import pandas as pd
import cobra

sys.path.insert(0, os.path.dirname(__file__))
MEMBERS = ['fastcore_loose', 'fastcore_default', 'fastcore_strict', 'gimme_default']
OUT = 'models'
cobra.Configuration().solver = 'glpk'


def load_reference():
    ref = pickle.load(open('cache/ref_medium_consistent.pkl', 'rb'))
    ref.id = 'FluxSprout_reference'
    ref.name = 'FluxSprout reference network (Zebrafish-GEM, flux-consistent under the embryonic medium)'
    return ref


def write_gz_sbml(model, path):
    tmp = path[:-3]
    cobra.io.write_sbml_model(model, tmp)
    with open(tmp, 'rb') as f, gzip.open(path, 'wb', compresslevel=9) as g:
        shutil.copyfileobj(f, g)
    os.remove(tmp)


def safe(name):
    return name.replace('/', '-').replace(' ', '_').replace('.', '')


def export_member(member):
    ref = load_reference()
    ctx = json.load(open(f'results/ctx_{member}.json'))
    wt = pd.read_csv(f'results/ko_{member}_baseline.csv', index_col=0).loc['__WT__']
    os.makedirs(f'{OUT}/context/{member}', exist_ok=True)
    rows = []
    for group, keep in ctx.items():
        path = f'{OUT}/context/{member}/{safe(group)}.xml.gz'
        if not os.path.exists(path):                    # resumable: existing files are only re-checked
            m = ref.copy()
            keep = set(keep)
            m.remove_reactions([r for r in m.reactions if r.id not in keep], remove_orphans=True)
            m.id = f'FluxSprout_{member}_{safe(group)}'
            m.name = f'FluxSprout context model, {member.replace("_", " ")}, group {group}'
            write_gz_sbml(m, path)
        back = cobra.io.read_sbml_model(path)          # read back exactly as a user would
        mu = back.slim_optimize()
        rows.append(dict(file=path, member=member, group=group, reactions=len(back.reactions), metabolites=len(back.metabolites),
                         genes=len(back.genes), growth_cobrapy=mu, growth_screen=float(wt[group]),
                         size_MB=round(os.path.getsize(path) / 1e6, 2)))
        print(member, group, len(back.reactions), round(mu, 5), round(float(wt[group]), 5), flush=True)
    return rows


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    ref = load_reference()
    write_gz_sbml(ref, f'{OUT}/FluxSprout_reference.xml.gz')
    cobra.io.save_json_model(ref, f'{OUT}/FluxSprout_reference.json')
    rows = []
    for ext, reader in [('xml.gz', cobra.io.read_sbml_model), ('json', cobra.io.load_json_model)]:
        p = f'{OUT}/FluxSprout_reference.{ext}'; b = reader(p)
        rows.append(dict(file=p, member='reference', group='all', reactions=len(b.reactions), metabolites=len(b.metabolites), genes=len(b.genes),
                         growth_cobrapy=b.slim_optimize(), growth_screen=ref.slim_optimize(), size_MB=round(os.path.getsize(p) / 1e6, 2)))
    for member in MEMBERS:                              # one process: GLPK can hang under multiprocessing
        rows += export_member(member)
    idx = pd.DataFrame(rows)
    idx['abs_diff'] = (idx.growth_cobrapy - idx.growth_screen).abs()
    idx.to_csv(f'{OUT}/model_index.csv', index=False)
    print('max |growth difference| between COBRApy and the screen:', idx.abs_diff.max())
