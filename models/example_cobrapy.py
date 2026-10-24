"""Load FluxSprout models with COBRApy and reproduce knockout results from the screen.

    pip install cobra
    python models/example_cobrapy.py
"""
import os
import cobra
from cobra.flux_analysis import single_gene_deletion

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))   # run from any folder

# Reference network: embryonic medium, NGAM >= 1, objective MAR13082 (generic biomass)
ref = cobra.io.read_sbml_model('models/FluxSprout_reference.xml.gz')      # or cobra.io.load_json_model('models/FluxSprout_reference.json')
print(ref.id, len(ref.reactions), 'reactions', len(ref.genes), 'genes, growth', round(ref.slim_optimize(), 4))

# Context models are gzipped SBML; COBRApy reads .xml.gz directly
for member, group in [('fastcore_default', 'ISV_tip'), ('fastcore_default', 'aDA'), ('gimme_default', 'ISV_tip'), ('gimme_default', 'aDA')]:
    m = cobra.io.read_sbml_model(f'models/context/{member}/{group}.xml.gz')
    wt = m.slim_optimize()
    genes = [g for g in ('lcat', 'idh2', 'sdha') if g in m.genes]
    ko = single_gene_deletion(m, gene_list=genes)
    ratio = {list(i)[0]: round(max(v / wt, 0.0), 3) for i, v in zip(ko.ids, ko.growth)}
    print(f'{member:17s} {group:8s} growth {wt:.4f}  knockout growth / wild type: {ratio}')

# Changing the medium: exchange reactions carry the uptake bounds (negative lower bound = uptake)
m = cobra.io.read_sbml_model('models/context/fastcore_default/ISV_tip.xml.gz')
with m:
    m.reactions.get_by_id('MAR09048').lower_bound = -1     # O2 uptake 1 mmol gDW-1 h-1 (hypoxia; baseline is 20)
    print('ISV_tip growth at low O2:', round(m.slim_optimize(), 4))
