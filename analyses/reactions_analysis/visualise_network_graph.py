import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from rdkit import Chem
from rdkit.Chem import inchi

# Load data
df_results = pd.read_csv('/Users/anasves/work/glycosphingotool/results_sphing-4-enine_hexadecanoate/sphingomapkey with reaction structures.tsv', sep='\t')
compound_table = pd.read_csv('/Users/anasves/work/glycosphingotool/analyses/chebi_lipidmaps_match/compounds_with_chebis.tsv', sep='\t')
compound_table.dropna(subset=['InChIKey'], inplace=True)

def smiles_to_inchikey(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return inchi.MolToInchiKey(mol) if mol else None
    except:
        return None

def is_ceramide_like(smiles):
    return 'CCCC' in smiles  # crude heuristic

def split_reaction_smiles(rxn_smiles):
    try:
        reactants, products = rxn_smiles.split(">>")
        return reactants.split("."), products.split(".")
    except ValueError:
        return [], []

def build_graph_from_rxnsmiles(df):
    G = nx.Graph()
    for _, row in df.iterrows():
        rxn = row['rxnSMILES'].replace('CHO', 'CO').replace('CH2O', 'CO')
        reactants, products = split_reaction_smiles(rxn)
        reactants = [r for r in reactants if is_ceramide_like(r)]
        products  = [p for p in products if is_ceramide_like(p)]
        reactant_keys = {smiles_to_inchikey(smi) for smi in reactants}
        product_keys  = {smiles_to_inchikey(smi) for smi in products}
        for rk in reactant_keys:
            for pk in product_keys:
                if rk and pk and rk != pk:
                    G.add_edge(rk, pk)
    return G

def filter_compounds(mode):
    df = compound_table.copy()
    if mode == 'only non-hypothetical':
        df.dropna(subset=['SMID'], inplace=True)
        df = df[df['See Other'] != 'hypothetical']
    elif mode == 'with SMID':
        df.dropna(subset=['SMID'], inplace=True)
    return df

def draw_panel(ax, G_full, compound_df, title, category_colors, default_color):
    inchikeys = compound_df['InChIKey'].tolist()
    G = G_full.subgraph(inchikeys)
    for _, row in compound_df.iterrows():
        ikey = row['InChIKey']
        category = row.get('category_assigned')
        if ikey in G.nodes:
            G.nodes[ikey]['category'] = category

    #pos = nx.forceatlas2_layout(G, max_iter=100, scaling_ratio=2.0)
    pos = nx.spring_layout(G, seed=42)
    node_colors = [
        category_colors.get(G.nodes[n].get('category'), default_color)
        for n in G.nodes
    ]
    nx.draw(
        G, pos, ax=ax, with_labels=False,
        node_color=node_colors, edge_color='lightgray',
        node_size=100
    )
    ax.set_title(title)
    ax.axis('off')

# Define category colors based on all data
all_categories = compound_table['category_assigned'].dropna().unique()
category_colors = {cat: plt.cm.tab20(i / len(all_categories)) for i, cat in enumerate(all_categories)}
default_color = 'black'

# Build full reaction graph once
G_full = build_graph_from_rxnsmiles(df_results)

# Plot panels
# Create figure with constrained layout to handle legends better
# Set global font scale (~3x default)
# Adjust global font settings (about 3x larger)
plt.rcParams.update({
    'font.size': 36,
    'axes.titlesize': 36,
    'legend.fontsize': 28,
    'legend.title_fontsize': 30,
})

# Create figure and subplots with extra spacing between them
fig, axes = plt.subplots(1, 3, figsize=(26, 10), gridspec_kw={'wspace': 0.3})
modes_titles = [
    ('only non-hypothetical', 'Non-hypothetical'),
    ('with SMID', 'With SMID'),
    ('all_compounds', 'All Compounds')
]
panel_labels = ['a.', 'b.', 'c.']

# Define vertical offset for label/title
label_y = 1.12
title_y = 1.05

for i, (ax, (mode, title)) in enumerate(zip(axes, modes_titles)):
    df_filtered = filter_compounds(mode)
    draw_panel(ax, G_full, df_filtered, "", category_colors, default_color)  # title will be added manually

    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Add Nature-style panel label
    ax.text(
        0, label_y, panel_labels[i],
        transform=ax.transAxes,
        fontsize=38, fontweight='bold',
        va='top', ha='left'
    )

    # Add subplot title above
    ax.text(
        0.5, title_y, title,
        transform=ax.transAxes,
        fontsize=36, ha='center', va='top'
    )

# Add vertical lines between panels (starting below title level)
for i in [1, 2]:
    x = (axes[i - 1].get_position().x1 + axes[i].get_position().x0) / 2
    line = Line2D(
        [x, x], [0.1, 0.88],  # start below the titles
        transform=fig.transFigure,
        color='gray', linewidth=2.5
    )
    fig.add_artist(line)

# Prepare legend
legend_elements = [
    mlines.Line2D([0], [0], marker='o', color='w',
                  label=cat, markerfacecolor=color, markersize=12)
    for cat, color in category_colors.items()
]
legend_elements.append(
    mlines.Line2D([0], [0], marker='o', color='w',
                  label='Uncategorized', markerfacecolor=default_color, markersize=12)
)

# Add legend below the panels, tighter and slightly smaller
fig.legend(
    handles=legend_elements,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.1),
    ncol=4,
    title='Categories'
)

# Tight layout leaving space above for titles and below for legend
plt.savefig("figure.png", dpi=600, bbox_inches='tight')
plt.show()