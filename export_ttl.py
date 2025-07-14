import re
import os
import pandas as pd
from rdkit import Chem
from rdflib import Graph, Namespace, RDF, RDFS, URIRef, Literal
import pandas as pd
import os
from datetime import datetime
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

df_chebi = pd.read_csv(os.path.join('analyses','chebi_lipidmaps_match','df_chebi.tsv'), sep='\t')

sphinganine_to_position_code = {"sphing-4-enine":"d18:1", "sphinganine":"d18:0"}
fa_to_position_code = { 
"_4Z_7Z_10Z_13Z_16Z_19Z_-docosahexaenoate":"4Z,7Z,10Z,13Z,16Z,19Z-22:6",
"_6Z_9Z_12Z_15Z_-octadecatetraenoate":"6Z,9Z,12Z,15Z-18:4",
"_15Z_-tetracosenoate":"15Z-24:1",
"_17Z_-hexacosenoate":"17Z-26:1",
"2-hydroxyarachidate":"20:0(2OH[R])",
"2-hydroxybehenate":"22:0(2OH[R])",
"2-hydroxyhexacosanoate":"26:0(2OH[R])",
"2-hydroxyhexadecanoate":"16:0(2OH[R])",
"2-hydroxynervonate":"15Z-24:1(2OH[R])",
"2-hydroxyoctadecanoate":"18:0(2OH[R])",
"2-hydroxytetracosanoate":"24:0(2OH[R])",
"all-cis-5_8_11_14_17-icosapentaenoate":"5Z,8Z,11Z,14Z,17Z-20:5",
"all-cis-icosa-8_11_14-trienoate":"8Z,11Z,14Z-20:3",
"arachidonate":"5Z,8Z,11Z,14Z-20:4",
"behenate":"22:0",
"cerotate":"26:0",
"hexadecanoate":"16:0",
"icosanoate":"20:0",
"linoleate":"9Z,12Z-18:2",
"octadecanoate":"18:0",
"oleate":"9Z-18:1",
"palmitoleate":"9Z-16:1",
"tetracosanoate":"24:0"
}

# RDF setup
g = Graph()
RH = Namespace("http://rdf.rhea-db.org/")
ANASVES_str = "http://rdf.example.org/anasves/"
ANASVES = Namespace(ANASVES_str)
GLYCOSPHINGO_str = "http://rdf.example.org/glycopshingolipids/"
GLYCOSPHINGO = Namespace(GLYCOSPHINGO_str)
GO = Namespace("http://purl.obolibrary.org/obo/GO_")
CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")

g.bind("rdf", RDF)
g.bind("rdfs", RDFS)
g.bind("rh", RH)
g.bind("anasves", ANASVES)

MNetIRI = ANASVES["glycopshingolipids"]
g.add((MNetIRI, RDFS.label, Literal("glycosphingolipids")))

compounds_total = [] # collect compounds for the .tsv for MetaNetX

def export_glyco_ttl(input_folder):
    

    # Load data
    df = pd.read_csv(os.path.join(input_folder, 'sphingomapkey with reaction structures.tsv'), sep="\t")
    df.dropna(subset=["RInChIKey"], inplace=True)
    df.drop_duplicates(subset=["RInChIKey"], inplace=True)

    def parse_equation(row):
        
        reactant_smiles = re.split(r'\.\s*|>>\s*', row['rxnSMILES'].split('>>')[0].replace('CHO', 'CO').replace('CH2O', 'CO'))
        product_smiles = re.split(r'\.\s*|>>\s*', row['rxnSMILES'].split('>>')[1].replace('CHO', 'CO').replace('CH2O', 'CO'))
        
        # try to get names from glyco_nomenclature if available
        reactant_names = re.split(r'\+\s*|=>\s*', row['glyco_nomenclature'].split('=>')[0])
        product_names = re.split(r'\+\s*|=>\s*', row['glyco_nomenclature'].split('=>')[1])

        # zip together names and SMILES
        def zip_side_with_inchikey(smiles_list, names_list):
            res = []
            for smiles, name in zip(smiles_list, names_list):
                mol = Chem.MolFromSmiles(smiles)
                inchikey = 'NA-NA-NA'
                if mol:
                    inchikey = Chem.MolToInchiKey(mol)
                res.append((smiles, name, inchikey))
            return res

        return zip_side_with_inchikey(reactant_smiles, reactant_names), zip_side_with_inchikey(product_smiles, product_names)

    g.add((RH[f"contains1"], RH.coefficient, Literal(1)))
    g.add((RH[f"contains1"], RDFS.subPropertyOf, RH.contains))

    for _, row in df.iterrows():
        accession = row["RInChIKey"].replace("Web-RInChIKey=", "")
        reaction_uri_str = f"{ANASVES}reaction_{accession.replace('-', '_')}"
        reaction_uri = URIRef(reaction_uri_str)
        left_uri = URIRef(f"{reaction_uri_str}_L")
        right_uri = URIRef(f"{reaction_uri_str}_R")
        reaction_ref = URIRef(reaction_uri)

        # Basic reaction info
        g.add((MNetIRI, ANASVES.includes, reaction_ref))
        g.add((reaction_ref, RH.accession, Literal(accession)))
        g.add((reaction_uri, RH.side, left_uri))
        g.add((reaction_uri, RH.side, right_uri))
        g.add((left_uri, RH.curatedOrder, Literal(1)))
        g.add((right_uri, RH.curatedOrder, Literal(2)))

        # Parse and add participants
        left_participants, right_participants = parse_equation(row)

        def add_compounds(side_uri, compounds, reaction_uri):
            for smiles, name, inchikey in compounds:
                if inchikey == "NA-NA-NA":
                    print('No INCHIKEY: can can balance error', smiles, name)
                    continue
                comp_id = inchikey.replace('-','_')
                part_uri = URIRef(f"{reaction_uri}_compound_{comp_id}")
                g.add((part_uri, RH.location, GO["0005575"]))
                comp_uri = URIRef(f"{GLYCOSPHINGO_str}Compound_{comp_id}")
                g.add((side_uri, RH[f"contains1"], part_uri))
                g.add((part_uri, RH.compound, comp_uri))

                chebi = df_chebi[df_chebi['inchikey_nostar']==inchikey]['COMPOUND_ID']
                if len(chebi)>0:
                    chebi = chebi.iloc[0]
                    chem_iri = CHEBI[str(chebi)]
                    g.add((comp_uri, RH.accession, Literal(f"CHEBI:{chebi}")))
                    g.add((comp_uri, RH.chebi, chem_iri))
                else:
                    sph = input_folder.replace('results_','').split('_')[0]
                    fa = input_folder.replace('results_sphing-4-enine_','').replace('results_sphinganine_','')
                    compounds_total.append({'#id':comp_id, "name":f"{name}({sphinganine_to_position_code[sph]}/{fa_to_position_code[fa]})", "SMILES": smiles})
                    #chem_iri = GLYCOSPHINGO[comp_id]
                    # g.add((comp_uri, RH.accession, Literal(comp_id)))
                    # g.add((comp_uri, RH.chebi, chem_iri))

        add_compounds(left_uri, left_participants, reaction_uri_str)
        add_compounds(right_uri, right_participants, reaction_uri_str)

output_dir='.'
# Output setup
if output_dir is None:
    output_dir = os.getcwd()
else:
    os.makedirs(output_dir, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_file = os.path.join(output_dir, f"glyco_reactions_{timestamp}.ttl")
all_results = [i for i in os.listdir('.') if i.startswith('results_')]
for res_folder in all_results:
    export_glyco_ttl(res_folder)

print(f"Saving RDF with {len(g)} triples to: {output_file}")
g.serialize(destination=output_file, format="ttl")
df = pd.DataFrame(compounds_total)
df.drop_duplicates(subset=['#id'], inplace=True)
df['ori_name']=''
df['synonym']=''
df['formula']=''
df['charge']=''
df['mass']=''
df['xref']=''
df['Mol_file/SDF']=''
df['InChI']=''
df['InChIKey']=''
df['alt_id']=''
df.to_csv('compounds_glycosphingolipids_for_metanetx.tsv', sep='\t', index=False, 
          columns=['#id','ori_name','name','synonym','formula','charge',
                   'mass','xref','Mol_file/SDF','SMILES','InChI','InChIKey','alt_id'])
