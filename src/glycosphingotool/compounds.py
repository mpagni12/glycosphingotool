import os
import pandas as pd
import re
from rdkit import Chem
import importlib.resources
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

"""
Based on what LIPIDMAPs has these are the sequences for the different series: 
Ganglio-series: GalNAcβ1-4Galβ1-4Glcβ-Cer 
Globo-series: Galα1-4Galβ1-4Glcβ-Cer 
Isoglobo-series: Galα1-3Galβ1-4Glcβ-Cer 
Lacto-series:  Galβ1-3GlcNAcβ1-3Galβ1-4Glcβ-Cer and 
Neolacto-series:  Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer
Mollu-series: Manα1-3Manβ1-4Glcβ-Cer 
Arthro-series:GlcNAcβ1-4Manβ1-4Glcβ-Cer
"""

dict_categories = {
"Galbeta1-4GlcNAcbeta1-3Galbeta1-4GlcCer":"Neolacto",
"GlcNAcbeta1-3Galbeta1-4GlcCer":"Lacto",
"Galalpha1-3Galbeta1-4GlcCer":"Isoglobo", 
"GalNAcbeta1-4Galbeta1-4GlcCer":"Ganglio",
"GalNacbeta1-4Galbeta1-4GlcCer":"Ganglio",
"Galalpha1-4Galbeta1-4GlcCer":"Globo",
"Manalpha1-3Manbeta1-4GlcCer":"Mollu",
"GlcNAcbeta1-4Manbeta1-4GlcCer":"Arthro",
"NeuAcalpha2-3Galbeta1-4GlcCer":"Ganglio"}

def process_participant_structures(output_folder, nacyl, sphingoid,
                                   input_xls=importlib.resources.files("glycosphingotool.assets").joinpath('SphingomapkeyV1.4.xls')):
    
    df = pd.read_csv(os.path.join(output_folder, 'sphingomapkey with reaction structures.tsv'), sep='\t')

    participant_structures = pd.DataFrame(columns=[
        'InChIKey', 'SMILES', 'InChI', 'name'
    ])

    for _, row in df.iterrows():
        # split reaction SMILES into participants
        participant_smiles = re.split(r'\.\s*|>>\s*', row['rxnSMILES'].replace('CHO', 'CO').replace('CH2O', 'CO'))

        # try to get names from glyco_nomenclature if available
        if 'glyco_nomenclature' in df.columns:
            participant_names = re.split(r'\+\s*|=>\s*', row['glyco_nomenclature'])

        else:
            participant_names = [''] * len(participant_smiles)

        # zip together names and SMILES
        for smiles, name in zip(participant_smiles, participant_names):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                participant_inchi = Chem.MolToInchi(mol)
                participant_inchikey = Chem.MolToInchiKey(mol)
                participant_structures = pd.concat([
                    participant_structures,
                    pd.DataFrame([{
                        'InChIKey': participant_inchikey,
                        'SMILES': smiles,
                        'InChI': participant_inchi,
                        'name': name.replace('α', 'alpha').replace('β', 'beta').strip()
                    }])
                ], ignore_index=True)


    participant_structures['is_cer'] = participant_structures['name'].apply(lambda x: x.endswith('Cer'))
    participant_structures = participant_structures[participant_structures['is_cer']==True]
    participant_structures.drop(columns=['is_cer'], inplace=True)

     # Merge with the input data
    df_sphingomapkey = pd.read_excel(input_xls, skiprows=2)

    participant_structures = participant_structures.merge(df_sphingomapkey, left_on='name', right_on='Formula', how='left')

    # Split into two parts
    with_value = participant_structures[participant_structures['Formula'].notna()]
    without_value = participant_structures[participant_structures['Formula'].isna()]
    participant_structures.drop(columns=['Formula'], inplace=True)

    # Drop duplicates in rows where 'value' is NaN
    without_value = without_value.drop_duplicates(subset=['InChIKey'])
    with_value = with_value.drop_duplicates(subset=['SMID', 'InChIKey'])
    without_value = without_value[~without_value['InChIKey'].isin(with_value['InChIKey'])]

    # Combine both
    participant_structures = pd.concat([with_value, without_value], ignore_index=True)

    participant_structures[['with sphingoid base', 'generic glycosphingolipid']] = participant_structures['SMILES'].apply(
    lambda s: pd.Series(generate_R_versions(s, nacyl, sphingoid))
)
    participant_structures.rename(columns={'name':'generated_glyconomenclature'}, inplace=True)
    participant_structures.drop(columns=['Formula'], inplace=True)
    
    # Assign category
    participant_structures['category_assigned'] = participant_structures['generated_glyconomenclature'].apply(assign_category)

    participant_structures.to_csv(os.path.join(output_folder, 'compounds.tsv'), sep='\t', index=False)

def generate_R_versions(smiles, nacyl, sphingoid):
    """
    "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCCC"
    "CCCCCCCCCCCCCCC"                                                       "[C@H](O)CCCCCCCCCCCCCCC"
    "*C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCCC"
    "*C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)*"
     
    Generates two R-derivatized SMILES by replacing:
    - nacyl group with *
    - both nacyl and sphingoid with *
    Tries to handle the two possible orders of amide linkage.
    """

    reverse_sphingoids = {"[C@H](O)CCCCCCCCCCCCCCC":"CCCCCCCCCCCCCCC[C@@H](O)",
                        "[C@H](O)/C=C/CCCCCCCCCCCCC":"CCCCCCCCCCCCC/C=C/[C@@H](O)",
                        "[C@H](O)[C@H](O)CCCCCCCCCCCCCC":"",
                        "[C@H](O)CCCCCCCCCCC(C)CC":"",
                        "C(=O)CCCCCCCCCCCCCCC":"",
                        "[C@H](O)CCCCCCCCCCCCC":"",
                        "C(=O)CCCCCCCCCCCCC":"",
                        "C(=O)CCCCCCCCCCCCCCCCC":"",
                        "C(=O)CCCCCCCCCCC":"",
                        "[C@H](O)/C=C/CCCCCCCCCCCC":"",
                        "[C@H](O)/C=C/CCCCCCCCC(C)CC":"",
                        "C(=O)CCCCCCCCCCC(C)CC":"",
                        "[C@H](O)CCCCCCCCCCCCCCCCC":"",
                        "[C@H](O)[C@H](O)CCCCCCCCCCCCCCCC":""}
    
    smiles_nacyl_removed = smiles
    
    # 1. try nacyl in front
    pattern1 = nacyl + "C(=O)N"
    if pattern1 in smiles:
        smiles_nacyl_removed = smiles.replace(pattern1, "*C(=O)N")
    else:
        # 2. try nacyl in back (amide reversed)
        pattern2 = "NC(=O)" + nacyl
        if pattern2 in smiles:
            smiles_nacyl_removed = smiles.replace(pattern2, "NC(=O)*")

    smiles_both_removed = smiles_nacyl_removed
    
    # sphingoid
    
    if sphingoid in smiles_nacyl_removed:
        smiles_both_removed = smiles_nacyl_removed.replace(sphingoid, "[C@H](O)*")
        return smiles_nacyl_removed, smiles_both_removed
    
    if reverse_sphingoids[sphingoid] in smiles_nacyl_removed:
        smiles_both_removed = smiles_nacyl_removed.replace(reverse_sphingoids[sphingoid], "*[C@@H](O)")

    return smiles_nacyl_removed, smiles_both_removed


def assign_category(nomenclature):
    for ending in dict_categories:
        if nomenclature.replace(')','').replace(']','').replace('}','').endswith(ending):
            return dict_categories[ending]
        cleaned = remove_branches_c(nomenclature).replace(')', '')
        if cleaned.endswith(ending):
            return dict_categories[ending]
        cleaned = remove_branches_square(nomenclature).replace(']', '')
        if cleaned.endswith(ending):
            return dict_categories[ending]

def remove_branches_c(s):
    # Remove all parentheses and their contents
    # Example: NeuAcalpha2-3(Galbeta1-4)Galbeta1-4GlcCer -> NeuAcalpha2-3Galbeta1-4GlcCer
    while True:
        new_s = re.sub(r'\([^()]*\)', '', s)
        if new_s == s:
            break
        s = new_s
    return s

def remove_branches_square(s):
    # This regex removes square brackets and everything inside, including nested parentheses
    # We remove all occurrences iteratively until none left
    pattern = r'\[[^\[\]]*\]'  # matches square brackets without nested square brackets

    while True:
        new_s = re.sub(pattern, '', s)
        if new_s == s:
            break
        s = new_s
    return s
