import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import os
from rdkit import RDLogger
import requests
from bs4 import BeautifulSoup
import argparse
import sys

RDLogger.DisableLog('rdApp.*')

tqdm.pandas()

##############################################################################################################
#             CHEBI            CHEBI            CHEBI            CHEBI            CHEBI            CHEBI
##############################################################################################################

def read_chebi():
    df_chebi = pd.read_csv('structures.csv')
    # ['ID', 'COMPOUND_ID', 'STRUCTURE', 'TYPE', 'DIMENSION', 'DEFAULT_STRUCTURE', 'AUTOGEN_STRUCTURE']

    df_chebi = df_chebi[df_chebi['TYPE']=='mol']
    print('Total number of molecules in ChEBI:', len(df_chebi['COMPOUND_ID'].unique()))
    df_chebi['SMILES'] = df_chebi['STRUCTURE'].progress_apply(MolToSmiles)
    df_chebi.dropna(subset=['SMILES'], inplace=True)
    df_chebi.drop(columns=['ID','STRUCTURE', 'TYPE','DIMENSION','DEFAULT_STRUCTURE','AUTOGEN_STRUCTURE'], inplace=True)
    # The structures 140708 and 140794 were corrected in ChEBI as. well, but new release with corrected structures include also structures submitted in this work
    df_chebi.loc[df_chebi['COMPOUND_ID'] == 140708, 'SMILES'] = '[1*][C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO[C@@H]5O[C@H](CO)[C@@H](O[C@@H]6O[C@H](CO)[C@H](O)[C@H](O)[C@H]6O)[C@H](O[C@@H]6O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]6O)[C@H]5NC(C)=O)[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O)[C@H]5O)[C@H]4NC(C)=O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)[2*]'
    df_chebi.loc[df_chebi['COMPOUND_ID'] == 140794, 'SMILES'] = '[1*][C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O[C@@]6(C(=O)[O-])O[C@@H]([C@@H]([C@@H](CO)O[C@@]7(C(=O)[O-])O[C@@H]([C@@H]([C@@H](CO)O)O)[C@H](NC(C)=O)[C@@H](O)C7)O)[C@H](NC(C)=O)[C@@H](O)C6)[C@H]5O)[C@H]4NC(C)=O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)[2*]'
    df_chebi['SMILES_nostar']=df_chebi['SMILES'].apply(lambda x: x.replace('[1*]','I').replace('[2*]','I').replace('*','I'))
    df_chebi['inchikey_nostar']=df_chebi['SMILES_nostar'].progress_apply(SmilesToInchikey)
    df_chebi.to_csv('df_chebi.tsv', sep='\t', index=False)
    return df_chebi

def SmilesToInchikey(smiles):
    try:
        return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
    except:
        pass

def MolToSmiles(mol):
    try:
        return Chem.MolToSmiles(Chem.MolFromMolBlock(mol))
    except:
        pass

def check_inchikey_in_chebi(smiles, df_chebi=pd.DataFrame()):
    res = df_chebi[df_chebi['inchikey_nostar']==smiles]
    return ', '.join(set([str(i) for i in res['COMPOUND_ID'].to_list()]))

def process_participants_structures_df(df_chebi, input_folder):

    # Example: Loop through files in the input folder
    file_path = os.path.join(input_folder, 'compounds.tsv')

    df_chebi.dropna(subset=['inchikey_nostar'], inplace=True)
    participants_structures_df = pd.read_csv(file_path, sep='\t')
    participants_structures_df['with sphingoid base_nostar']=participants_structures_df['with sphingoid base'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x.replace('*','I'))))
    participants_structures_df['generic glycosphingolipid_nostar']=participants_structures_df['generic glycosphingolipid'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x.replace('*','I'))))

    participants_structures_df['inchikey']=participants_structures_df['SMILES'].apply(lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x)))
    participants_structures_df['with sphingoid base inchikey_nostar']=participants_structures_df['with sphingoid base_nostar'].apply(lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x)))
    participants_structures_df['generic glycosphingolipid inchikey_nostar']=participants_structures_df['generic glycosphingolipid_nostar'].apply(lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x)))

    participants_structures_df['CHEBI_ID'] = participants_structures_df.progress_apply(lambda row: check_inchikey_in_chebi(row['inchikey'], df_chebi=df_chebi), axis=1)
    participants_structures_df['CHEBI_ID_R1'] = participants_structures_df.progress_apply(lambda row: check_inchikey_in_chebi(row['with sphingoid base inchikey_nostar'], df_chebi=df_chebi), axis=1)
    participants_structures_df['CHEBI_ID_R1_R2'] = participants_structures_df.progress_apply(lambda row: check_inchikey_in_chebi(row['generic glycosphingolipid inchikey_nostar'], df_chebi=df_chebi), axis=1)

    manual_mapping = pd.read_excel('Mapping to SphinGOMAP.xlsx')
    manual_mapping.dropna(subset=['CHEBI ID'], inplace=True)
    manual_mapping['chebiid']=manual_mapping['CHEBI ID'].apply(lambda x: int(x.split(':')[1]))
    participants_structures_df = participants_structures_df.merge(manual_mapping, how='outer', left_on = 'SMID', right_on='SphinGOMAP ID')
    participants_structures_df = participants_structures_df.merge(df_chebi, how='left', left_on='chebiid', right_on='COMPOUND_ID')
    participants_structures_df.drop_duplicates(subset=['generated_glyconomenclature'], inplace=True)
    participants_structures_df.to_csv('compounds_with_chebis.tsv', sep='\t', index=False)

    participants_structures_df.rename(columns={'CHEBI_ID_R1_R2':'Automatically found from structure ChEBI ID', 
                                               'chebiid': 'Manually mapped ChEBI ID',
                                               'generic glycosphingolipid':'generic glycosphingolipid SMILES',
                                               'Unnamed: 8': 'Lucila’s comment'}, inplace=True)
    participants_structures_df['match_automatic_and_manual_chebi'] = participants_structures_df.apply(lambda row: row['Manually mapped ChEBI ID'] in ([int(i) for i in str(row['Automatically found from structure ChEBI ID']).split(',') if i and i!='nan']) if row['Manually mapped ChEBI ID']>0 else '', axis=1)
    
    participants_structures_df.to_csv('compound_analysis.tsv', sep='\t', 
                                      columns=['SMID','SphinGOMAP ID',''
                                      'Automatically found from structure ChEBI ID',
                                      'Manually mapped ChEBI ID','Name',
                                      'generated_glyconomenclature',
                                      'Category','category_assigned','LBID','Common Name','References','Lucila’s comment',
                                      'generic glycosphingolipid SMILES', 'match_automatic_and_manual_chebi'], index=False)

def analyse_results():

    participants_structures_df = pd.read_csv('compounds_with_chebis.tsv', sep='\t')
    original_df = pd.read_excel(os.path.join('..','..','src','glycosphingotool','assets','SphingomapkeyV1.4.xls'), skiprows=2, usecols=['SMID', 'Formula'])
    original_df.dropna(subset=['SMID'], inplace=True)
    print('Unique SMID in original data', len(original_df['SMID'].unique()))
    print('Unique glyconomenclature formulas in original data', len(original_df['Formula'].unique()))

    participants_structures_df = participants_structures_df.merge(original_df, on='SMID', how='outer')

    print('Total sphingolipids (nomenclature) in df:', len(participants_structures_df['generated_glyconomenclature'].unique()))
    print('Total sphingolipids (InChIKey) in df:', len(participants_structures_df['InChIKey'].unique()))
    print()
    print('Not in SphinGOMAP, already in ChEBI')
    df_not_sphingomaps_chebi = participants_structures_df[participants_structures_df['SMID'].isna()]

    print("Structures not in SphinGOMAP (generated to fill the gaps)", len(df_not_sphingomaps_chebi['InChIKey'].unique()))
    print("CHEBI_ID", len(df_not_sphingomaps_chebi[df_not_sphingomaps_chebi["CHEBI_ID"].notna()]['InChIKey'].unique()))
    print("CHEBI_ID_R1", len(df_not_sphingomaps_chebi[df_not_sphingomaps_chebi["CHEBI_ID_R1"].notna()]['InChIKey'].unique()))
    print("CHEBI_ID_R1_R2", len(df_not_sphingomaps_chebi[df_not_sphingomaps_chebi["CHEBI_ID_R1_R2"].notna()]['InChIKey'].unique()))
    print()
    participants_structures_df = participants_structures_df[participants_structures_df['SMID'].notna()].copy()
    print('Total SMID df:', len(participants_structures_df['SMID'].unique()))

    # remove all generated compounds
    
    df_t = participants_structures_df[participants_structures_df['Formula']=='no molecule listed'].copy()
    print('No molecule listed:', len(df_t['SMID'].unique()))

    # drop the compounds with no formula listed from further stats
    participants_structures_df = participants_structures_df[~participants_structures_df['SMID'].isin(df_t['SMID'])]
    participants_structures_df = participants_structures_df[~participants_structures_df['SMID'].isna()]
    print('With molecule listed:', len(participants_structures_df['SMID'].unique()))

    # drop ceramide without sugar
    participants_structures_df = participants_structures_df[participants_structures_df['SMID']!=1]
    print('Without ceramide:', len(participants_structures_df['SMID'].unique()))
    vcf= pd.DataFrame(participants_structures_df['Formula'].value_counts())
    vcf.reset_index(inplace=True)
    print('Duplicated formulas', vcf[vcf['count']>1]['Formula'].to_list())
    participants_structures_df.drop_duplicates(subset=['Formula'], inplace=True)
    print('Without formula duplicates:', len(participants_structures_df['SMID'].unique()))

    participants_structures_df.drop_duplicates(subset=['InChIKey'], inplace=True)
    print('Without InChIKey duplicates:', len(participants_structures_df['SMID'].unique()))

    df_hypothetical = participants_structures_df[participants_structures_df['See Other']=='hypothetical']
    print("CHEBI_ID hypothetical", len(df_hypothetical[df_hypothetical["CHEBI_ID"].notna()]['InChIKey'].unique()))
    print("CHEBI_ID_R1 hypothetical", len(df_hypothetical[df_hypothetical["CHEBI_ID_R1"].notna()]['InChIKey'].unique()))
    print("CHEBI_ID_R1_R2 hypothetical", len(df_hypothetical[df_hypothetical["CHEBI_ID_R1_R2"].notna()]['InChIKey'].unique()))

    return participants_structures_df

def create_table_per_category_coverage_without_hypothetical(participants_structures_df, hypothetical=False):

    if hypothetical==True:
        participants_structures_df = participants_structures_df[participants_structures_df['See Other']=='hypothetical']
    elif hypothetical == False:
        participants_structures_df = participants_structures_df[participants_structures_df['See Other']!='hypothetical']
    elif hypothetical == 'all':
        pass

    print('non-hypothetical sphingomaps compounds in final result:', len(participants_structures_df[~participants_structures_df['SMID'].isna()]['InChIKey'].unique()))
    df_sphingomaps_chebi = participants_structures_df[participants_structures_df['SMID'].notna()].copy()
    df_sphingomaps_chebi.loc[df_sphingomaps_chebi['category_assigned'].isna(), 'category_assigned'] = 'No defined category'

    vc_spmk = pd.DataFrame(df_sphingomaps_chebi['category_assigned'].value_counts())
    vc_spmk.rename(columns={'count':'Number of structures'}, inplace=True)


    vc_r1 = pd.DataFrame(df_sphingomaps_chebi[df_sphingomaps_chebi["CHEBI_ID_R1"].notna()]['category_assigned'].value_counts())
    vc_r1.rename(columns={'count':'Structures in ChEBI (d18:1 / R1)'}, inplace=True)

    vc_r2 = pd.DataFrame(df_sphingomaps_chebi[df_sphingomaps_chebi["CHEBI_ID_R1_R2"].notna()]['category_assigned'].value_counts())
    vc_r2.rename(columns={'count':'Structures in ChEBI (R1 / R2)'}, inplace=True)

    vc = vc_spmk.merge(vc_r1, how='outer', on='category_assigned')
    vc = vc.merge(vc_r2, how='outer', on='category_assigned')
    vc.sort_values(by=['Number of structures'], inplace=True, ascending=False)
    vc = vc.fillna(0)

    vc['Structures in ChEBI (d18:1 / R1)'] = vc['Structures in ChEBI (d18:1 / R1)'].astype(int)
    vc['Structures in ChEBI (R1 / R2)'] = vc['Structures in ChEBI (R1 / R2)'].astype(int)
    vc.loc['Total'] = vc.sum(numeric_only=True)
    vc['r1 (%)'] = (vc['Structures in ChEBI (d18:1 / R1)'] / vc['Number of structures']) * 100
    vc['r1 (%)'] = vc['r1 (%)'].astype(int)
    vc['r2 (%)'] = (vc['Structures in ChEBI (R1 / R2)'] / vc['Number of structures']) * 100
    vc['r2 (%)'] = vc['r2 (%)'].astype(int)
    vc = vc.reset_index()
    vc.rename(columns={'category_assigned':'Category (Series)'}, inplace=True)
    vc.to_csv('sphingomapkey_compounds_in_chebi.tsv', sep='\t', index=False, 
              columns=['Category (Series)', 'Number of structures', 'Structures in ChEBI (d18:1 / R1)', 'r1 (%)', 'Structures in ChEBI (R1 / R2)', 'r2 (%)'])

def fetch_sphingomap_id(lm_id):
    # URL of the lipidmaps page
    url = f"https://lipidmaps.org/databases/lmsd/{lm_id}"

    # Make an HTTP GET request to fetch the page content
    response = requests.get(url)
    if response.status_code == 200:
        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Searching for the anchor tag that contains the ChEBI ID
        chebi_link = soup.find('a', href=lambda href: href and "chebiId" in href)

        # Extracting the ChEBI ID
        if chebi_link:
            chebi_id = chebi_link['href'].split('chebiId=')[1]
        else:
            chebi_id = None

        # Locate all bold labels and find 'SphingoMap ID' specifically
        labels = soup.find_all('div', class_='font-bold')
        for label in labels:
            if label.text.strip() == 'SphingoMap ID':
                # Assuming the SphingoMap ID is within the next div within the same larger grid structure
                sphingomap_id_div = label.find_next('div')
                if sphingomap_id_div:
                    return sphingomap_id_div.text.strip(), chebi_id
        return "SphingoMap ID not found.", chebi_id
    
    return "Failed to fetch page.", None

##############################################################################################################
#         LIPIDMAPS        LIPIDMAPS        LIPIDMAPS        LIPIDMAPS        LIPIDMAPS        LIPIDMAPS
##############################################################################################################

def match_lipidmap(participant_structures, df_chebi):

    if not os.path.exists('lmid2smid.tsv'):
        webscrap_lipidmaps_link_to_chebi()
    df_lipidmaps = pd.read_csv('lmid2smid.tsv', sep='\t')
    df_lipidmaps = df_lipidmaps[df_lipidmaps['SphingoMap ID'].notna()]
    df_lipidmaps = df_lipidmaps[df_lipidmaps['SphingoMap ID']!='SphingoMap ID not found.']
    df_lipidmaps = df_lipidmaps[df_lipidmaps['CHEBI'].notna()]

    print('ChebiID-Sphingomaps ID links in lipidmaps:', len(df_lipidmaps))
    df_lipidmaps['CHEBI_ID_lipidmaps'] = df_lipidmaps['CHEBI'].apply(lambda x: x.replace('CHEBI:',''))
    df_lipidmaps['Chebi_SMID'] = df_lipidmaps.apply(lambda row: row['CHEBI_ID_lipidmaps']+'-'+row['SphingoMap ID'], axis=1)
    
    participants_structures_df.dropna(subset=['SMID'], inplace=True)
    participants_structures_df['SMID'] = participants_structures_df['SMID'].astype(int)

    df_R1 = participant_structures[participant_structures['CHEBI_ID_R1'].notna()].copy()
    df_R1['CHEBI_ID_R1'] = df_R1['CHEBI_ID_R1'].astype(str).str.split(',')
    df_R1['CHEBI_ID_R1'] = df_R1['CHEBI_ID_R1'].apply(lambda lst: [x.strip() for x in lst])
    df_R1 = df_R1.explode('CHEBI_ID_R1')
    df_R1['CHEBI_ID_R1'] = df_R1['CHEBI_ID_R1'].astype(int)
    df_R1['Chebi_SMID'] = df_R1.apply(lambda row: f"{row['CHEBI_ID_R1']}-{row['SMID']}", axis=1)

    df_R1_R2 = participant_structures[participant_structures['CHEBI_ID_R1_R2'].notna()].copy()
    df_R1_R2['CHEBI_ID_R1_R2'] = df_R1_R2['CHEBI_ID_R1_R2'].astype(str).str.split(',')
    df_R1_R2['CHEBI_ID_R1_R2'] = df_R1_R2['CHEBI_ID_R1_R2'].apply(lambda lst: [x.strip() for x in lst])
    df_R1_R2 = df_R1_R2.explode('CHEBI_ID_R1_R2')
    df_R1_R2['CHEBI_ID_R1_R2'] = df_R1_R2['CHEBI_ID_R1_R2'].astype(int)
    df_R1_R2['Chebi_SMID'] = df_R1_R2.apply(lambda row: f"{row['CHEBI_ID_R1_R2']}-{row['SMID']}", axis=1)

    df_lipidmaps['direct_match_R1'] = df_lipidmaps['Chebi_SMID'].isin(df_R1['Chebi_SMID'])
    df_lipidmaps['direct_match_R2'] = df_lipidmaps['Chebi_SMID'].isin(df_R1_R2['Chebi_SMID'])

    df_sphingomapkey = participant_structures

    df_sphingomapkey['SMID'] = df_sphingomapkey['SMID'].astype(str)
    df_lipidmaps['SphingoMap ID'] = df_lipidmaps['SphingoMap ID'].astype(str)
    df_lipidmaps = df_lipidmaps.merge(df_sphingomapkey, left_on = 'SphingoMap ID', right_on='SMID', how='left')
    df_lipidmaps['same_compound_different_protonation_R1'] = df_lipidmaps.apply(lambda row: check_if_same_compound_but_different_protonation(row['CHEBI_ID_lipidmaps'], row['CHEBI_ID_R1'], df_chebi), axis=1)
    df_lipidmaps['same_compound_different_protonation_R2'] = df_lipidmaps.apply(lambda row: check_if_same_compound_but_different_protonation(row['CHEBI_ID_lipidmaps'], row['CHEBI_ID_R1_R2'], df_chebi), axis=1)
    df_lipidmaps['overall_match'] = df_lipidmaps.apply(lambda row: any([row['same_compound_different_protonation_R1'], row['same_compound_different_protonation_R2'],
                                                                       row['direct_match_R1'], row['direct_match_R2'] ]    
                                                                                        ), axis=1)
    df_lipidmaps.sort_values(by=['overall_match'], inplace=True)
    print('Total SMID in lipidsmaps', len(df_lipidmaps['SphingoMap ID'].unique()))
    print('Total match:', len(df_lipidmaps[df_lipidmaps['overall_match']==True]['SphingoMap ID'].unique()))
    df_lipidmaps.to_csv('df_lipidmaps_43.tsv', sep='\t', index=False, columns = [
        'lm_id',	'SphingoMap ID',	'overall_match', 'direct_match_R1',	'direct_match_R2',
          'same_compound_different_protonation_R1', 'same_compound_different_protonation_R2', 'Formula'	,
        'Category',	'Common Name','CHEBI_ID_lipidmaps', 'CHEBI_ID_R1',	'CHEBI_ID_R1_R2'	
    ])

def webscrap_lipidmaps_link_to_chebi():
    df_lipidmaps = pd.read_csv('lipidmaps_ids_cc0.tsv', sep='\t')
    df_no_structure = df_lipidmaps[df_lipidmaps['inchi_key'].isna()]
    df_no_structure['SP']=df_no_structure['lm_id'].apply(lambda x: x.startswith('LMSP'))
    df_no_structure = df_no_structure[df_no_structure['SP']==True]
    print('Total sphingolipids in LIPIDMAPS with *', len(df_no_structure['lm_id'].unique()))
    df_no_structure[['SphingoMap ID','CHEBI']] = df_no_structure.progress_apply(lambda row: fetch_sphingomap_id(row['lm_id']), axis=1, result_type='expand')
    df_no_structure.to_csv('lmid2smid.tsv', sep='\t', index=False)

def check_if_same_compound_but_different_protonation(chebi1, chebi2, df_chebi):
    if chebi1!='nan' and chebi2!='nan':
        inchikeys_1 = [i[:-2] for i in get_inchikeys_from_chebi(chebi1, df_chebi)] # crop the charge of inchikey
        inchikeys_2 = [i[:-2] for i in get_inchikeys_from_chebi(chebi2, df_chebi)]
        return len(set(inchikeys_1).intersection(set(inchikeys_2)))>0

def get_inchikeys_from_chebi(chebiid, df_chebi):
    if ',' in chebiid:
        return get_inchikeys_from_chebi(chebiid.split(',')[0], df_chebi)
    result = df_chebi[df_chebi['COMPOUND_ID']==int(chebiid)]
    return(result['inchikey_nostar'].to_list())

##############################################################################################################
#         Execute analysis
##############################################################################################################

parser = argparse.ArgumentParser(description="Analyze ChEBI matches from files in a folder.")
parser.add_argument(
    'input_folder',
    type=str,
    help='Path to the input folder containing files to analyze.'
)
args = parser.parse_args()

input_folder = args.input_folder

if not os.path.isdir(input_folder):
    print(f"Error: '{input_folder}' is not a valid directory.", file=sys.stderr)
    sys.exit(1)

if not os.path.exists('df_chebi.tsv'):
    df_chebi = read_chebi()
else:
    df_chebi = pd.read_csv('df_chebi.tsv', sep='\t')

df_chebi['COMPOUND_ID'] = df_chebi['COMPOUND_ID'].astype(int)

process_participants_structures_df(df_chebi, input_folder)
participants_structures_df = analyse_results()
create_table_per_category_coverage_without_hypothetical(participants_structures_df, hypothetical=False)
match_lipidmap(participants_structures_df, df_chebi)
