import pandas as pd
from rxnsmiles2rinchi import RInChI
from pyrheadb.RheaDB import RheaDB
import os


def rxnsmiles_to_rinchi(SMILES, rinchi):
    """Convert reaction SMILES to RInChI and RInChIKey."""
    try:
        return rinchi.rxn_smiles_to_rinchi_rinchikey(SMILES)
    except Exception:
        return None, None


def load_and_process_smiles(df, smiles_col='rxn_smiles_star'):
    """Replace wildcards in SMILES strings and generate RInChI keys."""
    df['rxnSMILES'] = df[smiles_col].str.replace(r'\[2\*]|\[1\*]|\*', 'CCCCCCCCCCCCCCC', regex=True)
    rinchi = RInChI()
    df[['RInChI', 'RInChIKey']] = df.apply(
        lambda row: rxnsmiles_to_rinchi(row['rxnSMILES'], rinchi),
        axis=1, result_type='expand'
    )
    return df


def compare_and_print(set_a, set_b, label_a, label_b):
    """Print overlap statistics between two sets."""
    print(f'Only {label_a}:', len(set_a - set_b))
    print(f'Only {label_b}:', len(set_b - set_a))
    print(f'{label_a} and {label_b}:', len(set_a & set_b))
    print(f'Total unique {label_b} rinchis:', len(set_b))


def check_overlap_manual_and_generated(df_generated_path):
    """Compare manually curated set of 55 Rhea reactions with generated ones that served as a test."""
    df_directions = pd.read_csv('rhea-directions.tsv', sep='\t')
    df_smiles = pd.read_csv('rhea-reaction-smiles.tsv', sep='\t', names=['RHEA_ID', 'rxn_smiles_star'])
    df_manual = pd.read_csv('Manual rhea reactions by Series.tsv', sep='\t', names=['RHEA:ID', 'description'])
    df_manual['rhea_id'] = df_manual['RHEA:ID'].str.replace('RHEA:', '').astype(int)

    df_manual = df_manual.merge(df_directions, left_on='rhea_id', right_on='RHEA_ID_MASTER')
    df_manual = df_manual.merge(df_smiles, left_on='RHEA_ID_LR', right_on='RHEA_ID')
    df_manual = load_and_process_smiles(df_manual)

    df_generated = pd.read_csv(df_generated_path, sep='\t')

    set_manual = set(df_manual['RInChIKey'].dropna())
    set_generated = set(df_generated['RInChIKey'].dropna())

    compare_and_print(set_manual, set_generated, 'manual', 'generated')


def check_overlap_all_rhea_and_generated(df_generated_path, output_file):
    """Compare all Rhea reactions with generated reactions from a specific file."""
    rdb = RheaDB()
    df_directions = pd.read_csv('rhea-directions.tsv', sep='\t')
    df_smiles = pd.read_csv('rhea-reaction-smiles.tsv', sep='\t', names=['RHEA_ID', 'rxn_smiles_star'])

    df_all = df_directions.merge(df_smiles, left_on='RHEA_ID_LR', right_on='RHEA_ID')
    df_all = load_and_process_smiles(df_all)

    df_generated = pd.read_csv(df_generated_path, sep='\t')

    set_rhea = set(df_all['RInChIKey'].dropna())
    set_generated = set(df_generated['RInChIKey'].dropna())

    print(f'Generated and rhea:', len(set_rhea & set_generated))
    print(f'Total unique generated rinchis:', len(set_generated))

    df_all = df_all.merge(df_generated, on='RInChIKey')
    df_all = df_all.merge(rdb.df_reactions, left_on='RHEA_ID_MASTER', right_on='MASTER_ID', how='left')
    df_all.drop_duplicates(subset=['MASTER_ID'], inplace=True)
    df_all.to_csv(output_file, sep='\t', index=False,
                  columns=['MASTER_ID', 'glyco_nomenclature', 'reaction_participant_names'])


######################
#        MAIN        #
######################
if __name__ == '__main__':
    check_overlap_manual_and_generated(
        df_generated_path=os.path.join('..','..','results_sphing-4-enine_hexadecanoate','sphingomapkey with reaction structures.tsv'))

    check_overlap_all_rhea_and_generated(
        df_generated_path=os.path.join('..','..','results_sphing-4-enine_hexadecanoate','sphingomapkey with reaction structures.tsv'),
        output_file='rhea_and_generated_r1.tsv'
    )

    check_overlap_all_rhea_and_generated(
        df_generated_path=os.path.join('..','..','results_sphinganine_hexadecanoate','sphingomapkey with reaction structures.tsv'),
        output_file='rhea_and_generated_r1_r2.tsv'
    )
