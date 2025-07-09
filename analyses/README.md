# Analyses presented in the paper
## Analyse matching between the generated data and ChEBI/LIPIDMAPS

The further assumes that you have run command to generate reactions for sphing-4-enine with hexadenoate and have corresponding results folder in your glycosphingotool folder
```bash
glycosphingotool process-all --output-folder 'results_sphing-4-enine_hexadecanoate' --nacyl 'CCCCCCCCCCCCCCC' --sphingoid '[C@H](O)/C=C/CCCCCCCCCCCCC'
```
* python scripts for analysis are supposed to be executed from within the directory
```bash
cd chebi_lipidmaps_match
```
* run script to analyse the compounds match
```bash
python analyse_chebi_match.py ../../results_sphing-4-enine_hexadecanoate/
```
* The table with percentages presented in the paper is in the 'sphingomapkey_compounds_in_chebi.tsv'.

output:
```bash
100%|████████████████████████████| 912/912 [00:03<00:00, 286.39it/s]
100%|████████████████████████████| 912/912 [00:03<00:00, 282.40it/s]
100%|████████████████████████████| 912/912 [00:03<00:00, 283.97it/s]
Unique SMID in original data 412
Unique glyconomenclature formulas in original data 368
Total sphingolipids (nomenclature) in df: 911
Total sphingolipids (InChIKey) in df: 911

Not in SphinGOMAP, already in ChEBI
Structures not in SphinGOMAP (generated to fill the gaps) 545
CHEBI_ID 0
CHEBI_ID_R1 13
CHEBI_ID_R1_R2 6

Total SMID df: 412
No molecule listed: 43
With molecule listed: 369
Without ceramide: 368
Duplicated formulas ['NeuAcalpha2-8NeuAcalpha2-3(Fucalpha1-2Galbeta1-3GalNAcbeta1-4)Galbeta1-4GlcCer', 'Galbeta1-3GalNAcbeta1-3Galalpha1-3Galbeta1-4GlcCer']
Without formula duplicates: 366
Without InChIKey duplicates: 366
CHEBI_ID hypothetical 0
CHEBI_ID_R1 hypothetical 2
CHEBI_ID_R1_R2 hypothetical 0
non-hypothetical sphingomaps compounds in final result: 287
ChebiID-Sphingomaps ID links in lipidmaps: 43
Total SMID in lipidsmaps 43
Total match: 43
```

## Analyse matching between Rhea and generated reactions
This analysis assumes you executed reaction generation for sphing-4-enine and sphinganine with hexadecanoate and have corresponding results folder in your glycosphingotool folder
```bash
glycosphingotool process-all --output-folder 'results_sphing-4-enine_hexadecanoate' --nacyl 'CCCCCCCCCCCCCCC' --sphingoid '[C@H](O)/C=C/CCCCCCCCCCCCC'
```
and
```bash
glycosphingotool process-all --output-folder 'results_sphinganine_hexadecanoate' --nacyl 'CCCCCCCCCCCCCCC' --sphingoid '[C@H](O)CCCCCCCCCCCCCCC'
```

* python scripts for analysis are supposed to be executed from within the directory
```bash
cd reactions_analysis
```
* run script to analyse the reactions match
```bash
python check_overlap_manual_and_generated.py
```
output
```bash
Only manual: 0
Only generated: 1932
manual and generated: 55
Total unique generated rinchis: 1987
Your Rhea DB version is 139
Using previously downloaded Rhea version
Generated and rhea: 85
Total unique generated rinchis: 1987
Your Rhea DB version is 139
Using previously downloaded Rhea version
Generated and rhea: 80
Total unique generated rinchis: 1987
```

* run script to generate network image
```bash
python visualise_network_graph.py
```
