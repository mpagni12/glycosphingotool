# Analyses presented in the paper
## Analyse matching between the generated data and ChEBI/LIPIDMAPS
```bash
cd chebi_lipidmaps_match
```

```bash
python analyse_chebi_match.py ../../results_sphing-4-enine_hexadecanoate/
```
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
