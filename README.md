
## Examples

* Convert glyco nomenclature of a sphingolipid into SMILES with defined n-acyl and sphingoid base

```bash
glycosphingotool convert "NeuAcalpha2-3Galbeta1-4GlcCer" --nacyl CCC --sphingoid "[C@H](O)/C=C/CC"
```
* Generate synthesis reactions for glyco nomenclature of a sphingolipid, generate reactions SMILES, RInChI and Web-RInChIKeys with defined n-acyl and sphingoid base

```bash
glycosphingotool generate "NeuAcalpha2-3Galbeta1-4GlcCer" --nacyl CCC --sphingoid "[C@H](O)/C=C/CC" --output-folder "NeuAcalpha2-3Galbeta1-4GlcCer"
```

* Process all the Excel file downloaded from SphingoMAP (the original link is currently broken at the original resource)

```bash
glycosphingotool process-all --output-folder 'results_SphingoMAP' --nacyl CCC --sphingoid "[C@H](O)/C=C/CCCCCC"
```

The source SphingomapkeyV1.4.xls can be found in src/glycosphingotool/assets

# Sphingoid bases
`[C@@H]([NH3+])CO` is removed since it is accounted in the code, the table summarises what should be used for --sphingoid option of the code

| Identifier                                                                                                           | Name                                   | SMILES for the code           |
|----------------------------------------------------------------------------------------------------------------------|----------------------------------------|-------------------|
| [SLM:000000035](https://www.swisslipids.org/#/entity/SLM:000000035/)                                                 | sphinganine                            | `[C@H](O)CCCCCCCCCCCCCCC`  |
| [SLM:000000412](https://www.swisslipids.org/#/entity/SLM:000000412/)                                                 | sphing-4-enine                         | `[C@H](O)/C=C/CCCCCCCCCCCCC`  |
| [SLM:000000386](https://www.swisslipids.org/#/entity/SLM:000000386/)                                                 | 4-hydroxysphinganine                   | `[C@H](O)[C@H](O)CCCCCCCCCCCCCC`  |
| [SLM:000000006](https://www.swisslipids.org/#/entity/SLM:000000006/)                                                 | 15-methylhexadecasphinganine           | `[C@H](O)CCCCCCCCCCC(C)CC`  |
| [SLM:000000384](https://www.swisslipids.org/#/entity/SLM:000000384/)                                                 | 3-dehydro-D-sphinganine                | `C(=O)CCCCCCCCCCCCCCC`  |
| [SLM:000389783](https://www.swisslipids.org/#/entity/SLM:000389783/)                                                 | hexadecasphinganine                    | `[C@H](O)CCCCCCCCCCCCC`  |
| [SLM:000000440](https://www.swisslipids.org/#/entity/SLM:000000440/)                                                 | 3-dehydrohexadecasphinganine           | `C(=O)CCCCCCCCCCCCC`  |
| [SLM:000000442](https://www.swisslipids.org/#/entity/SLM:000000442/)                                                 | 3-dehydroeicosasphinganine             | `C(=O)CCCCCCCCCCCCCCCCC`  |
| [SLM:000389965](https://www.swisslipids.org/#/entity/SLM:000389965/)                                                 | 3-dehydrotetradecasphinganine          | `C(=O)CCCCCCCCCCC`  |
| [SLM:000000623](https://www.swisslipids.org/#/entity/SLM:000000623/)                                                 | heptadecasphing-4-enine                | `[C@H](O)/C=C/CCCCCCCCCCCC`  |
| [SLM:000000003](https://www.swisslipids.org/#/entity/SLM:000000003/)                                                 | 15-methylhexadecasphing-4-enine        | `[C@H](O)/C=C/CCCCCCCCC(C)CC`  |
| [SLM:000000247](https://www.swisslipids.org/#/entity/SLM:000000247/)                                                 | 3-dehydro-15-methylhexadecasphinganine | `C(=O)CCCCCCCCCCC(C)CC` |
| [SLM:000000864](https://www.swisslipids.org/#/entity/SLM:000000864/)                                                 | eicosasphinganine                      | `[C@H](O)CCCCCCCCCCCCCCCCC`  |
| [SLM:000000865](https://www.swisslipids.org/#/entity/SLM:000000865/)                                                 | (4R)-hydroxyeicosasphinganine          | `[C@H](O)[C@H](O)CCCCCCCCCCCCCCCC`|

# N-acyls

| Identifier                                                                                   | Name                                   | SMILES for the code |
|----------------------------------------------------------------------------------------------|----------------------------------------|---------------------|
[CHEBI:7896](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:7896)                       | | |
[CHEBI:32372](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32372)                     | | |
[CHEBI:25629](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:25629)                     | | |
[CHEBI:30823](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:30823)                     | | |
[CHEBI:30245](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:30245)                     | | |
[CHEBI:77222](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77222)                     | | |
[CHEBI:32360](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32360)                     | | |
[CHEBI:71589](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:71589)                     | | |
[CHEBI:32395](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32395)                     | | |
[CHEBI:58562](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:58562)                     | | |
[CHEBI:23858](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:23858)                     | | |
[CHEBI:77016](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77016)                     | | |
[CHEBI:31014](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:31014)                     | | |
[CHEBI:32392](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32392)                     | | |
[CHEBI:31013](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:31013)                     | | |
[CHEBI:77221](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77221)                     | | |
[CHEBI:65097](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:65097)                     | | |
[CHEBI:76724](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76724)                     | | |
[CHEBI:76732](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76732)                     | | |
[CHEBI:76722](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76722)                     | | |
[CHEBI:76723](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76723)                     | | |
[CHEBI:84324](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:84324)                     | | |
[CHEBI:76728](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76728)                     | | |
