
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
[C@@H]([NH3+])CO is removed since it is accounted in the code, the table summarises what should be used for --sphingoid option of the code

| Identifier                                                                                                           | Name                                   | SMILES for the code           |
|----------------------------------------------------------------------------------------------------------------------|----------------------------------------|-------------------|
| [SLM:000000035](https://www.swisslipids.org/#/entity/SLM:000000035/)                                                 | sphinganine                            | "CCCCCCCCCCCCCCC[C@@H](O)"  |
| [SLM:000000412](https://www.swisslipids.org/#/entity/SLM:000000412/)                                                 | sphing-4-enine                         | "[C@H](O)/C=C/CCCCCCCCCCCCC"  |
| [SLM:000000386](https://www.swisslipids.org/#/entity/SLM:000000386/)                                                 | 4-hydroxysphinganine                   | "[C@H](O)[C@H](O)CCCCCCCCCCCCCC"  |
| [SLM:000000006](https://www.swisslipids.org/#/entity/SLM:000000006/)                                                 | 15-methylhexadecasphinganine           | "[C@H](O)CCCCCCCCCCC(C)CC"  |
| [SLM:000000384](https://www.swisslipids.org/#/entity/SLM:000000384/)                                                 | 3-dehydro-D-sphinganine                | "C(=O)CCCCCCCCCCCCCCC"  |
| [SLM:000389783](https://www.swisslipids.org/#/entity/SLM:000389783/)                                                 | hexadecasphinganine                    | "[C@H](O)CCCCCCCCCCCCC"  |
| [SLM:000000440](https://www.swisslipids.org/#/entity/SLM:000000440/)                                                 | 3-dehydrohexadecasphinganine           | "C(=O)CCCCCCCCCCCCC"  |
| [SLM:000000442](https://www.swisslipids.org/#/entity/SLM:000000442/)                                                 | 3-dehydroeicosasphinganine             | "C(=O)CCCCCCCCCCCCCCCCC"  |
| [SLM:000389965](https://www.swisslipids.org/#/entity/SLM:000389965/)                                                 | 3-dehydrotetradecasphinganine          | "C(=O)CCCCCCCCCCC"  |
| [SLM:000000623](https://www.swisslipids.org/#/entity/SLM:000000623/)                                                 | heptadecasphing-4-enine                | "[C@H](O)/C=C/CCCCCCCCCCCC"  |
| [SLM:000000003](https://www.swisslipids.org/#/entity/SLM:000000003/)                                                 | 15-methylhexadecasphing-4-enine        | "[C@H](O)/C=C/CCCCCCCCC(C)CC"  |
| [SLM:000000247](https://www.swisslipids.org/#/entity/SLM:000000247/)                                                 | 3-dehydro-15-methylhexadecasphinganine | "C(=O)CCCCCCCCCCC(C)CC" |
| [SLM:000000864](https://www.swisslipids.org/#/entity/SLM:000000864/)                                                 | eicosasphinganine                      | "[C@H](O)CCCCCCCCCCCCCCCCC"  |
| [SLM:000000865](https://www.swisslipids.org/#/entity/SLM:000000865/)                                                 | (4R)-hydroxyeicosasphinganine          | "CCCCCCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H]([NH3+])CO"|
