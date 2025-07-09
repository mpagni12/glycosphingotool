
## Examples

* Convert glyco nomenclature of a sphingolipid into SMILES with defined n-acyl and sphingoid base

```bash
glycosphingotool convert "NeuAcalpha2-3Galbeta1-4GlcCer" --nacyl CCC --sphingoid "[C@H](O)/C=C/CC"
```
* Generate synthesis reactions for glyco nomenclature of a sphingolipid, generate reactions SMILES, RInChI and Web-RInChIKeys with defined n-acyl and sphingoid base

```bash
glycosphingotool generate "NeuAcalpha2-3Galbeta1-4GlcCer" --nacyl CCC --sphingoid "[C@H](O)/C=C/CC" --output-folder "results_NeuAcalpha2-3Galbeta1-4GlcCer"
```

* Process all the Excel file downloaded from SphingoMAP (the original link is currently broken at the original resource)

```bash
glycosphingotool process-all --output-folder 'results_SphingoMAP' --nacyl CCC --sphingoid "[C@H](O)/C=C/CCCCCC"
```
The source SphingomapkeyV1.4.xls can be found in src/glycosphingotool/assets

* The systematic enumeration can be done as following:
```bash
bash generate_commands.sh > run_all.sh
```
```bash
bash run_all.sh
```

* Extract all the compounds that were generated
```bash
glycosphingotool extract-cmp --output-folder results_sphinganine_hexadecanoate --sphingoid "[C@H](O)CCCCCCCCCCCCCCC"
```

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
`C([O-])=O` is dropped for the SMILES since it is accounted by the code

| Identifier | Name | SMILES for the code |
|------------|------|---------------------|
| [CHEBI:7896](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:7896) | hexadecanoate | `CCCCCCCCCCCCCCC` |
| [CHEBI:32372](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32372) | palmitoleate | `CCCCCCC\C=C/CCCCCC` |
| [CHEBI:25629](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:25629) | octadecanoate | `CCCCCCCCCCCCCCCCC` |
| [CHEBI:30823](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:30823) | oleate | `CCCCCCC\C=C/CCCCCCCC` |
| [CHEBI:30245](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:30245) | linoleate | `CCCCCCC\C=C/C\C=C/CCCCC` |
| [CHEBI:77222](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77222) | (6Z,9Z,12Z,15Z)-octadecatetraenoate | `CCCC\C=C/C\C=C/C\C=C/C\C=C/CC` |
| [CHEBI:32360](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32360) | icosanoate | `CCCCCCCCCCCCCCCCCCC` |
| [CHEBI:71589](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:71589) | all-cis-icosa-8,11,14-trienoate | `CCCCCC\C=C/C\C=C/C\C=C/CCCCC` |
| [CHEBI:32395](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32395) | arachidonate | `CCCCC\C=C/C\C=C/C\C=C/C\C=C/CCC` |
| [CHEBI:58562](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:58562) | all-cis-5,8,11,14,17-icosapentaenoate | `CCC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CC` |
| [CHEBI:23858](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:23858) | behenate | `CCCCCCCCCCCCCCCCCCCC` |
| [CHEBI:77016](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77016) | (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoate | `CC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CC` |
| [CHEBI:31014](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:31014) | tetracosanoate | `CCCCCCCCCCCCCCCCCCCCCCC` |
| [CHEBI:32392](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:32392) | (15Z)-tetracosenoate | `CCCCCCCCCCCCC\C=C/CCCCCCCC` |
| [CHEBI:31013](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:31013) | cerotate | `CCCCCCCCCCCCCCCCCCCCCCCCC` |
| [CHEBI:77221](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:77221) | (17Z)-hexacosenoate | `CCCCCCCCCCCCCCC\C=C/CCCCCCCC` |
| [CHEBI:65097](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:65097) | 2-hydroxyhexadecanoate | `C(O)CCCCCCCCCCCCCC` |
| [CHEBI:76724](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76724) | 2-hydroxyoctadecanoate | `C(O)CCCCCCCCCCCCCCCC` |
| [CHEBI:76732](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76732) | 2-hydroxyarachidate | `C(O)CCCCCCCCCCCCCCCCCC` |
| [CHEBI:76722](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76722) | 2-hydroxybehenate | `C(O)CCCCCCCCCCCCCCCCCCCC` |
| [CHEBI:76723](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76723) | 2-hydroxytetracosanoate | `C(O)CCCCCCCCCCCCCCCCCCCCCC` |
| [CHEBI:84324](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:84324) | 2-hydroxynervonate | `C(O)CCCCCCCCCCCC\C=C/CCCCCCCC` |
| [CHEBI:76728](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:76728) | 2-hydroxyhexacosanoate | `C(O)CCCCCCCCCCCCCCCCCCCCCCCC` |

## Underlying algorithm

The identification of reactions necessary for synthesis is based on recursive graph backtracking

```python
import networkx as nx

def custom_dfs(graph, node, visited=None, result=None):
    """
    depth first search
    """
    if visited is None:
        visited = set()  # To keep track of visited nodes
    if result is None:
        result = []  # To store the traversal order

    # Mark the current node as visited and add to the result list
    visited.add(node)
    result.append(node)  # Process the node

    # Recur for all the neighbors (children in a tree context)
    for neighbor in graph.neighbors(node):
        if neighbor not in visited:
            custom_dfs(graph, neighbor, visited, result)
            result.append(node)  # Add the node again on return to represent backtracking

    return result

# Example graph creation
G = nx.DiGraph()
edges = [('A', 'B'), ('A', 'C'), ('B', 'D'), ('B', 'E'), ('C', 'F'), ('C', 'G')]
G.add_edges_from(edges)

# Start DFS from the root node 'A'
traversal_result = custom_dfs(G, 'A')
print("Traversal Result:", traversal_result)
```
