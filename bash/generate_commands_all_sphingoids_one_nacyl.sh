#!/usr/bin/env bash

# sphingoid base names
sphingoids_names=("sphinganine" "sphing-4-enine" "4-hydroxysphinganine" "15-methylhexadecasphinganine" "3-dehydro-D-sphinganine" "hexadecasphinganine" "3-dehydrohexadecasphinganine" "3-dehydroeicosasphinganine" "3-dehydrotetradecasphinganine" "heptadecasphing-4-enine" "15-methylhexadecasphing-4-enine" "3-dehydro-15-methylhexadecasphinganine" "eicosasphinganine" "(4R)-hydroxyeicosasphinganine")
sphingoids_smiles=("[C@H](O)CCCCCCCCCCCCCCC" "[C@H](O)/C=C/CCCCCCCCCCCCC" "[C@H](O)[C@H](O)CCCCCCCCCCCCCC" "[C@H](O)CCCCCCCCCCC(C)CC" "C(=O)CCCCCCCCCCCCCCC" "[C@H](O)CCCCCCCCCCCCC" "C(=O)CCCCCCCCCCCCC" "C(=O)CCCCCCCCCCCCCCCCC" "C(=O)CCCCCCCCCCC" "[C@H](O)/C=C/CCCCCCCCCCCC" "[C@H](O)/C=C/CCCCCCCCC(C)CC" "C(=O)CCCCCCCCCCC(C)CC" "[C@H](O)CCCCCCCCCCCCCCCCC" "[C@H](O)[C@H](O)CCCCCCCCCCCCCCCC")


# n-acyl names
nacyls_names=(
"hexadecanoate"
)

nacyls_smiles=(
"CCCCCCCCCCCCCCC"
)

# loop to generate commands
for s in "${!sphingoids_names[@]}"; do
  sph_name="${sphingoids_names[$s]}"
  sph_smiles="${sphingoids_smiles[$s]}"
  safe_sph=$(echo "$sph_name" | tr ' ,()' '_')

  for n in "${!nacyls_names[@]}"; do
    nacyl_name="${nacyls_names[$n]}"
    nacyl_smiles="${nacyls_smiles[$n]}"
    safe_nacyl=$(echo "$nacyl_name" | tr ' ,()' '_')

    folder="results_${safe_sph}_${safe_nacyl}"

    echo "glycosphingotool process-all --output-folder '$folder' --nacyl '$nacyl_smiles' --sphingoid '$sph_smiles'"
  done
done