#!/usr/bin/env bash

# sphingoid base names
sphingoids_names=("sphinganine" "sphing-4-enine")
sphingoids_smiles=("[C@H](O)CCCCCCCCCCCCCCC" "[C@H](O)/C=C/CCCCCCCCCCCCC")

# n-acyl names
nacyls_names=(
"hexadecanoate"
"palmitoleate"
"octadecanoate"
"oleate"
"linoleate"
"(6Z,9Z,12Z,15Z)-octadecatetraenoate"
"icosanoate"
"all-cis-icosa-8,11,14-trienoate"
"arachidonate"
"all-cis-5,8,11,14,17-icosapentaenoate"
"behenate"
"(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoate"
"tetracosanoate"
"(15Z)-tetracosenoate"
"cerotate"
"(17Z)-hexacosenoate"
"2-hydroxyhexadecanoate"
"2-hydroxyoctadecanoate"
"2-hydroxyarachidate"
"2-hydroxybehenate"
"2-hydroxytetracosanoate"
"2-hydroxynervonate"
"2-hydroxyhexacosanoate"
)

nacyls_smiles=(
"CCCCCCCCCCCCCCC"
"CCCCCCC\\C=C/CCCCCC"
"CCCCCCCCCCCCCCCCC"
"CCCCCCC\\C=C/CCCCCCCC"
"CCCCCCC\\C=C/C\\C=C/CCCCC"
"CCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CC"
"CCCCCCCCCCCCCCCCCCC"
"CCCCCC\\C=C/C\\C=C/C\\C=C/CCCCC"
"CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC"
"CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC"
"CCCCCCCCCCCCCCCCCCCC"
"CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC"
"CCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCC\\C=C/CCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCC\\C=C/CCCCCCCC"
"C(O)CCCCCCCCCCCCCC"
"C(O)CCCCCCCCCCCCCCCC"
"C(O)CCCCCCCCCCCCCCCCCC"
"C(O)CCCCCCCCCCCCCCCCCCCC"
"C(O)CCCCCCCCCCCCCCCCCCCCCC"
"C(O)CCCCCCCCCCCC\\C=C/CCCCCCCC"
"C(O)CCCCCCCCCCCCCCCCCCCCCCCC"
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