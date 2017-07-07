#!/usr/bin/env bash

#######################################################
# We define all clades we want to color or highlight
#######################################################
gotree subtree -i !{ncbi} -n "^Monotremata$" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_monotremata.txt
gotree subtree -i !{ncbi} -n "Metatheria$" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_metatheria.txt
gotree subtree -i !{ncbi} -n "Xenarthra" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_xenarthra.txt
gotree subtree -i !{ncbi} -n "Afrotheria" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_afrotheria.txt
gotree subtree -i !{ncbi} -n "Rodentia" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_rodentia.txt
gotree subtree -i !{ncbi} -n "Lagomorpha" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_lagomorpha.txt
gotree subtree -i !{ncbi} -n "Primata" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_primata.txt
gotree subtree -i !{ncbi} -n "Insectivora" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_insectivora.txt
gotree subtree -i !{ncbi} -n "Chiroptera" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_chiroptera.txt
gotree subtree -i !{ncbi} -n "Carnivora" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_carnivora.txt
gotree subtree -i !{ncbi} -n "Cetartiodactyla" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_cetartiodactyla.txt
gotree subtree -i !{ncbi} -n "Perissodactyla" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_perissodactyla.txt
gotree subtree -i !{ncbi} -n "Smutsia.Manis" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_manis.txt
gotree subtree -i !{ncbi} -n "^Tupaia$"  | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_scandentia.txt
# Sub cetartidactyla
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Cetacea" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_cetacea.txt
# Sub rodentia
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Sciurognathi" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_sciurognathi.txt
# Sub sciurognathi
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Muroidea" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_muroidea.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Dipodidae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_dipodidae.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Sciuridae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_sciuridae.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Heteromyidae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_heteromyidae.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Geomyidae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_geomyidae.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Perognathinae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_perognathinae.txt
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Dipodomys$" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_dipodomyinae.txt
# Sub carnivora
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Mustelinae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_mustelinae.txt
# Sub afrotheria
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Elephantidae" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_elephantidae.txt
#Simians
gotree subtree -i ../data/ncbitax/ncbi_labels.nw -n "^Anthropoidea.Simiiformes" | gotree stats tips | tail -n+2 | cut -f 4 > ncbi_simians.txt
# Dermoptera
echo "Galeopterus_variegatus" > ncbi_dermoptera.txt

#######################################################
# We define iTOL annotations / clade colors
#######################################################
declare -a CLADES
CLADES=(monotremata metatheria afrotheria rodentia primata insectivora chiroptera carnivora cetartiodactyla perissodactyla)
declare -a COLORS
COLORS=("#114477" "#4477AA" "#117755" "#44AA88" "#777711" "#AAAA44" "#DDDD77" "#771111" "#AA4444", "#DD7777")

cat > mammals_itol.txt <<EOF
TREE_COLORS
SEPARATOR SPACE
DATASET_LABEL Clades
COLOR #ff0000
DATA
EOF

LENGTH=${#CLADES[@]}
for ((i=0; i<$LENGTH; i++))
do
    F=ncbi_${CLADES[$i]}.txt
    C=${COLORS[$i]}
    for l in `cat $F`
    do
	echo "$l branch $C normal 5" >> mammals_itol.txt
    done
done

#######################################################
# We upload to iTOL
#######################################################
urltbe=$(gotree rename -i !{tbetree} -m !{namemap} -r | gotree reroot outgroup Ornithorhynchus_anatinus Zaglossus_bruijni Tachyglossus_aculeatus | gotree upload itol --name "tbe" mammals_itol.txt)
urlfbp=$(gotree rename -i !{fbptree} -m !{namemap} -r | gotree reroot outgroup Ornithorhynchus_anatinus Zaglossus_bruijni Tachyglossus_aculeatus | gotree upload itol --name "fbp" mammals_itol.txt)

#######################################################
# We define iTOL drawing options
#######################################################
cat > options_tbe.txt <<EOF
display_mode	2
label_display	1
align_labels	1
ignore_branch_length	1
bootstrap_display	1
bootstrap_type	1
bootstrap_symbol	1
bootstrap_slider_min	0.7
bootstrap_slider_max	1
bootstrap_symbol_min	20
bootstrap_symbol_max	20
bootstrap_symbol_color	#c8c7fc
current_font_size	20
line_width	2
inverted	0
EOF
cat > options_fbp.txt <<EOF
display_mode	2
label_display	1
align_labels	1
ignore_branch_length	1
bootstrap_display	1
bootstrap_type	1
bootstrap_symbol	1
bootstrap_slider_min	0.7
bootstrap_slider_max	1
bootstrap_symbol_min	20
bootstrap_symbol_max	20
bootstrap_symbol_color	#e9a640
current_font_size	20
line_width	2
inverted	0
EOF

#######################################################
# We download images from iTOL
#######################################################
gotree download itol -i $(basename $urlfbp) -f svg -o !{outprefix}_fbp.svg -c options_fbp.txt
gotree download itol -i $(basename $urltbe) -f svg -o !{outprefix}_tbe.svg -c options_tbe.txt
