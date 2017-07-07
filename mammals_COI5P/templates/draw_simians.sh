#!/usr/bin/env bash

# We create annotation file for simian clade using data from ncbi taxonomy
gotree subtree -i !{ncbi} -n "^Anthropoidea.Simiiformes" | gotree stats tips | tail -n+2 | cut -f 4 > simians.txt
echo "simians:$(cat simians.txt | tr '\n' ',')" > map.txt

# We rename the fbp tree using map file && annotate the simian clade && taking the subtree
gotree rename -i !{fbptree} -m !{namemap} -r | gotree reroot outgroup Ornithorhynchus_anatinus Zaglossus_bruijni Tachyglossus_aculeatus | gotree annotate -m map.txt | gotree subtree -n simians > simians_fbp.nw
gotree rename -i !{tbetree} -m !{namemap} -r | gotree reroot outgroup Ornithorhynchus_anatinus Zaglossus_bruijni Tachyglossus_aculeatus | gotree annotate -m map.txt | gotree subtree -n simians > simians_tbe.nw

# We create iTOL annotation file skeleton
echo -e "TREE_COLORS\nSEPARATOR SPACE\nDATASET_LABEL Simian_Subtree\nCOLOR #ff0000\nDATA"  > annot_itol.txt

# We take subclades of simians from ncbi and add them to itol annotation file
gotree subtree -i !{ncbi} -n "Hominoidea" | gotree stats tips | tail -n+2 | cut -f 4 | awk '{print $0 " range #0097ff Hominoidea"}' >> annot_itol.txt
gotree subtree -i !{ncbi} -n "Cercopithecinae"  | gotree stats tips | tail -n+2 | cut -f 4 | awk '{print $0 " range #7b1616 Cercopithecinae"}' >> annot_itol.txt
gotree subtree -i !{ncbi} -n "Colobinae" | gotree stats tips | tail -n+2 | cut -f 4 | awk '{print $0 " range #ff33ff Colobinae"}' >> annot_itol.txt
gotree subtree -i !{ncbi} -n "Platyrrhini" | gotree stats tips | tail -n+2 | cut -f 4 | awk '{print $0 " range #cf0e0e NewWorldMonkeys"}' >> annot_itol.txt

# outgroup for rooting the simian tree (new world monkeys)
gotree subtree -i !{ncbi} -n "Platyrrhini" | gotree stats tips | tail -n+2 | cut -f 4 > outgroup.txt

# We then reroot and upload fbp and tbe trees to itol
fbpurl=$(gotree reroot outgroup -l outgroup.txt -i simians_fbp.nw | gotree upload itol annot_itol.txt)
tbeurl=$(gotree reroot outgroup -l outgroup.txt -i simians_tbe.nw | gotree upload itol annot_itol.txt)

# We write iTOL tree drawing option files
cat > options_tbe.txt <<EOF
display_mode	1
label_display	1
align_labels	1
ignore_branch_length	1
vertical_shift_factor	0.5
horizontal_scale_factor	1.5
bootstrap_display	1
bootstrap_type	1
bootstrap_symbol	1
bootstrap_slider_min	0.7
bootstrap_slider_max	1
bootstrap_symbol_min	12
bootstrap_symbol_max	12
bootstrap_symbol_color	#c8c7fc
current_font_size	20
line_width	2
inverted	0
EOF

cat > options_fbp.txt <<EOF
display_mode	1
label_display	1
align_labels	1
ignore_branch_length	1
vertical_shift_factor	0.5
horizontal_scale_factor	1.5
bootstrap_display	1
bootstrap_type	1
bootstrap_symbol	1
bootstrap_slider_min	0.7
bootstrap_slider_max	1
bootstrap_symbol_min	12
bootstrap_symbol_max	12
bootstrap_symbol_color	#e9a640
current_font_size	20
line_width	2
EOF

# We finally download resulting images
gotree download itol -i $(basename $fbpurl) -f svg -o !{outprefix}_fbp.svg -c options_fbp.txt
gotree download itol -i $(basename $tbeurl) -f svg -o !{outprefix}_tbe.svg -c options_tbe.txt
