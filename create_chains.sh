#!/bin/bash
LANG=en_US
d=`pwd`

chains=$1 #f2-4_1.5A_1000poses_1000meanrk.chains
name=${chains%%.chains}

# first and last frag
f1=$2
f2=$3

######################################################
echo "merge fragments"
######################################################
# average the coordinates of the overlapping parts of fragments
$ssRNATTRACT/scripts/merge-npy.py --chainfile $name.chains --outp $name.npy \
  --npyfiles frag[$f1-$f2]-aa.npy \
  --preatoms frag[$f1-$f2]-aa.preatoms \
  --postatoms frag[$f1-$f2]-aa.postatoms

######################################################
echo "convert from reduced to all-atom coordinates"
######################################################
# first and last residues
g1=$f1
g2=$(($f2+2))

# get the name of residues in bound ligand
awk -v a=$g1 -v b=$g2 '$5>b{exit} $5>=a && $3=="P"{print $4}' \
    boundfrag/ligand-aa.pdb > resnames.txt

cat /dev/null > libraries.txt

# get the list of monomer libraries to be used for each residue
for res in `cat resnames.txt`; do
    echo $ssRNATTRACT/data/monomers_library/$res.npy >> libraries.txt
done

# Each merged residue is replaced by the closest monomer in library
$ssRNATTRACT/scripts/convert-merged-aa.py $name.npy $name-aa.npy \
    --lib `cat libraries.txt` --noext

######################################################
echo "Compute RMSD toward bound ligand"
######################################################
awk -v a=$g1 -v b=$g2 '$5<=b && $5>=a' boundfrag/ligand-aa.pdb > /tmp/lig.pdb
$ssRNATTRACT/scripts/rmsdnpy.py $name-aa.npy --ref /tmp/lig.pdb > $name.rmsd

echo "----------- best solutions ------------"
echo "rank rmsd"
echo "---------------------------------------"
sort -nk2 $name.rmsd|head

######################################################
echo "convert 100 best-rmsd solutions into PDB"
######################################################
n=100

sort -nk2 $name.rmsd|awk -v x=$n '{print $1}NR==x{exit}' > /tmp/sel
$ssRNATTRACT/scripts/select-npy.py $name-aa.npy $name-aa-best$n.npy --structures `cat /tmp/sel`
$ssRNATTRACT/scripts/npy2pdb.py $name-aa-best$n.npy /tmp/lig.pdb > $name-aa-best$n.pdb
