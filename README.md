# ssRNATTRACT
#################################
# requirements
#################################
python2 and numpy
python3 and numpy

add to your ~/.bashrc:
export ssRNATTRACT='full path to consero directory'

#################################
# Example of usage
#################################

# convert data from ATTRACT (RNA 3nt sequence AAUC)
$ssRNATTRACT/prep_assembly.sh AAU
$ssRNATTRACT/prep_assembly.sh AUC

for n in "-preatoms.npy" "-postatoms.npy" ".lrmsd"; do
    ln -s AAU$n frag1$n
    ln -s AUC$n frag2$n
done

# assemble fragments
# args: [first frag index] [last frag index] [Nb_poses] [Nb_chains] [meanrk] [min_overlap rmsd] [max_overlap rmsd]
$ssRNATTRACT/assemble.sh 1 2 10000 1000 10000 1.2 2.5

# convert chains of fragments into full all-atom models
# args: [chains] [first frag index] [last frag index]
$ssRNATTRACT/create_chains.sh f1-2_1.5A_10000poses.chains 1 2 