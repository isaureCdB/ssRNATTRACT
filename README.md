# ssRNATTRACT

## Requirements
python2 and numpy
python3 and numpy

add to your ~/.bashrc:
export ssRNATTRACT='full path to git directory'

have frag[i].rmsd for each fragment


## Example of usage

1. Convert data from ATTRACT (RNA 3nt sequence AAUC)\

$ssRNATTRACT/prep_assembly.sh AAU frag1\
$ssRNATTRACT/prep_assembly.sh AUC frag2

2. Assemble fragments\

args: [first frag index] [last frag index] [Nb_poses] [Nb_chains] [meanrk] [min_overlap rmsd] [max_overlap rmsd]\
exp: $ssRNATTRACT/assemble.sh 1 2 10000 1000 10000 1.2 2.5

3. Convert chains of fragments into full all-atom models\

args: [chains] [first frag index] [last frag index]\
expl: $ssRNATTRACT/create_chains.sh f1-2_1.5A_10000poses.chains 1 2 
