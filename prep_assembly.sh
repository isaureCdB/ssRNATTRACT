motif=$1 #sequence of the trinucleotide
frag=$2 #index of the fragment

listr=${motif}r.list
ligandr=`head -n 1 $listr`

#prepare npy arrays for assembling
if [ ! -f $motif-preatoms.npy ];then
  python2 $ATTRACTDIR/dump_coordinates.py $motif.dat $ligandr $motif-preatoms.npy `cat $ssRNATTRACT/data/prepostatoms/$motif.preatoms` --ens 2 $listr
  python2 $ATTRACTDIR/dump_coordinates.py $motif.dat $ligandr $motif-postatoms.npy `cat $ssRNATTRACT/data/prepostatoms/$motif.postatoms` --ens 2 $listr
fi

ln -s $motif-preatoms.npy frag$frag-preatoms.npy
ln -s $motif-postatoms.npy frag$frag-postatoms.npy
