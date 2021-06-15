motif=$1

listr=${motif}r.list
ligandr=`head -n 1 $listr`

#prepare npy arrays for assembling
python2 $ATTRACTDIR/dump_coordinates.py $motif.dat $ligandr $motif-preatoms.npy `cat $ssRNATTRACT/data/prepostatoms/$motif.preatoms` --ens 2 $listr
python2 $ATTRACTDIR/dump_coordinates.py $motif.dat $ligandr $motif-postatoms.npy `cat $ssRNATTRACT/data/prepostatoms/$motif.postatoms` --ens 2 $listr