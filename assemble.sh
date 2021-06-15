LANG=en_US

d=`pwd`

# first and last frag to assemble
f1=$1
f2=$2

# maximal number of poses to consider per fragment
maxstruc=$3

# maximal number of chains to sample
nchains=$4

# maximal geometric mean of poses ranks in each chain
# use meanrk=maxstruc if your poses ranking is not discriminant
meanrk=$5

# The maximal overlap RMSD between consecutive poses
# will be iteratively increased from $start to $end values
# until a number $nchains of chains are found.
start=$6  #recommanded 1.0
end=$7    #recommanded 3.0

a="frag["$f1"-"$f2"]-preatoms.npy"
b="frag["$f1"-"$f2"]-postatoms.npy"
nfrag=$(($f2-$f1+1))

# optional: give the lRMSD of each pose for each fragment
l="frag["$f1"-"$f2"].rmsd"

echo "overlap Nchains inf6  inf5  inf4 min_RMSD" >> frag$f1-${f2}.results
echo "-----------------------------------------" >> frag$f1-${f2}.results

# assemble chains with incremental overlap cutoff
# until $nchains chains are sampled

set -u -e
for rmsd in `seq $start 0.1 $end` ; do
    name=frag$f1-${f2}_incr-ov\_${rmsd}A
    echo ""
    echo "----------------------------- overlap $rmsd --------------------------------"

    # Create the graph of poses connectivity
    graph=f$f1-$f2\_$rmsd\A_$maxstruc\poses
    # !!! does not recompute the graph if file already exist !!!
    # if the file is wrong, delete it before re-runing the script

    #if [ ! -s $graph.json ];then
    #  $ssRNATTRACT/scripts/assemble.py --nfrag $nfrag --rmsd $rmsd --maxstruct $maxstruc \
    #  --preatoms frag[$f1-$f2]-preatoms.npy \
    #  --postatoms frag[$f1-$f2]-postatoms.npy \
     # > $graph.json
    #fi

    python2 $ssRNATTRACT/scripts/connect.py $nfrag $rmsd $maxstruc \
      --preatoms frag[$f1-$f2]-preatoms.npy \
      --postatoms frag[$f1-$f2]-postatoms.npy \
      > $graph.json
    
    # Walk the graph to extract chains of connected poses
    chains=$graph\_${meanrk}meanrk
    $ssRNATTRACT/scripts/make-chains.py --graph $graph.json --meanrank $meanrk \
      --preatoms frag[$f1-$f2]-preatoms.npy \
      --postatoms frag[$f1-$f2]-postatoms.npy \
      --rmsd frag[$f1-$f2].rmsd \
      > $chains.chains

    Nchains=`awk '$1=="#indices"{i+=1}END{print i}' $chains.chains`

    if [ $Nchains -gt $nchains ]  ;then
      mv $chains.chains $d/
      mv $graph.json $d/
      break
    fi
done

# count number of chains
# compute number of +/- correct chains (RMSD <6A, <5A, <4A )
# RMSD is approximated by averaging the RMSD of the fragments
awk -v r=$rmsd 'BEGIN{i=0;j=0; k=0;m=10}
    NR>1&&$2<6{i+=1}
    NR>1&&$2<5{j+=1}
    NR>1&&$2<4{k+=1}
    $2<m{m=$2}
    END{if (NR>1) printf "%.1f \t %i \t %3.0f \t %3.0f \t %3.0f \t %3.2f \n",
    r, NR-1, 100*i/(NR-1), 100*j/(NR-1),100*k/(NR-1), m }' $chains.chains >> frag$f1-${f2}.results