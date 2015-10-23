if [ "$#" -ne 11 ];
then
    echo "Usage: $0 N M alpha Pe/ep Da/ep^3 n_init t_max prefix interval last_only num_runs"
    exit
fi

N=$1
M=$2
alpha=$3
Pe_ep=$4
Da_ep3=$5
n_init=$6
t_max=$7
prefix=$8
interval=$9
last_only=${10}
num_runs=${11}

ep=$(echo "scale = 8; 1/($N + 1)" | bc)
Pe=$(echo "scale = 8; $Pe_ep*$ep" | bc)
Da=$(echo "scale = 8; $Da_ep3*$ep*$ep*$ep" | bc)
dx=$(echo "scale = 8; ($N + 1)/$M" | bc)
p=$(echo "scale = 8; (1/$dx + $Pe/2)/($Da*$dx)" | bc)
q=$(echo "scale = 8; (1/$dx - $Pe/2)/($Da*$dx)" | bc)

#n_init=$(echo "scale = 8; $alpha" | bc)

echo "N     = $N"
echo "M     = $M"
echo "alpha = $alpha"
echo "dx    = $dx"
echo "p     = $p"
echo "q     = $q"
echo "ep    = $ep"
echo "num_runs = $num_runs"

args="-n $M --n_init $n_init -a $alpha -b $p -p $p -q $q -t $t_max -f $prefix{1}.dat"

if [ "$last_only" -eq 1 ];
then
    args+=" -l"
else
    args+=" -i $interval"
fi

#~/Documents/phd/gillespie/main $args

parallel --progress ~/Documents/phd/gillespie/main "$args" ::: `seq -w 1 $num_runs`

if [ ! -z $prefix ];
then
    ~/Documents/phd/gillespie/stats $num_runs $prefix
else
    ~/Documents/phd/gillespie/stats $num_runs ""
fi
