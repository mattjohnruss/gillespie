if [ "$#" -ne 10 ];
then
    echo "Usage: $0 N L Pe/ep Da/ep^2 n_init t_max prefix interval last_only num_runs"
    exit
fi

N=$1
L=$2
Pe_ep=$3
Da_ep2=$4
n_init=$5
t_max=$6
prefix=$7
interval=$8
last_only=$9
num_runs=${10}

ep=$(echo "scale = 8; 1/($N + 1)" | bc)
Pe=$(echo "scale = 8; $Pe_ep*$ep" | bc)
Da=$(echo "scale = 8; $Da_ep2*$ep*$ep" | bc)
dx=$(echo "scale = 8; ($N + 1)/$L" | bc)
p=$(echo "scale = 8; (1/$dx + $Pe/2)/($Da*$dx)" | bc)
q=$(echo "scale = 8; (1/$dx - $Pe/2)/($Da*$dx)" | bc)

echo "N  = $N"
echo "L  = $L"
echo "dx = $dx"
echo "p  = $p"
echo "q  = $q"
echo "ep = $ep"
echo "num_runs = $num_runs"

args="-n $L -e $n_init -p $p -q $q -t $t_max -f $prefix{1}.dat"

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
