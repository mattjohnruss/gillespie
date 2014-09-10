rm -f mean_conv.dat
rm -f cov_conv.dat

for n in $@
do
    rm -rf $n
    mkdir $n
    cd $n
    ~/phd/gillespie/write_cmds $n
    parallel --progress ~/phd/gillespie/main -n 3 -a 1 -b 0.5 -c 100 -d 1.5 -t 100 -l -m 3 -v 2 -f {}.dat :::: cmds
    ~/phd/gillespie/stats $n ""
    rm cmds
    cd ..
    echo "$n $(tail $n/_mean.dat -n1 | cut -d' ' -f 2-)" >> mean_conv.dat
    echo -n "$n $(sed -n 6p $n/_covariance.dat | cut -d' ' -f 2-3)" >> cov_conv.dat
    echo " $(sed -n 7p $n/_covariance.dat | cut -d' ' -f 3)" >> cov_conv.dat
done
