#!/bin/bash
# $1 - n_urns
# $2 - every
# #$3 - datafile
if [ -z "$2" ]
then
    every=""
else
    every="every $2"
fi

echo "unset key" > all.plot
echo -n "replot file $every u 1:2 w l" >> all.plot

n=$(echo "$1 + 1" | bc)

for i in $(seq 3 $n)
do
    echo -n ", file $every u 1:$i w l" >> all.plot
done
