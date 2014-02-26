#!/bin/bash
echo "cat output*.dat | awk '"
echo "BEGIN {"
echo "for (i = 0; i <= $1; i++)"
echo "for (j = 1; j <= $2; j++)"
echo "total[i][j] = 0"

echo "} {"

echo "for (j = 1; j <= $2; j++)"
echo "total["'$1'"][j] += "'$(j+1) }'

echo "END {"
echo "for (i = 0; i <= $1; i++)"
echo -n "print i"
for j in $(seq 1 $2)
do
    echo -n ", total[i][$j]/$3"
done
echo ""
echo "}"
echo "'"
