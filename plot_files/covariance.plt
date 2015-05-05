set size ratio 1
#set palette gray negative
set xrange[0.5:n+0.5]
set yrange[0.5:n+0.5]
#set cbrange[-0.2:0.8]
set xtics 1,1,n
set ytics 1,1,n
set xtics offset -0.5,0

plot file u ($1+1):($2+1):3 index i matrix w image noti
