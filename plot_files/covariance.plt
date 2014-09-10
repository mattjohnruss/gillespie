set size ratio 1
#set palette gray negative
set xrange[-0.5:9.5]
set yrange[-0.5:9.5]
#set cbrange[-0.2:0.8]
set xtics 0,1,9
set ytics 0,1,9
set xtics offset -0.5,0

plot file index i matrix w image noti
