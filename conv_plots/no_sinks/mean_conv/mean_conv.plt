set terminal cairolatex pdf
unset title
unset key

set xlabel "$n$"
set format x "$\\scriptstyle{%.0t \\times 10^%T}$"

set output "mean_conv_1.tex"
set ylabel "$\\overline{n_1}$"
plot "mean_conv.dat" u 1:2 w lp, 140.740740741

set output "mean_conv_2.tex"
set ylabel "$\\overline{n_2}$"
plot "mean_conv.dat" u 1:3 w lp, 111.111111111

set output "mean_conv_3.tex"
set ylabel "$\\overline{n_3}$"
plot "mean_conv.dat" u 1:4 w lp, 66.6666666667

