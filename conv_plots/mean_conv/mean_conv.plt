set terminal cairolatex pdf
unset title
unset key

set xlabel "$n$"
set format x "$\\scriptstyle{%.0t \\times 10^%T}$"

set output "mean_conv_1.tex"
set ylabel "$\\overline{n_1}$"
plot "mean_conv_1.dat" w lp

set output "mean_conv_2.tex"
set ylabel "$\\overline{n_2}$"
plot "mean_conv_2.dat" w lp

set output "mean_conv_3.tex"
set ylabel "$\\overline{n_3}$"
plot "mean_conv_3.dat" w lp

