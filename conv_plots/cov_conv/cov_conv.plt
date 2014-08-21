set terminal cairolatex pdf
unset title
unset key

set xlabel "$n$"
set format x "$\\scriptstyle{%.0t \\times 10^%T}$"

set output "cov_conv_1_2.tex"
set ylabel "$\\overline{\\sigma_{12}}$"
plot "cov_conv_1_2.dat" w lp

set output "cov_conv_1_3.tex"
set ylabel "$\\overline{\\sigma_{13}}$"
plot "cov_conv_1_3.dat" w lp
