set terminal cairolatex pdf
unset title

set logscale xy

set xlabel "$n$"
set format x "$10^%T$"
set format y "$10^%T$"

f(x)=a*x+b

set output "cov_conv_fit_1_2.tex"
set ylabel "$| \\overline{\\sigma_{12}} - \\E{\\sigma_{12}} |$"
fit f(x) "cov_conv_1_2.dat" u (log($1)):(log(abs($2))) via a,b
plot [1000:] "cov_conv_1_2.dat" u 1:(abs($2)) w lp notitle,\
exp(b)*x**a w l title sprintf("${%.3f} \\times n^{%.3f}$", b, a)

set output "cov_conv_fit_1_3.tex"
set ylabel "$| \\overline{\\sigma_{13}} - \\E{\\sigma_{13}} |$"
fit f(x) "cov_conv_1_3.dat" u (log($1)):(log(abs($2))) via a,b
plot [1000:] "cov_conv_1_3.dat" u 1:(abs($2)) w lp notitle,\
exp(b)*x**a w l title sprintf("${%.3f} \\times n^{%.3f}$", b, a)
