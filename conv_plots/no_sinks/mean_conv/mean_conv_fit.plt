set terminal cairolatex pdf
unset title

set logscale xy

set xlabel "$n$"
set format x "$10^%T$"
set format y "$10^%T$"

f(x)=a*x+b

set output "mean_conv_fit_1.tex"
set ylabel "$| \\overline{n_1} - \\E{n_1} |$"
fit f(x) "mean_conv.dat" u (log($1)):(log(abs($2 - 140.740740741))) via a,b
plot [100:] "mean_conv.dat" u 1:(abs($2 - 140.740740741)) w lp notitle,\
exp(b)*x**a w l title sprintf("${%.3f} \\times n^{%.3f}$", b, a)

set output "mean_conv_fit_2.tex"
set ylabel "$| \\overline{n_2} - \\E{n_2} |$"
fit f(x) "mean_conv.dat" u (log($1)):(log(abs($3 - 111.111111111))) via a,b
plot [100:] "mean_conv.dat" u 1:(abs($3 - 111.111111111)) w lp notitle,\
exp(b)*x**a w l title sprintf("${%.3f} \\times n^{%.3f}$", b, a)

set output "mean_conv_fit_3.tex"
set ylabel "$| \\overline{n_3} - \\E{n_3} |$"
fit f(x) "mean_conv.dat" u (log($1)):(log(abs($4 - 66.6666666667))) via a,b
plot [100:] "mean_conv.dat" u 1:(abs($4 - 66.6666666667)) w lp notitle,\
exp(b)*x**a w l title sprintf("${%.3f} \\times n^{%.3f}$", b, a)
