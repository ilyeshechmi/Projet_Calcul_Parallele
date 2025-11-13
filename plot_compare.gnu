set term pngcairo size 1000,600
set output 'resultats/results_compare.png'
set xlabel 'x'
set ylabel 'u(x,t)'
set title 'Comparaison Rusanov vs solution exacte à t=10.0000'
set grid
plot 'resultats/results_compare.dat' using 1:2 with lines lw 2 lc 'blue' title 'u numérique', \
     'resultats/results_compare.dat' using 1:3 with lines lw 2 lc 'red' title 'u exacte'
