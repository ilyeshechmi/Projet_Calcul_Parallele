set term pngcairo size 1000,600
set output 'results.dat_compare.png'
set xlabel 'x'
set ylabel 'u(x,t)'
set title 'Comparaison Rusanov vs solution exacte à t=0.0000'
set grid
plot 'results.dat_compare.dat' using 1:2 with lines lw 2 lc 'blue' title 'u numérique', \
     'results.dat_compare.dat' using 1:3 with lines lw 2 lc 'red' title 'u exacte'
