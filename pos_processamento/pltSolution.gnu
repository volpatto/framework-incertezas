set terminal postscript eps enhanced color 
set encoding utf8
set output 'calibracao.eps'

#set terminal latex

#set key left top
set key rmargin right vertical
set xlabel 'Tempo (s)'
set ylabel 'V'
#set format x "%.3g
#set tics font ", 10"

set autoscale
set xrange [0:]
set yrange [0:]

set grid

filename = "caso1/solutionRK0.dat"
linha = 1
coluna = 2

filename2 = "caso1/data.dat"

y0 = 13.18
R = 440.0e3
C = 3.52051E-04
f(x) = y0*exp(-x/(R*C))

set style data points
plot f(x) title "Analítico" lw 2, filename using linha:coluna title "Numérico" pt 7 ps 0.5, filename2 using coluna:linha title "Experimental" pt 2 ps 1.2
