set terminal pdf
set output 'plot1d.pdf'

#set terminal latex

#set key left top
set key rmargin right vertical
set ylabel 'y(x)'
set xlabel 'x'
#set format x "%.3g
#set tics font ", 10"

set autoscale
#set xrange [0:]
#set yrange [0:]

set grid

filename = "caso1/alphaTriC.dat"
linha = 1
coluna = 2

set style data points
plot filename using linha:coluna title "y(x)" pt 7 ps 0.4
#plot filename using linha:coluna title "Numérico" pointtype 3 

