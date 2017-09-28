#set terminal epslatex color
set terminal postscript eps enhanced color 
set output 'alpha.eps'

#set terminal latex

#set key left top
#set key rmargin right vertical
set ylabel "{/Symbol a}(R)"
set xlabel 'R'
#set format x "%.3g
#set tics font ", 10"

set autoscale
#set xrange [0:]
set yrange [0:1.05]

set grid

filename = "caso1/alphaTriR.dat"
linha = 1
coluna = 2

set style data points
#plot filename using linha:coluna title "{/Symbol a}" pt 7 ps 0.4
plot filename using linha:coluna notitle pt 7 ps 0.4
#plot filename using linha:coluna title "Numérico" pointtype 3 

