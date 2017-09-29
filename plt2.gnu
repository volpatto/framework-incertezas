set terminal pdf
set output 'solution.pdf'

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

filename = "solution1.dat"
filename2 = "data2.dat"
linha = 1
coluna = 2

y0 = 13.18
R = 4.40e5
C = 3.5196e-4
#C = 3.5084000000000000E-004
f(x) = y0*exp(-x/(R*C))

set style data points
plot f(x) title "Analítico" lw 2, filename using linha:coluna title "Numérico" pointtype 3, filename2 using 2:1 title "Experimental" pointtype 1  
#plot filename using linha:coluna title "Numérico" pointtype 3 

