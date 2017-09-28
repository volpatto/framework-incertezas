set terminal pdf
set output 'fit.pdf'

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

filename = "data2.dat"
linha = 2
coluna = 1

a = 1
b = 1
f(x) = a*exp(-x/b)
FIT_LIMIT = 1e-10
fit f(x) filename using linha:coluna via a,b

set style data points
plot f(x) title "Analítico", filename using linha:coluna title "Experimental" pointtype 3 

