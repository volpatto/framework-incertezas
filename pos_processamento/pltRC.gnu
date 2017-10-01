set terminal pdf
set output 'RC.pdf'

#set terminal latex

#set key left top
set nokey
#set key rmargin right vertical
set xlabel 'R'
set ylabel 'C'
set format y "%1.3e
set tics font ", 9"

set autoscale
#set xrange [0:]
#set yrange [0:]

set grid

filename = "RC.dat"
linha = 1
coluna = 2

set style data points
#plot filename using linha:coluna title "RC" 
plot filename using linha:coluna
#plot filename using linha:coluna title "Numérico" pointtype 3 

