set terminal postscript eps enhanced color 
set encoding utf8
set output 'plotRandom.eps'
#set key left top
set key rmargin right vertical
set xlabel "t"
filename = 'caso1/solutionRK'
filename2 = "caso1/data.dat"

set autoscale

set ylabel 'V'

set grid

set style data points
set for [i=6:300:9] linetype i lc rgb "dark-orange" 

ii=1
ifi=5000
plot for [i=ii:ifi] filename.i.'.dat'  using 1:2 notitle pt 7 ps 0.5, filename2 using 2:1 title "Experimental" pt 2 ps 1.2
