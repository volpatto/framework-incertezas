set terminal postscript eps enhanced color 
set encoding utf8
set output 'solutionFuzzy.eps'
#set key left top
set key rmargin right vertical
set xlabel "t"
filename = 'caso1/fuzzy'
filename2 = "caso1/fuzzy"

set autoscale

set ylabel 'V'

set grid

set style data points
set for [i=6:300:9] linetype i lc rgb "dark-orange" 

ii=1
ifi=11
plot for [i=ii:ifi] filename.i.'.dat' using 1:5 notitle pt 7 ps 0.5, for [i=ii:ifi] filename2.i.'.dat' using 1:8 notitle pt 2 ps 0.5
