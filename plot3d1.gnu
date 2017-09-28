set terminal pdf
set output 'solucao3d.pdf'

set xlabel "t"
set ylabel "x"
set zlabel "N"

set autoscale

set grid
set ticslevel 0
#set samples 25, 25
#set isosamples 50, 50
set pm3d
set hidden3d
set palette model CMY rgbformulae 7,5,15
set view 60,45
splot "caso1/fuzzy.dat" using 2:3:1 w pm3d notitle
