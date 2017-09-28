set terminal postscript eps enhanced color 
set encoding utf8
set output 'alpha3d.eps'

set xlabel "R"
set ylabel "C"
set zlabel "{/Symbol a}"

set autoscale

set grid
set ticslevel 0
#set samples 25, 25
#set isosamples 50, 50
set hidden3d
set pm3d
set dgrid3d 30,30
#set palette model CMY rgbformulae 7,5,15
#set view 60,45
#splot "caso1/alpha3d.dat" using 1:2:3 w pm3d notitle
splot "caso1/alpha3d.dat" using 1:2:3 w lines notitle
