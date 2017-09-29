# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'errorbars.4.png'
set style data lines
set title "error represented by yerrorbars" 
set xlabel "Resistance [Ohm]" 
set ylabel "Power [W]" 
n(x)=1.53**2*x/(5.67+x)**2
GPFUN_n = "n(x)=1.53**2*x/(5.67+x)**2"
DEBUG_TERM_HTIC = 119
DEBUG_TERM_VTIC = 119
## Last datafile plotted: "battery.dat"
plot [0:50] "battery.dat" u 1:2:4 t "Power" w yerr, n(x) t "Theory" w lines
