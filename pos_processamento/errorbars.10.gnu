# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
# set output 'errorbars.10.png'
unset errorbars
set grid nopolar
set grid xtics mxtics ytics mytics noztics nomztics nortics nomrtics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set style data lines
set ytics  norangelimit logscale autofreq 
set title "Error on y represented by filledcurve shaded region" 
set xlabel "Time (sec)" 
set ylabel "Rate" 
set logscale y 10
n(x)=1.53**2*x/(5.67+x)**2
GPFUN_n = "n(x)=1.53**2*x/(5.67+x)**2"
DEBUG_TERM_HTIC = 119
DEBUG_TERM_VTIC = 119
Shadecolor = "#80E0A080"
## Last datafile plotted: "silver.dat"
plot 'silver.dat' using 1:($2+$3):($2-$3)       with filledcurve fc rgb Shadecolor title "Shaded error region",     '' using 1:2 smooth mcspline lw 2   title "Monotonic spline through data"
