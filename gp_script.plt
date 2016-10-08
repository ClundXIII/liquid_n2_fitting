# gnuplot script for liquid_n2_fitting
# (c) by Simon Michalke, 2016

set terminal postscript enhanced color
set output "fit_plot.eps"

set grid
plot "data.txt" with lines title "Measured Pressure", "fit_plot.txt" with lines title "Simulated Pressure"
