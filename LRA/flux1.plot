#!/usr/bin/env gnuplot

set terminal pdf
set output "LRA_flux1.pdf" 
set palette defined (0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
set xrange [0.0:165.0]
set yrange [0.0:165.0]
set xlabel "Core x dimeimsion [-]"
set ylabel "Core y dimension [-]"
set title "Group 1 Flux Distribution"
splot 'flux1.dat' matrix with image
