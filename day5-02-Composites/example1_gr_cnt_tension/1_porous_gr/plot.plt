set term png  enhanced font ',16' linewidth 1
set encoding iso_8859_1
#set title "Time  ps"
set output 'strain_stress.png'
set style fill solid
set multiplot layout 1,1 title ""

#set yrange [-600:400]
#set xlabel 'Time/ns'
set xtics  0.04
set ytics  10
set ylabel 'Stress (MPa)'
set xlabel 'Strain'
set key l 

#set ytics  200
p 'Stress-strain_porous_gr.txt' u ($1):($2) w l lc 1 lw 2 lt 1  title ''

