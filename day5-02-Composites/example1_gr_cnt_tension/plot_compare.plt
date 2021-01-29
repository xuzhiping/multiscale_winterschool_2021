set term png  enhanced font ',16' linewidth 1
set encoding iso_8859_1
set output 'strain_stress.png'
set style fill solid
set multiplot layout 1,1 title ""

set xtics  0.04
set ytics  10
set ylabel 'Stress (MPa)'
set xlabel 'Strain'
set key l 

p './2_gr_cnt/Stress-strain_gr_cnt.txt' u ($1):($2) w l lc 1 lw 2 lt 1  title '2D Gr-CNT','./1_porous_gr/Stress-strain_porous_gr.txt' u ($1):($2) w l lc 2 lw 2 lt 1  title 'Porous G',


