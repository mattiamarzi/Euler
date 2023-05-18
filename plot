set terminal gif animate delay 100 size 2560,1440
set output 'sol.gif'
stats 'sol.dat' nooutput

set xlabel font ",18"
set ylabel font ",18"
set xtics font ",18"
set ytics font ",18"
set xlabel 'Positions'
set grid
set key font ',18'

do for [i=1:int(STATS_blocks)] {
   set multiplot layout 2,2
   
   set size 0.5,0.5
   set xrange [-0.25:1.25]
   set yrange [-0.25:1.25]
   plot 'sol.dat' u 1:2 index (i-1) with l lw 3 lc 'red' title 'Density'

   set size 0.5,0.5
   set xrange [-0.25:1.25]
   set yrange [-0.25:1.5]
   plot 'sol.dat' u 1:3 index (i-1) with l lw 3 lc 'blue' title 'Velocity'

   set size 0.5,0.5
   set xrange [-0.25:1.25]
   set yrange [-0.25:1.25]
   plot 'sol.dat' u 1:4 index (i-1) with l lw 3 lc 'purple' title 'Pressure'

   set size 0.5,0.5
   set xrange [-0.25:1.25]
   set yrange [1.75:4.25]
   plot 'sol.dat' u 1:5 index (i-1) with l lw 3 lc 'orange' title 'Energy'

   unset multiplot

}