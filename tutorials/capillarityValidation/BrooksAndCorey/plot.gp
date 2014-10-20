set terminal postscript enhanced
set size square 0.65,0.65
set xtics 0.2 nomirror

set output "PC_BrooksAndCorey.eps"
set ytics 5 nomirror
set xlabel "Saturation"
set ylabel "dS/dx (m^{-1})"
set key at 0.8,21
plot [0.001:1.02] 9800.19/(1000*0.5*(x**(-1.5))) ti "Analytical",\
     'results/dSdx' using 1:2 every 5 ti "Numerical"

set output "Sb_BrooksAndCorey.eps"
set ytics 0.2 nomirror
set yrange [0:1.05]
set xlabel "Position (m)"
set ylabel "Saturation"
set key at 0.9,0.9
plot "results/Sb" using 1:2 w l lt 3 lw 2 ti "Numerical"
