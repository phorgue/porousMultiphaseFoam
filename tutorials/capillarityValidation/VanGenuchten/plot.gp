set terminal postscript enhanced
set size square 0.65,0.65
set xtics 0.2 nomirror

set output "PC_VanGenuchten.eps"
set ytics 10 nomirror
set xlabel "Saturation"
set ylabel "dS/dx (m^{-1})"
set key at 0.6,35
plot [0.001:1.02] 98.0019*(((1/(x**2))-1)**0.5) * (x**3) ti "Analytical",\
     "results/dSdx" using 1:2 ti "Numerical" 

set output "Sb_VanGenuchten.eps"
set ytics 0.2 nomirror
set yrange [0:1.05]
set xlabel "Position (m)"
set ylabel "Saturation"
set key at 0.9,0.9
plot "results/Sb" using 1:2 w l lt 3 lw 2 ti "Numerical"
