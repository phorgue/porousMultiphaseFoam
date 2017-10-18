set terminal postscript
set size square 0.7,0.7
set output "infiltration.eps"
set xlabel "Depth (m)"
set ylabel "Pressure head (m)"
set ytics 1 nomirror
set label "t=5 000s" at 0.2,-7 rotate by 90
set label "t=15 000s" at 0.38,-7 rotate by 90
set label "t=25 000s" at 0.57,-7 rotate by 90

set xtics nomirror
set key at 1,0.90
plot "sets/5000/acrossFlow_h.gplt" using (0.6-$1):2 every 10 w l title "", \
     "sets/15000/acrossFlow_h.gplt" using (0.6-$1):2 every 10 w l  title "", \
     "sets/25000/acrossFlow_h.gplt" using (0.6-$1):2 every 10 w l title "" 