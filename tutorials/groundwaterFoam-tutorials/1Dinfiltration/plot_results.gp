set terminal postscript
set size square 0.7,0.7
set output "infiltration.eps"
set xlabel "Depth (m)"
set ylabel "Pressure head (m)"
set ytics 1 nomirror
set label "t=17 000s" at 0.2,-7 rotate by 90
set label "t=50 000s" at 0.38,-7 rotate by 90
set label "t=95 000s" at 0.57,-7 rotate by 90

set xtics nomirror
set key at 1,0.90
plot "sets/17000/acrossFlow_h.gplt"  using (0.6-$1):2 every 10 title "" ps 0.75, \
     "sets/50000/acrossFlow_h.gplt" using (0.6-$1):2 every 10  title "" ps 0.75, \
     "sets/95000/acrossFlow_h.gplt" using (0.6-$1):2 every 10 title "" ps 0.75, \
     "reference/17ks.csv" with lines title "" , \
     "reference/50ks.csv" with lines title "" , \
     "reference/95ks.csv" with lines title "" ,