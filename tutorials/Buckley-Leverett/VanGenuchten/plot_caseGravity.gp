set terminal postscript
set size square 0.65,0.65
set output "BL_VanGenuchten_gravity.eps"
set xrange [0 : 1]
set xlabel "Position(m)"
set ylabel "Saturation"
set yrange [0:1]
set ytics 0.2

set key at 0.80,0.95
plot "postProcessing/sets/10000/acrossFlow_Sb.gplt"  using (1-$1):2 every 15 title "t = 10000s " ps 0.75, \
     "analytical/case2/10000s"  using 1:2 with lines title "",\
     "postProcessing/sets/20000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 20000s" ps 0.75, \
     "analytical/case2/20000s" using 1:2 with lines title "",\
     "postProcessing/sets/30000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 30000s" ps 0.75, \
     "analytical/case2/30000s" using 1:2 with lines title ""
