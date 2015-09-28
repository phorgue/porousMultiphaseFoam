set terminal postscript enhanced
set size square 0.65,0.65
set output "BL_Brooks_gravity.eps"
set xrange [0 : 1]
set xlabel "Position (m)"
set ylabel "Saturation"
set yrange [0:0.7]
set ytics 0.1
set key at 0.80,0.65
plot "postProcessing/sets/5000/acrossFlow_Sb.gplt"  using (1-$1):2 every 15 title "t = 5000s " ps 0.75, \
     "analytical/case2/5000s"  using 1:2 with lines title "",\
     "postProcessing/sets/10000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 10000s" ps 0.75, \
     "analytical/case2/10000s" using 1:2 with lines title "",\
     "postProcessing/sets/15000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 15000s" ps 0.75, \
     "analytical/case2/15000s" using 1:2 with lines title ""
