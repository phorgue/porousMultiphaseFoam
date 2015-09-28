set terminal postscript
set size square 0.65,0.65
set output "BL_VanGenuchten_noGravity.eps"
set xrange [0 : 1]
set xlabel "Position(m)"
set ylabel "Saturation"
set yrange [0:1]
set ytics 0.2 nomirror
set xtics nomirror
set key at 1,0.90
plot "postProcessing/sets/6000/acrossFlow_Sb.gplt"  using (1-$1):2 every 15 title "t = 6000s " ps 0.75, \
     "analytical/case1/6000s"  using 1:2 with lines title "",\
     "postProcessing/sets/12000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 12000s" ps 0.75, \
     "analytical/case1/12000s" using 1:2 with lines title "",\
     "postProcessing/sets/20000/acrossFlow_Sb.gplt" using (1-$1):2 every 15 title "t = 20000s" ps 0.75, \
     "analytical/case1/20000s" using 1:2 with lines title ""

