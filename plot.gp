set term png size 1900,1000 enhanced font "Terminal,10"

set grid

set auto x

set key left top

set title "Evaluation of recursive calls for graph named dh11_edges.csv"

set xlabel "Algorithm variants"
set ylabel "number of recursive calls

set style data histogram
set style fill solid border -1
set boxwidth 0.9

set xtic rotate by -45 scale 0

set multiplot layout 3, 3 rowsfirst

set yrange [0:3000]

plot "fichier.dat" u 2:xtic(1) t ""

