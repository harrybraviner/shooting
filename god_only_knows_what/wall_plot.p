set terminal pngcairo mono dashed size 800,480
set output 'z=6.5Li_AdS.png'
set title 'Shot from z=6.5 lower sign Lifshitz with delta_x=0.0001, init_e=0.005, [0.25,0.75]'
set xlabel '\rho'
set xrange [1:20]
set yrange [0:1]
plot './output.dat' using 1:4 w l ls 1 title '\zeta',\
'./output0.dat' using 1:4 w l ls 1 title '',\
'./output.dat' using 1:($2*0.5) w l ls 2 title '0.5*\hat{\beta}',\
'./output0.dat' using 1:($2*0.5) w l ls 2 title '',\
'./output.dat' using 1:($8*0.5) w l ls 3 title '0.5*e^{-2\hat{h}}',\
'./output0.dat' using 1:($8*0.5) w l ls 3 title ''
