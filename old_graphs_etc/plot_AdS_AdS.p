set term png size 800,600
set output 'AdS_AdS_plot.png'
set xrange[1:20]
set xlabel '\rho'
set yrange[-0.3:1]
set title 'Shot from \zeta=0.7 AdS space. delta_x=0.01, init_e=0.01, hit_e=0.005, [1.25,1.75]'
plot "./output.dat" using 1:4 w l lw 1 title '\zeta', "./output.dat" using 1:5 w l lw 2 title '\partial_{\rho} \zeta', "./output.dat" using 1:8 with l lw 3 title 'e^{-2\hat{h}}', "./output.dat" using 1:9 w l lw 4 title '\partial_{\rho} e^{-2\hat{h}}'
