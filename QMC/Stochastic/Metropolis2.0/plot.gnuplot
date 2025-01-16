set title 'Positions simulated in the atom'

set xlabel 'X axis'
set ylabel 'Y axis'
set zlabel 'Z axis'

set grid

splot 'Segment_walk' u 1:2:3 w p pt 6 ps 0.1 lc 6 t 'Trajectory', 'Origin' u 1:2:3 w p pt 7 ps 2 lc rgb 'black' t 'Nucleus'
pause -1
