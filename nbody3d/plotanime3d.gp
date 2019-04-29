set term gif animate delay 6
set output '3d_01.gif'
do for[i=0:100:1] {
  file = sprintf("nbody_%05d.dat", i)
  #time = sprintf("t=%d-%d[sec]", i*200, (i+1)*200)
  #set title time
  splot [-2:2] [-2:2] [-2:2] file u 1:2:3 pt 7 ps 1
}
