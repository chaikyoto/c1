set term gif animate delay 6
set output 'v_02.gif'
do for[i=0:100:1] {
  file = sprintf("nbody_%05d.dat", i)
  #time = sprintf("t=%d-%d[sec]", i*200, (i+1)*200)
  set xrange [-100:100]
  set yrange [-100:100]
  #set title time
  plot file u 1:2 pt 7 ps 1
}
