set term gif animate delay 6
set output 'disk.gif'
do for[i=0:200:1] {
  file = sprintf("dat/accd_n%05d.dat", i)
  #time = sprintf("t=%d-%d[sec]", i*200, (i+1)*200)
  set xrange [0:1]
  set yrange [0:1]
  #set title time
  plot file u 1:2 w l
}
