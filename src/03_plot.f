      program plot
c
c    Program to provide plots of Sin(x)
c
      implicit none
      character label*150
      real x
      integer i
      character xlabel*32,ylabel*32,title*32
      real fx
c
c   label   -   Character string 
c   xlabel  -   Contains a label for the x-axis
c   ylabel  -   Contains a label for the y-axis
c   title   -   Contains a title for the plot
c
c   Drive a separate true graphics program (gnuplot)
c
c   First set up the command file for gnuplot
c
      xlabel='''x'''
      ylabel='''y'''
      title="'sin(x)'"
      open (112,file='03_gnuxy')
c
c     write(112,*) 'set term wxt size 800, 800'
c
      label='set xlabel '//xlabel
      write(112,*)label
      write(112,*)'set xrange [0:6]'
      label='set ylabel '//ylabel
      write(112,*)label
      write(112,*)'set yrange [-1.2:1.2]'
      label='plot "03_dataxy" using 1:2 title '//title
      label=trim(label)//' with lines lt rgb "red"'
      write(112,*) label
      write (112,*) 'pause -1'
      close(112)
c
c   Generate x-y pairs for the graph
c
      open (112,file='03_dataxy')
      do 100 i=0,60
         x=.1*i
         fx=sin(x)
         write(112,*) x,fx
  100 continue
      close(112)
c
      print *, ' Hit the Return (Enter) key to continue'
c
c   Tell the system to run the program gnuplot
c   This call works on either IBM RS6000 or Sun, but is not part of
c   the Fortran standard.
c   Comment out the line if you aren't at a terminal with graphics
c
      call system ('gnuplot 03_gnuxy')
      stop
      end
