program plot

  implicit none
      
  real(4) :: fx, x
  integer :: i

  open (112,file='04_gnuxy')

  write(112,*) 'set grid'
  write(112,*) 'set xzeroaxis'
  write(112,*) 'set yzeroaxis'
  write(112,*) 'set border 0          # remove frame'
  write(112,*) 'set xtics axis        # place tics on axis rather than on border'
  write(112,*) 'set ytics axis'
  write(112,*) 'set ticscale 0        # [optional] labels only, no tics'
  write(112,*) 'set xtics add ("" 0)  # suppress origin label that lies on top of axis'
  write(112,*) 'set ytics add ("" 0)  # suppress origin label that lies on top of axis'
  write(112,*) ''
  write(112,*) '# if arrows are wanted only in the positive direction'
  write(112,*) 'set arrow 1 from 0,0 to graph 1, first 0 filled head'
  write(112,*) 'set arrow 2 from 0,0 to first 0, graph 1 filled head'
  write(112,*) ''
  write(112,*) '# if arrows in both directions from the origin are wanted'
  write(112,*) 'set arrow 3 from 0,0 to graph 0, first 0 filled head'
  write(112,*) 'set arrow 4 from 0,0 to first 0, graph 0 filled head'
  write(112,*) ''
  write(112,*) 'set nokey'
  write(112,*) 'set xrange [-4:4]'
  write(112,*) 'plot "04_dataxy_1" using 1:2 with lines lt rgb "blue"'
  write(112,*) 'pause -1'

  close(112)

  ! Generate x-y pairs for the graph

  open (112,file='04_dataxy_1')
  do i=-40,40
    x  = .1*i
    fx = x**3+x-10.0
    write(112,*) x, fx
  end do
  close(112)

  print *, ' Hit the Return (Enter) key to continue'

  call system ('gnuplot 04_gnuxy')
  stop

end
