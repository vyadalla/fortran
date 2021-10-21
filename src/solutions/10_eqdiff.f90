program eqdiff

  implicit none

  integer, parameter :: nt  = 200
  integer, parameter :: st = 1

  integer :: n
  real(8) :: pi, A, delta

  real(8), dimension(nt+st) :: U, U1, U2, UU
  real(8), dimension(nt+st) :: E, E1, E2

  open(100,file="10_results.txt")

  pi=2d0*dacos(0)
  A=1d0
  delta=pi/50d0

  U(:)=0d0
  E(:)=0d0
  U1(:)=0d0
  E1(:)=0d0
  U2(:)=0d0
  E2(:)=0d0
  UU(:)=0d0

  write(*,*) "Starting the loop..."
  write(100,6000) st,U(1),U1(1),U2(1),UU(1),E(1),E1(1),E2(1)

  do n=st,nt
    U(n+1)=U(n)+delta*A*dcos(n*delta)
    U1(n+1)=U1(n)+delta*A*dcos((n+0.5d0)*delta)
    U2(n+1)=U2(n)+delta*A*dcos((n+1d0)*delta)

    UU(n+1)=A*dsin(n*delta)

    E(n+1)=dabs(UU(n+1)-U(n+1))
    E1(n+1)=dabs(UU(n+1)-U1(n+1))
    E2(n+1)=dabs(UU(n+1)-U2(n+1))

    write(100,6000) (n+1),U(n+1),U1(n+1),U2(n+1),UU(n+1),E(n+1),E1(n+1),E2(n+1)
  end do

  write(*,5000) "End of the loop: (",nt,") iterations!"

  open (112,file='10_gnuxy')

  write(112,*) "set multiplot layout 2,1 rowsfirst"
  write(112,*) "set xlabel 'N'"
  write(112,*) "set xrange [1:201]"
  write(112,*)
  write(112,*) "set ylabel 'U'"
  write(112,*) "set yrange [-1.2:1.2]"
  write(112,*) "plot '10_results.txt' using 2 title 'U with N' with lines lt rgb 'blue', \"
  write(112,*) "     '10_results.txt' using 3 title 'U with N+0.5' with lines lt rgb 'red', \"
  write(112,*) "     '10_results.txt' using 4 title 'U with N+1' with lines lt rgb 'green', \"
  write(112,*) "     '10_results.txt' using 5 title 'U analytic' with lines lt rgb 'black'"
  write(112,*) 
  write(112,*) "set ylabel 'E'"
  write(112,*) "set yrange [-0.05:0.2]"
  write(112,*) "plot '10_results.txt' using 6 title 'E with N' with lines lt rgb 'blue', \"
  write(112,*) "     '10_results.txt' using 7 title 'E with N+0.5' with lines lt rgb 'red', \"
  write(112,*) "     '10_results.txt' using 8 title 'E with N+1' with lines lt rgb 'green'"
  write(112,*) 
  write(112,*) "pause -1"
  write(112,*) "unset multiplot"

  write(*,*) "See 10_results.txt"
  write(*,*)
  write(*,*) "Hit the Return (Enter) key to continue"

  close(112)
  close(100)

  call system('gnuplot 10_gnuxy')

  5000 format(1X,A,I3,A)
  6000 format(1X,I5,7F14.9)

end
