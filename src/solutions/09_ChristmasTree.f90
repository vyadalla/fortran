program ChristmasTree

  implicit none

  integer i, h, hmax, line, ball
  character(8) hmax_string
  character(1) sball

  ! get the command line argument
  call getarg(1,hmax_string)

  ! cast string to integer
  read(hmax_string,*) hmax

  ball=1 
  do h=1,hmax
    line=1
    ! write spaces to align the head of the tree
    write(*,'(a)',advance='no') repeat(' ',hmax-h)
    ! loop to decide when we have to create a new line
    do while(line.le.(2*h-1))
      ! modulo to decide when we have to put a ball or not
      sball='#'
      if(mod(ball,6).eq.0) sball='o'
      write(*,'(A)',advance='no') sball
      ! increment ball/line decision variables
      ball=ball+1
      line=line+1
    enddo
    ! create a new line
    write(*,*)
  enddo

end program
