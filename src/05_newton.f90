program newton

  implicit none

! Use a Newton iteration to solve the equation --> x**3+x-10
!
!  x     -    current approximation to the solution
!  f(x)  -    polynomial function
!  df    -    derivative of f with respect to x
!  xo    -    previous guess for solution
!  eps   -    convergence criterion
!  dx    -    change in solution approximation
!  it    -    number of iterations
!  itmax -    maximum number of iterations

  integer, parameter :: itmax = 1000 
  real(8), parameter :: eps   = 1.e-6

  integer :: it
  real(8) :: x,f,df,xo

  it = 0
  f  = 1d0

  write(*,*) 'Try to solve "x**3+x-10=0"'
  write(*,*) 'What is your initial guess for the solution?'
  read(*,*) x

  do while( dabs(f)>eps .and. it<=itmax )
    it=it+1
    xo=x
    call derivate(xo,f,df)
    x = xo -f/df
    write(*,*) 'it = ',it,', x = ',x,', f(x) = ',f
  end do

end program

! ******************************************************************************************

subroutine derivate(x,f,d)

! Evaluate the function f(x)=x**3+x-10
! also return the derivative of the function
                                             
   implicit none

   real(8), intent(in)  :: x
   real(8), intent(out) :: f, d

   d = 3*x**2 + 1d0
   f = x**3 + x - 10d0

end subroutine
