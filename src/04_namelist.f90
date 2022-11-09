PROGRAM test_namelist

  implicit none

  real*8 lon_min, lon_max, lat_min, lat_max

  NAMELIST/namlon/ lon_min, lon_max
  NAMELIST/namlat/ lat_min, lat_max

  write(*,*) 'Before:'
  call print_res(lon_min, lon_max, lat_min, lat_max)

  open(161,file='04_namelist.def',status='old',form='formatted')
  read(161,NML=namlon)

  write(*,*) 'Between:'
  call print_res(lon_min, lon_max, lat_min, lat_max)

  read(161,NML=namlat)
  close (161)

  write(*,*) 'After:'
  call print_res(lon_min, lon_max, lat_min, lat_max)

END

SUBROUTINE print_res(a,b,c,d)
  implicit none
  real*8, intent(in) :: a,b,c,d
  write(*,'(4(A12,F6.2))') '  lon_min = ',a,', lon_max = ',b, &
                          ', lat_min = ',c,', lat_max = ',d
  RETURN
END
