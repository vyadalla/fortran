MODULE arg
 implicit none
 integer :: a,b,c
 real(8) :: x
END MODULE arg

! * * * * * * * 

PROGRAM test_arg
 USE arg
 implicit none

 a = 2
 c = 1

 write(*,*) 'Before the call:'
 write(*,'(3(A5,I3))') ' a = ',a,', b = ',b,', c = ',c

 call sub

 write(*,*) 'After the call:'
 write(*,'(3(A5,I3))') 'a = ',a,', b = ',b,', c = ',c

END PROGRAM test_arg

! * * * * * * * 

SUBROUTINE sub
 USE arg, only : a,b,c    ! seuls a b et c sont utiles
 implicit none

 b = a + c
 c = c + 1

END SUBROUTINE sub
