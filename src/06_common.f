c * * * *
c syntaxe common
c common /nom de la zone commune/ liste des variables
c * * * * 

      PROGRAM test_arg

        implicit none
        integer a,b,c

        common /arg/ a,b,c

        a = 2
        c = 1

        print *, 'Before the call:'
        print *, 'a = ',a,', b = ',b,', c = ',c

        call sub

        print *, 'After the call:'
        print *, 'a = ',a,', b = ',b,', c = ',c

      END PROGRAM

      SUBROUTINE sub

        implicit none

        integer a,b,c
        common /arg/ a,b,c

        b = a + c
        c = c + 1
       
      END SUBROUTINE
