c
c   Program to demonstrate Arithmetic Assignments
c
      program arith
      implicit none
c
c     declare the data types for all Fortran variables
c
      real r2,r3,r4,r5,r6,ans1,ans2,ans3
      integer i2,i3,i4,i5,i6,ians1,ians2,ians3,ians4
c
c     r2 thru r6 take on the real values 2.0 thru 6.0
c
c     i2 thru i6 take on the integer values 2 thru 6
c
c     ans1, ans2, and ans3 will contain the answers from
c     real arithmetic
c
c     ians1 thru ians4 will contain the answers from
c     integer arithmetic
c
c
c     Set initial values of the variables with 2 valid forms
c     of data statements
      data r2/2./,r3/3./,r4/4.0/,r5/5.0/
      data i2,i3,i4,i5/2,3,4,5/
c
c     This ends the non-executable statements, nothing above
c     this point results in a machine instruction to perform
c     some operation.
c     Executable statements follow.
c
c     The result of any integer divide is truncated to the integer
c     value less than the correct decimal answer for the division
c     The result of this is that changing the order of operations
c     can make a big difference in the answers.  Notice how parentheses
c     force more expected results
c
      ians1=i2*i3/i5
      ians2=i3/i5*i2
      ians3=i2*(i3/i5)
      ians4=(i3/i5)*i2
      print *, '2*3/5 =', ians1, ', 3/5*2 =',ians2,
     &  ', 2*(3/5) =',ians3 ,', (3/5)*2 =',ians4
c
c     Real arithmetic behaves more uniformly
c
      ans1=r2*r3/r5
      ans2=r3/r5*r2
      ans3=(r3/r5)*r2
      print *, '2.0*3.0/5.0 =', ans1, ', 3.0/5.0*2.0 =',ans2,
     &  ', (3.0/5.0)*2.0 =',ans3
c
c     Watch how precedence of operations effects the following:
c
      ians1=i2+i5*i3**i2
      ians2=i5*i3**i2+i2
      ians3=i3**i2*i5+i2
      print *, '2+5*3**2 =',ians1,', 5*3**2+2 =',ians2,
     & ', 3**2*5+2 =',ians3
c
c     You can mix real and integers, but watch what happens
c
      ans1=r5+i3/i2
      print *, '5.0+3/2 =',ans1

c
c     You can do the same thing with constants in the expression
c
      ans2=5.0+3/2
      print *, '5.0+3/2 =',ans2
c
c     Look at what happens when I put a real in either the numerator
c     or denominator of the division term
      ans1=r5+i3/r2
      ans2=r5+r3/i2
      print *, '5.0+3/2.0 =',ans1, ', 5.0+3.0/2 =', ans2
c

c     Although Fortran normally works from left to right at a given
c     level of precedence (does all multiply and divide from left to
c     right before moving on to adds and subtracts).  It works
c     exponentiation from right to left when it hits 2 or more
c     sequential exponentiation operations
c
      ians1= i5**i3**i2
      ians2= (i5**i3)**i2
      ians3= i5**(i3**i2)
      print *, '5**3**2 =',ians1, ', (5**3)**2 =',ians2,
     &  ', 5**(3**2) =',ians3
c
c    When in doubt use parentheses to get the answer that you
c    really want.
c
      stop
      end
