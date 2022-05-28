      PROGRAM simpson_2d

!     Kuntal Ghosh
!     May 2022

!     Far superior to 2D trapezoidal rule!

      IMPLICIT NONE

      REAL*8 :: x1_start,x1_end,x2_start,x2_end,dx1,dx2,&
                integ,x1,x2
      REAL*8, EXTERNAL :: expression
      INTEGER*8 :: i,j,n,n1,n2
      INTEGER*8, ALLOCATABLE :: ss(:,:),sa(:),sb(:)
      REAL*8, ALLOCATABLE :: ff(:,:)

!     Enter the limits here      
      x1_start = 0.0d0
      x1_end = 1.0d0
      x2_start = 0.5d0
      x2_end = 1.0d0
!      n = 20000
      n1 = 2000
      n2 = 3000

      dx1 = (x1_end-x1_start)/DFLOAT(n1-1)
      dx2 = (x2_end-x2_start)/DFLOAT(n2-1)

!     Simpson's matrix is ss, function matrix is ff

      ALLOCATE (ss(n1,n2),ff(n1,n2),sa(n1),sb(n2))

      DO i = 1,n1
         IF (MOD(i,2)==0) THEN
            sa(i) = 4
         ELSE
            sa(i) = 2
         END IF
      END DO

      sa(1) = 1
      sa(n1) = 1

      DO i = 1,n2
         IF (MOD(i,2)==0) THEN
            sb(i) = 4
         ELSE
            sb(i) = 2
         END IF
      END DO

      sb(1) = 1
      sb(n2) = 1

      DO i = 1,n1
         DO j = 1,n2
            ss(i,j) = sa(i)*sb(j)
         END DO
      END DO

      DO i = 1,n1
         x1 = x1_start + (i-1)*dx1
         DO j = 1,n2
            x2 = x2_start + (j-1)*dx2
            ff(i,j) = expression(x1,x2)
         END DO
      END DO

      integ = 0.0d0
      DO i = 1,n1
         DO j = 1,n2
            integ = integ + ss(i,j)*ff(i,j)
         END DO
      END DO
                        
      integ = dx1*dx2*integ/9.0d0

      WRITE (*,*) "The integral is =", integ 
      WRITE (*,*) "Exact integral = ", (DATAN(1.0d0) - DATAN(0.5d0))

      DEALLOCATE (ss,ff,sa,sb)

      END PROGRAM simpson_2d


      REAL(KIND=8) FUNCTION expression(a,b) RESULT(c)

      IMPLICIT NONE

      REAL*8, INTENT(IN) :: a,b

!     Enter function here      
      c = (b**2-a**2)/(a**2 + b**2)**2

      END FUNCTION expression
