PROGRAM Smolarkiewicz
 IMPLICIT NONE
 INTEGER :: M, N, i, j, k, antidiffusion = 0, iterations 
 REAL :: dx, dt, eps, uv, sc
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid, A
 REAL, DIMENSION(:), ALLOCATABLE :: initial_x, psi_int, psi_tem, velocity_u, velocity_antidif

 PRINT*, "Reading dimensions from file"
 OPEN(10, file = "advec.dat")
 READ(10,*) M
 PRINT*, "Number of elements in space ", M
 ALLOCATE(initial_x(M))
 ALLOCATE(psi_int(M))
 ALLOCATE(psi_tem(M))
 ALLOCATE(velocity_u(M+1))
 ALLOCATE(velocity_antidif(M+1))
 ALLOCATE(A(M,M))
 READ(10,*) initial_x
 CLOSE(10)
 
 PRINT*, "How many timesteps?"
 READ*, N
 ALLOCATE(grid(N,M))
 PRINT*, "Give dx"
 READ*, dx
 PRINT*, "Give dt"
 READ*, dt
 PRINT*, "Give epsilon"
 READ*, eps
 PRINT*, "Give wind speed"
 READ*, uv
 PRINT*, "Antidiffusion step?"
 READ*, antidiffusion
 PRINT*, "Give Sc factor"
 READ*, sc
 PRINT*, "How many iterations?"
 READ*, iterations
 grid(1,:) = initial_x
 ! PRINT*, "Initial value"
 ! PRINT*, initial_x
 ! PRINT*, grid(1,:)
 velocity_u = uv
 ! PRINT*, velocity_u
 PRINT*, antidiffusion
 IF (antidiffusion == 1) THEN 
 OPEN(20, file = "wavea.dat")
 ELSE IF (antidiffusion == 0) THEN
 OPEN(20, file = "wave.dat")
 END IF
 DO j = 1, N-1
  ! Build multiplication matrix
  A =  MATRIX(velocity_u, dx, dt)
  
  ! If antidiffusion step
  IF(antidiffusion == 1) THEN
  psi_int = grid(j,:) 
  DO k = 1, iterations
  ! Multiply for first step
  CALL MVEC(A, psi_int, psi_tem)
  psi_int = psi_tem
  ! Get antidiffusion velocity vector
  velocity_antidif = ANTIDIF(velocity_u, psi_int, eps, dx, dt)
  velocity_antidif = velocity_antidif * sc
  ! Get matrix for second step
  A =  MATRIX(velocity_antidif, dx, dt) 
  ! Multiply again
  CALL MVEC(A, psi_int, psi_tem)
  psi_int = psi_tem
  ! If no antidiffusion step
  ENDDO
  grid(j+1,:) = psi_int
  ELSE IF (antidiffusion == 0) THEN 
  ! Find new vector directly with matrix multiplication
  CALL MVEC(A, grid(j,:), grid(j+1,:))
  END IF
  WRITE(20,*), grid(j+1,:)
 ENDDO
 CLOSE(20)
 DEALLOCATE(initial_x)
 DEALLOCATE(velocity_u)
 DEALLOCATE(grid)
 CONTAINS 
  ! This function builds the multiplication matrix based on the velocity vector u
  FUNCTION MATRIX(u, dx, dt)
   REAL :: dx, dt, div
   REAL, DIMENSION(:) :: u
   REAL, DIMENSION(size(u)-1,size(u)-1) :: MATRIX
   INTEGER :: i, dimen

   dimen = size(u) - 1
   div = dt/(2*dx)
   DO i=2, dimen-1, 1
    MATRIX(i,i-1) = div*(u(i)+abs(u(i))) 
    MATRIX(i,i) = 1 - div*(u(i+1)+abs(u(i+1))-u(i)+abs(u(i)))
    MATRIX(i,i+1) = -div*(u(i+1)-abs(u(i+1)))
   ENDDO
   MATRIX(1,1) = 1 - div*(u(2)+abs(u(2))-u(1)+abs(u(1)))
   MATRIX(dimen,dimen) = 1 - div*(u(dimen+1)+abs(u(dimen+1))-u(dimen)+abs(u(dimen))) 
   MATRIX(1,2) = -div*(u(2)-abs(u(2)))
   MATRIX(dimen,dimen-1) = div*(u(dimen)+abs(u(dimen)))
  END FUNCTION 

  FUNCTION ANTIDIF(u, psi, eps, dx, dt)
   REAL :: eps, dx, dt
   REAL, DIMENSION(:) :: u, psi
   REAL, DIMENSION(size(u)) :: ANTIDIF
   INTEGER :: i, dimen
   dimen = size(u) 
   DO i=2, dimen
    ANTIDIF(i) = ((abs(u(i))*dx-dt*u(i)*u(i))*(psi(i)-psi(i-1)))/((psi(i-1)+psi(i)+eps)*dx)
   ENDDO
  END FUNCTION 
 
  ! This subroutine multiplies a matrix A with a vector x and returns their product y
  SUBROUTINE MVEC(A,x,y)
   REAL, DIMENSION(:) :: x, y
   REAL, DIMENSION(:,:) :: A
   DO i=1, size(x)
    y(i) = sum(A(i,:) * x)
   ENDDO  
  END SUBROUTINE
END PROGRAM 
