PROGRAM Smolarkiewicz
 IMPLICIT NONE
 INTEGER :: M, N
 REAL :: dx, dy
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid 
 REAL, DIMENSION(:), ALLOCATABLE :: initial_x
 REAL, DIMENSION(:), ALLOCATABLE :: velocity_u
 

 PRINT*, "Reading from file"
 OPEN(10, file = "advec.dat")
 READ(10,*) M
 PRINT*, "Number of elements in space ", M
 ALLOCATE(initial_x(M))
 ALLOCATE(velocity_u(M+1))
 READ(10,*) initial_x
 READ(10,*) velocity_u
 
 PRINT*, "Give number of timesteps"
 READ*, N
 ALLOCATE(grid(N,M))

 PRINT*, "Give dx"
 READ*, dx
 PRINT*, "Give dy"
 READ*, dy
 
 grid(1,:) = initial_x

 PRINT*, grid(1,:)
 PRINT*, velocity_u
 CLOSE(10)
 
 DEALLOCATE(initial_x)
 DEALLOCATE(velocity_u)
 DEALLOCATE(grid)
END PROGRAM 
