! Implementation of Smolarkiewicz Positive Definite Advection Scheme with Small Implicit Diffusion (1983)
! G.A. Burger & J.M. Wolterink, October 2011
! 1 dimensional simulation
PROGRAM Smolarkiewicz
 IMPLICIT NONE

! Variables 
 INTEGER :: nargs, ierror, M, N, i, j, k, iterations
 REAL :: dx, dt, eps, uv, sc
 REAL, DIMENSION(:,:), ALLOCATABLE :: A
 REAL, DIMENSION(:), ALLOCATABLE :: initial_x, psi_int, psi_tem, velocity_u, velocity_antidif
 character(len=100) :: input_file

! Check whether there is an input file
 nargs = command_argument_count()
 if(nargs < 1) then
  print*, "Please provide the name of the input file"
  stop 1
 else
  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) then
   print*, "Please provide the name of the input file"
   stop 1
  end if
 end if

! Read input file
 call read_input_file(input_file)
 eps = 0.00000000000000000001
! Open output file
 OPEN(20, file = "wave.dat")

! Initialize and write to file
 write(20,*), initial_x

 psi_int = initial_x
! Initialize velocity vector
 velocity_u = uv

! Repeat for N time steps
 DO j = 1, N-1

! Build tridiagonal sparse multiplication matrix
  A =  MATRIX(velocity_u, dx, dt)

! Make copy of current situation
  psi_tem = psi_int  

! Multiply vector to get new situation
  CALL MVEC(A, psi_tem, psi_int)

! Apply antidiffusion step for some iterations (optional, iterations could be 0)
  DO k = 1, iterations

! Make local copy of current situation
	  psi_tem = psi_int

! Build antidiffusion vector u~
	  velocity_antidif = ANTIDIF(velocity_u, psi_tem, eps, dx, dt)
! Apply factor Sc (optional, sc could be 1)
	  velocity_antidif = velocity_antidif * sc
! Build new multiplication matrix with antidiffusion vector
	  A =  MATRIX(velocity_antidif, dx, dt) 
! Multiply
	  CALL MVEC(A, psi_tem, psi_int)
  ENDDO

! Write new situation to file
  WRITE(20,*), psi_int
 ENDDO
 CLOSE(20)

! Deallocation
 IF (ALLOCATED(initial_x)) DEALLOCATE(initial_x,STAT=ierror)
 IF (ierror /= 0) PRINT*, "initial_x : deallocation failed"
 IF (ALLOCATED(psi_int)) DEALLOCATE(psi_int,STAT=ierror)
 IF (ierror /= 0) PRINT*, "psi_int : deallocation failed"
 IF (ALLOCATED(velocity_u)) DEALLOCATE(velocity_u,STAT=ierror)
 IF (ierror /= 0) PRINT*, "velocity_u : deallocation failed"
 IF (ALLOCATED(velocity_antidif)) DEALLOCATE(velocity_antidif,STAT=ierror)
 IF (ierror /= 0) PRINT*, "velocity_antidif : deallocation failed"
 IF (ALLOCATED(A)) DEALLOCATE(A,STAT=ierror)
 IF (ierror /= 0) PRINT*, "A : deallocation failed"

 print*, "Program completed succesfully"

CONTAINS 
 
! This function builds the multiplication matrix based on the velocity vector u
! or the antidiffusion velocity vector u~. It returns a tridiagonal sparse matrix.
  FUNCTION MATRIX(u, dx, dt)
   REAL :: dx, dt, div
   REAL, DIMENSION(:) :: u
   REAL, DIMENSION(size(u)-1,3) :: MATRIX
   INTEGER :: i, dimen
   MATRIX = 0 	
   dimen = size(u) - 1
   div = dt/(2*dx)
   MATRIX(1,2) = 1 - div*(u(2)+abs(u(2))-u(1)+abs(u(1)))
   MATRIX(1,3) = -div*(u(2)-abs(u(2)))
   DO i=2, dimen-1, 1
    MATRIX(i,1) = div*(u(i)+abs(u(i))) 
    MATRIX(i,2) = 1 - div*(u(i+1)+abs(u(i+1))-u(i)+abs(u(i)))
    MATRIX(i,3) = -div*(u(i+1)-abs(u(i+1)))
   ENDDO
   MATRIX(dimen,2) = 1 - div*(u(dimen+1)+abs(u(dimen+1))-u(dimen)+abs(u(dimen))) 
   MATRIX(dimen,1) = div*(u(dimen)+abs(u(dimen)))
  END FUNCTION 

! This function returns the antidiffusion velocity vector based on the velocity vector and the 
! state psi of the system
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
 
! This subroutine multiplies a tridiagonal sparse matrix A with a vector x and returns their product vector y
  SUBROUTINE MVEC(A,x,y)
   REAL, DIMENSION(:) :: x, y
   REAL, DIMENSION(:,:) :: A
   !PRINT*, x
   M = size(x) 
   y(1) = A(1,2)*x(1) + A(1,3)*x(2) 
   DO i=2, M-1
    y(i) = sum(A(i,:) * x(i-1:i+1))
   ENDDO  
   y(M) = A(M,1)*x(m-1)+A(M,2)*x(2)
  END SUBROUTINE

  subroutine read_input_file(filename)
! Inputs all relevant parameters from the inputfile specified on the commandline
   character(len=*), intent(in) :: filename
   character(len=30) :: label
   character(len=20000) :: buffer
   integer, parameter :: i1fh = 11
   integer :: pos1, ios = 0, line = 0
   logical :: file_exists

   inquire(file=filename,exist=file_exists)
   if(.not. file_exists) then
    print*, "ERROR: File ", trim(filename), " does not exist"
    stop 1
   end if
   open(i1fh, file=filename)
   print*, "Reading parameters from file ", filename
   do while (ios == 0)
    read(i1fh,'(A)', iostat=ios) buffer
    if(ios == 0) then
     line = line + 1
     pos1 = scan(buffer,' ')
     label = buffer(1:pos1)
     buffer = buffer(pos1+1:)
     select case (label)
     case ('M')
      read(buffer,*, iostat=ios) M
      ALLOCATE(initial_x(M), STAT=ierror); IF (ierror /= 0) PRINT*, "initial_x : Allocation failed"
      ALLOCATE(psi_int(M), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_int : Allocation failed"
      ALLOCATE(psi_tem(M), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_tem : Allocation failed"
      ALLOCATE(velocity_u(M+1), STAT=ierror); IF (ierror /= 0) PRINT*, "velocity_u : Allocation failed"
      ALLOCATE(velocity_antidif(M+1), STAT=ierror); IF (ierror /= 0) PRINT*, "velocity_antidif : Allocation failed"
      ALLOCATE(A(M,3), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
     case ('dx')
      read(buffer,*, iostat=ios) dx
     case ('N')
      read(buffer,*, iostat=ios) N
     case ('dt')
      read(buffer,*, iostat=ios) dt
     case ('eps')
      read(buffer,*, iostat=ios) eps
     case ('iterations')
      read(buffer,*, iostat=ios) iterations
     case ('sc')
      read(buffer,*, iostat=ios) sc
     case ('uv')
      read(buffer,*, iostat=ios) uv
     case ('initial_x')
      read(buffer,*, iostat=ios) initial_x
     case ('velocity_u')
      read(buffer,*, iostat=ios) velocity_u
     case default
      if (.not.(label(1:1) == ' ' .or. label(1:1) == '!')) then
       print*, 'Skipping invalid label at line', line
      end if 
     end select
    end if
   end do
  end subroutine read_input_file
END PROGRAM 
