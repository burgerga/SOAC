PROGRAM Smolarkiewicz
 IMPLICIT NONE
<<<<<<< HEAD
 INTEGER :: M, N, i, j, antidiffusion = 0
 REAL :: dx, dt, eps, uv, sc
=======
 INTEGER :: count, ierror = 0, M, N, i, j, antidiffusion = 0
 REAL :: dx, dt, eps, uv
>>>>>>> 511e38922d89c9227b2792cfecafbcc450f56121
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid, A
 REAL, DIMENSION(:), ALLOCATABLE :: initial_x, psi_int, velocity_u, velocity_antidif
 character(len=100) :: input_file

 count = command_argument_count()
 if(count < 1) then
  print*, "Please provide the name of the input file"
  stop 1
 else
  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) then
   print*, "Please provide the name of the input file"
   stop 1
  end if
 end if

 call read_input_file(input_file)
 
<<<<<<< HEAD
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
=======
>>>>>>> 511e38922d89c9227b2792cfecafbcc450f56121
 grid(1,:) = initial_x
 velocity_u = uv
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
  PRINT*, "Antidiffusion"
  ! Multiply for first step
  CALL MVEC(A, grid(j,:), psi_int)
  ! Get antidiffusion velocity vector
  velocity_antidif = ANTIDIF(velocity_u, psi_int, eps, dx, dt)
  velocity_antidif = velocity_antidif * sc
  ! Get matrix for second step
  A =  MATRIX(velocity_antidif, dx, dt) 
  ! Multiply again
  CALL MVEC(A, psi_int, grid(j+1,:))
  ! If no antidiffusion step
  ELSE IF (antidiffusion == 0) THEN 
  ! Find new vector directly with matrix multiplication
  CALL MVEC(A, grid(j,:), grid(j+1,:))
  END IF
  WRITE(20,*), grid(j+1,:)
 ENDDO
 CLOSE(20)

 !deallocation
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
 IF (ALLOCATED(grid)) DEALLOCATE(grid,STAT=ierror)
 IF (ierror /= 0) PRINT*, "grid : deallocation failed"

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

  subroutine read_input_file(filename)
   character(len=*), intent(in) :: filename
   character(len=30) :: label
   character(len=200) :: buffer
   integer, parameter :: i1fh = 11
   integer :: pos, ios = 0, line = 0

   open(i1fh, file=filename)
   print*, "Reading parameters from file", filename
   do while (ios == 0)
    read(i1fh,'(A)', iostat=ios) buffer
    if(ios == 0) then
     line = line + 1
     pos = scan(buffer,' ')
     label = buffer(1:pos)
     buffer = buffer(pos+1:)
     select case (label)
     case ('M')
      read(buffer,*, iostat=ios) M
      ALLOCATE(initial_x(M), STAT=ierror); IF (ierror /= 0) PRINT*, "initial_x : Allocation failed"
      ALLOCATE(psi_int(M), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_int : Allocation failed"
      ALLOCATE(velocity_u(M+1), STAT=ierror); IF (ierror /= 0) PRINT*, "velocity_u : Allocation failed"
      ALLOCATE(velocity_antidif(M+1), STAT=ierror); IF (ierror /= 0) PRINT*, "velocity_antidif : Allocation failed"
      ALLOCATE(A(M,M), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
     case ('N')
      read(buffer,*, iostat=ios) N
      ALLOCATE(grid(N,M), STAT=ierror); IF (ierror /= 0) PRINT*, "grid : Allocation failed"
     case ('dx')
      read(buffer,*, iostat=ios) dx
     case ('dt')
      read(buffer,*, iostat=ios) dt
     case ('eps')
      read(buffer,*, iostat=ios) eps
     case ('uv')
      read(buffer,*, iostat=ios) uv
     case ('antidiffusion')
      read(buffer,*, iostat=ios) antidiffusion
     case ('initial_x')
      read(buffer,*, iostat=ios) initial_x
     case ('velocity_u')
      read(buffer,*, iostat=ios) velocity_u
     case default
      if (.not.(label(1:1) == ' ' .or. label(1:1) == '!')) then
       print*, 'Skipping invalid label at line ', line
      end if 
     end select
    end if
   end do
  end subroutine read_input_file
END PROGRAM 
