PROGRAM Smolarkiewicz2D
 IMPLICIT NONE
 INTEGER :: count, ierror = 0, MX, MY, N, i, j, k, iterations
 REAL :: dx, dy, dt, eps, uv, sc
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid, A, u, v, u_a, v_a
 REAL, DIMENSION(:), ALLOCATABLE :: initial, psi_int, psi_tem
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
 
 grid(1,:) = initial
 u = uv
 v = uv
 OPEN(20, file = "wave.dat")

 PRINT*, iterations
 DO j = 1, N-1
  ! Build multiplication matrix
  A =  MATRIX(u, v, dx, dy, dt, MX, MY)
  psi_tem = grid(j,:)  
  CALL MVEC(A, psi_tem, psi_int)
  DO k = 1, iterations
	  PRINT*, "Iterating"
	  psi_tem = psi_int
	  velocity_antidif = ANTIDIF(velocity_u, psi_tem, eps, dx, dy, dt)
	  velocity_antidif = velocity_antidif * sc
	  A =  MATRIX(velocity_antidif, dx, dy, dt, MX, MY) 
	  CALL MVEC(A, psi_tem, psi_int)
  ENDDO
  grid(j+1,:) = psi_int
  WRITE(20,*), grid(j+1,:)
 ENDDO
 CLOSE(20)

 !deallocation
 IF (ALLOCATED(initial)) DEALLOCATE(initial_x,STAT=ierror)
 IF (ierror /= 0) PRINT*, "initial : deallocation failed"
 IF (ALLOCATED(psi_int)) DEALLOCATE(psi_int,STAT=ierror)
 IF (ierror /= 0) PRINT*, "psi_int : deallocation failed"
 IF (ALLOCATED(u)) DEALLOCATE(u,STAT=ierror)
 IF (ierror /= 0) PRINT*, "u : deallocation failed"
 IF (ALLOCATED(v)) DEALLOCATE(v,STAT=ierror)
 IF (ierror /= 0) PRINT*, "v : deallocation failed"
 IF (ALLOCATED(u_a)) DEALLOCATE(u_a,STAT=ierror)
 IF (ierror /= 0) PRINT*, "u_a : deallocation failed"
 IF (ALLOCATED(v_a)) DEALLOCATE(v_a,STAT=ierror)
 IF (ierror /= 0) PRINT*, "v_a : deallocation failed"
 IF (ALLOCATED(A)) DEALLOCATE(A,STAT=ierror)
 IF (ierror /= 0) PRINT*, "A : deallocation failed"
 IF (ALLOCATED(grid)) DEALLOCATE(grid,STAT=ierror)
 IF (ierror /= 0) PRINT*, "grid : deallocation failed"

 CONTAINS 
  ! This function builds the multiplication matrix based on the velocity vector u
  FUNCTION MATRIX(u, v, dx, dy, dt, MX, MY)
   REAL :: dx, dt, div
   REAL, DIMENSION(:) :: u, v, uu, vv
   REAL, DIMENSION(MX*MY) :: MATRIX
   INTEGER :: i, dimen

   dimen = MX*MY
   alpha = dt/(2*dx)
   beta = dt/(2*dy)

   ! Build main diagonal
   DO i=1, dimen, 1
     uu = (u(:,i)+u(:,i+1))/2
     vv = (v(i,:)+v(i+1,:))/2
     MATRIX(i,i) = 1-alpha*(uu(i+1)+abs(uu(i+1)))+alpha*(uu(i)-abs(uu(i)))-beta*(vv(i+1)+abs(vv(i+1)))+beta*(vv(i)-abs(vv(i)))
   ENDDO

   ! Build other for diagonals
   DO i=1, dimen-1, 1
     vv = (v(i,:)+v(i+1,:))/2
     MATRIX(i,i-1) = beta*(vv(i)+abs(vv(i)))
   ENDDO
   
   DO i=1, dimen-1, 1
     vv = (v(i,:)+v(i+1,:))/2
     MATRIX(i,i+1) = -beta*(vv(i+1)+abs(vv(i+1)))
   ENDDO

   DO i=1, dimen-MY, 1
     uu = (u(:,i)+u(:,i+1))/2
     MATRIX(i,i+MY) = alpha*(uu(i)+abs(uu(i)))
   ENDDO

   DO i=MY,dimen, 1
     uu = (u(:,i)+u(:,i+1))/2
     MATRIX(i,i-MY) = alpha*(uu(i+1)+abs(uu(i+1)))
   ENDDO
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
     case ('MX')
      read(buffer,*, iostat=ios) MX
     case ('MY')
      read(buffer,*, iostat=ios) MY
      ALLOCATE(initial(MX*MY), STAT=ierror); IF (ierror /= 0) PRINT*, "initial : Allocation failed"
      ALLOCATE(psi_int(MX*MY), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_int : Allocation failed"
      ALLOCATE(psi_tem(MX*MY), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_tem : Allocation failed"
      ALLOCATE(u((MX+1),(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "u : Allocation failed"
      ALLOCATE(v((MX+1),(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v : Allocation failed"
      ALLOCATE(u_a((MX+1),(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "u_a : Allocation failed"
      ALLOCATE(v_a((MX+1),(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v_a : Allocation failed"
      ALLOCATE(A((MX*MY),(MX*MY)), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
     case ('N')
      read(buffer,*, iostat=ios) N
      ALLOCATE(grid(N,M), STAT=ierror); IF (ierror /= 0) PRINT*, "grid : Allocation failed"
     case ('dx')
      read(buffer,*, iostat=ios) dx
     case ('dx')
      read(buffer,*, iostat=ios) dx     case ('dt')
      read(buffer,*, iostat=ios) dt
     case ('eps')
      read(buffer,*, iostat=ios) eps
     case ('iterations')
      read(buffer,*, iostat=ios) iterations
     case ('sc')
      read(buffer,*, iostat=ios) sc
     case ('uv')
      read(buffer,*, iostat=ios) uv
     case ('initial')
      read(buffer,*, iostat=ios) initial_x
     case ('u')
      read(buffer,*, iostat=ios) u
     case ('v')
      read(buffer,*, iostat=ios) v
     case default
      if (.not.(label(1:1) == ' ' .or. label(1:1) == '!')) then
       print*, 'Skipping invalid label at line ', line
      end if 
     end select
    end if
   end do
  end subroutine read_input_file
END PROGRAM 
