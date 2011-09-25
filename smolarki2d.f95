PROGRAM Smolarkiewicz2D
 IMPLICIT NONE
 INTEGER :: count, ierror = 0, MX, MY, N, i, j, k, iterations, initial_pos
 REAL :: dx, dy, dt, eps, uv, sc, u, v
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid, A, m_u, m_v, m_u_a, m_v_a
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
 initial = 0 
 initial(initial_pos) = 1 
 PRINT*, size(grid,2) 
 PRINT*, size(initial)
 grid(1,:) = initial

 OPEN(20, file = "wave.dat")

 m_u = u
 m_v = v
 PRINT*, iterations
 DO j = 1, N-1

  ! Build multiplication matrix

  A =  MATRIX(m_u, m_v, dx, dy, dt, MX, MY)
  psi_tem = grid(j,:)   
  !PRINT*, psi_tem 
  !Do k = 1, 12

  !PRINT*, A(k,1:12)
  !ENDDO
  CALL MVEC(A, psi_tem, psi_int)

  DO k = 1, iterations
	  PRINT*, "Iterating"
	  psi_tem = psi_int

	  CALL ANTIDIF(m_u, m_v, m_u_a, m_v_a, psi_tem, eps, dx, dy, dt, MX, MY)
	  PRINT*, m_u(1,1:10), " en ", m_u_a(1,1:10)
	  m_v_a = m_v_a * sc
	  m_u_a = m_u_a * sc
	  A =  MATRIX(m_u_a, m_v_a, dx, dy, dt, MX, MY) 
	  CALL MVEC(A, psi_tem, psi_int)
  ENDDO
  !PRINT*, psi_int
  grid(j+1,:) = psi_int
  WRITE(20,*), grid(j+1,:)
 ENDDO


 CLOSE(20)
 !PRINT*, "V is ", v
 !deallocation
 IF (ALLOCATED(initial)) DEALLOCATE(initial,STAT=ierror)
 IF (ierror /= 0) PRINT*, "initial : deallocation failed"
 IF (ALLOCATED(psi_int)) DEALLOCATE(psi_int,STAT=ierror)
 IF (ierror /= 0) PRINT*, "psi_int : deallocation failed"
 IF (ALLOCATED(m_u)) DEALLOCATE(m_u,STAT=ierror)
 IF (ierror /= 0) PRINT*, "m_u : deallocation failed"
 IF (ALLOCATED(m_v)) DEALLOCATE(m_v,STAT=ierror)
 IF (ierror /= 0) PRINT*, "m_v : deallocation failed"
 IF (ALLOCATED(m_u_a)) DEALLOCATE(m_u_a,STAT=ierror)
 IF (ierror /= 0) PRINT*, "m_u_a : deallocation failed"
 IF (ALLOCATED(m_v_a)) DEALLOCATE(m_v_a,STAT=ierror)
 IF (ierror /= 0) PRINT*, "m_v_a : deallocation failed"
 IF (ALLOCATED(A)) DEALLOCATE(A,STAT=ierror)
 IF (ierror /= 0) PRINT*, "A : deallocation failed"
 IF (ALLOCATED(grid)) DEALLOCATE(grid,STAT=ierror)
 IF (ierror /= 0) PRINT*, "grid : deallocation failed"

 CONTAINS 
  ! This function builds the multiplication matrix based on the velocity vector u
FUNCTION MATRIX(m_u, m_v, dx, dy, dt, MX, MY)
   REAL :: dx, dy, dt, div, alpha, beta
   REAL, DIMENSION(:,:) :: m_u, m_v
   REAL, DIMENSION(MX*MY,MX*MY) :: MATRIX
   INTEGER :: i, dimen, MX, MY, ii, jj

   dimen = MX*MY
   alpha = dt/(2*dx)
   beta = dt/(2*dy)

  ! PRINT*, "Build main diagonal D1"
  ! PRINT*, size(m_u), " en ", size(m_v), " en ", size(A)
   DO i=1, dimen, 1
     ii = i/MY + 1
     jj = mod(i,MX)
     if(jj.EQ.0) jj = MX
    ! PRINT*, ii, " en ", jj, " en ", m_u(ii,jj), " en ", m_v(ii,jj)
     ! Interpolate u and v
     MATRIX(i,i) = 1-alpha*(m_u(ii+1,jj)-abs(m_u(ii+1,jj)))+alpha*(m_u(ii,jj)+abs(m_u(ii,jj)))
     MATRIX(i,i) = MATRIX(ii,jj) -beta*(m_v(ii,jj+1)-abs(m_v(ii,jj+1)))+beta*(m_v(ii,jj)+abs(m_v(ii,jj)))
   ENDDO

   ! Build other for diagonals
! D2
   DO i=2, dimen, 1
     ii = i/MY + 1
     jj = mod(i,MX) 
     if(jj.EQ.0) jj = MX
     MATRIX(i,i-1) = alpha*(m_u(ii,jj)+abs(m_u(ii,jj)))
   ENDDO
   
! D3
   DO i=1, dimen-1, 1
     ii = i/MY + 1
     jj = mod(i,MX) 
     if(jj.EQ.0) jj = MX
     MATRIX(i,i+1) = -alpha*(m_u(ii+1,jj)-abs(m_u(ii+1,jj)))
   ENDDO

! D4
   DO i=MY+1,dimen, 1
     ii = i/MY + 1
     jj = mod(i,MX) 
     if(jj.EQ.0) jj = MX

     MATRIX(i,i-MY) = beta*(m_v(ii,jj)+abs(m_v(ii,jj)))

   ENDDO

! D5
   DO i=1, dimen-MY, 1
     ii = i/MY + 1
     jj = mod(i,MX) 
     if(jj.EQ.0) jj = MX

     MATRIX(i,i+MY) = -beta*(m_v(ii,jj+1)-abs(m_v(ii,jj+1)))

   ENDDO
  END FUNCTION 


  ! Fixed for 2D
  SUBROUTINE ANTIDIF(u, v, u_a, v_a, psi, eps, dx, dy, dt, MX, MY)
   REAL :: eps, dx, dy, dt
   REAL, DIMENSION(:) :: psi
   REAL, DIMENSION(MX+1, MY) :: u, u_a
   REAL, DIMENSION(MX, MY+1) :: v, v_a 
   REAL, DIMENSION(MX, MY) :: psi_mat
   INTEGER :: i, j, MX, MY
   psi_mat = RESHAPE( psi, (/ MX, MY /))
   ! Build antidiffusion vector u_a
   DO j=1, MY
	   DO i=2, MX+1
	    u_a(i,j) = ((abs(u(i,j))*dx-dt*u(i,j)*u(i,j))*(psi_mat(i,j)-psi_mat(i-1,j)))/((psi_mat(i-1,j)+psi_mat(i,j)+eps)*dx)
if(psi_mat(i,j) > 0) then
PRINT*, u(i,j),  " naar ", u_a(i,j), " delen door ", ((psi_mat(i-1,j)+psi_mat(i,j)+eps)*dx), " met ", psi_mat(i,j)
	  endif
 ENDDO
   ENDDO
   ! Build antidiffusion vector v_a
   DO i=1, MX
	   DO j=2, MY+1
	    v_a(i,j) = ((abs(v(i,j))*dy-dt*v(i,j)*v(i,j))*(psi_mat(i,j)-psi_mat(i,j-1)))/((psi_mat(i,j-1)+psi_mat(i,j)+eps)*dy)
	   ENDDO
   ENDDO
  END SUBROUTINE
 
  ! This subroutine multiplies a matrix A with a vector x and returns their product y
  ! Same for 2D as for 1D
  SUBROUTINE MVEC(A,x,y)
   REAL, DIMENSION(:) :: x, y
   REAL, DIMENSION(:,:) :: A
   DO i=1, size(x)
    y(i) = sum(A(i,:) * x)
   ENDDO  
  END SUBROUTINE


  !! Read in files
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
      ALLOCATE(m_u((MX+1),MY), STAT=ierror); IF (ierror /= 0) PRINT*, "u : Allocation failed"
      ALLOCATE(m_v(MX,(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v : Allocation failed"
      ALLOCATE(m_u_a((MX+1),MY), STAT=ierror); IF (ierror /= 0) PRINT*, "u_a : Allocation failed"
      ALLOCATE(m_v_a(MX,(MY+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v_a : Allocation failed"
      ALLOCATE(A((MX*MY),(MX*MY)), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
     case ('N')
      read(buffer,*, iostat=ios) N
      ALLOCATE(grid(N,(MX*MY)), STAT=ierror); IF (ierror /= 0) PRINT*, "grid : Allocation failed"
     case ('dx')
      read(buffer,*, iostat=ios) dx
     case ('dy')
      read(buffer,*, iostat=ios) dy     
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
     case ('initial_pos')
      read(buffer,*, iostat=ios) initial_pos
     case ('m_u')
      read(buffer,*, iostat=ios) m_u
     case ('m_v')
      read(buffer,*, iostat=ios) m_v
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
