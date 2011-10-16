! Same as Smolarkiewicz2D except with sparse matrices
PROGRAM Smolarkiewicz2DS
 IMPLICIT NONE
 INTEGER :: count, ierror = 0, MX, MY, N, i, j, k, q, iterations, initial_pos
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
  PRINT*, "Got matrix"
  psi_tem = grid(j,:)   
  CALL MVEC(A, psi_tem, psi_int,MX,MY)
  !if(j.EQ.2) THEN
!	OPEN(30, file = "matrix.dat")!
!	DO q=1, MX*MY
!	 WRITE(30,*), A(q,:)
 ! 	ENDDO
!	CLOSE(30)
 ! ENDIF


  DO q = 1, iterations
	  PRINT*, "Iterating"
	  psi_tem = psi_int
	  CALL ANTIDIF(m_u, m_v, m_u_a, m_v_a, psi_tem, eps, dx, dy, dt, MX, MY)
	  A =  MATRIX(m_u_a, m_v_a, dx, dy, dt, MX, MY) 
	  PRINT*, A(1,:)
	  PRINT*, A(2,:)
	  PRINT*, A(3,:)
	  PRINT*, A(4,:)
	  PRINT*, A(5,:)
	  CALL MVEC(A, psi_tem, psi_int,MX,MY)
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
   REAL, DIMENSION(MX*MY,5) :: MATRIX
   INTEGER :: dimen, MX, MY, ii, jj, di

   MATRIX = 0 
   dimen = MX*MY
   alpha = dt/(2*dx)
   beta = dt/(2*dy)
  ! ii and jj is location in grid, rows are ii columns are jj. eg. (4,3) is fourth row, third column
  ! PRINT*, "Build main diagonal D1"

   DO di=1, dimen, 1
     ! Get location of psi element on normal grid  
    CALL GETGRIDLOCATION(di, MX, MY, ii, jj)
    ! PRINT*, ii, " en ", jj, " en ", m_u(ii,jj), " en ", m_v(ii,jj)
     ! Interpolate u and v
     MATRIX(di,3) = 1-alpha*(m_u(ii,jj+1)+abs(m_u(ii,jj+1))-(m_u(ii,jj)-abs(m_u(ii,jj))))
     MATRIX(di,3) = MATRIX(di,3) -beta*(m_v(ii+1,jj)+abs(m_v(ii+1,jj))-(m_v(ii,jj)-abs(m_v(ii,jj))))
   ENDDO

   ! Build other for diagonals
! D2
   DO di=2, dimen, 1
     CALL GETGRIDLOCATION(di, MX, MY, ii, jj)
     MATRIX(di,2) = alpha*(m_u(ii,jj)+abs(m_u(ii,jj)))
   ENDDO
   
! D3
   DO di=1, dimen-1, 1
     CALL GETGRIDLOCATION(di, MX, MY, ii, jj) 
     MATRIX(di,4) = -alpha*(m_u(ii,jj+1)-abs(m_u(ii,jj+1)))
   ENDDO

! D4
   DO di=MY+1,dimen, 1
     CALL GETGRIDLOCATION(di, MX, MY, ii, jj) 
     MATRIX(di,1) = beta*(m_v(ii,jj)+abs(m_v(ii,jj)))
   ENDDO

! D5
   DO di=1, dimen-MY, 1
     CALL GETGRIDLOCATION(di, MX, MY, ii, jj)
     MATRIX(di,5) = -beta*(m_v(ii+1,jj)-abs(m_v(ii+1,jj)))
   ENDDO

  END FUNCTION 


  ! Fixed for 2D
  SUBROUTINE ANTIDIF(u, v, u_a, v_a, psi, eps, dx, dy, dt, MX, MY)
   REAL :: eps, dx, dy, dt
   REAL, DIMENSION(:) :: psi
   REAL, DIMENSION(MX+1, MY) :: u, u_a
   REAL, DIMENSION(MX, MY+1) :: v, v_a 
   INTEGER :: ii, jj, MX, MY, fi, se
	
   ! Build antidiffusion vector u_a
   DO jj=1, MY
     DO ii=1, MX-1
	CALL GETVECTORLOCATION(fi, MX, MY, ii+1, jj) 
	CALL GETVECTORLOCATION(se, MX, MY, ii, jj) 
	u_a(ii+1,jj) = ((abs(u(ii+1,jj))*dx-dt*u(ii+1,jj)*u(ii+1,jj))*(psi(fi)-psi(se)))/((psi(fi)+psi(se)+eps)*dx)
	!IF(abs(u_a(ii,jj)-u(ii,jj)) .GT. u(ii,jj)) THEN	
	!PRINT*, u_a(ii,jj), " ", psi(fi), " ", psi(se), " ", u(ii,jj), " ", ii, " ", jj, " ", fi, " ", se, " "
	!ENDIF
     ENDDO
   ENDDO
   ! Build antidiffusion vector v_a
   DO ii=1, MX
     DO jj=1, MY-1 
	CALL GETVECTORLOCATION(fi, MX, MY, ii, jj+1) 
	CALL GETVECTORLOCATION(se, MX, MY, ii, jj) 
	v_a(ii,jj+1) =  ((abs(v(ii,jj+1))*dy-dt*v(ii,jj+1)*v(ii,jj+1))*(psi(fi)-psi(se)))/((psi(fi)+psi(se)+eps)*dy)
	!PRINT*, v_a(ii,jj), " ", psi(fi), " ", psi(se), " ", v(ii,jj), " ", ii, " ", jj, " ", fi, " ", se, " "
     ENDDO
  ENDDO
  END SUBROUTINE
 
  ! This subroutine multiplies a matrix A with a vector x and returns their product y
  ! Same for 2D as for 1D
  SUBROUTINE MVEC(A,x,y,MX,MY)
   REAL, DIMENSION(:) :: x, y
   REAL, DIMENSION(:,:) :: A
   INTEGER :: MX, MY
   y(1) = A(1,3) * x(1) + A(1,4) * x(2) + A(1,5) * x(MY) 
   DO i=2, MY
    y(i) = A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1) + A(i,5) * x(i+MY)
   ENDDO  
   DO i = MY+1, (MX-1)*MY
    y(i) = A(i,1) * x(i-MY) + A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1) + A(i,5) * x(i+MY)
   ENDDO
   DO i = (MX-1)*MY+1, MX*MY
    y(i) = A(i,1) * x(i-MY) + A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1)
   ENDDO
  END SUBROUTINE

  SUBROUTINE GETGRIDLOCATION(di, MX, MY, ii, jj)
   !ii, jj are locations in u and v matrices, di is location in gridvector
   integer :: di, MX, MY, ii, jj
   ii = (di-1) / MX + 1
   jj = mod(di, MX) 
   if(jj.EQ.0) jj = MX
  END SUBROUTINE

  SUBROUTINE GETVECTORLOCATION(di, MX, MY, ii, jj)
   !ii, jj are locations in u and v matrices, di is location in gridvector
   integer :: di, MX, MY, ii, jj
   di = (ii-1)*MX+jj
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
      ALLOCATE(A((MX*MY),5), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
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
