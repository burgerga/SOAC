! Implementation of Smolarkiewicz Positive Definite Advection Scheme with Small Implicit Diffusion (1983)
! G.A. Burger & J.M. Wolterink, October 2011
! 2 dimensional simulation
PROGRAM Smolarkiewicz2D
 IMPLICIT NONE

! Variables 
 INTEGER :: count, ierror = 0, M, N, i, j, k, q, cone = 0, iterations, initial_pos,cloudsize, hc, xor, yor 
 REAL :: dx, dy, dt, eps, uv, sc, u, v, angvel = 0.1
 REAL, DIMENSION(:,:), ALLOCATABLE :: grid, A, m_u, m_v, m_u_a, m_v_a, initialmat
 REAL, DIMENSION(:), ALLOCATABLE :: initial, psi_int, psi_tem
 character(len=100) :: input_file

! Check whether there is an input file
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

! Read input file
 call read_input_file(input_file)
 
! Open output file
 if(iterations .EQ. 0) then
 OPEN(20, file = "wave0.dat")
 elseif(iterations .EQ. 1) then
 OPEN(20, file = "wave1.dat")
 elseif(iterations .EQ. 2) then
 OPEN(20, file = "wave2.dat")
 elseif(iterations .EQ. 3) then
 OPEN(20, file = "wave3.dat")
 endif
 OPEN(80, file = "mua.dat")
 OPEN(90, file = "mva.dat")
! Initialize system. Initial is an M*M vector which represents the 2D space.
 initial = 0 
 hc = cloudsize/2
 DO j = 1 , cloudsize 
 
	initial(initial_pos-(j-1)*M-hc:initial_pos-(j-1)*M+hc) = 1
 ENDDO
 !initial(initial_pos-(cloudsize/2):initial_pos+(cloudsize/2)) = 1 


! Set velocity matrices
 m_u = u
 m_v = v
 PRINT*, cone
 IF(cone .EQ. 1) THEN
  xor = 25*dx
  yor = 50*dy
  !dx = 1 
  !dy = 1
  !dt = 0.1
  PRINT*, 'Cone!' 
  initial = 0 
  ALLOCATE(initialmat(M,M), STAT=ierror); IF (ierror /= 0) PRINT*, "initial : Allocation failed"
  DO j = 1, M 
	DO k = 1, M
		initialmat(j,k) = 3.87-(0.3/dx)*(sqrt(abs(xor-REAL(j*dx))**2+abs(yor-REAL(k*dx))**2))
		if(initialmat(j,k)<0) then
			initialmat(j,k) = 0 		
		endif
        ENDDO
  ENDDO
  DO j = 1, M+1 
	DO k = 1, M
		m_u(j,k) = -angvel*(j*dx-50*dx)
        ENDDO
  ENDDO
  DO j = 1, M 
	DO k = 1, M+1
		m_v(j,k) = angvel*(k*dx-50*dy)
        ENDDO
  ENDDO
  DO j = 1, M
	initial((j-1)*M+1:j*M) = initialmat(j,:)
  ENDDO
  DEALLOCATE(initialmat,STAT=ierror)
 ENDIF
  psi_int = initial
! Repeat for N time steps 
 DO j = 1, N-1

! Build five-diagonal multiplication matrix
  A =  MATRIX(m_u, m_v, dx, dy, dt, M)

! Make copy of current situation
  psi_tem = psi_int  

! Multiply vector to get new situation
  CALL MVEC(A, psi_tem, psi_int,M)

! Apply antidiffusion step for some iterations (optional, iterations could be 0)
  DO k = 1, iterations

! Make local copy of current situation
	  psi_tem = psi_int
! Build antidiffusion velocity matrices u~ and v~
	  CALL ANTIDIF(m_u, m_v, m_u_a, m_v_a, psi_tem, eps, dx, dy, dt, M)
	m_u_a = m_u_a * sc
        m_v_a = m_v_a * sc
  IF(MOD(j, 10).EQ.0) THEN
	  WRITE(80,*), reshape(m_u_a,(/1, (M+1)*M/))
          WRITE(90,*), reshape(m_v_a,(/1, M*(M+1)/))
  ENDIF
!print*, maxval(m_v_a,2)
	  !PRINT*, m_u_a
! Build new multiplication matrix with antidiffusion vector
	  A =  MATRIX(m_u_a, m_v_a, dx, dy, dt, M) 
! Multiply
	  CALL MVEC(A, psi_tem, psi_int,M)

  ENDDO
! Write new situation to file
  IF(MOD(j, 10).EQ.0) THEN
	  WRITE(20,*), psi_int
  ENDIF
 ENDDO

! Close output file
 CLOSE(20)
 CLOSE(80)
 CLOSE(90)

! Deallocation
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

! This function builds the multiplication matrix based on the velocity matrices u, v
! or the antidiffusion velocity matrices u~, v~. It returns a five-diagonal sparse matrix.
FUNCTION MATRIX(m_u, m_v, dx, dy, dt, M)
   REAL :: dx, dy, dt, div, alpha, beta
   REAL, DIMENSION(:,:) :: m_u, m_v
   REAL, DIMENSION(M*M,5) :: MATRIX
   INTEGER :: dimen, M, ii, jj, di

   MATRIX = 0 
   dimen = M*M
   alpha = dt/(2*dx)
   beta = dt/(2*dy)
! ii and jj is location in grid, rows are ii columns are jj. eg. (4,3) is fourth row, third column

   DO di=1, dimen, 1
! 1D coordinate to 2D
     CALL GETGRIDLOCATION(di, M, ii, jj)
     MATRIX(di,3) = 1-alpha*(m_u(ii,jj+1)+abs(m_u(ii,jj+1))-(m_u(ii,jj)-abs(m_u(ii,jj))))
     MATRIX(di,3) = MATRIX(di,3) -beta*(m_v(ii+1,jj)+abs(m_v(ii+1,jj))-(m_v(ii,jj)-abs(m_v(ii,jj))))
   ENDDO

! Build other four diagonals
   DO di=2, dimen, 1
     CALL GETGRIDLOCATION(di, M, ii, jj)
     MATRIX(di,2) = alpha*(m_u(ii,jj)+abs(m_u(ii,jj)))
   ENDDO
   
   DO di=1, dimen-1, 1
     CALL GETGRIDLOCATION(di, M, ii, jj) 
     MATRIX(di,4) = -alpha*(m_u(ii,jj+1)-abs(m_u(ii,jj+1)))
   ENDDO

   DO di=M+1+1,dimen, 1
     CALL GETGRIDLOCATION(di, M, ii, jj) 
     MATRIX(di,1) = beta*(m_v(ii,jj)+abs(m_v(ii,jj)))
   ENDDO

   DO di=1, dimen-M-1, 1
     CALL GETGRIDLOCATION(di, M, ii, jj)
     MATRIX(di,5) = -beta*(m_v(ii+1,jj)-abs(m_v(ii+1,jj)))
   ENDDO
   MATRIX(1,:) = 0
   MATRIX(M*M,:) = 0
  END FUNCTION 



! Build antidiffusion matrices
  SUBROUTINE ANTIDIF(u, v, u_a, v_a, psi, eps, dx, dy, dt, M)
   REAL :: eps, dx, dy, dt, up, vp, psipsum
   REAL, DIMENSION(:) :: psi
   REAL, DIMENSION(M+1, M) :: u, u_a
   REAL, DIMENSION(M, M+1) :: v, v_a 
   INTEGER :: ii, jj, M, fi, se
   REAL, DIMENSION(:), ALLOCATABLE :: psips
   ALLOCATE(psips(M), STAT=ierror); IF (ierror /= 0) PRINT*, "psips : Allocation failed"
   ! Build antidiffusion matrix u_a
   u_a = 0 ; v_a = 0 ;
   DO ii = 1, M+1
	!PRINT*, ii
	psips = psi((ii-1)*M+1:ii*M)

        DO jj = 1, M-1
		up  = u(ii+1,jj)
		psipsum = psips(jj)+psips(jj+1)
		if(psipsum>0.0001) then
		u_a(ii+1,jj) = ((abs(up)*dx-dt*up*up)*(psips(jj+1)-psips(jj)))/(psipsum*dx)
		else
		u_a(ii+1,jj) = 0
		endif
  	ENDDO
   ENDDO
   DO ii = 1, M+1
	!PRINT*, ii
	psips = psi(ii:(M-1)*M+ii:M)

	DO jj = 1, M
		vp = v(ii, jj+1)
		psipsum = psips(jj)+psips(jj+1) 
		if(psipsum>0.0001) then
		v_a(ii,jj+1) = ((abs(vp)*dy-dt*vp*vp)*(psips(jj+1)-psips(jj)))/(psipsum*dy)
		else
		v_a(ii,jj+1) = 0
		endif
	ENDDO
   ENDDO
  END SUBROUTINE
 
! This subroutine multiplies a sparse matrix A with a vector x and returns their product y
  SUBROUTINE MVEC(A,x,y,M)
   REAL, DIMENSION(:) :: x, y
   REAL, DIMENSION(:,:) :: A
   INTEGER :: M
   y(1) = A(1,3) * x(1) + A(1,4) * x(2) + A(1,5) * x(M) 
   DO i=2, M
    y(i) = A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1) + A(i,5) * x(i+M)
   ENDDO  
   DO i = M+1, (M-1)*M
    y(i) = A(i,1) * x(i-M) + A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1) + A(i,5) * x(i+M)
   ENDDO
   DO i = (M-1)*M+1, M*M
    y(i) = A(i,1) * x(i-M) + A(i,2) * x(i-1) + A(i,3) * x(i) + A(i,4) * x(i+1)
   ENDDO
   y(1:M) = 0
   y((M-1)*M+1:M*M) = 0
   DO i = 2 , M-1 
	y((i-1)*M+1) = 0
	y(i*m) = 0
   ENDDO
  END SUBROUTINE

! Translate vector to matrix coordinate
  SUBROUTINE GETGRIDLOCATION(di, M, ii, jj)
   !ii, jj are locations in u and v matrices, di is location in gridvector
   integer :: di, M, ii, jj
   ii = (di-1) / M + 1
   jj = mod(di, M) 
   if(jj.EQ.0) jj = M
  END SUBROUTINE

! Translate matrix to vector coordinate
  SUBROUTINE GETVECTORLOCATION(di, M, ii, jj)
   !ii, jj are locations in u and v matrices, di is location in gridvector
   integer :: di, M, ii, jj
   di = (ii-1)*M+jj
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
     case ('M')
      read(buffer,*, iostat=ios) M
      ALLOCATE(initial(M*M), STAT=ierror); IF (ierror /= 0) PRINT*, "initial : Allocation failed"
      ALLOCATE(psi_int(M*M), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_int : Allocation failed"
      ALLOCATE(psi_tem(M*M), STAT=ierror); IF (ierror /= 0) PRINT*, "psi_tem : Allocation failed"
      ALLOCATE(m_u((M+1),M), STAT=ierror); IF (ierror /= 0) PRINT*, "u : Allocation failed"
      ALLOCATE(m_v(M,(M+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v : Allocation failed"
      ALLOCATE(m_u_a((M+1),M), STAT=ierror); IF (ierror /= 0) PRINT*, "u_a : Allocation failed"
      ALLOCATE(m_v_a(M,(M+1)), STAT=ierror); IF (ierror /= 0) PRINT*, "v_a : Allocation failed"
      ALLOCATE(A((M*M),5), STAT=ierror); IF (ierror /= 0) PRINT*, "A : Allocation failed"
     case ('N')
      read(buffer,*, iostat=ios) N
      ALLOCATE(grid(N,(M*M)), STAT=ierror); IF (ierror /= 0) PRINT*, "grid : Allocation failed"
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
     case ('cone')
      read(buffer,*, iostat=ios) cone
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
     case ('cloudsize')
      read(buffer,*, iostat=ios) cloudsize
     case default
      if (.not.(label(1:1) == ' ' .or. label(1:1) == '!')) then
       print*, 'Skipping invalid label at line ', line
      end if 
     end select
    end if
   end do
  end subroutine read_input_file
END PROGRAM 
