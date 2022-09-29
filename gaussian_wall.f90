module lapackMod

contains

   subroutine sym_3diag_eigen(d, e, vl, vu, m, w, z)
      implicit none (type, external)

      external :: DSTEVX
      double precision, external :: DLAMCH

      double precision, intent(in) :: d(:), e(:), vl, vu
      double precision, intent(out) :: w(:), z(:,:)
      integer, intent(out) :: m

      character :: jobz, range
      integer :: info, ldz, n, il, iu
      double precision :: abstol

      integer :: ifail(size(d)), iwork(5*size(d))
      double precision :: work(size(iwork))

      jobz = 'V'  ! Compute eigenvalues and eigenvectors
      range = 'V'  ! All eigenvalues in the half-open interval (VL,VU] will be found

      n = size(d)

      abstol = 2*DLAMCH('s')

      ldz = size(z, 1)

      write(*,*)
      write(*,*) 'Eigenmodes calculation started ...'

      call DSTEVX(jobz, range, n, d, e, vl, vu, il, iu, &
                  abstol, m, w, z, ldz, work, iwork, ifail, info)

      if(info.ne.0) then

         write(*,'(a)') 'SGEEV LAPACK - Failure!'
         write(*,'(a)') 'The eigenmodes could not be calculated.'
         write(*,'(a)') 'LAPACK routine SGEEV returned a nonzero'
         write(*,'(a,i8)') 'Value of the error flag, INFO=',info
         stop

      else

         write(*,*) 'Eigenmodes calculation finished successfully.'

      end if

   end subroutine sym_3diag_eigen

end module lapackMod

program main
   use lapackMod
	implicit none

   integer :: i

   integer, parameter :: M = 1500  ! Amount of points
   double precision, parameter :: R = 5.0 ! Grid width
   double precision, parameter :: delta = 2 * R / M  ! Grid spacing

   double precision, parameter :: x_0 = -R  ! The left boundary point
   double precision, parameter :: x_M = R  ! The rigth boundary point
   double precision, parameter, dimension(M-1) :: x = [(x_0 + i * delta, i = 1, M-1)] ! Grid points

   double precision, parameter :: v_0 = -10.76  ! Parameter of the potential function (v_0 < 0)

   double precision, parameter, dimension(M-1) :: diag_elem = [(2 / delta**2 + v_0 * exp(-x(i)**2), i = 1, M-1)]
   ! The M-1 diagonal elements of the tridiagonal matrix

   double precision, parameter, dimension(M-2) :: subdiag_elem = [(-1 / delta**2, i = 1, M-2)]
   ! The M-2 subdiagonal elements of the tridiagonal matrix

   double precision :: d(M-1), e(M-2), w(M-1), z(M-1, M-1)
   double precision :: vl, vu
   integer :: num_eigen

   d = diag_elem; e = subdiag_elem
   vl = v_0; vu = 0.0

   call sym_3diag_eigen(d, e, vl, vu, num_eigen, w, z)

   write (*, *)
   write (*, *) "Number of eigenvalues found =", num_eigen
   
   if (num_eigen .gt. 0) then
      write (*, *) 'Eigenvalue number  Eigenvalue'
      
      do i = 1, num_eigen
         write (*, *) i, w(i)
      end do
   end if

end program
