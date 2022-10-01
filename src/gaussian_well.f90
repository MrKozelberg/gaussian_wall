module lapackMod

   implicit none

   public

   integer :: i

   integer :: M  ! Amount of points
   double precision :: R  ! Grid width
   double precision :: delta  ! Grid spacing

   double precision :: x_0  ! The left boundary point
   double precision :: x_M  ! The rigth boundary point
   double precision, allocatable :: x ! Grid points

   double precision, allocatable :: diag_elem != [(2 / delta**2 + v_0 * exp(-x(i)**2), i = 1, M-1)]
   ! The M-1 diagonal elements of the tridiagonal matrix

   double precision, allocatable :: subdiag_elem
   ! The M-2 subdiagonal elements of the tridiagonal matrix

   double precision, allocatable :: w(:), z(:, :)
   double precision :: v_0
   integer :: num_eigen

contains

   subroutine init(v0, m_, r_)

      double precision, intent(in) :: v0, r_
      integer, intent(in) :: m_

      M = m_
      R = r_; x_0 = -R; x_M = R
      delta = 2 * R / M
      x = [(x_0 + i * delta, i = 1, M-1)]

      v_0 = v0
      diag_elem = [(2 / delta**2 + v_0 * exp(-x(i)**2), i = 1, M-1)]
      subdiag_elem = [(-1 / delta**2, i = 1, M-2)]

   end subroutine init

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

      jobz = 'N'  ! Compute eigenvalues only
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

   subroutine solver()
      double precision :: d(M-1), e(M-2)
      double precision :: vl, vu

      d = diag_elem; e = subdiag_elem
      vl = v_0; vu = 0.0

      call sym_3diag_eigen(d, e, vl, vu, num_eigen, w, z)
   end subroutine solver

   subroutine num_eigen_vs_v0()
      double precision, dimension(70) :: v0_ = [(-i/3.5, i=1, 70)]

      double precision :: rr = 15.0
      integer :: io

      open(newunit=io, file="../output/num_eigen_vs_v0.txt")

      write(io, *) '# v0 # num_eigen' 
      
      do i = 1, size(v0_)
         call init(v0_(i), 10000, rr)

         call solver()

         write(io, *) v0_(i), num_eigen

      end do

      close(io)

   end subroutine num_eigen_vs_v0

   subroutine time_vs_M()
      integer, dimension(6) :: mm = [(100*2**i, i=1, 6)]

      double precision :: rr = 6.0, v0 = -1.0
      real(4) :: start_time, finish_time
      integer :: io

      open(newunit=io, file="../output/time_vs_M.txt")

      write(io, *) '# M # T, sec' 

      do i = 1, size(mm)

         call cpu_time(start_time)

         call init(v0, mm(i), rr)

         call solver()

         call cpu_time(finish_time)

         write(io, *) M, finish_time - start_time

      end do

      close(io)

   end subroutine time_vs_M

end module lapackMod

program main
   use lapackMod
	implicit none

   call num_eigen_vs_v0()

end program
