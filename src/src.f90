module task
	implicit none

	integer :: i, j, M, N
	real(8) :: R, delta, x_0, x_M, v_0
	real(8), allocatable :: x(:), diag(:), subdiag(:), w(:), z(:,:)

contains
	
	subroutine init(v0, rr, mm)
		real(8), intent(in) :: v0, rr
		integer, intent(in) :: mm

		v_0 = v0; R = rr; M = mm

		delta = 2*R/M
		x_0 = -R; x_M = R

		if (.not. allocated(x)) allocate(x(M-1)) 
		x = [(x_0+i*delta, i = 1, M-1)]

		if (.not. allocated(diag)) allocate(diag(M-1))
		diag = [(2/delta**2+v_0*exp(-x(i)**2), i = 1, M-1)]

		if (.not. allocated(subdiag)) allocate(subdiag(M-2))
		subdiag(:) = -1/delta**2

		if (.not. allocated(w)) allocate(w(M-1))
		if (.not. allocated(z)) allocate(z(M-1,M-1))
	end subroutine init

	subroutine dealloc()
		deallocate(x)
		deallocate(diag)
		deallocate(subdiag)
		deallocate(w)
		deallocate(z)
	end subroutine dealloc

	subroutine eigen_stdm(d, e, vl, vu, n, w, z)
		implicit none (type, external)

      	external :: DSTEVX
    	real(8), external :: DLAMCH

    	real(8), intent(in) :: d(:), e(:), vl, vu
      	real(8), intent(out) :: w(:), z(:,:)
      	integer, intent(out) :: n  ! Number of eigenvalues have been found

      	character :: jobz, range
      	integer :: info, ldz, il, iu
      	real(8) :: abstol

      	integer :: ifail(size(d)), iwork(5*size(d))
      	real(8) :: work(size(iwork))

      	jobz = 'N'  ! Compute eigenvalues only
      	range = 'V'  ! All eigenvalues in the half-open interval (VL,VU] will be found

      	abstol = 2*DLAMCH('s')

      	ldz = size(z, 1)

      	write(*,*)
      	write(*,*) 'Eigenmodes calculation started ...'

      	call DSTEVX(jobz, range, size(d), d, e, vl, vu, il, iu, &
                  	abstol, n, w, z, ldz, work, iwork, ifail, info)

      	if(info.ne.0) then

        	write(*,'(a)') 'SGEEV LAPACK - Failure!'
         	write(*,'(a)') 'The eigenmodes could not be calculated.'
         	write(*,'(a)') 'LAPACK routine SGEEV returned a nonzero'
         	write(*,'(a,i8)') 'Value of the error flag, INFO=',info
         	stop

      	else

         	write(*,*) 'Eigenmodes calculation finished successfully.'

      	end if
	end subroutine eigen_stdm

	subroutine solver()
      	real(8) :: d(M-1), e(M-2), vl, vu

      	d = diag; e = subdiag
      	vl = v_0; vu = 0.0

      	call eigen_stdm(d, e, vl, vu, N, w, z)
   	end subroutine solver

   	subroutine N_vs_v0()
   		real(8), dimension(70) :: v0_ = [(-i/3.5, i=1, 70)]

      	real(8) :: rr = 15.0
     	integer :: io, mm

     	mm = 10000

      	open(newunit=io, file="../output/N_vs_v0.txt")

      	write(io, *) '# v0 # num_eigen' 
      
      	do i = 1, size(v0_)

         	call init(v0_(i), rr, mm)

         	call solver()

         	write(io, *) v0_(i), N

      	end do

      	close(io)
   	end subroutine N_vs_v0

   	subroutine T_vs_M()
   		integer, dimension(300) :: mm = [(100*i, i=1, 300)]

      	real(8) :: rr = 6.0, v0 = -1.0
      	real(4) :: start_time, finish_time
      	integer :: io

      	open(newunit=io, file="../output/T_vs_M.txt")

      	write(io, *) '# M # T, sec' 

      	do i = 1, size(mm)

         	call cpu_time(start_time)

         	call init(v0, rr, mm(i))

         	call solver()

         	call cpu_time(finish_time)

         	call dealloc()

         	write(io, *) M, finish_time - start_time

      	end do

      	close(io)
   	end subroutine T_vs_M

   	subroutine E_vs_M_vs_R(v0, ch)
   		integer, dimension(150) :: mm = [(200*i, i=1, 150)]
   		real(8), dimension(150) :: rr = [(2+real(i,8)/15.0, i=1, 150)]

   		real(8), intent(in) :: v0
   		character, intent(in) :: ch

   		! We consider only first three levels
   		real(8), dimension(size(mm)) :: e1, e2, e3
   		integer :: io1, io2, io3, ior, iom

   		open(newunit=ior, file="../output/r_"//ch//".txt")

   		do i = 1, size(rr)
   			write(ior, *) rr(i)
   		end do

   		close(ior)

   		open(newunit=iom, file="../output/m_"//ch//".txt")

   		do i = 1, size(mm)
   			write(iom, *) mm(i)
   		end do

   		close(iom)

   		open(newunit=io1, file="../output/e1_"//ch//".txt")
   		open(newunit=io2, file="../output/e2_"//ch//".txt")
   		open(newunit=io3, file="../output/e3_"//ch//".txt")

   		do i = 1, size(rr)
   			e1(:) = 0.0; e2(:) = 0.0; e3(:) = 0.0

   			do j = 1, size(mm)
   				call init(v0, rr(i), mm(j))

   				call solver()

   				if (N .eq. 0) then

   					e1(j) = 1
   					e2(j) = 1
   					e3(j) = 1

   				else if (N .eq. 1) then

   					e1(j) = w(1)
   					e2(j) = 1
   					e3(j) = 1

   				else if (N .eq. 2) then

   					e1(j) = w(1)
   					e2(j) = w(2)
   					e3(j) = 1

   				else

   					e1(j) = w(1)
   					e2(j) = w(2)
   					e3(j) = w(3)

   				end if

   				call dealloc()

   			end do

   			write(io1, *) e1(:)
   			write(io2, *) e2(:)
   			write(io3, *) e3(:)

   		end do
      	
      	close(io3)
      	close(io2)
      	close(io1)

   	end subroutine E_vs_M_vs_R

end module task

program src
	use task
	implicit none

	! test
	! real(8) :: v0, rr
	! integer :: mm

	! v0 = -5.0; rr = 15.0; mm = 30000

	! call init(v0, rr, mm)

	! call solver()

	! write(*, *) 'N = ', N

	! do i = 1, N
	! 	write(*, *) i, w(i)
	! end do

	! First part of the task
	! call N_vs_v0()

	! Second part of the task
	! call T_vs_M()

	! Third part of the task
	call E_vs_M_vs_R(real(-1.0,8),"-1")
end program src