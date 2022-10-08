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

	subroutine iterations_method()
		implicit none
		real(8) :: d(M-1), e(M-2), H(M-1,M-1), H_(M-1,M-2)
		real(8) :: E_max, phi_1(M-1), phi_2(M-1)

		do i = 1, ! <---


	end subroutine iterations_method
	
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