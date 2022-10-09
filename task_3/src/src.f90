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

	subroutine iterations_method(iter_num)
		implicit none
		integer, intent(in) :: iter_num
		real(8) :: H(M-1,M-1), H_(M-1,M-1), id(M-1,M-1)
		real(8) :: E_max, E_min, phi_1(M-1), phi_2(M-1)

		! Hamiltonian
		H(:,:) = 0.0
		do i = 1, M-1
			H(i,i) = diag(i)
			if (i .eq. 1) then
				H(i,i+1) = subdiag(i)
			else if (i .eq. M-1) then
				H(i-1,i) = subdiag(i)
			else
				H(i,i+1) = subdiag(i)
				H(i-1,i) = subdiag(i)
			end if
		end do

		! phi_1 = [(exp(-x(i)**2), i=1, M-1)]
		! phi_1 = phi_1 / norm2(phi_1)
		phi_1(:) = 1/delta
		do i = 1, 1000
			phi_2 = matmul(H, phi_1)
			phi_2 = phi_2 / norm2(phi_2)
			phi_1 = phi_2
		end do

		E_max = dot_product(matmul(H, phi_1), phi_1)

		print *, 'E_max is ', E_max

		id(:,:) = 0.0
		do i = 1, M-1
			id(i,i) = 1.0
		end do

		H_ = H - E_max/2.0 * id

		phi_1 = [(1/cosh(x(i))**(2), i=1, M-1)]
		phi_1 = phi_1 / norm2(phi_1)
		do i = 1, iter_num
			phi_2 = matmul(H_, phi_1)
			phi_2 = phi_2 / norm2(phi_2)
			phi_1 = phi_2
		end do

		E_min = dot_product(matmul(H_, phi_1), phi_1)

		print *, 'E_1 is ', E_min + E_max/2.0


	end subroutine iterations_method
	
   	! subroutine T_vs_M()
   	! 	integer, dimension(300) :: mm = [(100*i, i=1, 300)]

    !   	real(8) :: rr = 6.0, v0 = -1.0
    !   	real(4) :: start_time, finish_time
    !   	integer :: io

    !   	open(newunit=io, file="../output/T_vs_M.txt")

    !   	write(io, *) '# M # T, sec' 

    !   	do i = 1, size(mm)

    !      	call cpu_time(start_time)

    !      	call init(v0, rr, mm(i))

    !      	call solver()

    !      	call cpu_time(finish_time)

    !      	call dealloc()

    !      	write(io, *) M, finish_time - start_time

    !   	end do

    !   	close(io)
   	! end subroutine T_vs_M

end module task

program src
	use task
	implicit none

	! test
	real(8) :: v0, rr
	integer :: mm

	v0 = -1.0; rr = 6.0; mm = 1000

	call init(v0, rr, mm)

	call iterations_method(10000)

	! write(*, *) 'N = ', N

	! do i = 1, N
	! 	write(*, *) i, w(i)
	! end do

	! First part of the task
	! call N_vs_v0()

	! Second part of the task
	! call T_vs_M()

	! Third part of the task
	! call E_vs_M_vs_R(real(-1.0,8),"-1")
end program src