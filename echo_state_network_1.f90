program echo_state_network
	implicit none
	
	!general constants and variables
	real (8) :: pi
	integer :: i, j, k, l
	!lorenz system
	real (8) :: sigma, r, beta
	real (8) :: dt, total_t, t
	real (8) :: x_0, y_0, z_0
	real (8), dimension (:), allocatable :: x, y, z
	integer :: steps
	!downsampling 
	real (8), dimension(:,:), allocatable :: u_k
	integer :: downsample_length
	!echo state network variables, matrices
	real (8), dimension (:,:), allocatable :: w, v, x_i, Res_x, Res_y, Res_z,  y_i
	real (8) :: alpha, c, a
	real (8) :: w_sum, v_sum
	integer :: n
	
	!define general constants 
	pi = 4.d0*datan(1.d0)
	
	!generate Lorenz system trajectory
	
	sigma = 10.d0
	r = 28.d0
	beta = 8.d0/3.d0
	
	x_0 = 1
	y_0 = 1
	z_0 = 1
	
	dt = 0.005 !seconds
	total_t = 1000 !seconds
	steps = int(total_t/dt) 
	
	allocate(x(steps), y(steps), z(steps))
	
	call iterate_heun(x,y,z)
	
	!downsample to dt = 0.02, and combine trajectories to form u_k 
	
	downsample_length = int(steps/4)
	
	allocate(u_k(3,downsample_length))
		
	do i = 1, downsample_length
		u_k(1,i) = x(i*4)
		u_k(2,i) = y(i*4)
		u_k(3,i) = z(i*4)
	end do
	
	deallocate(x, y, z)
	
	!generate reservoir (w) and input (v) matrices
	
	n = 50 !reservoir nodes
	
	allocate(w(n,n), v(3,n))
	
	call generate_reservoir_matrix(w)
	call normalize_spectral_radius(w)
	call generate_input_matrix(v)
	
	!define other esn components, variable, starting values, iterate esn. 
	
	allocate(x_i(downsample_length, N))
		
	alpha = 0.005
	c = 1
	a = 1
	
	x_i = 1
	
	w_sum = 0
	v_sum = 0
	
	do l = 1, downsample_length - 1 
		do i = 1, N
			do j = 1, N
				w_sum = w(j,i)*x_i(l, j)
			end do
		
			do k = 1, 3
				v_sum = v(k,i)*u_k(k, l)
			end do
		
			x_i(l+1, i)= (1-alpha)*x_i(l,i)+ &
				alpha*tanh(c*w_sum + a*v_sum)
		
			w_sum = 0
			v_sum = 0
		end do
	end do
	
	!esn readout
	
	open (10,file='y_i.dat') 
	
	allocate(y_i(3, downsample_length), Res_x(downsample_length,N), &
				Res_y(downsample_length,N), Res_z(downsample_length,N))
	
	y_i = 0
	
	do i = 3000, downsample_length - 5
		do j = 1, N
			Res_x(i,j) = u_k(1, i )/x_i(i,j)	
					!write(*,*) u_k(1,i+5), Res(i,j)*x_i(i,j), x_i(i,j)
		end do

	end do
	
	do i = 3000, downsample_length - 5
		do j = 1, N
			Res_x(i,j) = u_k(1, i + 5)/x_i(i,j)	
					!write(*,*) u_k(1,i+5), Res(i,j)*x_i(i,j), x_i(i,j)
		end do

	end do
	
	do i = 3000, downsample_length - 5
		do j = 1, N
			Res_y(i,j) = u_k(2, i + 5)/x_i(i,j)	
					!write(*,*) u_k(1,i+5), Res(i,j)*x_i(i,j), x_i(i,j)
		end do

	end do
	
	do i = 3000, downsample_length - 5
		do j = 1, N
			Res_z(i,j) = u_k(3, i + 5)/x_i(i,j)	
					!write(*,*) u_k(1,i+5), Res(i,j)*x_i(i,j), x_i(i,j)
		end do

	end do

	
	do i = 3000, downsample_length
		do j = 1, N
			y_i(1,i) = y_i(1,i) + x_i(i, j)*Res_x(i,j)
			y_i(2,i) = y_i(2,i) + x_i(i, j)*Res_y(i,j)
			y_i(3,i) = y_i(3,i) + x_i(i, j)*Res_z(i,j)	
		end do
		write(10,*) y_i(1,i)/50, y_i(2,i)/50, y_i(3,i)/50
	end do
	
	
		
	contains
	
	!start: lorenz and heun
	
	real function dot_u_1(u_1, u_2)
		real (8), intent (in) ::u_1, u_2
		dot_u_1 = sigma*(u_2 - u_1)
	end function dot_u_1
	
	real function dot_u_2(u_1, u_2, u_3)
		real (8) :: u_1, u_2, u_3
		dot_u_2 = u_1*(r-u_3)-u_2
	end function dot_u_2
		
	real function dot_u_3(u_1, u_2, u_3)
		real (8) :: u_1, u_2, u_3
		dot_u_3 = u_1*u_2 - beta*u_3
	end function dot_u_3
	
	subroutine iterate_heun(x, y, z)
		real (8), dimension (:), intent (inout) :: x, y, z
		real (8) :: p_kx, p_ky, p_kz, h_kx, h_ky, h_kz
		integer :: i
		
		h_kx = x_0
		h_ky = y_0
		h_kz = z_0	
		
		open (10,file='trajectory.dat') 
		
		do i = 1, steps 
					
			p_kx = h_kx + dot_u_1(p_kx, p_ky)*dt
			p_ky = h_ky + dot_u_2(p_kx, p_ky, p_kz)*dt
			p_kz = h_kz + dot_u_3(p_kx, p_ky, p_kz)*dt
			
			h_kx = h_kx + (dot_u_1(h_kx, h_ky) + dot_u_1(p_kx, p_ky))*dt/2
			h_ky = h_ky + (dot_u_2(h_kx, h_ky, h_kz) + dot_u_2(p_kx, p_ky, p_kz))*dt/2
			h_kz = h_kz + (dot_u_3(h_kx, h_ky, h_kz)+ dot_u_3(p_kx, p_ky, p_kz))*dt/2
			
			x(i) = h_kx
			y(i) = h_ky
			z(i) = h_kz
			
			write(10,*) h_kx, h_ky, h_kz
			
		end do		
		close (10)
	end subroutine iterate_heun
	
	!end: lorenz and huen 
	
	!start generate reservoir (w) and input (v) matrices
	
	subroutine generate_reservoir_matrix(W)
		real (8), dimension (:,:), intent (inout) :: W
		real (8) :: p, weight
		
		w = 0
		weight = 2
		
		do i = 1, N
			do j = 1, N
				call random_number(p)
				if (p.lt.0.05) then
					do while (abs(weight).gt.1)
						weight =  norm_rand(0.5, 1.0)/sqrt(2*pi)
						
							
					end do									
					W(i,j) = abs(weight)
					weight = 2			
				end if
			end do
			
		end do	
	end subroutine generate_reservoir_matrix
	
		subroutine normalize_spectral_radius(w)
		
		real (8), dimension(:,:), intent(inout) :: w
		character :: jobvl, jobvr
		real (8), dimension(:), allocatable :: work, wr, wi, magnitudes
		real (8), dimension(:,:), allocatable :: vr, vl, a
		integer :: n, lda, ldvl, ldvr, lwork, info, i
		
	
		a = w
		jobvl = 'N'
		jobvr = 'N'
		n = size(a,1)
		lda = n
		ldvl = n
		ldvr = n
		lwork = 10000
		
		allocate(wr(n), wi(n), vl(ldvl, n), &
				vr(ldvr, n), work(lwork))
			
		call dgeev(jobvl, jobvr, n, a, lda, wr, wi, &
						vl, ldvl, vr, ldvr, work, lwork, info)		
		do i = 1, n
			magnitudes = sqrt(wr**2 + wi**2)
		end do	
					
		w = w/maxval(magnitudes)
		
		end subroutine normalize_spectral_radius
		
		subroutine generate_input_matrix(v)
			real (8), dimension(:,:), intent (inout) :: v
			real (8) :: p, weight, num
			integer :: i, j
			
			v = 0
			
			call random_number(num)
			v(1,int(num*N)) = norm_rand(0.5, 1.0)/sqrt(2*pi)	
			call random_number(num)
			v(2,int(num*N)) = norm_rand(0.5, 1.0)/sqrt(2*pi)	
			call random_number(num)
			v(3,int(num*N)) = norm_rand(0.5, 1.0)/sqrt(2*pi)		
			
			
			do i = 1, 3
				do j = 1, N
					call random_number(p)
					if (p.lt.0.05) then
						do while (abs(weight).gt.1)
							weight =  norm_rand(0.5, 1.0)/sqrt(2*pi)	
						end do									
						v(i,j) = abs(weight)
						weight = 2			
					end if
				end do
			end do	
	end subroutine generate_input_matrix
	
	! uses the marsaglia polar method to create pseudo-random number pairs that
	! have a normal distribution. The function saves one of the random numbers from
	! the generated pair as a spare to be returned for the next call to the function.
	! Returns a real scalar
	function norm_rand(mean, std_dev)
		real :: norm_rand
		real, intent(in) :: mean, std_dev
		real :: x, y, r
		real, save :: spare
		logical, save :: has_spare
		! use a spare saved from a previous run if one exists
		if (has_spare) then
			has_spare = .FALSE.
			norm_rand = mean + (std_dev * spare)
			return
		else
			r = 1.0
			do while ( r >= 1.0 )
				! generate random number pair between 0 and 1
				call random_number(x)
				call random_number(y)
				! normalise random numbers to be in square of side-length = R
				x = (x * 2.0) - 1.0
				y = (y * 2.0) - 1.0
				r = x*x + y*y
			end do

			! calculate the co-efficient to multiply random numbers x and y
			! by to achieve normal distribution
			r = sqrt((-2.0 * log(r)) / r)

			norm_rand = mean + (std_dev * x * r)
			spare = y * r
			has_spare = .TRUE.
			return
		end if
	end function norm_rand	
	
		subroutine write_matrix(Matrix)
		real (8), dimension(:,:), intent(inout) :: Matrix
		integer :: i
		do i=1, size(Matrix, 1)
			write(*,'(20G12.4)') Matrix(i,:) 
		end do	
	end subroutine

end program
