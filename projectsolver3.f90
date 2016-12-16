program solvermain
    implicit none
    
    interface temp
        subroutine temp(dr, dz, dt, k1, Rho, Cp, heat_gen_peak, Nr, Nz, n)
            double precision, intent(in) :: dr, dz(4), dt, k1(4), Rho(4), Cp(4), heat_gen_peak
			integer, intent(in) :: n, Nr(4), Nz(4)
        end subroutine temp
    end interface
    
 	integer :: n, Nr(4), Nz(4)
	double precision :: dr, dz(4), dt, k1(4), Rho(4), Cp(4), heat_gen_peak
	

	!test conditions to run from command line
	dr = 0.0001
	dz(1) = 0.000075
	dz(2) = 0.00000025
	dz(3) = 0.000000025
	dz(4) = 0.00000015
	dt = 0.02000240009600384
	k1(1) = 35.0
	k1(2) = 73.0
	k1(3) = 21.0
	k1(4) = 5.3
	Rho(1) = 4000.0
	Rho(2) = 2140.0
	Rho(3) = 4500.0
	Rho(4) = 11000.0
	Cp(1) = 880.0
	Cp(2) = 125.6
	Cp(3) = 544.0
	Cp(4) = 250.0
	heat_gen_peak = 30000000000.0
	Nr(1) = 47
	Nr(2) = 20
	Nr(3) = 20
	Nr(4) = 20
	Nz(1) = 20
	Nz(2) = 20
	Nz(3) = 20
	Nz(4) = 20
	n = 500

	call temp(dr, dz, dt, k1, Rho, Cp, heat_gen_peak, Nr, Nz, n)
end program solvermain


subroutine temp(dr, dz, dt, k1, Rho, Cp, heat_gen_peak, Nr, Nz, n)
	implicit none

	integer, intent(in) :: n, Nr(4), Nz(4)
	double precision, intent(in) :: dr, dz(4), dt, k1(4), Rho(4), Cp(4), heat_gen_peak
	double precision :: T_new(Nr(1) * Nz(1) + Nr(2) * Nz(2) + Nr(3) * Nz(3) + Nr(4) * Nz(4)), &
						T(Nr(1) * Nz(1) + Nr(2) * Nz(2) + Nr(3) * Nz(3) + Nr(4) * Nz(4)), &
						T_ic(Nr(1) * Nz(1) + Nr(2) * Nz(2) + Nr(3) * Nz(3) + Nr(4) * Nz(4)), &
						b(Nr(1) * Nz(1) + Nr(2) * Nz(2) + Nr(3) * Nz(3) + Nr(4) * Nz(4))
						

	integer :: i, j, m, p, c, c_inner, c_outer, c_upper, c_lower
	double precision :: r_i, alpha(4), beta(4), gamma1(4), mu(4), tau, omega, max_tollerance, current_tollerance, q
	

	open (unit = 1, file = "temp.txt")
	read (1, *) (T_ic(i), i=1, (Nr(1) * Nz(1) + Nr(2) * Nz(2) + Nr(3) * Nz(3) + Nr(4) * Nz(4)))
	close (1)

	do i = 1, 4
		alpha(i) = (k1(i) * dt) / (2*Rho(i) * Cp(i)*dr)
		beta(i) = (k1(i) * dt) / (Rho(i) * Cp(i) * (dr**2))
		gamma1(i) = (k1(i) * dt) / (Rho(i) * Cp(i) * dz(i)**2)
		mu(i) = 1.0 + (2 * beta(i)) + (2 * gamma1(i))
	end do	


	!inital first guess set to zero
	do c = 1, (Nr(1)*Nz(1)+Nr(2)*Nz(2)+Nr(3)*Nz(3)+Nr(4)*Nz(4))
		T(c) = 0.0
	end do
	
	
	!Jacobi Iteration
	
	!starting loop for each time step, n = # of time step
	do p = 1, n

		!for p > 1, the initial guess of T(c) will be already set to the solution of the previous time step

	
		!defining right side vector b(c) - adding in effects of heat generation		
		!layer 1 elements		
		do j = 1, Nz(1)
			do i = 1, Nr(1)
				!setting cardinal index			
				c = ((j-1)*Nr(1)) + i

				!no heat generation
				b(c) = T_ic(c)
			end do
		end do
		
		!layer 2 elements
		do j = 1, Nz(2)
			q = heat_gen_peak*( (j*dz(2)) / (dz(2)*Nz(2)+dz(3)*Nz(3)+dz(4)*Nz(4)) )		
			do i = 1, Nr(2)
				!setting cardinal index
				c = ((j-1)*Nr(2)) + i + (Nr(1)*Nz(1))
				
				!adding term for heat generation to right side vector for layer 2 elements
				b(c) = T_ic(c) + ((q*dt) / (Rho(2)*Cp(2)))

			end do
		end do

		!layer 3 elements
		do j = 1, Nz(3)
			q = heat_gen_peak*(((j*dz(3)) + (dz(2)*Nz(2))) / (dz(2)*Nz(2)+dz(3)*Nz(3)+dz(4)*Nz(4)))
			do i = 1, Nr(3)
				!setting cardinal index
				c = ((j-1)*Nr(3)) + i + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2))
				
				!adding term for heat generation to right side vector for layer 4 elements
				b(c) = T_ic(c) + ((q*dt) / (Rho(3)*Cp(3)))
			end do
		end do

		!layer 4 elements
		do j = 1, Nz(4)
			q = heat_gen_peak*(( (j*dz(4)) + (dz(2)*Nz(2)) + (dz(3)*Nz(3)) ) / ( dz(2)*Nz(2)+dz(3)*Nz(3)+dz(4)*Nz(4) ))
			do i = 1, Nr(4)
				!setting cardinal index
				c = ((j-1)*Nr(4)) + i + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2)) + (Nr(3)*Nz(3))
				
				!adding term for heat generation to right side vector for layer 4 elements
				b(c) = T_ic(c) + ((q*dt) / (Rho(4)*Cp(4)))
			end do
		end do

		!setting max tollerance to non-zero value
		max_tollerance = 1.0
		!initializing counter for number of iterations		
		m = 0

		!while loop for Jacobi Iteration
		do while (max_tollerance .gt. 0.000001)
			
			!incrementing iteration counter		
			m = m + 1
		
			!layer 1 iteration step	
			
			!setting T_new for layer 1 for r not equal to zero
			do j = 1, Nz(1)
				do i = 2, Nr(1)
					!setting cardinal index			
					c = ((j-1)*Nr(1)) + i

					c_inner = c - 1
			
					if (i .ne. Nr(1)) then
						c_outer = c + 1
					else
						!insulated boundary						
						c_outer = c -1
					end if

					if (j .ne. 1) then			
						c_lower = c - Nr(1)
					else
						!insulated boundary
						c_lower = c + Nr(1)
					end if
	
					if (j .ne. Nz(1)) then			
						c_upper = c + Nr(1)
					else if ((j .eq. Nz(1)) .and. (i .gt. Nr(2))) then
						!insulated boundary						
						c_upper = c - Nr(1) 
					else
						!borders layer 2					
						c_upper = (Nr(1)*Nz(1)) + i
					end if

					!setting radius dependent variables
					r_i = dr*(i-1)
					tau = (alpha(1) / r_i) + beta(1)
					omega = (alpha(1) / r_i) - beta(1)

					T_new(c) = (1 / mu(1)) * (b(c) - (omega*T(c_inner) - tau*T(c_outer) - gamma1(1)*T(c_upper) - gamma1(1)*T(c_lower)))
					
				end do
			end do

			!setting T_new for layer 1 for r = 0
			do j = 1, Nz(1)
				c = ((j-1)*Nr(1)) + 1
		
				c_outer = c + 1
			
				if (j .ne. 1) then			
					c_lower = c - Nr(1)
				else
					c_lower = c + Nr(1)
				end if
	
				if (j .ne. Nz(1)) then			
					c_upper = c + Nr(1)
				else
					!j at Nz(1) - top of layer	
					c_upper = (Nr(1)*Nz(1)) + 1				
				end if

				T_new(c) = (1 / mu(1)) * (b(c) - (-2*beta(1)*T(c_outer) - gamma1(1)*T(c_upper) - gamma1(1)*T(c_lower)))
			end do

			!layer 2 iteration step	
			
			!setting T_new for layer 2 for r not equal to zero
			do j = 1, Nz(2)
				do i = 2, Nr(2)
					!setting cardinal index			
					c = ((j-1)*Nr(2)) + i + (Nr(1)*Nz(1))

					c_inner = c - 1
			
					if (i .ne. Nr(2)) then
						c_outer = c + 1
					else
						!insulated boundary						
						c_outer = c -1
					end if

					if (j .ne. 1) then			
						c_lower = c - Nr(2)
					else
						!j equal 1 - bottom of layer
						c_lower = c - Nr(1)
					end if
	
					if (j .ne. Nz(2)) then			
						c_upper = c + Nr(2)
					else
						!top of layer 2 - same formula
						c_upper = c + Nr(2)						
					end if

					!setting radius dependent variables
					r_i = dr*(i-1)
					tau = (alpha(2) / r_i) + beta(2)
					omega = (alpha(2) / r_i) - beta(2)

					T_new(c) = (1 / mu(2)) * (b(c) - (omega*T(c_inner) - tau*T(c_outer) - gamma1(2)*T(c_upper) - gamma1(2)*T(c_lower)))
				end do
			end do

			
			!setting T_new for layer 2 for r = 0
			do j = 1, Nz(2)
				c = ((j-1)*Nr(2)) + 1 + (Nr(1)*Nz(1))
		
				c_outer = c + 1
			
				if (j .ne. 1) then			
					c_lower = c - Nr(2)
				else
					c_lower = c - Nr(1)
				end if
	
				if (j .ne. Nz(2)) then			
					c_upper = c + Nr(2)
				else
					!top of layer 2 - same formula
					c_upper = c + Nr(2)				
				end if

				T_new(c) = (1 / mu(2)) * (b(c) - (-2*beta(2)*T(c_outer) - gamma1(2)*T(c_upper) - gamma1(2)*T(c_lower)))
			end do

			!layer 3 iteration step	
			
			!setting T_new for layer 3 for r not equal to zero
			do j = 1, Nz(3)
				do i = 2, Nr(3)
					!setting cardinal index			
					c = ((j-1)*Nr(3)) + i + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2))

					c_inner = c - 1
			
					if (i .ne. Nr(3)) then
						c_outer = c + 1
					else
						!insulated boundary						
						c_outer = c -1
					end if

					if (j .ne. 1) then			
						c_lower = c - Nr(3)
					else
						!j equal 1 - bottom of layer
						c_lower = c - Nr(2)
					end if
	
					if (j .ne. Nz(3)) then			
						c_upper = c + Nr(3)
					else
						!top of layer 3 - same formula
						c_upper = c + Nr(3)						
					end if

					!setting radius dependent variables
					r_i = dr*(i-1)
					tau = (alpha(3) / r_i) + beta(3)
					omega = (alpha(3) / r_i) - beta(3)

					T_new(c) = (1 / mu(3)) * (b(c) - (omega*T(c_inner) - tau*T(c_outer) - gamma1(3)*T(c_upper) - gamma1(3)*T(c_lower)))
				end do
			end do

			!setting T_new for layer 3 for r = 0
			do j = 1, Nz(3)
				c = ((j-1)*Nr(3)) + 1 + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2))
		
				c_outer = c + 1
			
				if (j .ne. 1) then			
					c_lower = c - Nr(3)
				else
					c_lower = c - Nr(2)
				end if
	
				if (j .ne. Nz(3)) then			
					c_upper = c + Nr(3)
				else
					!top of layer 3 - same formula
					c_upper = c + Nr(3)				
				end if

				T_new(c) = (1 / mu(3)) * (b(c) - (-2*beta(3)*T(c_outer) - gamma1(3)*T(c_upper) - gamma1(3)*T(c_lower)))
			end do

			!layer 4 iteration step	
			
			!setting T_new for layer 4 for r not equal to zero
			do j = 1, Nz(4)
				do i = 2, Nr(4)
					!setting cardinal index			
					c = ((j-1)*Nr(4)) + i + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2)) + (Nr(3)*Nz(3))

					c_inner = c - 1
			
					if (i .ne. Nr(4)) then
						c_outer = c + 1
					else
						!insulated boundary						
						c_outer = c -1
					end if

					if (j .ne. 1) then			
						c_lower = c - Nr(4)
					else
						!j equal 1 - bottom of layer
						c_lower = c - Nr(3)
					end if
	
					if (j .ne. Nz(4)) then			
						c_upper = c + Nr(4)
					else
						!top of layer 4 - insulated
						c_upper = c - Nr(4)						
					end if

					!setting radius dependent variables
					r_i = dr*(i-1)
					tau = (alpha(4) / r_i) + beta(4)
					omega = (alpha(4) / r_i) - beta(4)

					T_new(c) = (1 / mu(4)) * (b(c) - (omega*T(c_inner) - tau*T(c_outer) - gamma1(4)*T(c_upper) - gamma1(4)*T(c_lower)))
				end do
			end do

			!setting T_new for layer 4 for r = 0
			do j = 1, Nz(4)
				c = ((j-1)*Nr(4)) + 1 + (Nr(1)*Nz(1)) + (Nr(2)*Nz(2)) + (Nr(3)*Nz(3))
		
				c_outer = c + 1
			
				if (j .ne. 1) then			
					c_lower = c - Nr(4)
				else
					c_lower = c - Nr(3)
				end if
	
				if (j .ne. Nz(4)) then			
					c_upper = c + Nr(4)
				else
					!top of layer 4 - insulated
					c_upper = c - Nr(4)				
				end if

				T_new(c) = (1 / mu(4)) * (b(c) - (-2*beta(4)*T(c_outer) - gamma1(4)*T(c_upper) - gamma1(4)*T(c_lower)))
			end do
		
			!updating the max tollerance value every 200 iterations	
			if (mod(m,200) .eq. 0) then
				max_tollerance = abs(T_new(1) - T(1))
				
				do c = 2, (Nr(1)*Nz(1)+Nr(2)*Nz(2)+Nr(3)*Nz(3)+Nr(4)*Nz(4))
					current_tollerance = abs(T_new(c) - T(c))
					if (current_tollerance .gt. max_tollerance) then
						max_tollerance = current_tollerance
					end if
				end do
			end if
		
			!copying for next iteration
			do c = 1, (Nr(1)*Nz(1)+Nr(2)*Nz(2)+Nr(3)*Nz(3)+Nr(4)*Nz(4))
				T(c) = T_new(c)
			end do
	
		end do
		!ended while loop - solution for time step determined


		!copying solved temperature to the initial temperature for next time step
		do c = 1, (Nr(1)*Nz(1)+Nr(2)*Nz(2)+Nr(3)*Nz(3)+Nr(3)*Nz(4))
			T_ic(c) = T(c)
		end do
	
	end do
	!completed n # of time steps; T = solution at time dt * n

	

	!output temp_out to text file

	open (unit = 1, file = "temp.txt", action = "write", status = "replace")
	do c = 1, (Nr(1)*Nz(1)+Nr(2)*Nz(2)+Nr(3)*Nz(3)+Nr(4)*Nz(4))
		write (1,*) (T(c))
	end do
	close (1)

end subroutine temp
