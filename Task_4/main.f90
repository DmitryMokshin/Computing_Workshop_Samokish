program main
    use :: init_data
    use :: SIEI
implicit none
	integer, parameter :: num_of_alphas = 8
    integer :: i, n, k
    real(mp), dimension(0:8, 1:num_of_coef) :: c
    real(mp) :: h, x
    real(mp), dimension(0:num_of_alphas) :: alphas
    
    alphas = (/( 0.01_mp ** i, i = 0, num_of_alphas )/)

	do i = 0, 8
		c(i, :) = coefficients_of_the_series(kernel_T, f_1, alphas(i))
    end do 

    n = 20
    h = (b - a) / n

    open(15, file='result.dat', status='replace')
    
    write(15, '(20x, a, $)') '|'
    do i = 0, num_of_alphas
		write(15, '(a, e9.2, 5x, a, $)') 'alpha = ', alphas(i), '|'
	end do
	
	write(15, *) 

    do i = 0, n
		x = a + h * i
		write(15, '(e20.12, a, $)') x, '|'  
		do k = 0, num_of_alphas
			write(15, '(e20.12, 2x, a, $)') result_fun(x, c(k, :)), '|'
		end do
		write(15,*)
    end do

    close(15)

    open(17, file='result_coef.dat', status='replace')

    do i = 0, num_of_alphas
		write(17, '(a, e9.2, 5x, a, $)') 'alpha = ', alphas(i), '|'
	end do
	
	write(17, *) 
    
	do k = 1, num_of_coef
		do i = 0, num_of_alphas
			write(17, '(e20.12, 2x, a, $)') c(i, k), '|'
		end do
		write(17, *)
    end do

    close(17)

    !write(*, *) 'alpha = ', alphas
    !write(*, *) error_method_collocations(c, kernal_integral_equation)

end program main
