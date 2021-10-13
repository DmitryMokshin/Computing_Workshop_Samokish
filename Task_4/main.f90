program main
    use :: init_data
    use :: SIEI
implicit none

    integer :: i, n
    real(mp), dimension(1:num_of_coef) :: c
    real(mp) :: h

    c = coefficients_of_the_series(kernel_T, f_1)
    
    n = 100
    h = (b - a) / n

    open(15, file='result.dat', status='replace')

    do i = 0, n
        write(15, *) a + h * i, result_fun(a + h * i, c)
    end do

    close(15)

    open(17, file='result_coef.dat', status='replace')

    do i = 1, num_of_coef
        write(17, *) c(i)
    end do

    close(17)

    write(*, *) 'alpha = ', alpha
    write(*, *) error_method_collocations(c, kernal_integral_equation)

end program main