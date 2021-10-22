program main
    use :: init_data
    use :: SIEI
implicit none
    integer, parameter :: num_of_alphas = 7
    integer :: i, n, k
    real(mp), dimension(0:num_of_alphas, 1:num_of_coef) :: c
    real(mp) :: h, x, init_alpha, step_alpha
    real(mp), dimension(0:num_of_alphas) :: alphas
    character(2) :: num_of_coef_char

    init_alpha = 10.0_mp ** (-0.0_mp)
    step_alpha = 0.01_mp
    alphas = (/( init_alpha * step_alpha ** i, i = 0, num_of_alphas )/)

    num_of_coef_char = char(ichar('0')+num_of_coef/10)//char(ichar('0')+mod(num_of_coef,10))

    do i = 0, num_of_alphas
        c(i, :) = coefficients_of_the_series(kernel_T, f_1, alphas(i))
    end do 

    n = 20
    h = (b - a) / n

    open(15, file='result'//num_of_coef_char//'.csv', status='replace')

    write(15, '(a, $)') 'x'
    do i = 0, num_of_alphas
        write(15, '(a, a, e9.2, $)') '|', 'y, \alpha = ', alphas(i)
    end do

    write(15, *) 

    do i = 0, n
        x = a + h * i
        write(15, '(e20.12, $)') x
        do k = 0, num_of_alphas
            write(15, '(a, e20.12, $)') '|', result_fun(x, c(k, :))
        end do
        write(15,*)
    end do

    close(15)

    open(17, file='result_coef'//num_of_coef_char//'.csv', status='replace')

    write(17, '(a, e9.2, $)') 'c_i, \alpha = ', alphas(0)
    do i = 1, num_of_alphas
        write(17, '(a, e9.2, $)') '| c_i, \alpha = ', alphas(i)
    end do

    write(17, *) 

    do k = 1, num_of_coef
        write(17, '(e20.12,  $)') c(0, k)
        do i = 1, num_of_alphas
            write(17, '(a, e20.12,  $)') '|', c(i, k)
        end do
        write(17, *)
    end do

    close(17)

    do i = 0, n
        x = a + h * i
        write(*, '(e20.12, $)') x
        do k = 0, num_of_alphas
            write(*, '(e20.12, $)') error_of_results(x, c(k, :))
        end do
        write(*,*)
    end do

    n = 100
    h = (b - a) / n

    open(20, file='resultpolynomlegendre'//num_of_coef_char//'.csv', status='replace')

    write(20, '(a, $)') 'x'
    do i = 0, num_of_coef - 1
        write(20, '(a, a, i2, $)') '|', 'P, N = ', i
    end do

    write(20, *) 

    do i = 0, n
        x = a + h * i
        write(20, '(e20.12, $)') x
        do k = 0, num_of_coef - 1
            write(20, '(a, e20.12, $)') '|', legendre_polynom_rec(x, k)
        end do
        write(20,*)
    end do

end program main
