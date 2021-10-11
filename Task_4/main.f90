program main
    use :: init_data
    use :: coord_grid_nodes
    use :: SIEI
implicit none

    real(mp), dimension(1:num_of_coef, 0:num_of_coef) :: P_matrix
    integer :: i, n
    real(mp), dimension(1:num_of_coef) :: c
    real(mp) :: h

    c = coefficients_of_the_series(kernal_integral_equation, f)
    n = 100
    h = (b - a) / n

    open(15, file='resut.dat', status='replace')

    do i = 0, 100
        write(15, *) a + h * i, result_fun(a + h * i, c)
    end do

    close(15)

end program main