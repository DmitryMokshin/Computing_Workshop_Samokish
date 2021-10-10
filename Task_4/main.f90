program main
    use :: init_data
    use :: coord_grid_nodes
implicit none

    real(mp), dimension(1:num_of_coef, 0:num_of_coef) :: P_matrix
    integer :: i

    P_matrix = legendre_polynomial_coefficients()

    do i = 1, num_of_coef
        write(*, *) i, sum(P_matrix(i, :)), legendre_polynom(0.5_mp, P_matrix(i, :))
    end do

end program main