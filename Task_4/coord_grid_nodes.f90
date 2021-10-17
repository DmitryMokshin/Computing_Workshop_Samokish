module coord_grid_nodes
    use :: init_data
    implicit none
contains

    function solution_cheb() result(vector_roots_cheb)
        real(mp), dimension(1:num_of_coef) :: vector_roots_cheb
        integer :: k

        do k = 1, num_of_coef
            vector_roots_cheb(k) = 0.5_mp * (a + b) + 0.5_mp * (b - a) * cos((2.0_mp * k - 1) / 2.0_mp / num_of_coef * pi)
        end do
        
    end function solution_cheb

    function legendre_polynomial_coefficients() result(result_coefficient)
        real(mp), dimension(1:num_of_coef, 0:num_of_coef) :: result_coefficient
        integer :: i
        real(mp), dimension(0:num_of_coef, 0:num_of_coef) :: P_n_matrix
        
        P_n_matrix = 0.0_mp
        P_n_matrix(0, num_of_coef) = 1.0_mp; P_n_matrix(1, num_of_coef - 1) = 1.0_mp

        do i = 2, num_of_coef
            P_n_matrix(i, 0:num_of_coef)=(2.0_mp * i - 1.0_mp) / i * &
            & eoshift(P_n_matrix(i - 1,0:num_of_coef), 1, 0.0_mp, 1) - (i - 1.0_mp) / i * P_n_matrix(i - 2,0:num_of_coef)
        end do

        result_coefficient=P_n_matrix(0:num_of_coef - 1, 0:num_of_coef)

    end function legendre_polynomial_coefficients

    function legendre_polynom(y, leg_pol_coef) result(value_polynom)
        real(mp) :: y, value_polynom, x
        real(mp), dimension(0:num_of_coef) :: leg_pol_coef, vector_x
        integer :: i

        x = 2.0_mp / (b - a) * y - (b + a) / (b - a)

        vector_x = (/(x ** i, i = num_of_coef, 0, -1)/)
        value_polynom = dot_product(leg_pol_coef, vector_x)

    end function legendre_polynom

    function legendre_polynom_rec(y, k)
        real(mp) :: x, legendre_polynom_rec, y
        real(mp), dimension(-1:1) :: P
        integer :: k, i

        x = 2.0_mp / (b - a) * y - (b + a) / (b - a)

        P(-1) = 1.0_mp
        P(0) = x

        if (k == 0) then
            legendre_polynom_rec = P(-1)
        end if

        if (k == 1) then
            legendre_polynom_rec = P(0)
        end if

        if (k > 1) then
            i = 1
            do 
                P(1) = (2.0_mp * i + 1.0_mp) / (i + 1.0_mp) * x * P(0) - i / (i + 1.0_mp) * P(-1)
                if (i + 1 == k) then
                    legendre_polynom_rec = P(1)
                    exit
                end if
                P(-1) = P(0)
                P(0) = P(1)
                i = i + 1
            end do
        end if

    end function legendre_polynom_rec

end module coord_grid_nodes