module SIEI
    use :: init_data
    use :: gauss_quad
    use :: coord_grid_nodes
    use :: SSE
    implicit none
contains

    function kernel_T(x, t)
        real(mp) :: kernel_T, x, t

        kernel_T = gauss_quad_integral(a, b, 27, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernal_integral_equation(ksi, x) * kernal_integral_equation(ksi, t)

            end function dot_fun

    end function kernel_T

    function f_1(x)
        real(mp) :: f_1, x, t

        f_1 = gauss_quad_integral(a, b, 27, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernal_integral_equation(ksi, x) * f(ksi)

            end function dot_fun

    end function f_1

    function coefficients_of_the_series(kernel_K, fun_f) result(c)
        real(mp), dimension(1:num_of_coef) :: c
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K

            function fun_f(x)
                use :: init_data
            implicit none
                    real(mp) :: x, fun_f
            end function fun_f

        end interface
        integer :: k, i
        real(mp), dimension(1:num_of_coef, 0:num_of_coef) :: pol_matrix
        real(mp), dimension(1:num_of_coef) :: roots_cheb, vector_b
        real(mp), dimension(1:num_of_coef, 1:num_of_coef) :: matrix_A

        roots_cheb = solution_cheb()

        pol_matrix = legendre_polynomial_coefficients()

        do i = 1, num_of_coef
            vector_b(i) = fun_f(roots_cheb(i))
        end do

        do i = 1, num_of_coef
            do k = 1, num_of_coef
                matrix_A(i, k) = alpha * legendre_polynom(roots_cheb(i), pol_matrix(k, :)) &
                & - integral(pol_matrix(k, :), roots_cheb(i), kernel_K)
            end do
        end do

        call SSE_Mod_Gauss(matrix_A, vector_b, c)

    end function coefficients_of_the_series

    function integral(pol_coef, root_cheb, kernel_K)
        real(mp), dimension(0:) :: pol_coef
        real(mp) :: integral, root_cheb
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K
        end interface

        integral = gauss_quad_integral(a, b, 27, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernel_K(root_cheb, ksi) * legendre_polynom(ksi, pol_coef)

            end function dot_fun

    end function integral

    function result_fun(x, c)
        real(mp) :: x, result_fun
        real(mp), dimension(1:num_of_coef, 0:num_of_coef) :: pol_matrix
        real(mp), dimension(1:) :: c
        integer :: k

        pol_matrix = legendre_polynomial_coefficients()

        result_fun = 0.0_mp

        do k = 1, num_of_coef
            result_fun = result_fun + c(k) * legendre_polynom(x, pol_matrix(k, :))
        end do

    end function result_fun

    function scalar_product_l2(psi, phi)
        real(mp) :: scalar_product_l2
        interface
            function psi(t)
                use :: init_data
            implicit none
                    real(mp) :: t, psi
            end function psi

            function phi(t)
                use :: init_data
            implicit none
                    real(mp) :: t, phi
            end function phi

        end interface

        scalar_product_l2 = gauss_quad_integral(a, b, 27, dot_fun)

        contains

        function dot_fun(ksi)
            real(mp) :: ksi, dot_fun

            dot_fun = psi(ksi) * phi(ksi) 

        end function dot_fun

    end function scalar_product_l2

    function error_method_collocations(c, kernel_K)
        real(mp), dimension(1:) :: c
        real(mp) :: delta, eps, error_method_collocations, h
        integer :: i, n
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K

        end interface

        eps = sqrt(scalar_product_l2(int_for_dot, int_for_dot))
        delta = sqrt(scalar_product_l2(difference_f, difference_f))

        write(*, *) 'eps = ', eps
        write(*, *) 'delta = ', delta

        error_method_collocations = delta ** 2 / eps

        open(16, file='error_result.dat', status='replace')

        n = 100
        h = (b - a) / n

        do i = 0, n
            write(16, *) a + h * i, int_for_dot(a + h * i)
        end do

        close(16)

        contains

        function u(t)
            real(mp) :: u, t

            u = result_fun(t, c)

        end function u

        function int_for_dot(t)
            real(mp) :: t, int_for_dot

            int_for_dot = integral_Ku(t, kernel_K, u) - f(t)

        end function int_for_dot

        function difference_f(t)
            real(mp) :: difference_f, t

            difference_f = f(t) - f_1(t)

        end function difference_f

    end function error_method_collocations

    function integral_Ku(x, kernel_K, u)
        real(mp) :: integral_Ku, x
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K

            function u(x)
                use :: init_data
            implicit none
                    real(mp) :: x, u
            end function u

        end interface

        integral_Ku = gauss_quad_integral(a, b, 27, dot_fun)

        contains

        function dot_fun(ksi)
            real(mp) :: ksi, dot_fun

            dot_fun = kernel_K(x, ksi) * u(ksi)

        end function dot_fun

    end function integral_Ku

end module SIEI