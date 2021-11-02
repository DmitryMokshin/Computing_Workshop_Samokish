module SIEI
    use :: init_data
    use :: gauss_quad
    use :: coord_grid_nodes
    use :: SSE
    implicit none
contains

    function kernel_T(x, t)
        real(mp) :: kernel_T, x, t

        kernel_T = gauss_quad_integral(a, b, 15, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernal_integral_equation(ksi, x) * kernal_integral_equation(ksi, t)

            end function dot_fun

    end function kernel_T

    function f_1(x)
        real(mp) :: f_1, x, t

        f_1 = gauss_quad_integral(a, b, 15, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernal_integral_equation(ksi, x) * f(ksi)

            end function dot_fun

    end function f_1

    function coefficients_of_the_series(kernel_K, fun_f, alpha) result(c)
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
        real(mp) :: alpha
        real(mp), dimension(1:num_of_coef) :: roots_cheb, vector_b
        real(mp), dimension(1:num_of_coef, 1:num_of_coef) :: matrix_A

        roots_cheb = solution_cheb()

        do i = 1, num_of_coef
            vector_b(i) = fun_f(roots_cheb(i))
        end do

        do i = 1, num_of_coef
            do k = 1, num_of_coef
                matrix_A(i, k) = alpha * legendre_polynom_rec(roots_cheb(i), k - 1) &
                & + integral(roots_cheb(i), k, kernel_K)
            end do
        end do

        call SSE_Mod_Gauss(matrix_A, vector_b, c)
        
    end function coefficients_of_the_series

    function integral(root_cheb, k, kernel_K)
        real(mp) :: integral, root_cheb
        integer :: k
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K
        end interface

        integral = gauss_quad_integral(a, b, 15, dot_fun)

        contains

            function dot_fun(ksi)
                real(mp) :: ksi, dot_fun

                dot_fun = kernel_K(root_cheb, ksi) * legendre_polynom_rec(ksi, k - 1)

            end function dot_fun

    end function integral

    function result_fun(x, c)
        real(mp) :: x, result_fun
        real(mp), dimension(1:) :: c
        integer :: k

        result_fun = 0.0_mp

        do k = 1, num_of_coef
            result_fun = result_fun + c(k) * legendre_polynom_rec(x, k - 1)
        end do

    end function result_fun

    function operator_A(x, u, alpha)
        real(mp) :: x, alpha, operator_A
        interface
            function u(x)
                use :: init_data
            implicit none
                    real(mp) :: x, u
            end function u
        end interface

        operator_A = alpha * u(x) + gauss_quad_integral(a, b, 15, dot_T_u)

        contains

        function dot_T_u(ksi)
            real(mp) :: ksi, dot_T_u

                dot_T_u = kernel_T(x, ksi) * u(ksi)

        end function dot_T_u

    end function operator_A

    function error_of_results(x, result_c)
        real(mp), dimension(:) :: result_c
        real(mp) :: error_of_results, x

        error_of_results = gauss_quad_integral(a, b, 15, dot_fun) - f(x)

        contains

        function u(t)
            real(mp) :: t, u

            u = result_fun(t, result_c)

        end function u

        function dot_fun(ksi)
            real(mp) :: ksi, dot_fun

            dot_fun = kernal_integral_equation(x, ksi) * u(ksi)

        end function dot_fun


    end function error_of_results

    function error_of_results_method(x, result_c, alpha)
        real(mp), dimension(:) :: result_c
        real(mp) :: error_of_results_method, x, alpha

        error_of_results_method = operator_A(x, u, alpha) - f_1(x)

        contains

        function u(t)
            real(mp) :: t, u

            u = result_fun(t, result_c)

        end function u

    end function error_of_results_method

    function error_of_results_true(x)
        real(mp) :: error_of_results_true, x

        error_of_results_true = gauss_quad_integral(a, b, 15, dot_fun) - f(x)

        contains

        function dot_fun(ksi)
            real(mp) :: ksi, dot_fun

            dot_fun = kernal_integral_equation(x, ksi) * u_result(ksi)

        end function dot_fun


    end function error_of_results_true

end module SIEI
