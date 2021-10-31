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

        do i = 0, num_of_coef - 1
            do k = 0, num_of_coef - 1
                matrix_A(i + 1, k + 1) = matrix_A_i_k(i, k, kernel_K, alpha)
            end do
        end do

        do i = 0, num_of_coef - 1
            vector_b(i + 1) = vector_b_i(i, fun_f)
        end do

        ! do i = 1, num_of_coef
        !     write(*, *) matrix_A(i, :), vector_b(i)
        ! end do
        
        write(*, *) 'Det matrix: ', FindDet(matrix_A, num_of_coef)

        call SSE_Mod_Gauss(matrix_A, vector_b, c)

        write(*, *) 'discrepancy :', norm2(matmul(matrix_A, c) - vector_b)

    end function coefficients_of_the_series

    function matrix_A_i_k(i, k, kernel_K, alpha)
        real(mp) :: matrix_A_i_k, alpha
        integer :: i, k
        interface
            function kernel_K(x, t)
                use :: init_data
            implicit none
                    real(mp) :: x, t, kernel_K
            end function kernel_K

        end interface

        matrix_A_i_k = scalar_dot(dot_A_phi, phi_k)

        contains 

        function dot_A_phi(t)
            real(mp) :: dot_A_phi, t

            dot_A_phi = operator_A(t, phi_i, kernel_K, alpha)

        end function dot_A_phi

        function phi_k(x)
            real(mp) phi_k, x

            phi_k = legendre_polynom_rec(x, k)

        end function phi_k

        function phi_i(x)
            real(mp) phi_i, x

            phi_i = legendre_polynom_rec(x, i)

        end function phi_i

    end function matrix_A_i_k

    function vector_b_i(i, fun_f)
        real(mp) :: vector_b_i
        integer :: i
        interface
            function fun_f(x)
                use :: init_data
            implicit none
                    real(mp) :: x, fun_f
            end function fun_f

        end interface

        vector_b_i = scalar_dot(fun_f, phi_i)

        contains

            function phi_i(t)
                real(mp) :: t, phi_i

                phi_i = legendre_polynom_rec(t, i)

            end function phi_i

    end function vector_b_i

    function operator_A(x, u, kernel_K, alpha)
        real(mp) :: x, alpha, operator_A
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

        operator_A = alpha * u(x) + gauss_quad_integral(a, b, 15, dot_K_u)

        contains

            function dot_K_u(t)
                real(mp) :: t, dot_K_u

                dot_K_u = kernel_K(x, t) * u(t)

            end function dot_K_u

    end function operator_A

    function result_fun(x, c)
        real(mp) :: x, result_fun
        real(mp), dimension(1:) :: c
        integer :: k

        result_fun = 0.0_mp

        do k = 1, num_of_coef
            result_fun = result_fun + c(k) * legendre_polynom_rec(x, k - 1)
        end do

    end function result_fun

    function scalar_dot(phi, psi)
        real(mp) :: scalar_dot
        interface
            function phi(x)
                use :: init_data
            implicit none
                    real(mp) :: x, phi
            end function phi

            function psi(x)
                use :: init_data
            implicit none
                    real(mp) :: x, psi
            end function psi

        end interface

        scalar_dot = gauss_quad_integral(a, b, 15, dot_phi_psi)

        contains

        function dot_phi_psi(t)
            real(mp) :: t, dot_phi_psi

            dot_phi_psi = phi(t) * psi(t)

        end function dot_phi_psi

    end function scalar_dot

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

        error_of_results_method = operator_A(x, u, kernal_integral_equation, alpha) - f(x)

        contains

        function u(t)
            real(mp) :: t, u

            u = result_fun(t, result_c)

        end function u

    end function error_of_results_method

end module SIEI
