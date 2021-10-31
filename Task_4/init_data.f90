module init_data
    implicit none
    integer, parameter :: mp = 8, num_of_coef = 12
    real(mp), parameter :: pi = acos(-1.0_mp), a = 0.0, b = 1.0
contains

    function kernal_integral_equation(x, t)
        real(mp) :: x, t
        real(mp) :: kernal_integral_equation

            kernal_integral_equation = 1.0_mp / (1.0_mp + (x + 1.0_mp) * (t + 1.0_mp)) ** 2

    end function kernal_integral_equation

    function f(x) result(result_fun)
        real(mp) :: x, result_fun

            result_fun = log(2.0_mp * (x + 2.0_mp) / (2.0_mp * x + 3.0_mp)) / (x + 1.0_mp) &
            & - log(2.0_mp) / (x + 1.0_mp) / (2.0_mp * x + 3.0_mp)

    end function f

    function u_result(x)
        real(mp) :: u_result, x

        u_result = log(1.0_mp + x)

    end function u_result

end module init_data
