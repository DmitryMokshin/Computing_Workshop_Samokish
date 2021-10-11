module init_data
    implicit none
    integer, parameter :: mp = 8, num_of_coef = 5
    real(mp), parameter :: a = -1.0, b = 1.0, pi = acos(-1.0_mp), alpha=1.0
contains

    function kernal_integral_equation(x, t)
        real(mp) :: x, t
        real(mp) :: kernal_integral_equation

            ! kernal_integral_equation = 1.0_mp / (x + 1.0_mp) ** 2 / (t + 1.0_mp) ** 2
            kernal_integral_equation = x * t ** 2.0_mp - x

    end function kernal_integral_equation

    function f(x) result(result_fun)
        real(mp) :: x, result_fun

            ! result_fun = log(2.0_mp * (x + 2.0_mp) / (2.0_mp * x + 3.0_mp)) - log(2.0_mp) / (x + 1.0_mp) / (2.0_mp * x + 3.0_mp)
            result_fun = 1.0_mp + 4.0_mp / 3.0_mp * x
        
    end function f

end module init_data