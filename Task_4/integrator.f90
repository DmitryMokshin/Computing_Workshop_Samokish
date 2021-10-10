module gauss_quad
    use init_data
    implicit none
contains

    function gauss_quad_integral(left_int, right_int, n, fun_for_int) result(integral)
        real(mp) :: integral, left_int, right_int
        real(mp), dimension(1:n) :: solution_legendre, coefficient_gauss_quad
        integer :: n, i
        character(2) :: arg_num
        interface
            function fun_for_int(x) result(f_res)
                use init_data
            implicit none
                real(mp) :: x
                real(mp) :: f_res
            end function fun_for_int
        end interface
        
        arg_num = char(ichar('0')+n/10)//char(ichar('0')+mod(n,10))
        open(10, file='quad'//arg_num//'.dat', status='old')

        do i = 1, n
            read(10,*) coefficient_gauss_quad(i), solution_legendre(i)
        end do

        integral = 0.0_mp
        do i = 1, n
            integral = integral + coefficient_gauss_quad(i) * fun_for_int((right_int + left_int) &
            & / 2.0_mp + (right_int - left_int) * solution_legendre(i) / 2.0_mp)
        end do
        integral = integral * (right_int - left_int) / 2.0_mp

        close(10)
        
    end function gauss_quad_integral

end module gauss_quad