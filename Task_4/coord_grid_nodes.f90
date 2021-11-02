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