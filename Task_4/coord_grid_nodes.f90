module coord_grid_nodes
    use :: init_data
    implicit none
contains

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

    REAL(mp) FUNCTION FindDet(matrix, n)
        IMPLICIT NONE
        REAL(mp), DIMENSION(n,n) :: matrix
        INTEGER, INTENT(IN) :: n
        REAL(mp) :: m, temp
        INTEGER :: i, j, k, l
        LOGICAL :: DetExists = .TRUE.
        l = 1
        !Convert to upper triangular form
        DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
        DetExists = .FALSE.
        DO i = k+1, n
        IF (matrix(i,k) /= 0) THEN
        DO j = 1, n
        temp = matrix(i,j)
        matrix(i,j)= matrix(k,j)
        matrix(k,j) = temp
        END DO
        DetExists = .TRUE.
        l=-l
        EXIT
        ENDIF
        END DO
        IF (DetExists .EQV. .FALSE.) THEN
        FindDet = 0
        return
        END IF
        ENDIF
        DO j = k+1, n
        m = matrix(j,k)/matrix(k,k)
        DO i = k+1, n
        matrix(j,i) = matrix(j,i) - m*matrix(k,i)
        END DO
        END DO
        END DO

        !Calculate determinant by finding product of diagonal elements
        FindDet = l
        DO i = 1, n
            FindDet = FindDet * matrix(i,i)
        END DO

    END FUNCTION FindDet

end module coord_grid_nodes