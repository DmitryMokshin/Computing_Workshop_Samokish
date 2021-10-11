module SSE
    use init_data
implicit none
contains

    subroutine SSE_Mod_Gauss(A,B,x)
        real(mp) :: s, tmp
        real(mp), dimension(:, :) :: A
        real(mp), dimension(:) :: B,x
        real(mp), dimension(size(x, dim=1)+1) :: tr
        real(mp), dimension(size(x, dim=1), size(x, dim=1) + 1) :: inter_mat
        real(mp), parameter :: eps=1.0_mp / 10000.0_mp
        integer :: n, i, j, k, p, q, coord(2), cdr(size(x,dim=1))

        n = size(x,dim=1)
        s = 0;
        do i = 1, n
            inter_mat(i, 1:n) = A(i, :)
        enddo
        inter_mat(:, n + 1)=B
    
        do k = 1, n		   
            coord=maxloc(abs(inter_mat(k:n, k:n)))
            p = coord(1) + k - 1
            q = coord(2) + k - 1			
            cdr(k) = q				
            tr(1:n) = inter_mat(:, q)
            inter_mat(:, q) = inter_mat(:, k)
            inter_mat(:, k) = tr(1:n)			
            tr(:) = inter_mat(p, :)
            inter_mat(p, :) = inter_mat(k, :)
            inter_mat(k, :) = tr(:)	
            forall(i = k+1:n, j = k:n+1)
                inter_mat(i, j) = inter_mat(i, j) - inter_mat(i, k) * inter_mat(k, j)&
                &/ inter_mat(k, k)
            end forall   
        end do
			
        forall(i = 1:n)
            forall(j = i:n+1)
                inter_mat(i, j) = inter_mat(i, j) / inter_mat(i, i)
            end forall
        end forall
		
        x(n) = inter_mat(n, n + 1)
        do i = n-1, 1, -1
            s = 0;
            do j = i + 1,n
                s = s + inter_mat(i, j) / inter_mat(i, i) * x(j)
            end do
            x(i) = inter_mat(i, n + 1) - s
        end do

        do i = n, 1, -1	
            tmp = x(i)
            x(i) = x(cdr(i))
            x(cdr(i)) = tmp	
        end do		
	
    end subroutine SSE_Mod_Gauss

end module SSE