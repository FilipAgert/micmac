module test_utils
    use constants
    implicit none
    private
    public :: eq, eq_r, eq_mat, print_matrix, eq_mat_comm

    real(r_kind), parameter :: epsilon = 1e-6




    contains

    logical function eq(n1,n2)
        integer, intent(in) :: n1,n2
        eq = n1 == n2
    end function

    logical function eq_r(r1,r2)
        real(r_kind), intent(in) :: r1, r2
        if(abs(r1-r2) < epsilon) then
            eq_r = .true.
        else
            eq_r = .false.
        endif
    end function

    logical function eq_mat(mat1,mat2)
        real(r_kind), intent(in), dimension(:,:) :: mat1, mat2
        integer :: row, col
        do row = 1, size(mat1,1)
            do col = 1,size(mat2,2)
                if(abs(mat1(row,col)-mat2(row,col)) < epsilon) then
                    eq_mat = .true.
                else
                    eq_mat = .false.
                    return
                endif
            end do
        end do
    end function


    logical function eq_mat_comm(mat1,mat2,discardN) !!In truncated basis, last rows and cols wont commute. dont compare with them
        real(r_kind), intent(in), dimension(:,:) :: mat1, mat2
        integer :: row, col, discardN
        do row = 1, size(mat1,1) - discardN
            do col = 1,size(mat2,2) - discardN
                if(abs(mat1(row,col)-mat2(row,col)) < epsilon) then
                    eq_mat_comm = .true.
                else
                    eq_mat_comm = .false.
                    return
                endif
            end do
        end do
    end function
    subroutine print_matrix(A)
        implicit none
        real(r_kind), dimension(:,:), intent(in) :: A
        integer :: i, j, n, m
    
        n = size(A, 1)
        m = size(A, 2)
        write(*,'(A,f7.1)') "Sum abs matrix:", sum(abs(A))
        ! do i = 1, n
        !     write(*,'(100(f5.1,1x))') (A(i,j), j=1,m)
        ! end do
    end subroutine print_matrix
    

end module test_utils