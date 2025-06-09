module Hamiltonian
    use constants
    implicit none


    private
    public :: diagonalize


    contains

    subroutine diagonalize(E, V, H) !!Diagnoalize H and return V and E as eigenvectors sorted with lowest energy first.
        real(r_kind), intent(in) :: H(:,:)
        real(r_kind), intent(out) :: E(size(H,1)), V(size(H,1),size(H,1))
        external dsyev


    end subroutine
end module