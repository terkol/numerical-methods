program polynomial_roots
    implicit none
    integer, parameter :: N = 4
    real(8), dimension(N) :: p
    complex(8), dimension(N-1) :: roots

    ! Coefficients of the polynomial P(x) = p(1) + p(2)*x + p(3)*x^2 + ... + p(N)*x^(N-1)
    p = (/ -6.0d0, 11.0d0, -6.0d0, 1.0d0 /)  ! Example polynomial x^3 - 6x^2 + 11x - 6

    call myroots(N, p, roots)

    print *, "Roots of the polynomial:"
    print *, roots

contains

    subroutine myroots(N, p, roots)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: p(N)
        complex(8), intent(out) :: roots(N-1)
        real(8), dimension(N-1, N-1) :: A
        real(8), dimension(N-1) :: wr, wi
        real(8), dimension(N-1, N-1) :: vl, vr
        integer :: info, lwork
        real(8), allocatable :: work(:)

        ! Construct the companion matrix
        A = 0.0d0
        A(2:N-1, 1:N-2) = 1.0d0
        A(:, N-1) = -p(1:N-1) / p(N)

        ! Allocate workspace for LAPACK
        lwork = -1
        allocate(work(lwork))

        ! Compute eigenvalues using LAPACK DGEEV
        call DGEEV('N', 'N', N-1, A, N-1, wr, wi, vl, N-1, vr, N-1, work, lwork, info)

        if (info /= 0) then
            print *, "Error computing eigenvalues, info =", info
            stop
        end if

        ! Store roots as complex numbers
        roots = cmplx(wr, wi, kind=8)

        ! Deallocate workspace
        deallocate(work)
    end subroutine myroots

end program polynomial_roots
