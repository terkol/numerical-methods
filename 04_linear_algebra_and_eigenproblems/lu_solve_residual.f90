program problem4
    implicit none
    integer :: N, i, j, m
    double precision, allocatable :: A(:,:), A_copy(:,:), x(:), b(:), b_copy(:)
    ! Read given file to get N, A, b
    read(5,*) N
    print "(a, i0)", "N = ", N
    allocate(A(N, N), b(N), x(N))
    do i = 1, N
        read(5,*) (A(i, j), j = 1, N)
    end do
    read(5,*) (b(i), i = 1, N)
    allocate(A_copy(N,N), b_copy(N))
    A_copy = A
    b_copy = b
    call solvex(N, A_copy, b_copy, x)   ! Solve x
    ! Set m
    m = 0

    if (m > 30 .or. m == 0) then
        print "(a,a,a,e12.6)", "Norm-", "inf", " of residual: ", residual(N, A, x, b, huge(m))  ! Print this if m is 0 or very big
    else 
        print "(a,i0,a,e12.6)", "Norm-", m, " of residual: ", residual(N, A, x, b, m) ! Else print this
    end if
contains

    function residual(N, A, x, b, m) result(norm)   ! Return norm of residue of given matrix equation
        integer :: N, m
        double precision :: A(N, N), x(N), b(N), norm, Ax(N), res(N)
        Ax = matmul(A, x)   ! Do matrix multiplication
        res = Ax - b    ! Calculate residue
        norm = (sum(abs(res)**m))**(1.0d0/m)    ! Calculate norm for given m
    end function residual

    subroutine solvex(N, A, b, x)   ! Solve x
        integer :: N
        double precision :: A(N,N), b(N), x(N)
        integer :: nrhs=1, pivot(N), ok
        call dgesv(N, nrhs, A, N, pivot, b, N, ok)  ! Use lapack routine to calculate new b
        if (ok /= 0) then
            stop    ! If something went wrong, stop
        end if
        x = b   ! Make x equal to new b
    end subroutine solvex 

end program problem4
