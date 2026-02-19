program roots
    implicit none
    double precision :: a, b, x0, BB, b_root, n_root, g_root, start_time, end_time, tol
    ! Initial values
    a = 0.0d0  
    b = 0.5d0
    x0 = 1.0d0
    BB = 20.0d0
    tol = 1.0d-10

    ! Get roots and calculate execution times
    call CPU_TIME(start_time)
    b_root = bisect_f(a, b)
    call CPU_TIME(end_time)
    print*, "Bisection root: ", b_root, "Execution time: ", end_time - start_time

    call CPU_TIME(start_time)
    n_root = newton_f(x0)
    call CPU_TIME(end_time)
    print*, "Newton root:    ", n_root, "Execution time: ", end_time - start_time

    g_root = newton_g(x0, BB)
    print*, "Newton_g root:  ", g_root

contains
    function bisect_f(a, b) result(c)
        double precision :: a, b, c, fa, fb, fc
        integer :: n, nmax = 1000
        fa = f(a)
        fb = f(b)
        if (f(a)*f(b) < 0) then
            do n = 1, nmax          ! Loop
                if (abs(fc) > tol) then
                    c = (a + b)/2.0d0   ! Calculate middle point between a and b
                    fc = f(c)           ! Calculate f(x) at middle point
                    if (fa * fc < 0.0d0) then   ! If f(a) and f(c) have the same sign
                        b = c                   ! Make b equal to the middle point
                        fb = fc
                    else if (fb * fc < 0.0d0) then  ! If f(b) and f(c) have the same sign
                        a = c                       ! Make a equal to the middle point
                        fa = fc
                    endif
                endif
            enddo
        else
            print*,"Bisection doesn not bracket a root"
        endif
    end function bisect_f

    function newton_f(x0) result(x)
        double precision :: x0, x, fx, dfx, h = 1e-10
        integer :: n, nmax = 1000
        x = x0          ! Make x equal to given starting point
        do n = 1, nmax  ! Loop
            if (abs(fx) > tol) then
                fx = f(x)   ! Calculate f(x)
                dfx = (f(x+h)-f(x-h))/(2.0d0*h)     ! Calculate f'(x)
                if (abs(dfx) < h) then              ! If derivative is too small, exit
                    print*, "Derivative is too small"
                    exit
                endif
                x = x - fx/dfx  ! New value for x
            endif
        enddo
        end function newton_f

        function newton_g(x0, B) result(x)
            double precision :: x0, B, x, fx, dfx, h = 1.0d-8
            integer :: n, nmax = 1000
            x = x0          ! Make x equal to given starting point
            do n = 1, nmax  ! Loop
                fx = x + exp(-B*x**2)*cos(x)                    ! Calculate f(x)
                dfx = 1-exp(-B*x**2)*(2*B*x*cos(x)+sin(x))      ! Calculate f'(x)
                if (abs(dfx) < h) then                          ! If derivative is too small, exit
                    print*, "Derivative is too small"
                    exit
                endif
                x = x - fx/dfx  ! New value for x
            enddo
        end function newton_g
    
        function f(x) result(y)
            double precision :: x, y, pi
            pi = 4.0d0*datan(1.0d0)
            y = sin(3.0d0*pi*((x**3)/((x**2)-1.0d0))) + 0.5d0   ! Given function
        end function f
    
end program roots