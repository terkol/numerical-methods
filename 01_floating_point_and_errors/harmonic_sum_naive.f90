program harmonic_sum
    implicit none
    integer, parameter :: rk64 = selected_real_kind(p=15, r=300)
    print*, "harmonic()", harmonic()
    print*, "harmonic_bunch(50)", harmonic_bunch(50)
    print*, "harmonic_bunch(100)", harmonic_bunch(100)
    print*, "harmonic_bunch(500)", harmonic_bunch(500)
contains
    function harmonic() result(sum)
        integer :: n
        real(rk64) :: sum
        sum = 0.0_rk64
        do n = 1, 1000000
            sum = sum + 1.0_rk64/real(n, rk64)     ! divide 1 by the current iteration and add the quotient to sum
        end do
    end function harmonic
    
    function harmonic_bunch(bunch) result(sum)
        integer :: bunch
        integer :: n
        integer :: m
        real(rk64) :: sum_bunch
        real(rk64) :: sum
        sum = 0.0_rk64
        do n = 1, 1000000/bunch ! million divided by the amount of bunches so the next loop can loop through the bunches themselves for a total of million iterations
            sum_bunch = 0.0_rk64
            do m = 1, bunch     ! Fill the bunch
                sum_bunch = sum_bunch + 1.0_rk64 / real((n - 1.0_rk64) * bunch + m, rk64) ! divide one by the (current iteration of the outer loop - 1) * size of the bunch + iteration of inner loop
            end do  
            sum = sum + sum_bunch ! after inner loop is over and the sum_bunch contains all the quotients of the bunch, the bunch is added to the total sum
        end do
    end function
end program harmonic_sum