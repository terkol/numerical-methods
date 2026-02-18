program kahan
    implicit none
    integer, parameter :: rk64 = selected_real_kind(p=15, r=300)
    integer :: n
    real(rk64) :: last, current, lastdelta
    last = 0.0_rk64
    lastdelta = 0.0_rk64
    do n = 1, 9
        current = harmonic_kahan(10**n)
        print*, "N: ",10**n, "S:",  current, "ΔS:", current-last, "Δ(ΔS):", lastdelta - (current - last)
        lastdelta = current - last 
        last = current 
    enddo
contains
    function harmonic_kahan(max) result(s)
        integer, intent(in) :: max
        integer :: m
        real(rk64) :: s, x, y, t, e
        s = 0.0_rk64 
        e = 0.0_rk64 
        do m = 1, max 
            x = 1.0_rk64/real(m, rk64)  
            y = x - e
            t = s + y
            e = (t - s) - y
            s = t
        enddo
    end function harmonic_kahan
end program kahan