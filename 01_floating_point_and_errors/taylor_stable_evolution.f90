program funk
    implicit none
    integer :: n
    double precision :: half
    double precision :: f1_dir, f2_dir, f1_tay, f2_tay
    ! Call funkx subroutine 100 times, 50 times approaching zero and 50 times increasing the distance to zero
    do n = 1, 100
    ! Dividing n with 10000 to get closer to zero and then subtracting half to get the range of the loop to both sides of zero:
        half = dble(n)/10000.0d0 - 0.005d0

        call funkx(half, f1_dir, f2_dir)
        call funkx_taylor(half, f1_tay, f2_tay)

        print*, "x:", half, &
                " f1_direct:", f1_dir, " f1_taylor:", f1_tay, " Δf1:", f1_dir - f1_tay, &
                " f2_direct:", f2_dir, " f2_taylor:", f2_tay, " Δf2:", f2_dir - f2_tay
    end do

contains
    subroutine funkx(x, f1, f2)
        implicit none
        double precision, intent(in)  :: x
        double precision, intent(out) :: f1, f2

        if (x == 0.0d0) then
            f1 = 0.0d0
            f2 = 0.0d0
        else
            f1 = (cos(x) - 1.0d0) / (x ** 2.0d0)
            f2 = (exp(x) - exp(-x)) / (2.0d0 * x)
        end if
    end subroutine funkx

    subroutine funkx_taylor(x, f1, f2)
        implicit none
        double precision, intent(in)  :: x
        double precision, intent(out) :: f1, f2

        f1 = -(1.0d0/2.0d0) + (x**2.0d0)/24.0d0
        f2 =  1.0d0 + (x**2.0d0)/6.0d0
    end subroutine funkx_taylor
end program funk
