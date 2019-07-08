!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 01/05/2019
!  For computing direct and inverse Fourier transforms.
module fourier
    use types, only: p
    use params

    implicit none

    private
    public initialize_fourier, fourier_inv, fourier_dir

    real(p) :: work(ix) !! Work array required by FFTPACK. Contains trigonometric functions etc.
    integer :: ifac(15) !! Work array required by FFTPACK. Contains prime factors

contains
    !> Initializes the Fourier transforms.
    subroutine initialize_fourier
        call rffti1(ix, work, ifac)
    end subroutine

    !> Transforms Fourier coefficients to grid-point data.
    function fourier_inv(input, kcos) result(output)
        use geometry, only: cosgr

        real(p), intent(in) :: input(2*mx,il) !! Input field
        integer, intent(in) :: kcos           !! Scale output by cos(lat) (1) or not (0)
        real(p)             :: output(ix,il)  !! Output field

        integer :: j, m
        real(p) :: fvar(ix), ch(ix)

        do j = 1,il
            fvar(1) = input(1,j)

            do m = 3, 2*mx
                fvar(m-1) = input(m,j)
            end do
            do m = 2*mx, ix
                fvar(m) = 0.0
            end do

            ! Inverse FFT
            call rfftb1(ix, fvar, ch, work, ifac)

            ! Copy output into grid-point field, scaling by cos(lat) if needed
            if (kcos == 1) then
                output(:,j) = fvar
            else
                output(:,j) = fvar * cosgr(j)
            end if
        end do
    end function

    !> Transforms grid-point data to Fourier coefficients.
    function fourier_dir(input) result(output)
        real(p), intent(in) :: input(ix,il)    !! Input field
        real(p)             :: output(2*mx,il) !! Output field

        integer :: j, m
        real(p) :: fvar(ix), scale
        real(p) :: ch(ix)

        ! Copy grid-point data into working array
        do j = 1, il
            fvar = input(:,j)

            ! Direct FFT
            call rfftf1(ix, fvar, ch, work, ifac)

            ! Copy output into spectral field, dividing by no. of long.
            scale = 1.0/float(ix)

            ! Mean value (a(0))
            output(1,j) = fvar(1)*scale
            output(2,j) = 0.0

            do m = 3, 2*mx
                output(m,j) = fvar(m-1)*scale
            end do
        end do
    end function
end module
