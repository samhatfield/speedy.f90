module fourier
    use params

    implicit none

    private
    public initialize_fourier, fourier_inv, fourier_dir

    real :: wsave(kx,2*ix+15)

contains
    subroutine initialize_fourier
        integer k
        call rffti1(ix, wsave(1,ix+1), wsave(1,2*ix+1))

        do k = 2, kx
            wsave(k,2*ix+15) = wsave(1,2*ix+15)
        end do
    end subroutine

    ! From Fourier coefficients to grid-point data
    function fourier_inv(input, kcos, k) result(output)
    	use geometry, only: cosgr

        real, intent(in) :: input(2*mx,il)
        real :: output(ix,il)
        integer, intent(in) :: kcos, k
        integer :: j, m
        real :: fvar(ix)

        do j = 1,il
            fvar(1) = input(1,j)

            do m = 3, 2*mx
                fvar(m-1) = input(m,j)
            end do
            do m = 2*mx, ix
                fvar(m) = 0.0
            end do

            ! Inverse FFT
            call rfftb1(ix, fvar, wsave(k,:), wsave(k,ix+1), wsave(k,2*ix+1))

            ! Copy output into grid-point field, scaling by cos(lat) if needed
            if (kcos == 1) then
                output(:,j) = fvar
            else
                output(:,j) = fvar * cosgr(j)
            end if
        end do
    end function

    ! From grid-point data to Fourier coefficients
    function fourier_dir(input, k) result(output)
        real, intent(in) :: input(ix,il)
        integer, intent(in) :: k
        real :: output(2*mx,il)
        integer :: j, m
        real :: fvar(ix), scale

        ! Copy grid-point data into working array
        do j = 1, il
            fvar = input(:,j)

            ! Direct FFT
            call rfftf1(ix, fvar, wsave(k,:), wsave(k,ix+1), wsave(k,2*ix+1))

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
