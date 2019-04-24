module fourier
    use mod_atparam

    implicit none

    private
    public initialize_fourier, gridx, specx

    real, dimension(2*ix+15) :: wsave(2*ix+15)

contains
    subroutine initialize_fourier
        call rffti(ix,wsave)
    end subroutine

    ! From Fourier coefficients to grid-point data
    subroutine gridx(varm,vorg,kcos)
    	use geometry, only: cosgr

        real, intent(in) :: varm(mx2,il)
        real, intent(inout) :: vorg(ix,il)
        integer, intent(in) :: kcos
        integer :: j, m
        real :: fvar(ix)

    	do j = 1,il
    		fvar(1) = varm(1,j)

            do m=3,mx2
              fvar(m-1)=varm(m,j)
            end do
            do m=mx2,ix
              fvar(m)=0.0
            end do

            ! Inverse FFT
            call rfftb(ix,fvar,wsave)

            ! Copy output into grid-point field, scaling by cos(lat) if needed
            if (kcos.eq.1) then
                vorg(:,j) = fvar
            else
                vorg(:,j) = fvar * cosgr(j)
            end if
        end do
    end

    ! From grid-point data to Fourier coefficients
    subroutine specx(vorg,varm)
        real, intent(in) :: vorg(ix,il)
        real, intent(inout) :: varm(mx2,il)
        integer :: j, m
        real :: fvar(ix), scale

        ! Copy grid-point data into working array
        do j=1,il
            fvar = vorg(:,j)

            ! Direct FFT
            CALL RFFTF (IX,FVAR,WSAVE)
            !CALL DFFTF (IX,FVAR,WSAVE)

            ! Copy output into spectral field, dividing by no. of long.
            scale=1./float(ix)

            ! Mean value (a(0))
            varm(1,j)=fvar(1)*scale
            varm(2,j)=0.0

            do m=3,mx2
                varm(m,j)=fvar(m-1)*scale
            end do
        end do
    end
end module
