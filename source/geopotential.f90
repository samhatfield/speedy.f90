module geopotential
    use params

    implicit none

    private
    public initialize_geopotential, get_geopotential

    ! Constants for hydrostatic equation
    real :: xgeop1(kx), xgeop2(kx)

contains
    subroutine initialize_geopotential
        use physical_constants, only: rgas
        use geometry, only: hsg, fsg

        integer :: k

        ! Coefficients to compute geopotential
        do k = 1, kx
          xgeop1(k) = rgas*log(hsg(k+1)/fsg(k))
          if (k /= kx) xgeop2(k+1) = rgas*log(fsg(k+1)/hsg(k+1))
        end do
    end subroutine

    ! Compute spectral geopotential from spectral temperature T and spectral
    ! topography phis, as in GFDL Climate Group GCM
    function get_geopotential(t, phis) result(phi)
        use geometry, only: hsg, fsg

        complex, intent(in) :: t(mx,nx,kx), phis(mx,nx)
        complex :: phi(mx,nx,kx)
        integer :: k
        real :: corf

        ! 1. Bottom layer (integration over half a layer)
        phi(:,:,kx) = phis + xgeop1(kx) * t(:,:,kx)

        ! 2. Other layers (integration two half-layers)
        do k = kx-1,1,-1
            phi(:,:,k) = phi(:,:,k+1) + xgeop2(k+1)*t(:,:,k+1)&
                & + xgeop1(k)*t(:,:,k)
        end do

        ! 3. lapse-rate correction in the free troposphere
        do k = 2,kx-1
            corf=xgeop1(k)*0.5*log(hsg(k+1)/fsg(k))/log(fsg(k+1)/fsg(k-1))
            phi(1,:,k) = phi(1,:,k) + corf*(t(1,:,k+1) - t(1,:,k-1))
        end do
    end function
end module
