subroutine geop(jj)
    ! subroutine geop (jj)
    !
    ! Purpose : compute spectral geopotential from spectral temperature T
    !           and spectral topography PHIS, as in GFDL Climate Group GCM
    ! Input :   jj = time level index (1 or 2)
    ! Modified common blocks : DYNSP2

    use mod_atparam
    use mod_dynvar
    use mod_dyncon1, only: xgeop1, xgeop2, hsg, fsg

    implicit none

    integer, intent(in) :: jj
    integer :: k
    real :: corf

    ! 1. Bottom layer (integration over half a layer)
    phi(:,:,kx) = phis + xgeop1(kx) * t(:,:,kx,jj)

    ! 2. Other layers (integration two half-layers)
    do k = kx-1,1,-1
        phi(:,:,k) = phi(:,:,k+1) + xgeop2(k+1)*t(:,:,k+1,jj)&
            & + xgeop1(k)*t(:,:,k,jj)
    end do

    ! 3. lapse-rate correction in the free troposphere
    do k = 2,kx-1
        corf=xgeop1(k)*0.5*log(hsg(k+1)/fsg(k))/log(fsg(k+1)/fsg(k-1))
        phi(1,:,k) = phi(1,:,k) + corf*(t(1,:,k+1,jj) - t(1,:,k-1,jj))
    end do
end
