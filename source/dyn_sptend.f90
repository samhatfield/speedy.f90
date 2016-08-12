subroutine sptend (divdt,tdt,psdt,j4)
    !   subroutine sptend (divdt,tdt,psdt,j4)
    !
    !   Purpose : compute spectral tendencies of divergence, temperature
    !             and log_surf.pressure)
    !   Input/output : divdt = divergence tendency (spec.)
    !                  tdt   = temperature tendency (spec.)
    !                  psdt  = tendency of log_surf.pressure (spec.)
    !                  j4    = time level index (1 or 2)
    
    use mod_atparam
    use mod_dynvar
    use mod_dyncon1, only: rgas, dhs, dhsr
    use mod_dyncon2, only: tref, tref2, tref3

    implicit none

    complex, intent(inout) :: psdt(mx,nx), divdt(mx,nx,kx), tdt(mx,nx,kx)
    integer, intent(in) :: j4

    complex :: dumk(mx,nx,kxp), dmeanc(mx,nx), sigdtc(mx,nx,kxp)
    complex :: tempc(mx,nx,3)
    complex :: dumc(mx,nx,2), zero

    integer :: k

    zero = (0.,0.)

    ! Vertical mean div and pressure tendency
    dmeanc(:,:) = zero
    do k=1,kx
        dmeanc = dmeanc + div(:,:,k,j4) * dhs(k)
    end do

    psdt = psdt - dmeanc
    psdt(1,1) = zero

    ! Sigma-dot "velocity" and temperature tendency
    sigdtc(:,:,1) = zero
    sigdtc(:,:,kxp) = zero

    do k=1,kxm
        sigdtc(:,:,k+1) = sigdtc(:,:,k) - dhs(k)*(div(:,:,k,j4) - dmeanc)
    end do

    dumk(:,:,1) = zero
    dumk(:,:,kxp) = zero

    do k=2,kx
        dumk(:,:,k) = sigdtc(:,:,k) * (tref(k) - tref(k-1))
    end do

    do k=1,kx
        tdt(:,:,k) = tdt(:,:,k) - (dumk(:,:,k+1) + dumk(:,:,k)) * dhsr(k)&
            & + tref3(k) * (sigdtc(:,:,k+1) + sigdtc(:,:,k))&
            & - tref2(k) * dmeanc
    end do

    ! Geopotential and divergence tendency
    call geop(j4)

    do k=1,kx
        dumc(:,:,1) = phi(:,:,k) + rgas*tref(k)*ps(:,:,j4)
        call lap(dumc(1,1,1),dumc(1,1,2))
        divdt(:,:,k) = divdt(:,:,k) - dumc(:,:,2)
    end do
end   
