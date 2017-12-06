subroutine grtend(vordt,divdt,tdt,psdt,trdt,j1,j2)
    !   subroutine grtend (vordt,divdt,tdt,psdt,trdt,j1,j2)
    !
    !   Purpose: compute non-linear tendencies in grid-point space
    !            from dynamics and physical parametrizations,
    !            and convert them to spectral tendencies
    !
    !   dF/dt = T_dyn(F(J2)) + T_phy(F(J1))
    !
    !   Input:  j1 = time level index for physical tendencies 
    !           j2 = time level index for dynamical tendencies 
    !   Output: vordt = spectral tendency of vorticity
    !           divdt = spectral tendency of divergence
    !           tdt   = spectral tendency of temperature
    !           psdt  = spectral tendency of log(p_s)
    !           trdt  = spectral tendency of tracers

    USE mod_atparam
    USE mod_dynvar
    use mod_dyncon1, only: akap, rgas, dhs, fsg, dhsr, fsgr, coriol
    use mod_dyncon2, only: tref, tref3

    implicit none

    !** notes ****
    ! -- TG does not have to be computed at both time levels every time step,
    !     I have left it this way to retain parallel structure with subroutine
    !     using latitude loop
    ! -- memory can be reduced considerably eliminating TGG, computing VORG
    !     only when needed, etc -- I have not optimized this subroutine for
    !     routine use on the YMP
    ! -- results from grtend1.F should duplicate results from grtend.F
    !                              -- Isaac
    !************

    complex, dimension(mx,nx,kx), intent(inout) :: vordt, divdt, tdt
    complex, intent(inout) :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
    integer, intent(in) :: j1, j2

    complex :: dumc(mx,nx,2)

    real, dimension(ix,il,kx) :: utend, vtend, ttend
    real :: trtend(ix,il,kx,ntr)

    real, dimension(ix,il,kx) :: ug, vg, tg, vorg, divg, tgg, puv
    real, dimension(ix,il) :: px, py, umean, vmean, dmean
    real :: trg(ix,il,kx,ntr), sigdt(ix,il,kxp)
    real :: temp(ix,il,kxp), sigm(ix,il,kxp)

    integer :: iitest = 0, k, i, itr, j

    if (iitest.eq.1) print*,'inside GRTEND'

    ! -------------
    ! Grid converts 
    do k=1,kx
        call grid(vor(:,:,k,j2),vorg(:,:,k),1)
        call grid(div(:,:,k,j2),divg(:,:,k),1)
        call grid(  t(:,:,k,j2),  tg(:,:,k),1)

        do itr=1,ntr
          call grid(tr(:,:,k,j2,itr),trg(:,:,k,itr),1)
        end do

        call uvspec(vor(:,:,k,j2),div(:,:,k,j2),dumc(:,:,1),dumc(:,:,2))
        call grid(dumc(:,:,2),vg(:,:,k),2)
        call grid(dumc(:,:,1),ug(:,:,k),2)

        do j=1,il
            do i=1,ix
                vorg(i,j,k)=vorg(i,j,k)+coriol(j)
            end do
        end do
    end do

    umean(:,:) = 0.0
    vmean(:,:) = 0.0
    dmean(:,:) = 0.0

    do k=1,kx
        umean(:,:) = umean(:,:) + ug(:,:,k) * dhs(k)
        vmean(:,:) = vmean(:,:) + vg(:,:,k) * dhs(k)
        dmean(:,:) = dmean(:,:) + divg(:,:,k) * dhs(k)
    end do

    ! Compute tendency of log(surface pressure)
    ! ps(1,1,j2)=zero
    call grad(ps(:,:,j2),dumc(:,:,1),dumc(:,:,2))
    call grid(dumc(:,:,1),px,2)
    call grid(dumc(:,:,2),py,2)

    call spec(-umean * px - vmean * py,psdt)
    psdt(1,1) = (0.0, 0.0)

    ! Compute "vertical" velocity  
    sigdt(:,:,1) = 0.0
    sigdt(:,:,kxp) = 0.0
    sigm(:,:,1) = 0.0
    sigm(:,:,kxp) = 0.0

    ! (The following combination of terms is utilized later in the 
    !     temperature equation)
    do k=1,kx
        puv(:,:,k) = (ug(:,:,k) - umean) * px + (vg(:,:,k) - vmean) * py
    end do

    do k=1,kx
        !cspj sigdt is the vertical velocity (in sigma coords)
        sigdt(:,:,k+1) = sigdt(:,:,k) - dhs(k)*(puv(:,:,k)+divg(:,:,k)-dmean)
        sigm(:,:,k+1) = sigm(:,:,k) - dhs(k)*puv(:,:,k)
    end do
 
    ! Subtract part of temperature field that is used as reference for 
    ! implicit terms
    do k=1,kx
        tgg(:,:,k) = tg(:,:,k) - tref(k)
    end do

    ! Zonal wind tendency
    temp(:,:,1) = 0.0
    temp(:,:,kxp) = 0.0

    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k) * (ug(:,:,k) - ug(:,:,k-1))
    end do

    do k=1,kx
        utend(:,:,k) = vg(:,:,k) * vorg(:,:,k) - tgg(:,:,k)*rgas*px - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
    end do

    ! Meridional wind tendency
    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k) * (vg(:,:,k) - vg(:,:,k-1))
    end do

    do k=1,kx
        vtend(:,:,k) = -ug(:,:,k)*vorg(:,:,k) - tgg(:,:,k)*rgas*py - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
    end do

    ! Temperature tendency
    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k)*(tgg(:,:,k) - tgg(:,:,k-1)) + sigm(:,:,k)*(tref(k) - tref(k-1))
    end do

    do k=1,kx
        ttend(:,:,k) = tgg(:,:,k)*divg(:,:,k) - (temp(:,:,k+1)+temp(:,:,k))*dhsr(k) &
            & + fsgr(k)*tgg(:,:,k)*(sigdt(:,:,k+1) + sigdt(:,:,k)) + tref3(k)*(sigm(:,:,k+1) + sigm(:,:,k)) &
            & + akap*(tg(:,:,k)*puv(:,:,k) - tgg(:,:,k)*dmean(:,:))
    end do

    ! Tracer tendency
    do itr=1,ntr 
        do k=2,kx
            temp(:,:,k) = sigdt(:,:,k)*(trg(:,:,k,itr) - trg(:,:,k-1,itr)) 
        end do

        !spj for moisture, vertical advection is not possible between top 
        !spj two layers
        !kuch three layers
        !if(iinewtrace.eq.1)then
        temp(:,:,2:3) = 0.0
        !endif

        do k=1,kx
            trtend(:,:,k,itr) = trg(:,:,k,itr)*divg(:,:,k)-(temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
        end do
    end do

    !******************** physics ****************************

    call geop(j1)

    call phypar(vor(:,:,:,j1),div(:,:,:,j1),t(:,:,:,j1),tr(:,:,:,j1,1),phi,&
        & ps(:,:,j1),utend,vtend,ttend,trtend)

    !*********************************************************

    do k=1,kx
        !  convert u and v tendencies to vor and div spectral tendencies
        !  vdspec takes a grid u and a grid v and converts them to 
        !  spectral vor and div
        call vdspec(utend(:,:,k), vtend(:,:,k), vordt(:,:,k), divdt(:,:,k), 2)

        !  divergence tendency
        !  add -lapl(0.5*(u**2+v**2)) to div tendency,
        call spec(0.5 * (ug(:,:,k)*ug(:,:,k) + vg(:,:,k)*vg(:,:,k)), dumc(:,:,1))
        call lap(dumc(:,:,1), dumc(:,:,2))
        divdt(:,:,k) = divdt(:,:,k) - dumc(:,:,2)

        !  temperature tendency
        !  and add div(vT) to spectral t tendency
        call vdspec(-ug(:,:,k)*tgg(:,:,k), -vg(:,:,k)*tgg(:,:,k), dumc(:,:,1), tdt(:,:,k), 2)
        call spec(ttend(:,:,k), dumc(:,:,2))
        tdt(:,:,k) = tdt(:,:,k) + dumc(:,:,2)

        ! tracer tendency
        do itr=1,ntr
            call vdspec(-ug(:,:,k)*trg(:,:,k,itr), -vg(:,:,k)*trg(:,:,k,itr), dumc(:,:,1), trdt(:,:,k,itr), 2)
            call spec(trtend(:,:,k,itr), dumc(:,:,2))
            trdt(:,:,k,itr) = trdt(:,:,k,itr) + dumc(:,:,2)
        end do
    end do
end 
