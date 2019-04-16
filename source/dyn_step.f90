subroutine step(j1,j2,dt,alph)
    !   subroutine step (j1,j2,dt,alph)
    !
    !   Purpose: perform one time step starting from F(1) and F(2)
    !            and using the following scheme:
    !
    !   Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
    !   F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
    !   F(2) = Fnew
    !
    !   Input:
    !   If J1=1, J2=1 : forward time step (eps=0)
    !   If J1=1, J2=2 : initial leapfrog time step (eps=0)
    !   If J1=2, J2=2 : leapfrog time step with time filter (eps=ROB)
    !   DT = time step (if DT < or = 0, tendencies are computed but
    !                   no time stepping is performed)
    !   alph = 0   : forward step for gravity wave terms
    !   alph = 1   : backward implicit step for g.w.
    !   alph = 0.5 : centered implicit step for g.w.

    use mod_dyncon0, only: tdrs
    use mod_atparam
    use mod_dynvar
    use mod_hdifcon
    use mod_tsteps, only: rob, wil

    implicit none

    integer, intent(in) :: j1, j2
    real, intent(in) :: dt, alph
    complex, dimension(mx,nx,kx) ::  ordt, divdt, tdt, vordt
    complex :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
    real :: eps, sdrag

    complex :: ctmp(mx,nx,kx)

    integer :: n, itr, k, m

    ! 1. Computation of grid-point tendencies
    ! (converted to spectral at the end of GRTEND)
    call grtend(vordt,divdt,tdt,psdt,trdt,1,j2)

    ! 2. Computation of spectral tendencies
    if (alph.eq.0.) then
        call sptend(divdt,tdt,psdt,j2)
    else
        call sptend(divdt,tdt,psdt,1)

        ! implicit correction
        call implic(divdt,tdt,psdt)
    endif

    ! 3. Horizontal diffusion

    ! 3.1 Diffusion of wind and temperature
    call hordif(kx,vor,vordt,dmp, dmp1)
    call hordif(kx,div,divdt,dmpd,dmp1d)

    do k=1,kx
        do m=1,mx
            do n=1,nx
                ctmp(m,n,k) = t(m,n,k,1)+tcorh(m,n)*tcorv(k)
            enddo
        enddo
    enddo

    call hordif(kx,ctmp,tdt,dmp,dmp1)

    ! 3.2 Stratospheric diffusion and zonal wind damping
    sdrag = 1./(tdrs*3600.)
    do n = 1,nx
        vordt(1,n,1) = vordt(1,n,1)-sdrag*vor(1,n,1,1)
        divdt(1,n,1) = divdt(1,n,1)-sdrag*div(1,n,1,1)
    enddo

    call hordif(1,vor, vordt,dmps,dmp1s)
    call hordif(1,div, divdt,dmps,dmp1s)
    call hordif(1,ctmp,tdt,  dmps,dmp1s)

    ! 3.4 Diffusion of tracers
    do k=1,kx
        do m=1,mx
            do n=1,nx
                ctmp(m,n,k) = tr(m,n,k,1,1)+qcorh(m,n)*qcorv(k)
            enddo
        enddo
    enddo

    call hordif(kx,ctmp,trdt,dmpd,dmp1d)

    if (ntr.gt.1) then
        do itr=2,ntr
            call hordif(kx,tr(1,1,1,1,itr),trdt(1,1,1,itr),dmp,dmp1)
        enddo
    endif

    ! 4. Time integration with Robert filter
    if (dt.le.0.) return

    if (j1.eq.1) then
        eps = 0.
    else
        eps = rob
    endif

    call timint(j1,dt,eps,wil,1,ps,psdt)

    call timint(j1,dt,eps,wil,kx,vor,vordt)
    call timint(j1,dt,eps,wil,kx,div,divdt)
    call timint(j1,dt,eps,wil,kx,t,  tdt)

    do itr=1,ntr
        call timint(j1,dt,eps,wil,kx,tr(1,1,1,1,itr),trdt(1,1,1,itr))
    enddo
end

subroutine hordif(nlev,field,fdt,dmp,dmp1)
    !   Aux. subr. HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
    !   Purpose : Add horizontal diffusion tendency of FIELD
    !             to spectral tendency FDT at NLEV levels
    !             using damping coefficients DMP and DMP1

    USE mod_atparam

    implicit none

    integer, intent(in) :: nlev
    complex, intent(in) :: field(mxnx,kx)
    complex, intent(inout) :: fdt(mxnx,kx)
    real, intent(in) :: dmp(mxnx), dmp1(mxnx)
    integer :: k, m

    do k=1,nlev
        do m=1,mxnx
            fdt(m,k)=(fdt(m,k)-dmp(m)*field(m,k))*dmp1(m)
        enddo
    enddo
end

subroutine timint(j1,dt,eps,wil,nlev,field,fdt)
    !  Aux. subr. timint (j1,dt,eps,wil,nlev,field,fdt)
    !  Purpose : Perform time integration of field at nlev levels
    !            using tendency fdt

    use mod_atparam

    implicit none

    integer, intent(in) :: j1, nlev
    real, intent(in) :: dt, eps, wil
    complex, intent(in) :: fdt(mxnx,nlev)
    complex, intent(inout) :: field(mxnx,nlev,2)
    real :: eps2
    complex :: fnew(mxnx)
    integer :: k, m

    eps2 = 1.-2.*eps

    if (ix.eq.iy*4) then
        do k=1,nlev
            call trunct(fdt(1,k))
        enddo
    endif

    ! The actual leap frog with the robert filter
    do k=1,nlev
        do m=1,mxnx
            fnew (m)     = field(m,k,1) + dt*fdt(m,k)
            field(m,k,1) = field(m,k,j1) +  wil*eps*(field(m,k,1)&
                & -2*field(m,k,j1)+fnew(m))

            ! and here comes Williams' innovation to the filter
            field(m,k,2) = fnew(m)-(1-wil)*eps*(field(m,k,1)&
                &-2*field(m,k,j1)+fnew(m))
        enddo
    enddo
end