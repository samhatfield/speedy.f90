! Read topography and climatological boundary conditions
subroutine inbcon
    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps, only: isst0
    use mod_atparam
    use mod_surfcon
    use mod_cli_land
    use mod_cli_sea
    use mod_dyncon1, only: grav, radang

    implicit none

    real*4 :: r4inp(ix,il)
    real   :: inp(ix,il)
    real   :: veg(ix,il), swl1(ix,il), swl2(ix,il)

    integer :: i, idep2, irec, irecl, it, j, jrec
    real :: rad2deg, rsw, sdep1, sdep2, swroot, swwil2, thrsh

    ! Set threshold for land-sea mask definition
    ! (ie minimum fraction of either land or sea)
    thrsh = 0.1

    ! Read surface geopotential (i.e. orography)
    call load_boundary_file(1,20,inp,0)
    phi0 = grav*inp

    ! Also store spectrally truncated surface geopotential
    call truncg (ntrun,phi0,phis0)

    ! Read land-sea mask
    call load_boundary_file(1,20,fmask,1)

    ! Initialize land-sfc boundary conditions

    ! Fractional and binary land masks
    fmask_l = fmask
    do j=1,il
        do i=1,ix
            if (fmask_l(i,j) >= thrsh) then
                bmask_l(i,j) = 1.0
                if (fmask(i,j) > (1.0 - thrsh)) fmask_l(i,j) = 1.0
            else
                bmask_l(i,j) = 0.0
                fmask_l(i,j) = 0.0
            end if
        end do
    end do

    ! Annual-mean surface albedo
    call load_boundary_file(1,20,alb0,2)

    ! Land-surface temp.
    do it = 1,12
        call load_boundary_file(1,23,inp,it-1)

        call fillsf(inp,ix,il,0.)

       stl12(1:ix,1:il,it) = inp
    end do

    call forchk(bmask_l, 12, 0.0, 400.0, 273.0, stl12)

    ! Snow depth
    do it = 1,12
        call load_boundary_file(1,24,inp,it-1)

        snowd12(1:ix,1:il,it) = inp
    end do

    CALL forchk(bmask_l, 12, 0.0, 20000.0, 0.0, snowd12)

    ! Read soil moisture and compute soil water availability using vegetation fraction
    ! Read vegetation fraction
    call load_boundary_file(1,20,veg,3)
    call load_boundary_file(1,20,inp,4)

    ! Combine high and low vegetation fractions
    veg = max(0.,veg+0.8*inp)

    ! Read soil moisture
    sdep1 = 70.
    idep2 = 3
    sdep2 = idep2*sdep1

    swwil2= idep2*swwil
    rsw   = 1./(swcap+idep2*(swcap-swwil))

    do it = 1,12
        call load_boundary_file(1,26,swl1,3*it-3)
        call load_boundary_file(1,26,swl2,3*it-2)

        ! Combine soil water content from two top layers
        do j = 1,il
            do i = 1,ix
                swroot = idep2*swl2(i,j)
                inp(i,j) = min(1.,rsw*(swl1(i,j)+veg(i,j)*max(0.,swroot-swwil2)))
            end do
        end do

        soilw12(1:ix,1:il,it) = inp
    end do

    call forchk(bmask_l, 12, 0.0, 10.0, 0.0, soilw12)

    ! Initialize sea-surface boundary conditions

    ! Fractional and binary sea masks
    do j=1,il
        do i=1,ix
            fmask_s(i,j) = 1.0 - fmask(i,j)

            if (fmask_s(i,j) >= thrsh) then
                bmask_s(i,j) = 1.0
                if (fmask_s(i,j) > (1.0 - thrsh)) fmask_s(i,j) = 1.0
            else
                bmask_s(i,j) = 0.0
                fmask_s(i,j) = 0.0
            end if
        end do
    end do

    ! Grid latitudes for sea-surface variables
    rad2deg = 90.0/asin(1.)
    deglat_s = rad2deg*radang

    ! SST
    do it = 1,12
        call load_boundary_file(1,21,inp,it-1)

        call fillsf(inp,ix,il,0.)

        sst12(1:ix,1:il,it) = inp
    end do

    call forchk(bmask_s, 12, 100.0, 400.0, 273.0, sst12)

    ! Sea ice concentration
    do it = 1,12
        call load_boundary_file(1,22,inp,it-1)

        inp = max(inp,0.)

        sice12(1:ix,1:il,it) = inp
    end do

    call forchk(bmask_s, 12, 0.0, 1.0, 0.0, sice12)

    ! SST anomalies for initial and prec./following months
    if (isstan > 0) then
        print *, 'isst0 = ', isst0
        do it=1,3
            if ((isst0 <= 1 .and. it /= 2) .or. isst0 > 1) then
                call load_boundary_file(1,30,inp,isst0-2+it-1)
            end if

            sstan3(1:ix,1:il,it) = inp
        end do

        call forchk(bmask_s, 3, -50.0, 50.0, 0.0, sstan3)
    end if

    ! Climatological fields for the ocean model (TO BE RECODED)
    ! Annual-mean heat flux into sea-surface`
    hfseacl = 0.0

    if (icsea >= 1) then
        irecl = 4*ix*il
        irec = 0

        open ( unit=31, file='fort.31', status='old',&
            & form='unformatted', access='direct', recl=irecl )

        do it = 1,12
            irec=irec+2
            read (31,rec=irec) r4inp

            do j = 1,il
                do i = 1,ix
                    hfseacl(i,j) = hfseacl(i,j)+r4inp(i,j)
                end do
            end do
        end do

        do j = 1,il
            do i = 1,ix
                if (bmask_s(i,j).gt.0.) then
                    hfseacl(i,j) = hfseacl(i,j)/(12.*fmask_s(i,j))
                else
                    hfseacl(i,j) = 0.
                end if
            end do
        end do

        call forchk (bmask_s, 1, -1000.0, 1000.0, 0.0, hfseacl)
    end if

    ! Ocean model SST climatology:
    ! defined by adding SST model bias to obs. climatology
    ! (bias may be defined in a different period from climatology)

    if (icsea >= 3) then
        do it = 1,12
            read (32) r4inp

            do j = 1,il
                do i = 1,ix
                    sstom12(i,j,it) = sst12(i,j,it)+r4inp(i,j)
                end do
            end do
        end do

        call forchk (bmask_s, 12, 100.0, 400.0, 273.0, sstom12)
    end if
end

! Check consistency of sfc fields with land-sea mask and set undefined values to a constant
! (to avoid over/underflow)
subroutine forchk (fmask,nf,fmin,fmax,fset,field)
    use mod_atparam, only: ix, il

    implicit none

    real, intent(in) :: fmask(ix,il)
    integer, intent(in) :: nf
    real, intent(in) :: fmin, fmax, fset
    real, intent(inout) :: field(ix,il,nf)

    integer :: i, j, jf, nfault

    do jf = 1, nf
        nfault = 0

        do i = 1, ix
            do j = 1, il
                if (fmask(ix,il) > 0.0) then
                    if (field(ix,il,jf) < fmin .or. field(ix,il,jf) > fmax) then
                        nfault = nfault + 1
                    end if
                else
                    field(ix,il,jf) = fset
                end if
            end do
        end do

        print *, 'Number of faulty points for field: ', jf, ' = ',nfault
    end do

    print *, 'Undefined values set to ', fset
end

! Compute a spectrally-filtered grid-point field
! Input   : itr : spectral truncation (triangular)
!         : fg1 : original grid-point field
! Output  : fg2 : filtered grid-point field
subroutine truncg (itr,fg1,fg2)
    use mod_atparam

    implicit none

    integer, intent(in) :: itr

    real, dimension(ix,il), intent(inout) :: fg1 (ix,il), fg2(ix,il)
    complex :: fsp(mx,nx), zero
    integer :: n, m, itwn

    print *, 'Filter applied at wavenumber ', itr

    zero = (0.,0.)

    call spec (fg1,fsp)

    do n=1,nx
        do m=1,mx
            itwn=isc*(m-1)+n-1
            if (itwn.gt.itr) fsp(m,n)=zero
        end do
    end do

    call grid (fsp,fg2,1)
end

! Replace missing values in surface fields
! NB: it is assumed that non-missing values exist near the Equator
subroutine fillsf(sf,ix,il,fmis)
    implicit none

    real :: sf(ix,il), sf2(0:ix+1)
    integer, intent(in) :: ix, il
    real, intent(in) :: fmis

    integer :: khem, j, j1, j2, j3, i, nmis
    real :: fmean

    do khem = 1,2
       if (khem == 1) then
            j1 = il/2
            j2 = 1
            j3 = -1
        else
            j1 = j1+1
            j2 = il
            j3 = 1
        end if

        do j=j1,j2,j3
            sf2(1:ix) = sf(1:ix,j)

            nmis = 0
            do i=1,ix
                if (sf(i,j) < fmis) then
                    nmis = nmis+1
                    sf2(i) = 0.
                end if
            end do

            if (nmis < ix) fmean = sum(sf2(1:ix))/float(ix-nmis)

            do i=1,ix
                if (sf(i,j).lt.fmis) sf2(i) = fmean
            end do

            sf2(0)      = sf2(ix)
            sf2(ix+1) = sf2(1)
            do i=1,ix
                if (sf(i,j).lt.fmis) sf(i,j) = 0.5*(sf2(i-1)+sf2(i+1))
            end do
        end do
    end do
end

subroutine load_boundary_file(ioflag,iunit,fld,offset)
    ! if ioflag = 1 : read  field on from unit=iunit
    ! if ioflag = 2 : write field on from unit=iunit

    use mod_atparam

    implicit none

    integer, intent(in) :: ioflag, iunit, offset
    real     :: fld(ix,il)
    real(4) :: inp(ix,il)
    integer :: i

    open(unit=iunit, form='unformatted', access='direct', recl=ix*4, convert='little_endian')
    if (ioflag <= 1) then
        do i = 1, il
            read(iunit,rec=offset*il+i) inp(:,il+1-i)
        end do

        fld = inp

        ! Fix undefined values
        where (fld <= -999) fld = 0.0
    else
        inp = fld
        do i = il*offset+1, il*offset+il
            write(iunit,rec=i) inp(:,i)
        end do
    endif
    close(unit=iunit)
end
