! Read topography and climatological boundary conditions
subroutine inbcon
    use mod_tsteps, only: isst0
    use mod_atparam
    use mod_surfcon
    use land_model, only: fmask_l, bmask_l, stl12, snowd12, soilw12
    use sea_model, only: fmask_s, bmask_s, deglat_s, sst12, sice12, sstan3, hfseacl, sstom12, &
        & sea_coupling_flag, sst_anomaly_coupling_flag
    use mod_dyncon1, only: grav, radang
    use input_output, only: load_boundary_file, load_boundary_file_old

    implicit none

    real, dimension(ix,il) :: veg_low, veg_high, veg, swl1, swl2
    integer :: i, idep2, month, j, jrec
    real :: rsw, sdep1, sdep2, swroot, swwil2, thrsh

    ! Set threshold for land-sea mask definition
    ! (ie minimum fraction of either land or sea)
    thrsh = 0.1

    ! Read surface geopotential (i.e. orography)
    phi0 = grav*load_boundary_file("surface.nc", "orog")

    ! Also store spectrally truncated surface geopotential
    call truncg (ntrun,phi0,phis0)

    ! Read land-sea mask
    fmask = load_boundary_file("surface.nc", "lsm")

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
    alb0 = load_boundary_file("surface.nc", "alb")

    ! Land-surface temp.
    do month = 1,12
        stl12(:,:,month) = load_boundary_file("land.nc", "stl", month)

        call fillsf(stl12(:,:,month), ix, il, 0.0)
    end do

    call forchk(bmask_l, 12, 0.0, 400.0, 273.0, stl12)

    ! Snow depth
    do month = 1,12
        snowd12(:,:,month) = load_boundary_file("snow.nc", "snowd", month)
    end do

    call forchk(bmask_l, 12, 0.0, 20000.0, 0.0, snowd12)

    ! Read soil moisture and compute soil water availability using vegetation fraction
    ! Read vegetation fraction
    veg_high = load_boundary_file("surface.nc", "vegh")
    veg_low  = load_boundary_file("surface.nc", "vegl")

    ! Combine high and low vegetation fractions
    veg = max(0.,veg_high+0.8*veg_low)

    ! Read soil moisture
    sdep1 = 70.
    idep2 = 3
    sdep2 = idep2*sdep1

    swwil2 = idep2*swwil
    rsw    = 1.0/(swcap+idep2*(swcap-swwil))

    do month = 1,12
        ! Combine soil water content from two top layers
        swl1 = load_boundary_file("soil.nc", "swl1", month)
        swl2 = load_boundary_file("soil.nc", "swl2", month)

        do j = 1,il
            do i = 1,ix
                swroot = idep2*swl2(i,j)
                soilw12(i,j,month) = min(1.0, rsw*(swl1(i,j) + veg(i,j)*max(0.0, swroot - swwil2)))
            end do
        end do
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
    deglat_s =  radang*90.0/asin(1.0)

    ! SST
    do month = 1,12
        sst12(:,:,month) = load_boundary_file("sea_surface_temperature.nc", "sst", month)

        call fillsf(sst12(:,:,month), ix, il, 0.0)
    end do

    call forchk(bmask_s, 12, 100.0, 400.0, 273.0, sst12)

    ! Sea ice concentration
    do month = 1,12
        sice12(:,:,month) = max(load_boundary_file("sea_ice.nc", "icec", month), 0.0)
    end do

    call forchk(bmask_s, 12, 0.0, 1.0, 0.0, sice12)

    ! SST anomalies for initial and prec./following months
    if (sst_anomaly_coupling_flag > 0) then
        print *, 'isst0 = ', isst0
        do month=1,3
            if ((isst0 <= 1 .and. month /= 2) .or. isst0 > 1) then
                sstan3(:,:,month) = load_boundary_file("sea_surface_temperature_anomaly.nc", &
                    & "ssta", isst0-2+month-1, 420)
            end if
        end do

        call forchk(bmask_s, 3, -50.0, 50.0, 0.0, sstan3)
    end if

    ! Climatological fields for the ocean model (TO BE RECODED)
    ! Annual-mean heat flux into sea-surface`
    hfseacl = 0.0

    if (sea_coupling_flag >= 1) then
        stop "Model behaviour when sea_coupling_flag >= 1 not implemented yet"
    end if

    ! Ocean model SST climatology:
    ! defined by adding SST model bias to obs. climatology
    ! (bias may be defined in a different period from climatology)

    if (sea_coupling_flag >= 3) then
        stop "Model behaviour when sea_coupling_flag >= 3 not implemented yet"
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
                if (fmask(i,j) > 0.0) then
                    if (field(i,j,jf) < fmin .or. field(i,j,jf) > fmax) then
                        nfault = nfault + 1
                    end if
                else
                    field(i,j,jf) = fset
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
