subroutine inbcon(grav0,radlat)
    !
    ! subroutine inbcon (grav0,radlat)
    !
    ! Purpose : Read topography and climatological boundary conditions 
    ! Input :   grav0  = gravity accel.
    !           radlat = grid latitudes in radiants

    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps, only: isst0
    use mod_atparam
    use mod_surfcon
    use mod_cli_land
    use mod_cli_sea
 	  							
    implicit none

    real, intent(in) :: grav0, radlat(il)
    integer, parameter :: nlon = ix, nlat = il, ngp = ix*il

    real*4 :: r4inp(nlon,nlat), dummy4
    real   :: inp(nlon,nlat), phis1(nlon,nlat)
    real   :: veg(nlon,nlat), swl1(nlon,nlat), swl2(nlon,nlat)

    integer :: iitest=1, i, idep2, irec, irecl, it, j, jrec
    real :: rad2deg, rsw, sdep1, sdep2, swroot, swwil2, thrsh

    ! Set threshold for land-sea mask definition
    ! (ie minimum fraction of either land or sea)

    thrsh = 0.1

    ! 1. Read topographical fields (orography, land-sea mask)
    if (iitest >= 1) print *,' read orography' 

    call load_boundary_file(1,20,inp,0)

    phi0 = grav0*inp

    call truncg (ntrun,phi0,phis0)

    if (iitest >= 1) print *,' read fractional land-sea mask'  

    call load_boundary_file(1,20,inp,1)

    fmask = inp

    ! 2. Initialize land-sfc boundary conditions

    ! 2.1 Fractional and binary land masks
    do j=1,il
        do i=1,ix
            fmask_l(i,j) = fmask(i,j)

            if (fmask_l(i,j).ge.thrsh) then
              bmask_l(i,j) = 1.
              if (fmask(i,j).gt.(1.-thrsh)) fmask_l(i,j) = 1.
            else
              bmask_l(i,j) = 0.
              fmask_l(i,j) = 0.
            end if

            fmask1(i,j) = fmask_l(i,j)
        end do
    end do

    ! 2.2 Annual-mean surface albedo
    if (iitest >= 1) print *,' read surface albedo' 
 
    call load_boundary_file(1,20,inp,2)

    alb0 = inp

    ! 2.3 Land-surface temp.
    if (iitest >= 1) print *,' reading land-surface temp.'

    do it = 1,12
        call load_boundary_file(1,23,inp,it-1)

        call fillsf(inp,nlon,nlat,0.)

       stl12(1:nlon,1:nlat,it) = inp
    end do

    if (iitest == 1) print *,' checking land-surface temp.'

    call forchk(bmask_l,stl12,ngp,12,0.,400.,273.)

    ! 2.4 Snow depth
    if (iitest >= 1) print *,' reading snow depth'  

    do it = 1,12
        call load_boundary_file(1,24,inp,it-1)
  
        snowd12(1:nlon,1:nlat,it) = inp
    end do

    if (iitest >= 1) print *,' checking snow depth'

    CALL FORCHK (bmask_l,snowd12,ngp,12,0.,20000.,0.)
       
    ! 2.5 Read soil moisture and compute soil water availability 
    !     using vegetation fraction

    if (iitest >= 1) print *,' reading soil moisture'  

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
        do j = 1,nlat
            do i = 1,nlon
                swroot = idep2*swl2(i,j)
                inp(i,j) = min(1.,rsw*(swl1(i,j)+veg(i,j)*max(0.,swroot-swwil2)))		
            end do
        end do

        soilw12(1:nlon,1:nlat,it) = inp
    end do

    if (iitest >= 1) print *,' checking soil moisture'

    call forchk(bmask_l,soilw12,ngp,12,0.,10.,0.)

    ! 3. Initialize sea-sfc boundary conditions
    
    ! 3.1 Fractional and binary sea masks
    do j=1,il
        do i=1,ix
            fmask_s(i,j) = 1.-fmask(i,j)

            if (fmask_s(i,j).ge.thrsh) then
                bmask_s(i,j) = 1.
                if (fmask_s(i,j).gt.(1.-thrsh)) fmask_s(i,j) = 1.
            else
                bmask_s(i,j) = 0.
                fmask_s(i,j) = 0.
            end if
        end do
    end do

    ! Grid latitudes for sea-sfc. variables
    rad2deg = 90.0/asin(1.)
    deglat_s = rad2deg*radlat
         
    ! 3.2 SST 
    if (iitest >= 1) print *,' reading sst' 

    do it = 1,12
        call load_boundary_file(1,21,inp,it-1)

        call fillsf(inp,nlon,nlat,0.)

        sst12(1:nlon,1:nlat,it) = inp
    end do

    if (iitest >= 1) print *,' checking sst'

    call forchk(bmask_s,sst12,ngp,12,100.,400.,273.)

    ! 3.3 Sea ice concentration
    if (iitest >= 1) print *,' reading sea ice'  

    do it = 1,12
        call load_boundary_file(1,22,inp,it-1)

        inp = max(inp,0.)

        sice12(1:nlon,1:nlat,it) = inp
    end do

    if (iitest >= 1) print *,' checking sea ice'

    call forchk(bmask_s,sice12,ngp,12,0.,1.,0.)

    ! 3.4 SST anomalies for initial and prec./following months
    if (isstan > 0) then
        if (iitest >= 1) print *,' reading sst anomalies' 

        print *, 'isst0 = ', isst0
        do it=1,3
            if ((isst0 <= 1 .and. it /= 2) .or. isst0 > 1) then
                call load_boundary_file(1,30,inp,isst0-2+it-1)
            end if

            sstan3(1:nlon,1:nlat,it) = inp
        end do

        if (iitest >= 1) print *,' checking sst anomalies'

        call forchk(bmask_s,sstan3,ngp,3,-50.,50.,0.)
    end if

    ! 4. Climatological fields for the ocean model (TO BE RECODED)
    ! 4.1. Annual-mean heat flux into sea-surface

    hfseacl = 0.0

    if (icsea >= 1) then
        if (iitest >= 1) print *,' reading sfc heat fluxes' 

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

        if (iitest >= 1) print *,' checking sfc heat fluxes'

        call forchk (bmask_s,hfseacl,ix*il,1,-1000.,1000.,0.)
    end if

    ! 4.2. Ocean model SST climatology:
    !      defined by adding SST model bias to obs. climatology
    !      (bias may be defined in a different period from climatology)

    if (icsea >= 3) then
        if (iitest >= 1) print *,' reading ocean model SST bias' 

        !irecl = 4*ix*il
        !irec = 0

        !open ( unit=32, file='fort.32', status='old',&
        !   & form='unformatted', access='direct', recl=irecl )

        do it = 1,12
            ! irec=irec+1
            ! read (32,rec=irec) r4inp
            read (32) r4inp

            do j = 1,il
                do i = 1,ix
                    sstom12(i,j,it) = sst12(i,j,it)+r4inp(i,j)
                end do
            end do  
        end do

        if (iitest >= 1) print *,' checking ocean model SST'

        call forchk (bmask_s,sstom12,ix*il,12,100.,400.,273.)
    end if
end

subroutine forchk (fmask,field,ngp,nf,fmin,fmax,fset)
    ! Aux. routine forchk: Check consistency of sfc fields with land-sea mask 
    ! and set undefined values to a constant (to avoid over/underflow)

    implicit none

    real, intent(in) :: fmask(ngp)
    real, intent(inout) :: field(ngp,nf)
    integer, intent(in) :: ngp, nf
    real, intent(in) :: fmin, fmax, fset

    integer :: jf, jgp, nfault

    do jf = 1,nf
        nfault=0

        do jgp = 1,ngp
            if (fmask(jgp).gt.0.0) then
                if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax) then
                    nfault = nfault+1
                end if
            else
                field(jgp,jf) = fset
            end if
        end do

        print *, ' field: ', jf, '   no. of faulty points:', nfault
    end do

    print *, ' undefined values set to', fset
end

subroutine ftland (stl,phi0,phis0,fmaskl)
    use mod_dyncon0, only: gamma
    use mod_dyncon1, only: gcos, grav
    use mod_atparam

    implicit none

    integer, parameter :: nlon = ix, nlat = il

    real, dimension(nlon, nlat), intent(inout) :: stl, phi0, phis0, fmaskl
    real :: stl2(nlon,nlat), sumt, sumw
    integer :: nl8, nlat1, nlat2, i, idtr, itr, j, jband, jfil
    real :: gam

    nl8 = nlat/8
    gam = 0.001*gamma/grav

    nlat1 = 1
    nlat2 = nl8

    do jband=1,8
        sumt=0.
        sumw=0.

        do j=nlat1,nlat2
            do i=1,nlon
                stl(i,j)=stl(i,j)+gam*phi0(i,j)
                sumt=sumt+gcos(j)*fmaskl(i,j)*stl(i,j)
                sumw=sumw+gcos(j)*fmaskl(i,j)
            end do
        end do

        SUMT=SUMT/SUMW

        do j=nlat1,nlat2
            do i=1,nlon
                if (fmaskl(i,j).eq.0.) stl(i,j)=sumt
            end do
        end do
  
        nlat1=nlat1+nl8
        nlat2=nlat2+nl8
    end do

    itr=7
    idtr=(ntrun-6)/3

    do jfil=1,4
        call truncg (itr,stl,stl2)

        do j=1,nlat
            do i=1,nlon
                if (fmaskl(i,j).eq.0.) stl(i,j)=stl2(i,j)
            end do
        end do

        itr=min(itr+idtr,ntrun)
    end do

    call truncg (itr,stl,stl2)

    stl = stl2 - gam * phis0
end

subroutine truncg (itr,fg1,fg2)
    ! subroutine truncg (itr,fg1,fg2)
    ! Purpose : compute a spectrally-filtered grid-point field
    ! Input   : itr : spectral truncation (triangular)
    !         : fg1 : original grid-point field
    ! Output  : fg2 : filtered grid-point field

    USE mod_atparam

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

subroutine fillsf(sf,nlon,nlat,fmis)
    ! subroutine fillsf (sf,nlon,nlat)
    ! Purpose: replace missing values in surface fields
    ! NB: it is assumed that non-missing values exist near the Equator

    implicit none

    real :: sf(nlon,nlat), sf2(0:nlon+1)
    integer, intent(in) :: nlon, nlat
    real, intent(in) :: fmis

    integer :: khem, j, j1, j2, j3, i, nmis
    real :: fmean

    do khem = 1,2
       if (khem == 1) then
            j1 = nlat/2
            j2 = 1
            j3 = -1
        else
            j1 = j1+1
            j2 = nlat
            j3 = 1
        end if

        do j=j1,j2,j3
            sf2(1:nlon) = sf(1:nlon,j)

            nmis = 0
            do i=1,nlon
                if (sf(i,j) < fmis) then
                    nmis = nmis+1
                    sf2(i) = 0.
                end if
            end do

            if (nmis < nlon) fmean = sum(sf2(1:nlon))/float(nlon-nmis) 

            do i=1,nlon
                if (sf(i,j).lt.fmis) sf2(i) = fmean
            end do

            sf2(0)      = sf2(nlon)
            sf2(nlon+1) = sf2(1)
            do i=1,nlon
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

    integer, parameter :: nlon = ix, nlat = il, ngp = ix*il
    integer, intent(in) :: ioflag, iunit, offset
    real     :: fld(nlon,nlat)
    real(4) :: inp(nlon,nlat)
    integer :: i

    open(unit=iunit, form='unformatted', access='direct', recl=nlon*4, convert='little_endian')
    if (ioflag <= 1) then
        do i = 1, nlat
            read(iunit,rec=offset*nlat+i) inp(:,nlat+1-i)
        end do

        fld = inp

        ! Fix undefined values
        where (fld <= -999) fld = 0.0
    else
        inp = fld
        do i = nlat*offset+1, nlat*offset+nlat
            write(iunit,rec=i) inp(:,i)
        end do
    endif
    close(unit=iunit)
end
