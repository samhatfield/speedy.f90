subroutine ini_sea(istart)
    ! subroutine ini_sea(istart)

    ! Input : istart = restart flag ( 0 = no, 1 = yes)

    use mod_cpl_flags, only: icsea
    use mod_atparam
    use mod_cli_sea, only: deglat_s
    use mod_var_sea

    implicit none

    integer, intent(in) :: istart

    ! 1. Compute climatological fields for initial date
    call atm2sea(0)

    ! 2. Initialize prognostic variables of ocean/ice model
    !    in case of no restart or no coupling
    sst_om(:)  = sstcl_ob(:)      ! SST 
    tice_om(:) = ticecl_ob(:)     ! sea ice temperature
    sice_om(:) = sicecl_ob(:)     ! sea ice fraction

    if (icsea.le.0) sst_om(:) = 0.

    ! 3. Compute additional sea/ice variables
    wsst_ob(:) = 0.
    if (icsea.ge.4) call sea_domain('elnino',deglat_s,wsst_ob)

    call sea2atm(0)
end

subroutine atm2sea(jday)
    ! subroutine atm2sea(jday)

    use mod_cpl_flags, only: icsea, icice, isstan
    use mod_atparam
    use mod_cplvar_sea, only: vsea_input
    use mod_date, only: iday, imont1, tmonth
    use mod_flx_sea, only: hflux_s, hflux_i
    use mod_cli_sea, only: fmask_s, sst12, sice12, sstan3, hfseacl, sstom12
    use mod_var_sea, only: sstcl_ob, sicecl_ob, ticecl_ob, sstan_ob, sstcl_om,&
        & sst_om, tice_om

    implicit none

    integer, intent(in) :: jday
    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    real :: fmasks(ngp)                  ! sea fraction
    real :: hfyearm(ngp)                 ! annual mean heat flux into the ocean
    integer :: j
    real :: sstcl0, sstfr

    ! 1. Interpolate climatological fields and obs. SST anomaly
    !    to actual date

    ! Climatological SST
    call forin5(ngp,imont1,tmonth,sst12,sstcl_ob)

    ! Climatological sea ice fraction
    call forint(ngp,imont1,tmonth,sice12,sicecl_ob)

    ! SST anomaly
    if (isstan.gt.0) then 
        if (iday.eq.1.and.jday.gt.0) call OBS_SSTA
        call forint (ngp,2,tmonth,sstan3,sstan_ob)
    end if

    ! Ocean model climatological SST
    if (icsea.ge.3) then
        call forin5 (ngp,imont1,tmonth,sstom12,sstcl_om)
    end if

    ! Adjust climatological fields over sea ice

    ! SST at freezing point
    sstfr = 273.2-1.8

    do j=1,ngp
        sstcl0 = sstcl_ob(j)

        if (sstcl_ob(j).gt.sstfr) then
            sicecl_ob(j) = min(0.5,sicecl_ob(j))
            ticecl_ob(j) = sstfr
            if (sicecl_ob(j).gt.0.) then
                sstcl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/(1.-sicecl_ob(j))
            end if
        else
            sicecl_ob(j) = max(0.5,sicecl_ob(j))
            ticecl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/sicecl_ob(j)
            !ticecl_ob(j) = sstcl_ob(j)
            sstcl_ob(j)  = sstfr
        end if

        if (icsea.ge.3) sstcl_om(j) = sstcl_om(j)+(sstcl_ob(j)-sstcl0)
    end do

    hfyearm = reshape(hfseacl, (/ngp/))
    fmasks = reshape(fmask_s, (/ngp/))

    if (jday.le.0) return
        ! 2. Set input variables for mixed-layer/ocean model
        if (icsea.gt.0.or.icice.gt.0) then
            vsea_input(:,1) = sst_om(:)
            vsea_input(:,2) = tice_om(:)
            vsea_input(:,3) = sicecl_ob(:)
            !vsea_input(:,4) = hflux_s(:)*fmasks(:)
            !vsea_input(:,5) = hflux_i(:)*fmasks(:)
            vsea_input(:,4) = hflux_s(:)
            vsea_input(:,5) = hflux_i(:)
            vsea_input(:,6) = sstcl_ob(:)
            vsea_input(:,7) = ticecl_ob(:)
            !vsea_input(:,8) = hfyearm(:)*fmasks(:)
            vsea_input(:,8) = hfyearm(:)
        end if

        ! 3. Call message-passing routines to send data (if needed)
end

subroutine sea2atm(jday)
    ! subroutine sea2atm(jday)

    use mod_cpl_flags, only: icsea, icice, isstan
    use mod_atparam
    use mod_cplvar_sea, only: vsea_output
    use mod_var_sea

    implicit none

    integer, intent(in) :: jday

    if (jday.gt.0.and.(icsea.gt.0.or.icice.gt.0)) then
        ! 1. Run ocean mixed layer or 
        !    call message-passing routines to receive data from ocean model
        call sea_model 

        ! 2. Get updated variables for mixed-layer/ocean model
        sst_om(:)   = vsea_output(:,1)      ! sst
        tice_om(:)  = vsea_output(:,2)      ! sea ice temperature 
        sice_om(:)  = vsea_output(:,3)      ! sea ice fraction

        !sice_om(:)  = sicecl_ob(:)
    end if

    ! 3. Compute sea-sfc. anomalies and full fields for atm. model
    ! 3.1 SST
    sstan_am(:) = 0.

    if (icsea.le.1) then
        if (isstan.gt.0) sstan_am(:) = sstan_ob(:)

        ! Use observed SST (climatological or full field)
        sst_am(:) = sstcl_ob(:) + sstan_am(:)
    else if (icsea.eq.2) then
        ! Use full ocean model SST
        sst_am(:) = sst_om(:)
    else if (icsea.ge.3) then
        ! Define SST anomaly from ocean model ouput and climatology 
        sstan_am(:) = sst_om(:) - sstcl_om(:)

        ! Merge with observed SST anomaly in selected area
        if (icsea.ge.4) then
            sstan_am(:) = sstan_am(:) + wsst_ob(:)*(sstan_ob(:)-sstan_am(:))
        end if

        ! Add observed SST climatology to model SST anomaly 
        sst_am(:) = sstcl_ob(:) + sstan_am(:)
    end if

    ! 3.2 Sea ice fraction and temperature
    if (icice.gt.0) then
        sice_am(:) = sice_om(:)
        tice_am(:) = tice_om(:)
    else
        sice_am(:) = sicecl_ob(:)
        tice_am(:) = ticecl_ob(:)
    end if

    sst_am(:)  = sst_am(:)+sice_am(:)*(tice_am(:)-sst_am(:))
    ssti_om(:) = sst_om(:)+sice_am(:)*(tice_am(:)-sst_om(:))
end

subroutine rest_sea(imode)
    ! subroutine rest_sea(imode)

    ! Purpose : read/write sea variables from/to a restart file
    ! Input :   IMODE = 0 : read model variables from a restart file
    !                 = 1 : write model variables  to a restart file

    use mod_cpl_flags, only: icsea, icice
    use mod_atparam
    use mod_var_sea, only: sst_om, tice_om, sice_om, sst_am, tice_am, sice_am

    implicit none

    integer, intent(in) :: imode
    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    real :: sst_c(ngp)              ! sst corrected for sea-ice values
    real :: sstfr

    if (imode.eq.0) then
        read (3)  sst_om(:)       ! sst 
        read (3) tice_om(:)       ! sea ice temperature
        read (3) sice_om(:)       ! sea ice fraction
    else
        !    write sea/ice model variables from coupled runs,
        !    otherwise write fields used by atmospheric model
        sstfr = 273.2-1.8

        if (icsea.gt.0) then
            write (10) sst_om(:) 
        else
            sst_c(:) = max(sst_am(:),sstfr)
            write (10) sst_c(:)
        end if

        if (icice.gt.0) then
            write (10) tice_om(:) 
            write (10) sice_om(:) 
        else
            write (10) tice_am(:)
            write (10) sice_am(:) 
        end if
    end if
end

subroutine obs_ssta 
    ! subroutine obs_ssta 

    ! Purpose : update observed SST anomaly array

    use mod_atparam
    use mod_cli_sea, only: sstan3, bmask_s
    use mod_date, only: imonth
    use mod_tsteps, only: iyear0, issty0

    implicit none
 
    integer :: i, j, next_month
    integer, parameter :: nlon = ix, nlat = il, ngp = ix*il
    real   :: inp(nlon,nlat)

    sstan3(:,:,1) = sstan3(:,:,2)
    sstan3(:,:,2) = sstan3(:,:,3)

    ! Compute next month given initial SST year
    next_month = (iyear0 - issty0) * 12 + imonth

    ! Read next month SST anomalies
    call load_boundary_file(1,30,inp,next_month-1)

    sstan3(1:nlon,1:nlat,3)   = inp

    call forchk(bmask_s,sstan3(1,1,3),nlon*nlat,1,-50.,50.,0.)

 100  continue

    print *, ' warning: end-of-file reached on ssT anomaly file'
    print *, ' sst anomaly will be kept constant'
end  
