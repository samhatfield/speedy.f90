subroutine ini_land(istart)
    ! subroutine ini_land (istart)
    !
    ! Input : istart = restart flag ( 0 = no, 1 = yes)

    use mod_atparam
    use mod_var_land, only: stlcl_ob, stl_lm

    implicit none

    integer, intent(in) :: istart

    ! 1. Compute climatological fields for initial date
    call atm2land(0)

    ! 2. Initialize prognostic variables of land model
    !    in case of no restart or no coupling
    if (istart.le.0 .or. istart == 2) then
        stl_lm(:)  = stlcl_ob(:)      ! land sfc. temperature 
    end if

    ! 3. Compute additional land variables
    call land2atm(0)
end

subroutine atm2land(jday)
    use mod_cpl_flags, only: icland
    use mod_atparam
    use mod_cpl_land_model, only: vland_input
    use mod_flx_land, only: hflux_l
    use mod_cli_land, only: stl12, snowd12, soilw12
    use mod_date, only: imont1, tmonth
    use mod_var_land, only: stlcl_ob, snowdcl_ob, soilwcl_ob, stl_lm

    implicit none

    integer, intent(in) :: jday
    integer, parameter :: nlon=ix, nlat=il, ngp=nlon*nlat

    ! 1. Interpolate climatological fields to actual date

    ! Climatological land sfc. temperature
    call forin5(ngp,imont1,tmonth,stl12,stlcl_ob)

    ! Climatological snow depth
    call forint(ngp,imont1,tmonth,snowd12,snowdcl_ob)

    ! Climatological soil water availability
    call forint(ngp,imont1,tmonth,soilw12,soilwcl_ob)

    if (jday.le.0) return

    ! 2. Set input variables for mixed-layer/ocean model
    if (icland.gt.0) then
        vland_input(:,1) = stl_lm(:)
        vland_input(:,2) = hflux_l(:)
        vland_input(:,3) = stlcl_ob(:)
    end if

    ! 3. Call message-passing routines to send data (if needed)
end

subroutine land2atm(jday)
    use mod_cpl_flags, only: icland
    use mod_atparam
    use mod_cpl_land_model, only: land_model, vland_output
    use mod_var_land

    implicit none

    integer, intent(in) :: jday

    if (jday.gt.0.and.icland.gt.0) then
        ! 1. Run ocean mixed layer or 
        !    call message-passing routines to receive data from ocean model
        call land_model 

        ! 2. Get updated variables for mixed-layer/ocean model
        stl_lm(:) = vland_output(:,1)      ! land sfc. temperature
    end if

    ! 3. Compute land-sfc. fields for atm. model
    ! 3.1 Land sfc. temperature
    if (icland.le.0) then
        ! Use observed climatological field
        stl_am(:) = stlcl_ob(:)
    else
        ! Use land model sfc. temperature
        stl_am(:) = stl_lm(:)
    end if

    ! 3.2 Snow depth and soil water availability
    snowd_am(:) = snowdcl_ob(:)
    soilw_am(:) = soilwcl_ob(:)
end

subroutine rest_land(imode)
    ! subroutine rest_land (imode)

    ! Purpose : read/write land variables from/to a restart file
    ! Input :   IMODE = 0 : read model variables from a restart file
    !                 = 1 : write model variables  to a restart file

    use mod_cpl_flags, only: icland
    use mod_atparam
    use mod_var_land, only: stl_am, stl_lm

    implicit none

    integer, intent(in) :: imode

    if (imode.eq.0) then
        read (3)  stl_lm(:)       ! Land sfc. temperature 
    else
        ! Write land model variables from coupled runs,
        ! otherwise write fields used by atmospheric model
        if (icland.gt.0) then
            write (10) stl_lm(:) 
        else
            write (10) stl_am(:)
        end if
    end if
end
