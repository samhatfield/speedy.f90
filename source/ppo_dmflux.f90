! Add up fluxes to provide daily averages used in sea/land models
subroutine dmflux(iadd)
    ! Input: iadd = 0 to initialize storage arrays to 0
    !             > 0 to increment arrays with current flux values

    use mod_tsteps, only: nsteps
    use mod_atparam
    use mod_physcon, only: alhc, sbc
    use mod_surfcon, only: fmask
    use mod_cpl_land_model, only: fmask_l, hflux_l
    use mod_cpl_sea_model, only: tice_am, sice_am, hflux_s, hflux_i
    use mod_physvar
    use mod_radcon, only: albsea, albice, emisfc
    use mod_date, only: model_datetime

    implicit none

    integer, parameter :: ngp=ix*il

    integer, intent(in) :: iadd
    integer :: j

    real :: difice(ix,il)

    real :: fland(ngp), esbc, rsteps, sstfr, sstfr4

    fland = reshape(fmask_l,(/ngp/))

    ! 1. Initialization
    if (iadd <= 0) then
        ! Set all daily-mean arrays to zero
        hflux_l = 0.0
        hflux_s = 0.0
        hflux_i = 0.0
        return
    end if

    rsteps = 1./real(nsteps)

    ! SST at freezing point
    sstfr  = 273.2-1.8

    sstfr4 = sstfr**4
    esbc   = emisfc*sbc

    ! Store fluxes over land (SI units, all heat fluxes downw.)
    hflux_l(:) = hflux_l(:) + hfluxn(:,1)*rsteps

    ! Difference in net (downw.) heat flux between ice and sea surface
    difice = (albsea - albice)*reshape(ssrd, (/ ix,il /)) + esbc*(sstfr4 - tice_am**4)&
        & + reshape(shf(:,2), (/ ix, il /)) + reshape(evap(:,2), (/ ix, il /))*alhc

    hflux_s = hflux_s + rsteps*reshape(hfluxn(:,2), (/ ix,il /))
    hflux_i = hflux_i + rsteps*(reshape(hfluxn(:,2), (/ ix, il /)) + difice*(1.0 - sice_am))

    ! Store fluxes for daily-mean output

    ! Multiply net heat fluxes by land or sea fractions
    hfluxn(:,1) = hfluxn(:,1)*fland(:)
    hfluxn(:,2) = hfluxn(:,2)*(1.-fland(:))
end
