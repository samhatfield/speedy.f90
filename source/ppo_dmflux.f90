! Add up fluxes to provide daily averages used in sea/land models
subroutine dmflux(iadd)
    ! Input: iadd = 0 to initialize storage arrays to 0
    !             > 0 to increment arrays with current flux values

    use params
    use physical_constants, only: alhc, sbc
    use land_model, only: fmask_l, hflux_l
    use sea_model, only: tice_am, sice_am, hflux_s, hflux_i
    use physics, only: hfluxn, shf, evap, ssrd
    use mod_radcon, only: albsea, albice, emisfc

    implicit none

    integer, intent(in) :: iadd

    real :: difice(ix,il)

    real :: esbc, rsteps, sstfr, sstfr4

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
    hflux_l = hflux_l + hfluxn(:,:,1)*rsteps

    ! Difference in net (downw.) heat flux between ice and sea surface
    difice = (albsea - albice)*ssrd + esbc*(sstfr4 - tice_am**4) + shf(:,:,2) + evap(:,:,2)*alhc

    hflux_s = hflux_s + rsteps*hfluxn(:,:,2)
    hflux_i = hflux_i + rsteps*(hfluxn(:,:,2) + difice*(1.0 - sice_am))

    ! Store fluxes for daily-mean output

    ! Multiply net heat fluxes by land or sea fractions
    hfluxn(:,:,1) = hfluxn(:,:,1)*fmask_l
    hfluxn(:,:,2) = hfluxn(:,:,2)*(1.0 - fmask_l)
end
