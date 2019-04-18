subroutine dmflux(iadd)
    ! subroutine dmflux (iadd)
    !
    ! Purpose: Add up fluxes to provide daily averages
    !          used in sea/land models
    ! Input: IADD = 0 to initialize storage arrays to 0
    !             > 0 to increment arrays with current flux values

    use mod_tsteps, only: nsteps
    use mod_atparam
    use mod_flx_sea
    use mod_physcon, only: alhc, sbc
    use mod_surfcon, only: fmask
    use mod_cpl_land_model, only: fmask_l, hflux_l
    use mod_cpl_sea_model, only: tice_am, sice_am
    use mod_physvar
    use mod_radcon, only: albsea, albice, emisfc
    use mod_date, only: model_datetime

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer, intent(in) :: iadd
    integer :: j

    real :: prec(ngp), difice(ngp)

    real :: fland(ngp), esbc, rstep1, rstep2, rsteps, sstfr, sstfr4

    fland = reshape(fmask_l,(/ngp/))

    ! 1. Initialization
    if (iadd.le.0) then
        if (model_datetime%hour /= 0) then
            ! Read from flux file
            open (100,file='fluxes.grd',form='unformatted',access='direct',recl=8*ngp)
            read (100,rec=4) (hflux_l(j),j=1,ngp)

            read (100,rec=5) (prec_s(j),j=1,ngp)
            read (100,rec=6) (snowf_s(j),j=1,ngp)
            read (100,rec=7) (evap_s(j),j=1,ngp)
            read (100,rec=8) (ustr_s(j),j=1,ngp)
            read (100,rec=9) (vstr_s(j),j=1,ngp)
            read (100,rec=10) (ssr_s(j),j=1,ngp)
            read (100,rec=11) (slr_s(j),j=1,ngp)
            read (100,rec=12) (shf_s(j),j=1,ngp)
            read (100,rec=13) (ehf_s(j),j=1,ngp)
            read (100,rec=14) (hflux_s(j),j=1,ngp)
            read (100,rec=15) (hflux_i(j),j=1,ngp)
            close (100)
        else
            ! Set all daily-mean arrays to zero
            hflux_l(:) = 0.

            prec_s(:)  = 0.
            snowf_s(:) = 0.
            evap_s(:)  = 0.
            ustr_s(:)  = 0.
            vstr_s(:)  = 0.
            ssr_s(:)   = 0.
            slr_s(:)   = 0.
            shf_s(:)   = 0.
            ehf_s(:)   = 0.
            hflux_s(:) = 0.
            hflux_i(:) = 0.
        end if
        return
    end if

    rsteps = 1./real(nsteps)
    rstep1 = rsteps*0.001
    rstep2 = rsteps*alhc

    ! SST at freezing point
    sstfr  = 273.2-1.8

    sstfr4 = sstfr**4
    esbc   = emisfc*sbc

    ! Total precipitation
    prec(:) = precls(:)+precnv(:)

    ! 2. Store fluxes over land (SI units, all heat fluxes downw.)
    hflux_l(:) = hflux_l(:) + hfluxn(:,1)*rsteps

    ! 3. Store fluxes over sea (SI units, all heat fluxes downw.)
    prec_s(:) = prec_s(:) + prec(:)  *rstep1
    evap_s(:) = evap_s(:) + evap(:,2)*rstep1

    ustr_s(:) = ustr_s(:) - ustr(:,2)*rsteps
    vstr_s(:) = vstr_s(:) - vstr(:,2)*rsteps

    ssr_s(:) = ssr_s(:) + ssr(:)   *rsteps
    slr_s(:) = slr_s(:) - slr(:)   *rsteps
    shf_s(:) = shf_s(:) - shf(:,2) *rsteps
    ehf_s(:) = ehf_s(:) - evap(:,2)*rstep2

    ! Difference in net (downw.) heat flux between ice and sea surface
    difice(:) = (albsea-albice)*ssrd(:)+ esbc*(sstfr4-tice_am(:)**4)&
        & + shf(:,2)+evap(:,2)*alhc

    hflux_s(:) = hflux_s(:) + rsteps* hfluxn(:,2)
    hflux_i(:) = hflux_i(:) + rsteps*(hfluxn(:,2)+difice(:)*(1.-sice_am(:)))

    ! 4.1 Store fluxes for daily-mean output

    ! Multiply net heat fluxes by land or sea fractions
    hfluxn(:,1) = hfluxn(:,1)*fland(:)
    hfluxn(:,2) = hfluxn(:,2)*(1.-fland(:))

    ! End of flux increment
end
