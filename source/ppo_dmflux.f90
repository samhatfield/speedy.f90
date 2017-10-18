subroutine dmflux(iadd)
    ! subroutine dmflux (iadd)
    !
    ! Purpose: Add up fluxes to provide daily averages 
    !          used in sea/land models and daily/time-mean output
    ! Input: IADD = 0 to initialize storage arrays to 0
    !             > 0 to increment arrays with current flux values  

    use mod_tsteps, only: nsteps, istart
    use mod_atparam
    use mod_tmean, only: save2d_2, save2d_d2
    use mod_flx_land
    use mod_flx_sea
    use mod_physcon, only: alhc, sbc
    use mod_surfcon, only: fmask, fmask1
    use mod_var_sea, only: tice_am, sice_am
    use mod_physvar
    use mod_radcon, only: albsea, albice, emisfc
    use mod_date, only: ihour

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer, intent(in) :: iadd
    integer :: j

    real :: prec(ngp), difice(ngp)

    real :: fland(ngp), esbc, rstep1, rstep2, rsteps, sstfr, sstfr4

    fland = reshape(fmask1,(/ngp/))

    ! 1. Initialization
    if (iadd.le.0) then
        if (ihour /= 0 .and. istart /= 0) then
            ! Read from flux file
            open (100,file='fluxes.grd',form='unformatted',access='direct',recl=8*ngp)
            read (100,rec=1) (prec_l(j),j=1,ngp)
            read (100,rec=2) (snowf_l(j),j=1,ngp)
            read (100,rec=3) (evap_l(j),j=1,ngp)
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
            ! NB storage arrays for time-mean output are re-initialized 
            ! by subroutines TMOUT and DMOUT after write-up
    
            prec_l(:)  = 0.
            snowf_l(:) = 0.
            evap_l(:)  = 0.
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
    prec_l(:) = prec_l(:) + prec(:)  *rstep1
    evap_l(:) = evap_l(:) + evap(:,1)*rstep1

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

    ! Surface water budget (in mm/day)
    save2d_d2(:,1) = save2d_d2(:,1) + prec(:)  *86.400
    save2d_d2(:,2) = save2d_d2(:,2) + evap(:,3)*86.400

    ! Surface momentum budget 
    save2d_d2(:,3) = save2d_d2(:,3) - ustr(:,3)
    save2d_d2(:,4) = save2d_d2(:,4) - vstr(:,3)

    ! OLR
    save2d_d2(:,5) = save2d_d2(:,5) + olr(:)

    ! Surface energy budget
    save2d_d2(:,6) = save2d_d2(:,6) + hfluxn(:,1)
    save2d_d2(:,7) = save2d_d2(:,7) + hfluxn(:,2)

    ! 4.2 Store fluxes for time-mean output
    ! Surface water budget (in mm/day)
    save2d_2(:,1) = save2d_2(:,1) + precls(:)*86.400
    save2d_2(:,2) = save2d_2(:,2) + precnv(:)*86.400
    save2d_2(:,3) = save2d_2(:,3) + evap(:,3)*86.400

    ! Surface momentum budget 
    save2d_2(:,4) = save2d_2(:,4) - ustr(:,3)
    save2d_2(:,5) = save2d_2(:,5) - vstr(:,3)

    ! Top-of-atmosphere energy budget
    save2d_2(:,6) = save2d_2(:,6) + tsr(:)
    save2d_2(:,7) = save2d_2(:,7) + olr(:)

    ! Surface energy budget
    save2d_2(:,8)  = save2d_2(:,8)  + ssr(:)
    save2d_2(:,9)  = save2d_2(:,9)  + slr(:)
    save2d_2(:,10) = save2d_2(:,10) + shf(:,3)
    save2d_2(:,11) = save2d_2(:,11) + hfluxn(:,1)
    save2d_2(:,12) = save2d_2(:,12) + hfluxn(:,2)

    ! End of flux increment
end
