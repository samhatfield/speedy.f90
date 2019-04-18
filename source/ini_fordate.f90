subroutine fordate(imode)
    !
    !   subroutine fordate (imode)
    !
    !   purpose :	compute forcing fields for the current date
    !             and correction terms for horiz. diffusion
    !
    !   input : imode : 0 = initialization step, 1 = daily update

    use mod_lflags, only: lco2
    use mod_dyncon0, only: refrh1
    use mod_atparam
    use mod_hdifcon, only: tcorh, qcorh
    use mod_physcon, only: rd
    use mod_surfcon, only: phis0, alb0, sd2sc
    use mod_date, only: model_datetime, tyear
    use mod_cpl_land_model, only: stl_am, snowd_am, fmask_l
    use mod_cpl_sea_model, only: fmask_s, sstcl_ob, sst_am, sice_am
    use mod_radcon, only: ablco2, ablco2_ref, albsea, albice, snowc, albsn,&
        & alb_l, alb_s, albsfc

    implicit none

    integer, parameter :: nlon = ix, nlat = il, nlev = kx, ngp = nlon * nlat

    integer, intent(in) :: imode
    real, dimension(nlon, nlat) :: corh, tsfc, tref, psfc, qsfc, qref
    real :: gamlat(nlat)

    real :: del_co2, dummy, pexp
    integer :: i, j, iyear_ref

    ! time variables for interpolation are set by newdate

    ! 1. time-independent parts of physical parametrizations
    if (imode == 0) then
        call radset
        call sflset(phis0)

        ablco2_ref = ablco2
    end if

    ! 2. daily-mean radiative forcing
    ! incoming solar radiation
    call sol_oz(tyear)

    ! total surface albedo

    do i = 1, ix
        do j = 1, il
            snowc(i,j)  = min(1.0, snowd_am(i,j)/sd2sc)
            alb_l(i,j)  = alb0(i,j) + snowc(i,j) * (albsn - alb0(i,j))
            alb_s(i,j)  = albsea + sice_am(i,j) * (albice - albsea)
            albsfc(i,j) = alb_s(i,j) + fmask_l(i,j) * (alb_l(i,j) - alb_s(i,j))
        end do
    end do

    ! linear trend of co2 absorptivity (del_co2: rate of change per year)
    iyear_ref = 1950
    del_co2   = 0.005
    ! del_co2   = 0.0033

    if (lco2) then
        ablco2 = ablco2_ref * exp(del_co2 * (model_datetime%year + tyear - iyear_ref))
    end if

    ! 3. temperature correction term for horizontal diffusion
    call setgam(tyear,gamlat)

    do j = 1, nlat
        do i = 1, nlon
            corh(i,j) = gamlat(j) * phis0(i,j)
        end do
    end do

    call spec(corh,tcorh)

    ! 4. humidity correction term for horizontal diffusion
    do j = 1, nlat
        pexp = 1./(rd * gamlat(j))
        do i = 1, nlon
            tsfc(i,j) = fmask_l(i,j) * stl_am(i,j) + fmask_s(i,j) * sst_am(i,j)
            tref(i,j) = tsfc(i,j) + corh(i,j)
            psfc(i,j) = (tsfc(i,j)/tref(i,j))**pexp
        end do
    end do

    call shtorh(0, ngp, tref,   1., -1., dummy, dummy, qref)
    call shtorh(0, ngp, tsfc, psfc,  1., dummy, dummy, qsfc)

    corh = refrh1 * (qref - qsfc)

    call spec(corh,qcorh)
end

subroutine setgam(tyear,gamlat)
    ! aux. routine gamlat : compute reference lapse rate
    !                       as a function of latitude and date

    use mod_dyncon0, only: gamma
    use mod_atparam
    use mod_dyncon1, only: grav

    implicit none

    real, intent(in) :: tyear
    integer, parameter :: nlon = ix, nlat = il, nlev = kx, ngp = nlon * nlat
    integer :: j

    real, intent(inout) :: gamlat(nlat)

    gamlat(1) = gamma/(1000. * grav)
    do j = 2, nlat
        gamlat(j) = gamlat(1)
    end do
end
