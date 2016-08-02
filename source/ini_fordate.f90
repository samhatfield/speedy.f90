subroutine fordate(imode)
    !
    !   subroutine fordate (imode)
    !   
    !   purpose :	compute forcing fields for the current date
    !             and correction terms for horiz. diffusion
    !
    !   input : imode : 0 = initialization step, 1 = daily update

    use lflags, only: lco2

    implicit none

    include "atparam.h"
    include "atparam1.h"

    integer, parameter :: nlon = ix, nlat = il, nlev = kx, ngp = nlon * nlat

    include "com_date.h"

    include "com_dyncon0.h"
    include "com_physcon.h"

    include "com_radcon.h"
    include "com_hdifcon.h"

    include "com_surfcon.h"

    include "com_cli_sea.h"
    include "com_cli_land.h"

    include "com_var_sea.h"
    include "com_var_land.h"

    integer, intent(in) :: imode
    real, dimension(nlon, nlat) :: corh, tsfc, tref, psfc, qsfc, qref
    real :: gamlat(nlat)

    real :: fland(ngp), alb_0(ngp)
    equivalence (fland,fmask_l), (alb_0,alb0)

    real :: del_co2, dummy, pexp
    integer :: i, j, ij, iitest = 0, iyear_ref

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

    do j = 1, ngp
        snowc(j)  = min(1., snowd_am(j)/sd2sc)
        alb_l(j)  = alb_0(j) + snowc(j) * (albsn - alb_0(j))
        alb_s(j)  = albsea + sice_am(j) * (albice - albsea)
        albsfc(j) = alb_s(j) + fland(j) * (alb_l(j) - alb_s(j))
    end do

    ! linear trend of co2 absorptivity (del_co2: rate of change per year)
    iyear_ref = 1950
    del_co2   = 0.005
    ! del_co2   = 0.0033

    if (lco2) then
        ablco2 = ablco2_ref * exp(del_co2 * (iyear + tyear - iyear_ref))
    end if

    ! 3. temperature correction term for horizontal diffusion
    call setgam(tyear,gamlat)

    do j = 1, nlat
        do i = 1, nlon
            corh(i,j) = gamlat(j) * phis0(i,j)
        end do
    end do

    if (iitest > 1.and.imode == 0) then
        call outest(19,phis0)
        call outest(19,corh)
    end if

    call spec(corh,tcorh)

!   4. humidity correction term for horizontal diffusion
    ij = 0
    do j = 1, nlat
        pexp = 1./(rd * gamlat(j))
        do i = 1, nlon
            ij = ij + 1
!            tsfc(i,j) = fmask_l(i,j)*stlcl_ob(ij)
!     &               +fmask_s(i,j)*sstcl_ob(ij)
            tsfc(i,j) = fmask_l(i,j) * stl_am(ij)&
                & + fmask_s(i,j) * sst_am(ij)
            tref(i,j) = tsfc(i,j) + corh(i,j)
            psfc(i,j) = (tsfc(i,j)/tref(i,j))**pexp
        end do
    end do

    call shtorh(0, ngp, tref,   1., -1., dummy, dummy, qref)
    call shtorh(0, ngp, tsfc, psfc,  1., dummy, dummy, qsfc)

    do j = 1, nlat
        do i = 1, nlon
            corh(i,j) = refrh1 * (qref(i,j) - qsfc(i,j))
        end do
    end do

    if (iitest > 1.and.imode == 0) call outest(19,corh)

    call spec(corh,qcorh)
end

subroutine setgam(tyear,gamlat)  
    ! aux. routine gamlat : compute reference lapse rate 
    !                       as a function of latitude and date

    implicit none

    include "atparam.h"
    include "atparam1.h"

    integer, intent(in) :: tyear
    integer, parameter :: nlon = ix, nlat = il, nlev = kx, ngp = nlon * nlat
    integer :: j
                                            
    include "com_dyncon0.h"                 
    include "com_physcon.h"

    real, intent(inout) :: gamlat(nlat)

    gamlat(1) = gamma/(1000. * gg)
    do j = 2, nlat
        gamlat(j) = gamlat(1)
    end do
end

subroutine outest(iunit,fout)
    ! aux. routine outest : write one field on a test output file 

    implicit none

    include "atparam.h"

    integer, intent(in) :: iunit
    real, intent(in) :: fout(ix, il)
    integer :: i, j

    real*4 :: r4out(ix,il)

    do j = 1, il
        do i = 1, ix
            r4out(i,j) = fout(i,j)
        end do
    end do

    write (iunit) r4out
end
