module mod_spectral
    use mod_atparam

    implicit none

    private
    public el2, elm2, el4, trfilt, l2, ll, mm, nsh2, sia, coa, wt, wght, cosg,&
        & cosgr, cosgr2, gradx, gradym, gradyp, sqrhlf, consq, epsi, repsi,&
        & emm, ell, poly, cpol, uvdx, uvdym, uvdyp, vddym, vddyp

    ! Initial. in parmtr
    real, dimension(mx,nx) :: el2, elm2, el4, trfilt
    integer :: l2(mx,nx), ll(mx,nx), mm(mx), nsh2(nx)

    ! Initial. in parmtr
    real, dimension(iy) :: sia, coa, wt, wght
    real, dimension(il) :: cosg, cosgr, cosgr2

    ! Initial. in parmtr
    real :: gradx(mx), gradym(mx,nx), gradyp(mx,nx)

    ! Initial. in parmtr
    real :: sqrhlf, consq(mxp), epsi(mxp,nxp), repsi(mxp,nxp), emm(mxp), ell(mxp,nxp)

    ! Initial. in parmtr
    real :: poly(mx,nx)

    ! Initial. in parmtr
    real :: cpol(mx2,nx,iy)

    ! Initial. in parmtr
    real, dimension(mx,nx) :: uvdx, uvdym, uvdyp

    ! Initial. in parmtr
    real, dimension(mx,nx) :: vddym, vddyp
end module
