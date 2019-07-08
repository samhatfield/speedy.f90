!> author: Sam Hatfield, Fred Kucharski, Franco Molteni
!  date: 29/04/2019
!  For storing variables used by multiple physics schemes.
module auxiliaries
    use types, only: p
    use params

    implicit none

    private
    public precnv, precls, snowcv, snowls, cbmf, tsr, ssrd, ssr, slrd, slr, olr, slru
    public ustr, vstr, shf, evap, hfluxn

    ! Physical variables shared among all physics schemes
    real(p), dimension(ix,il)   :: precnv !! Convective precipitation  [g/(m^2 s)], total
    real(p), dimension(ix,il)   :: precls !! Large-scale precipitation [g/(m^2 s)], total
    real(p), dimension(ix,il)   :: snowcv !! Convective precipitation  [g/(m^2 s)], snow only
    real(p), dimension(ix,il)   :: snowls !! Large-scale precipitation [g/(m^2 s)], snow only
    real(p), dimension(ix,il)   :: cbmf   !! Cloud-base mass flux
    real(p), dimension(ix,il)   :: tsr    !! Top-of-atmosphere shortwave radiation (downward)
    real(p), dimension(ix,il)   :: ssrd   !! Surface shortwave radiation (downward-only)
    real(p), dimension(ix,il)   :: ssr    !! Surface shortwave radiation (net downward)
    real(p), dimension(ix,il)   :: slrd   !! Surface longwave radiation (downward-only)
    real(p), dimension(ix,il)   :: slr    !! Surface longwave radiation (net upward)
    real(p), dimension(ix,il)   :: olr    !! Outgoing longwave radiation (upward)
    real(p), dimension(ix,il,3) :: slru   !! Surface longwave emission (upward)

    ! Third dimension -> 1:land, 2:sea, 3: weighted average
    real(p), dimension(ix,il,3) :: ustr   !! U-stress
    real(p), dimension(ix,il,3) :: vstr   !! V-stress
    real(p), dimension(ix,il,3) :: shf    !! Sensible heat flux
    real(p), dimension(ix,il,3) :: evap   !! Evaporation [g/(m^2 s)]
    real(p), dimension(ix,il,3) :: hfluxn !! Net heat flux into surface
end module
