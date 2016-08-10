subroutine forint(ngp,imon,fmon,for12,for1)  
    ! Aux. routine FORINT : linear interpolation of monthly-mean forcing

    implicit none

    integer, intent(in) :: ngp, imon
    real, intent(in) :: fmon, for12(ngp,*)
    real, intent(inout) :: for1(ngp)
    integer :: imon2
    real :: wmon

    if (fmon.le.0.5) then
        imon2 = imon-1
        if (imon.eq.1) imon2 = 12
        wmon = 0.5-fmon
    else
        imon2 = imon+1
        if (imon.eq.12) imon2 = 1
        wmon = fmon-0.5
    end if

    for1 = for12(:,imon) + wmon*(for12(:,imon2) - for12(:,imon))
end

subroutine forin5(ngp,imon,fmon,for12,for1)
    ! Aux. routine FORIN5 : non-linear, mean-conserving interpolation 
    !                       of monthly-mean forcing fields

    implicit none

    integer, intent(in) :: ngp, imon
    real, intent(in) :: fmon, for12(ngp,12)
    real, intent(inout) :: for1(ngp)
    integer :: im1, im2, ip1, ip2
    real :: c0, t0, t1, t2, t3, wm1, wm2, w0, wp1, wp2

    im2 = imon-2
    im1 = imon-1
    ip1 = imon+1
    ip2 = imon+2

    if (im2.lt.1)  im2 = im2+12
    if (im1.lt.1)  im1 = im1+12
    if (ip1.gt.12) ip1 = ip1-12
    if (ip2.gt.12) ip2 = ip2-12
 
    c0 = 1./12.
    t0 = c0*fmon
    t1 = c0*(1.-fmon)
    t2 = 0.25*fmon*(1-fmon)

    wm2 =        -t1   +t2
    wm1 =  -c0 +8*t1 -6*t2
    w0  = 7*c0      +10*t2     
    wp1 =  -c0 +8*t0 -6*t2
    wp2 =        -t0   +t2 

    for1 = wm2*for12(:,im2) + wm1*for12(:,im1) + w0*for12(:,imon) +&
        & wp1*for12(:,ip1) + wp2*for12(:,ip2)
end
