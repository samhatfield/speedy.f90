subroutine impint(dt,alph)
    ! subroutine impint(dt,alph)
    !
    ! Purpose : initialize constants for implicit computation of
    !           horizontal diffusion and gravity waves
    ! Input :   dt   = time step
    !           alph = stepping coefficient for gravity wave scheme
    !                  (0.0 = forward, 0.5 = centred, 1.0 = backward)
    ! Initialized common blocks : dync5, dync6, hdifc2

    ! IMPINT initializes constants for the implicit gravity wave computation.
    ! It is assumed that that all implicit steps are of length DELT2 and use 
    ! the forward/backward parameter ALPH.  IMPINT has to be re-called 
    ! whenever either of these two parameters is changed. IMPINT should
    ! be called even if the explicit option is chosen for the gravity wave
    ! terms (the reference state temperature TREF is subtracted from some
    ! terms anyway to reduce roundoff error; also the constants needed for
    ! the biharmonic diffusion, which is assumed always to be backwards 
    ! implicit, are defined in IMPINT)

    use mod_dyncon0, only: gamma
    use mod_atparam
    use mod_dyncon1, only: akap, rgas, hsg, dhs, fsg, fsgr, a, grav
    use mod_dyncon2
    use mod_hdifcon, only: dmp, dmpd, dmps, dmp1, dmp1d, dmp1s

    implicit none
	  								
    real, intent(in) :: dt, alph
    real :: dsum(kx), ya(kx,kx)
    integer :: indx(kx), m, n, k, k1, k2, l, ll, mm
    real :: rgam, xi, xxi, xxx

    ! 1. Constants for backwards implicit biharmonic diffusion
    do m=1,mx
        do n=1,nx 
            dmp1 (m,n)=1./(1.+dmp (m,n)*dt)
            dmp1d(m,n)=1./(1.+dmpd(m,n)*dt)
            dmp1s(m,n)=1./(1.+dmps(m,n)*dt)
        end do
    end do

    ! 1. Constants for implicit gravity wave computation
    ! reference atmosphere, function of sigma only
    rgam = rgas*gamma/(1000.*grav)

    do k=1,kx
        tref(k)=288.*max(0.2,fsg(k))**rgam
        print *, '  tref = ', tref(k)
        tref1(k)=rgas*tref(k)
        tref2(k)=akap*tref(k)
        tref3(k)=fsgr(k)*tref(k)
    end do

    ! Other constants 
    xi=dt*alph
    xxi = xi/(a*a)

    dhsx = xi * dhs

    do n=1,nx
        do m=1,mx      
            mm=isc*(m-1)+1
            ll=mm+n-2
            elz(m,n)=float(ll)*float(ll+1)*xxi
        end do
    end do
 
    !T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')

    xa(:kx,:kxm) = 0.0

    do k=1,kx
        do k1=1,kx
            ya(k,k1)=-akap*tref(k)*dhs(k1)
        end do
    end do

    do k=2,kx
        xa(k,k-1)=0.5*(akap*tref(k)/fsg(k)-(tref(k)-tref(k-1))/dhs(k))
    end do

    do k=1,kxm
        xa(k,k)=0.5*(akap*tref(k)/fsg(k)-(tref(k+1)-tref(k))/dhs(k))
    end do

    !sig(k)=xb(k,k')*d(k')
    dsum(1)=dhs(1)
    do k=2,kx
        dsum(k)=dsum(k-1)+dhs(k)
    end do

    do k=1,kxm
        do k1=1,kx
            xb(k,k1)=dhs(k1)*dsum(k)
            if(k1.le.k) xb(k,k1)=xb(k,k1)-dhs(k1)
        end do
    end do

    !t(k)=tex(k)+xc(k,k')*d(k')
    do k=1,kx
        do k1=1,kx
            xc(k,k1)=ya(k,k1)
            do k2=1,kxm
                xc(k,k1)=xc(k,k1)+xa(k,k2)*xb(k2,k1)
            end do
        end do
    end do

    !P(K)=XD(K,K')*T(K') 
    xd = 0.0

    do k=1,kx
        do k1=k+1,kx
            xd(k,k1)=rgas*log(hsg(k1+1)/hsg(k1))
        end do
    end do
    do k=1,kx
        xd(k,k)=rgas*log(hsg(k+1)/fsg(k))
    end do

    !P(K)=YE(K)+XE(K,K')*D(K')
    do k=1,kx
        do k1=1,kx
            xe(k,k1)=0.
            do k2=1,kx
                xe(k,k1)=xe(k,k1)+xd(k,k2)*xc(k2,k1)
            end do
        end do
    end do

    do l=1,lmax
        xxx=(float(l)*float(l+1))/(a*a)
        do k=1,kx
            do k1=1,kx
                xf(k,k1,l)=xi*xi*xxx*(rgas*tref(k)*dhs(k1)-xe(k,k1))
            end do
        end do
        do k=1,kx
            xf(k,k,l)=xf(k,k,l)+1.
        end do
    end do

    do l=1,lmax
        call inv(xf(1,1,l),xj(1,1,l),indx,kx)
    end do

    do k=1,kx
        do k1=1,kx
            xc(k,k1)=xc(k,k1)*xi
        end do
    end do
end
