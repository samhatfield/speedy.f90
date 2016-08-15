subroutine inirdf(indrdf)
    !  subroutine inirdf (indrdf)
    !
    !  Purpose : Initialize random diabatic forcing 
    !  Input :   inirdf = index of forcing perturbation

    use mod_atparam
    use mod_physcon, only: slat
    use mod_randfor, only: randfh

    implicit none

    real, external :: ran1

    integer, intent(in) :: indrdf
    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    real :: redgrd(0:36,0:18), randf2(nlon,nlat), rnlon(0:18), colat(nlat)
    real :: ampl, flat1, flat2, fran, freq0, rdeg, rlon
    integer :: i, iseed, j, jlat, jlat1, jlat2, jlon, jlon1, nf, ntrfor

    integer :: nlonrg(0:18) = (/ 1,  6, 12, 18, 24, 28, 32, 34, 36, 36,&
        & 36, 34, 32, 28, 24, 18, 12,  6,  1 /)

    ! RMS aplitude of (non-null) horizontal perturbation
    ampl = 0.5

    ! Frequency of grid points with null perturbation
    freq0 = 0.

    ! ntrfor = spectral truncation of random forcing
    ntrfor = 18

    ! 1. Initialization 
    iseed = -abs(indrdf)
    if (indrdf.lt.0.) ampl = -ampl

    do i=1,20
       fran = ran1(iseed)
    end do

    do jlat=0,18
       rnlon(jlat) = float(nlonrg(jlat))/float(nlon)
    end do

    rdeg = 9./asin(1.)
    do j=1,nlat
       colat(j)=rdeg*asin(slat(j))+9.
    end do

    do nf=1,2
        ! 2. Fill reduced grid with normally-distributed random numbers
        do jlat=0,18
            call gausts(nlonrg(jlat),0.,ampl,0.,0,iseed,redgrd(1,jlat))

            if (freq0.gt.0.) then
                do jlon=1,nlonrg(jlat)
                    fran = RAN1(iseed)
                    if (fran.lt.freq0) redgrd(jlon,jlat) = 0.
                end do
            end if

           redgrd(0,jlat) = redgrd(nlonrg(jlat),jlat)
        end do

        ! 3. Interpolate random field to gaussian grid
        do j=1,nlat
            jlat1 = int(colat(j))
            jlat2 = jlat1+1

            do i=1,nlon
                rlon  = (i-1)*rnlon(jlat1)
                jlon  = int(rlon)
                flat1 = redgrd(jlon,jlat1) + (rlon-jlon) *&
                    & (redgrd(jlon+1,jlat1)-redgrd(jlon,jlat1))
 
                rlon  = (i-1)*rnlon(jlat2)
                jlon1 = int(rlon)
                flat2 = redgrd(jlon,jlat2) + (rlon-jlon) *&
                    & (redgrd(jlon+1,jlat2)-redgrd(jlon,jlat2))
 
                randf2(i,j) = flat1+(colat(j)-jlat1)*(flat2-flat1)
           end do
        end do

        ! 4. Spectral filter of gaussian-grid field
        call truncg(ntrfor,randf2,randfh(1,1,nf))
      end do
end

subroutine gausts(nt,av,sd,ac,ndis,iseed,ts)
    !  subroutine gausts (nt,av,sd,ac,ndis,iseed,ts)
    !  Computes a gaussian-dist. time series (ts) of (nt) random values, with
    !  assigned average (av), stand. dev. (sd), and lag-1 autocorrelation (ac). 
    !  Autocor. may be discontinued at the limits between (ndis) sub-series. 
    !  Uses function ran1 to generate uniform deviates from seed (iseed)
    !  Adapted from Numerical Recipes, Chapter 7.2 

    implicit none

    real, external :: ran1

    integer, intent(in) :: nt, ndis, iseed
    real, intent(in) :: av, sd, ac
    real, intent(inout) :: ts(nt)
    integer :: j, nt2, j2, jd
    real :: v1, v2, rsq, fact, sd2

    rsq = 2.0

    ! 1. Generate a time series of (nt) gaussian deviates
    do j=2,nt,2
        do while(rsq.gt.1.or.rsq.eq.0.)
            v1=2.*ran1(iseed)-1.
            v2=2.*ran1(iseed)-1.
            rsq=v1*v1+v2*v2
        end do

        fact=sqrt(-2.*log(rsq)/rsq)
        ts(j-1)=v1*fact
        ts(j)  =v2*fact
    end do

    ! 2. Introduce autocorrelation (if requested)
    if (ac.ne.0.) then
        nt2=nt/max(1,ndis)
        sd2=sqrt(1.-ac*ac)

        j=0
        do jd=1,ndis
            j=j+1
            do j2=2,nt2
                j=j+1
                ts(j)=ac*ts(j-1)+sd2*ts(j)
            end do
        end do
    end if

    ! 3. Set assigned average and standard deviation
    do j=1,nt
        ts(j)=sd*ts(j)+av
    end do
end

function ran1(idum)
    !   function ran1 (idum)
    !   Returns a uniform random deviate between 0.0 and 1.0
    !   Set IDUM to any negative value to (re)initialize the sequence
    !   From Numerical Recipes, Chapter 7.1 

    implicit none

    integer :: idum
    integer, parameter :: im=714025, ia=1366, ic=150889
    real, parameter :: rm=1./im
    integer, save :: iy = -1
    integer :: j
    real :: ran1

    integer, save :: ir(97) = 0.0

    if (idum.lt.0.or.iy.lt.0) then
        ! Initialize the shuffle array
        idum=mod(ic+abs(idum),im)

        do j=1,97
            idum=mod(ia*idum+ic,im)
            ir(j)=idum
        end do

        idum=mod(ia*idum+ic,im)
        iy=idum
    end if

    ! Get one integer number from the shuffle table 
    j=1+(97*iy)/im
    !if (j.gt.97.or.j.lt.1) stop ' error in random no. generator'
    iy=ir(j)

    ! Turn the selected integer into a real no. between 0 and 1
    ran1=iy*rm

    ! Replace the selected integer with another random integer
    idum=mod(ia*idum+ic,im)
    ir(j)=idum
end
