subroutine iogrid(imode)
    !  SUBROUTINE IOGRID (IMODE)
    !  Created by Takemasa Miyoshi
    !  Converted to FORTRAN 90 by Sam Hatfield
    !
    !  Purpose : read or write a gridded file in sigma coordinate
    !  Input :   IMODE = 1 : read model variables from a gridded file (sigma)
    !                  = 2 : write model variables  to a gridded file (p)
    !                  = 3 : write a GrADS control file (for p)
    !                  = 4 : write model variables  to a gridded file (sigma)
    !                  = 5 : write a GrADS control file (for sigma)
    !  Initialized common blocks (if IMODE = 1) : DYNSP1, SFCANOM 
    !

    use mod_atparam, only: ix, iy, nx, mx, il, kx
    use mod_physcon, only: p0, gg, rd, sig, sigl, pout
    use mod_dynvar
    use mod_dyncon1
    use mod_date
    use mod_tsteps
    use mod_tmean
    use mod_flx_land
    use mod_flx_sea

    implicit none

    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat
    integer, intent(in) :: imode

    ! Check what this file corresponds to in my version
!    include "com_anomvar.h"

    complex, dimension(mx,nx) :: ucostmp, vcostmp
    real, dimension(ngp,kx) :: ugr, vgr, tgr, qgr, phigr
    real :: psgr(ngp)
    real, dimension(ngp,kx) :: ugr1, vgr1, tgr1, qgr1, phigr1
    real :: rrgr1(ngp), aref, phi1, phi2, textr, tref
    real(4), dimension(ngp,kx) :: ugr4, vgr4, tgr4, qgr4, phigr4
    real(4), dimension(ngp) :: psgr4(ngp), rrgr4(ngp)

    ! For vertical interpolation !adapted from ppo_tminc.f
    integer :: k0(ngp)
    real :: w0(ngp), zout(ngp), zinp(nlev), rdzinp(nlev)

    ! File names etc.
    character(len=14) :: filename='yyyymmddhh.grd'
    character(len=14) :: ctlname='yyyymmddhh.ctl'
    character(len=16) :: filenamep='yyyymmddhh_p.grd'
    character(len=16) :: ctlnamep='yyyymmddhh_p.ctl'
    character(len=3) :: cmon3='JAN'
    integer :: irec
    integer :: iitest=0
    integer :: j, k

    if (imode.eq.1) then
        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
            & 'Read gridded dataset for year/month/date/hour: ',&
            & iyear,'/',imonth,'/',iday,'/',ihour

        open (90,form='unformatted',access='direct',recl=4*ngp)
        irec=1
        do k=kx,1,-1
            read (90,rec=irec) (ugr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            read (90,rec=irec) (vgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            read (90,rec=irec) (tgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            read (90,rec=irec) (qgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        read (90,rec=irec) (psgr4(j),j=1,ngp)
        close (90)

        ugr = ugr4
        vgr = vgr4
        tgr = tgr4
        qgr = qgr4 *1.0d3
        psgr = psgr4
        psgr = log(psgr/p0)
        if(iitest==1) print *,' UGR  :',minval(ugr),maxval(ugr)
        if(iitest==1) print *,' VGR  :',minval(vgr),maxval(vgr)
        if(iitest==1) print *,' TGR  :',minval(tgr),maxval(tgr)
        if(iitest==1) print *,' QGR  :',minval(qgr),maxval(qgr)
        if(iitest==1) print *,' PSGR :',minval(psgr),maxval(psgr)

        ! Conversion from gridded variable to spectral variable
        do k=1,kx
            call vdspec(ugr(1,k),vgr(1,k),vor(1,1,k,1),div(1,1,k,1),2)
            call spec(tgr(1,k),t(1,1,k,1))
            call spec(qgr(1,k),tr(1,1,k,1,1))
            if(ix.eq.iy*4) then
                call trunct(vor(1,1,k,1))
                call trunct(div(1,1,k,1))
                call trunct(t(1,1,k,1))
                call trunct(tr(1,1,k,1,1))
            end if
        end do
        call spec(psgr(1),ps(1,1,1))
        if (ix.eq.iy*4) call trunct(ps(1,1,1))
    else if (imode.eq.2.or.imode.eq.4) then
        ! 2. Write date and model variables to the gridded file (2:P,4:sigma)

        ! Conversion from spectral model variable to gridded variable
        do k=1,kx
           call uvspec(vor(1,1,k,1),div(1,1,k,1),ucostmp,vcostmp)
           call grid(ucostmp,ugr(1,k),2)
           call grid(vcostmp,vgr(1,k),2)
        end do

        do k=1,kx
           call grid(t(1,1,k,1),tgr(1,k),1)
           call grid(tr(1,1,k,1,1),qgr(1,k),1)
           call grid(phi(1,1,k),phigr(1,k),1)
        end do

        call grid(ps(1,1,1),psgr(1),1)

        ! Vertical interpolation from sigma level to pressure level (ppo_tminc.f)
        if (imode.eq.2) then ! p-level output
            zinp(1) = -sigl(1)
            do k=2,nlev
               zinp(k) = -sigl(k)
               rdzinp(k) = 1.0d0/(zinp(k-1)-zinp(k))
            end do

            do k=1,kx
                do j=1,ngp
                    zout(j) = psgr(j) - log(pout(k))
                end do

                call setvin(zinp,rdzinp,zout,ngp,kx,k0,w0)

                call verint(tgr1(1,k),tgr,ngp,kx,k0,w0)

                do j=1,ngp
                    if(zout(j).lt.zinp(nlev)) then
                        textr = max(tgr1(j,k),tgr(j,nlev))
                        aref = rd*0.006d0/gg * (zinp(nlev)-zout(j))
                        tref = tgr(j,nlev)*(1.0d0+aref+0.5*aref*aref)
                        tgr1(j,k) = textr + 0.7d0*(tref-textr)
                    end if
                end do

                do j=1,ngp
                    w0(j) = max(w0(j),0.0d0)
                end do

                call verint(ugr1(1,k),ugr,ngp,kx,k0,w0)
                call verint(vgr1(1,k),vgr,ngp,kx,k0,w0)
                call verint(qgr1(1,k),qgr,ngp,kx,k0,w0)

                do j=1,ngp
                   phi1 = phigr(j,k0(j))&
                       & +0.5*rd*(tgr1(j,k)+tgr(j,k0(j)))*(zout(j)-zinp(k0(j)))
                   phi2 = phigr(j,k0(j)-1)&
                       & +0.5*rd*(tgr1(j,k)+tgr(j,k0(j)-1))*(zout(j)-zinp(k0(j)-1))
                   phigr1(j,k) = phi1 + w0(j)*(phi2-phi1)
                end do

                do j=1,ngp
                  rrgr1(j) = save2d_d2(j,1)& ! g/m^2/s
                      &*3.6*4.0/real(nsteps)*6.0 ! mm/6hr
                end do
            end do
        else  ! sigma-level output
            ugr1 = ugr
            vgr1 = vgr
            tgr1 = tgr
            qgr1 = qgr
            phigr1 = phigr
        end if

        ! Output
        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
            & 'Write gridded dataset for year/month/date/hour: ',&
            & iyear,'/',imonth,'/',iday,'/',ihour

        ugr4 = ugr1
        vgr4 = vgr1
        tgr4 = tgr1
        qgr4 = qgr1*1.0d-3 ! kg/kg
        phigr4 = phigr1/gg   ! m
        psgr4 = p0*exp(psgr)! Pa
        rrgr4 = rrgr1

        if (imode.eq.2) then
            write (filenamep(1:4),'(i4.4)') iyear
            write (filenamep(5:6),'(i2.2)') imonth
            write (filenamep(7:8),'(i2.2)') iday
            write (filenamep(9:10),'(i2.2)') ihour
            open (99,file=filenamep,form='unformatted',access='direct',&
                & recl=4*ix*il)
        else
            write (filename(1:4),'(i4.4)') iyear
            write (filename(5:6),'(i2.2)') imonth
            write (filename(7:8),'(i2.2)') iday
            write (filename(9:10),'(i2.2)') ihour
            open (99,file=filename,form='unformatted',access='direct',&
                & recl=4*ix*il)
        end if
        irec=1
        do k=kx,1,-1
            write (99,rec=irec) (ugr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            write (99,rec=irec) (vgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            write (99,rec=irec) (tgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        do k=kx,1,-1
            write (99,rec=irec) (qgr4(j,k),j=1,ngp)
            irec=irec+1
        end do
        if (imode.eq.2) then !Z output is only for p-level
            do k=kx,1,-1
                write (99,rec=irec) (phigr4(j,k),j=1,ngp)
                irec=irec+1
            end do
        end if
        write (99,rec=irec) (psgr4(j),j=1,ngp)
        irec=irec+1
        write (99,rec=irec) (rrgr4(j),j=1,ngp)
        close (99)
        if(iitest==1) print *,' UGR  :',minval(ugr4),maxval(ugr4)
        if(iitest==1) print *,' VGR  :',minval(vgr4),maxval(vgr4)
        if(iitest==1) print *,' TGR  :',minval(tgr4),maxval(tgr4)
        if(iitest==1) print *,' QGR  :',minval(qgr4),maxval(qgr4)
        if(iitest==1) print *,' PHIGR:',minval(phigr4),maxval(phigr4)
        if(iitest==1) print *,' PSGR :',minval(psgr4),maxval(psgr4)
        if(iitest==1) print *,' RRGR :',minval(rrgr4),maxval(rrgr4)

        open (100,file='fluxes.grd',form='unformatted',access='direct',recl=8*ix*il)
        write (100,rec=1) (prec_l(j),j=1,ngp)
        write (100,rec=2) (snowf_l(j),j=1,ngp)
        write (100,rec=3) (evap_l(j),j=1,ngp)
        write (100,rec=4) (hflux_l(j),j=1,ngp)

        write (100,rec=5) (prec_s(j),j=1,ngp)
        write (100,rec=6) (snowf_s(j),j=1,ngp)
        write (100,rec=7) (evap_s(j),j=1,ngp)
        write (100,rec=8) (ustr_s(j),j=1,ngp)
        write (100,rec=9) (vstr_s(j),j=1,ngp)
        write (100,rec=10) (ssr_s(j),j=1,ngp)
        write (100,rec=11) (slr_s(j),j=1,ngp)
        write (100,rec=12) (shf_s(j),j=1,ngp)
        write (100,rec=13) (ehf_s(j),j=1,ngp)
        write (100,rec=14) (hflux_s(j),j=1,ngp)
        write (100,rec=15) (hflux_i(j),j=1,ngp)
        close (100)
    else if (imode.eq.3.or.imode.eq.5) then
        ! 3. Write a GrADS control file (3:p,5:sigma)
        if (imonth.eq.1) then
            cmon3='JAN'
        else if (imonth.eq.2) then
            cmon3='FEB'
        else if (imonth.eq.3) then
            cmon3='MAR'
        else if (imonth.eq.4) then
            cmon3='APR'
        else if (imonth.eq.5) then
            cmon3='MAY'
        else if (imonth.eq.6) then
            cmon3='JUN'
        else if (imonth.eq.7) then
            cmon3='JUL'
        else if (imonth.eq.8) then
            cmon3='AUG'
        else if (imonth.eq.9) then
            cmon3='SEP'
        else if (imonth.eq.10) then
            cmon3='OCT'
        else if (imonth.eq.11) then
            cmon3='NOV'
        else if (imonth.eq.12) then
            cmon3='DEC'
        end if

        if (imode.eq.3) then !p-level
            write (ctlnamep(1:4),'(I4.4)') iyear
            write (ctlnamep(5:6),'(I2.2)') imonth
            write (ctlnamep(7:8),'(I2.2)') iday
            write (ctlnamep(9:10),'(I2.2)') ihour
            open (11,file=ctlnamep,form='formatted')
            write (11,'(A)') 'DSET ^%y4%m2%d2%h2_p.grd'
        else !sigma-level
            write (ctlname(1:4),'(I4.4)') iyear
            write (ctlname(5:6),'(I2.2)') imonth
            write (ctlname(7:8),'(I2.2)') iday
            write (ctlname(9:10),'(I2.2)') ihour
            open (11,file=ctlname,form='formatted')
            write (11,'(A)') 'DSET ^%y4%m2%d2%h2.grd'
        end if
        write (11,'(A)') 'TITLE SPEEDY MODEL OUTPUT'
        write (11,'(A)') 'UNDEF -9.99E33'
        write (11,'(A)') 'OPTIONS template big_endian'
        write (11,'(A)') 'XDEF 96 LINEAR 0.0 3.75'
        write (11,'(A,48F8.3)') 'YDEF 48 LEVELS ',&
            & (RADANG(J)*90.0d0/ASIN(1.0d0),J=1,48)
        if (imode.eq.3) then
            write (11,'(A)') 'ZDEF 8 LEVELS 925 850 700 500 300 200 100 30'
        else
            write (11,'(A,7F6.3)') 'ZDEF 8 LEVELS ',(sig(k),k=8,1,-1)
        end if
        if (ndaysl.ne.0) then
            write (11,'(A,I4,A,I2.2,A,I2.2,A,I4.4,A)') 'TDEF ',&
               & ndaysl*4+1,' LINEAR ',ihour,'Z',iday,cmon3,iyear,' 6HR'
        else
            write (11,'(A,I4,A,I2.2,A,I2.2,A,I4.4,A)') 'TDEF ',&
                & 2,' LINEAR ',ihour,'Z',iday,cmon3,iyear,' 6HR'
        end if
        if (imode.eq.3) then !p-level
            write (11,'(A)') 'VARS 7'
        else !sigma-level
           write (11,'(A)') 'VARS 6'
        end if
        write (11,'(A)') 'U 8 99 U-wind [m/s]'
        write (11,'(A)') 'V 8 99 V-wind [m/s]'
        write (11,'(A)') 'T 8 99 Temperature [K]'
        write (11,'(A)') 'Q 8 99 Specific Humidity [kg/kg]'
        if (imode.eq.3) then
          write (11,'(A)') 'Z 8 99 Geopotential Height [m]'
        end if
        write (11,'(A)') 'PS 0 99 Surface Pressure [Pa]'
        write (11,'(A)') 'RAIN 0 99 Precipitation [mm/6hr]'
        write (11,'(A)') 'ENDVARS'
        close (11)
    else
        print *,'Hey, look at the usage! (IOGRID)'
        stop
    end if

    return

    ! 4. Stop integration if gridded file is not found
    200 continue

    print*, ' Hey, what are you doing?',&
        & ' fort.2 should contain time setting'

    stop 'invalid gridded data input'
end
