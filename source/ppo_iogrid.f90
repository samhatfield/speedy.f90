!  Read or write a gridded file in sigma coordinate
!  Input :   imode = 1 : read model variables from a gridded file (sigma)
!                  = 4 : write model variables  to a gridded file (sigma)
!                  = 5 : write a GrADS control file (for sigma)
!  Created by Takemasa Miyoshi
!  Converted to FORTRAN 90 by Sam Hatfield
subroutine iogrid(imode)
    use mod_atparam, only: ix, iy, nx, mx, il, kx
    use mod_physcon, only: p0, sig
    use mod_dynvar
    use mod_dyncon1, only: radang, grav
    use mod_date, only: model_datetime
    use mod_tsteps
    use mod_flx_land
    use mod_flx_sea

    implicit none

    integer, parameter :: ngp=ix*il
    integer, intent(in) :: imode

    complex, dimension(mx,nx) :: ucostmp, vcostmp
    real, dimension(ngp,kx) :: ugr, vgr, tgr, qgr, phigr
    real, dimension(ngp) :: psgr, rrgr
    real(4), dimension(ngp,kx) :: ugr4, vgr4, tgr4, qgr4, phigr4
    real(4), dimension(ngp) :: psgr4(ngp), rrgr4(ngp)

    ! File names etc.
    character(len=16) :: filename='yyyymmddhhmm.grd'
    character(len=16) :: ctlname='yyyymmddhhmm.ctl'
    character(len=3) :: cmon3='JAN'
    integer :: irec
    integer :: j, k

    if (imode.eq.1) then
        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
            & 'Read gridded dataset for year/month/date/hour: ',&
            & model_datetime%year,'/',model_datetime%month,'/',model_datetime%day,'/',model_datetime%hour

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
    else if (imode.eq.4) then
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

        ! Output
        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)',&
            & 'Write gridded dataset for year/month/date/hour/minute: ', &
            & model_datetime%year,'/',model_datetime%month,'/',model_datetime%day,'/', &
            & model_datetime%hour,'/',model_datetime%minute

        ugr4 = ugr
        vgr4 = vgr
        tgr4 = tgr
        qgr4 = qgr*1.0d-3 ! kg/kg
        phigr4 = phigr/grav   ! m
        psgr4 = p0*exp(psgr)! Pa
        rrgr4 = rrgr

        write (filename(1:4),'(i4.4)') model_datetime%year
        write (filename(5:6),'(i2.2)') model_datetime%month
        write (filename(7:8),'(i2.2)') model_datetime%day
        write (filename(9:10),'(i2.2)') model_datetime%hour
        write (filename(11:12),'(i2.2)') model_datetime%minute
        open (99,file=filename,form='unformatted',access='direct',&
            & recl=4*ix*il)

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
        write (99,rec=irec) (psgr4(j),j=1,ngp)
        irec=irec+1
        write (99,rec=irec) (rrgr4(j),j=1,ngp)
        close (99)

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
    else if (imode.eq.5) then
        ! 3. Write a GrADS control file (3:p,5:sigma)
        if (model_datetime%month.eq.1) then
            cmon3='JAN'
        else if (model_datetime%month.eq.2) then
            cmon3='FEB'
        else if (model_datetime%month.eq.3) then
            cmon3='MAR'
        else if (model_datetime%month.eq.4) then
            cmon3='APR'
        else if (model_datetime%month.eq.5) then
            cmon3='MAY'
        else if (model_datetime%month.eq.6) then
            cmon3='JUN'
        else if (model_datetime%month.eq.7) then
            cmon3='JUL'
        else if (model_datetime%month.eq.8) then
            cmon3='AUG'
        else if (model_datetime%month.eq.9) then
            cmon3='SEP'
        else if (model_datetime%month.eq.10) then
            cmon3='OCT'
        else if (model_datetime%month.eq.11) then
            cmon3='NOV'
        else if (model_datetime%month.eq.12) then
            cmon3='DEC'
        end if

        write (ctlname(1:4),'(I4.4)') model_datetime%year
        write (ctlname(5:6),'(I2.2)') model_datetime%month
        write (ctlname(7:8),'(I2.2)') model_datetime%day
        write (ctlname(9:10),'(I2.2)') model_datetime%hour
        write (ctlname(11:12), '(I2.2)') model_datetime%minute
        open (11,file=ctlname,form='formatted')
        write (11,'(A)') 'DSET ^%y4%m2%d2%h2%m2.grd'

        write (11,'(A)') 'TITLE SPEEDY MODEL OUTPUT'
        write (11,'(A)') 'UNDEF -9.99E33'
        write (11,'(A)') 'OPTIONS template big_endian'
        write (11,'(A)') 'XDEF 96 LINEAR 0.0 3.75'
        write (11,'(A,48F8.3)') 'YDEF 48 LEVELS ',&
            & (RADANG(J)*90.0d0/ASIN(1.0d0),J=1,48)

        write (11,'(A,8F6.3)') 'ZDEF 8 LEVELS ',(sig(k),k=8,1,-1)

        write (11,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I4.4,A)') 'TDEF ',&
            & 2,' LINEAR ',model_datetime%hour,':',model_datetime%minute,'Z',&
            & model_datetime%day,cmon3,model_datetime%year,' 6HR'
        write (11,'(A)') 'VARS 6'

        write (11,'(A)') 'U 8 99 U-wind [m/s]'
        write (11,'(A)') 'V 8 99 V-wind [m/s]'
        write (11,'(A)') 'T 8 99 Temperature [K]'
        write (11,'(A)') 'Q 8 99 Specific Humidity [kg/kg]'
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
