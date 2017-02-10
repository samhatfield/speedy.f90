subroutine diagns(jj,istep)
    ! subroutine diagns(jj,istep)

    ! Purpose: print global means of eddy kinetic energy and temperature 
    ! Input : jj    = time level index (1 or 2)
    !         istep = time step index


    use mod_tsteps, only: nstdia, nstppr, nstout, ihout
    use mod_atparam
    use mod_dynvar

    implicit none

    integer, intent(in) :: jj, istep
    integer, parameter :: nlon=ix, nlat=il, nlev=kx, ngp=nlon*nlat

    integer :: k, m, n, kk
    complex :: temp(mx,nx)
    real :: diag(kx,3), sqhalf

    ! 1. Get global-mean temperature and compute eddy kinetic energy 
    sqhalf = sqrt(0.5)

    do k=1,kx
        diag(k,1)=0.
        diag(k,2)=0.
        diag(k,3)=sqhalf*real(t(1,1,k,jj))

        call invlap(vor(1,1,k,jj),temp)

        do m=2,mx
            do n=1,nx
                diag(k,1)=diag(k,1)-real(temp(m,n)*conjg(vor(m,n,k,jj)))
            end do
        end do

        call invlap(div(1,1,k,jj),temp)

        do m=2,mx
            do n=1,nx
                diag(k,2)=diag(k,2)-real(temp(m,n)*conjg(div(m,n,k,jj)))
            end do
        end do
    end do

    ! 2. Print results to screen
    if (mod(istep,nstdia).eq.0) then
        print 2001, istep, (diag(k,1),k=1,kx)
        print 2002,        (diag(k,2),k=1,kx)
        print 2003,        (diag(k,3),k=1,kx)
    end if

    ! 3. Stop integration if model variables are out of range
    do k=1,kx
        if (diag(k,1).gt.500.or.diag(k,2).gt.500.or.diag(k,3).lt.180.or.&
            & diag(k,3).gt.320.) then

            print 2001, istep, (diag(kk,1),kk=1,kx)
            print 2002,        (diag(kk,2),kk=1,kx)
            print 2003,        (diag(kk,3),kk=1,kx)

            ! Write model fields at t-1 on output file 
            if (ihout .eqv. .false.) then !Only when no hourly output
                call tmout(0)
                call tminc
                nstout=nstppr
                call tmout(1)
            end if

            stop '*** model variables out of accepted range ***'
        end if
    end do

    2001 format(' step =',i6,' reke =',(10f8.2))
    2002 format         (13x,' deke =',(10f8.2))
    2003 format         (13x,' temp =',(10f8.2))
end
