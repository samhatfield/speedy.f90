subroutine ludcmp(a,n,np,indx,d)
    implicit none

    real, intent(inout) :: a(np,np), d
    integer, intent(inout) :: indx(n)
    integer, intent(in) :: n, np
    integer, parameter :: nmax = 100, tiny = 1.0e-20
    integer :: i, j, k, imax
    real :: vv(nmax), aamax, dum, sum

    d = 1.0

    do i=1,n
        aamax=0.
        do j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j)) 
        end do
        if(aamax.eq.0.) stop 'singular'
        vv(i)=1./aamax
    end do

    do j=1,n
        if(j.gt.1) then
            do i=1,j-1
                sum=a(i,j)
                if(i.gt.1) then
                    do k=1,i-1
                        sum=sum-a(i,k)*a(k,j)
                    end do
                    a(i,j)=sum
                end if
            end do
        end if

        aamax=0.
        do i=j,n
            sum=a(i,j)
            if(j.gt.1) then
                do k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
                end do
                a(i,j)=sum
            end if
            dum=vv(i)*abs(sum)
            if(dum.ge.aamax) then
                imax=i
                aamax=dum
            end if
        end do

        if(j.ne.imax) then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
        end if

        indx(j)=imax
        if(j.ne.n) then
            if(a(j,j).eq.0) a(j,j)=tiny
            dum=1./a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            end do
        end if
    end do

    if(a(n,n).eq.0.) a(n,n)=tiny
end

subroutine lubksb(a,n,np,indx,b)
    implicit none

    real, intent(inout) :: a(np,np), b(n)
    integer, intent(in) :: n, np, indx(n)
    integer :: ii, i, ll, j
    real :: sum

    ii=0

    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0) then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            end do
        else if(sum.ne.0) then
            ii=i
        end if
        b(i)=sum
    end do

    do i=n,1,-1
        sum=b(i)
        if(i.lt.n) then
          do j=i+1,n
            sum=sum-a(i,j)*b(j)
          end do
        end if
        b(i)=sum/a(i,i)
    end do
end

subroutine inv(a,y,indx,n)
    implicit none

    real, intent(inout) :: a(n,n), y(n,n)
    integer, intent(inout) :: indx(n)
    integer, intent(in) :: n
    integer :: i
    real :: d

    y = 0.0

    do i=1,n
        y(i,i)=1.
    end do

    call ludcmp(a,n,n,indx,d)

    do i=1,n
        call lubksb(a,n,n,indx,y(1,i))
    end do
end
