module mod_input
    implicit none

    private
    public load_boundary_file

contains
    function load_boundary_file(iunit,offset) result(field)
        use mod_atparam

        integer, intent(in) :: iunit, offset
        real     :: field(ix,il)
        real(4) :: inp(ix,il)
        integer :: i

        open(unit=iunit, form='unformatted', access='direct', recl=ix*4, convert='little_endian')

        do i = 1, il
            read(iunit,rec=offset*il+i) inp(:,il+1-i)
        end do

        field = inp

        ! Fix undefined values
        where (field <= -999) field = 0.0

        close(unit=iunit)
    end
end module
