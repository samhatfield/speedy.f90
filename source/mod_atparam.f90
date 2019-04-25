module mod_atparam
    implicit none

    private
    public trunc, ix, iy, il
    public nx, mx, mxnx
    public kx, ntr

    integer, parameter :: trunc = 30, ix = 96, iy = 24, il = 2*iy
    integer, parameter :: nx = trunc+2, mx = trunc+1, mxnx = mx*nx
    integer, parameter :: kx = 8, ntr=1
end module
