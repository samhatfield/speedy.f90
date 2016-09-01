module mod_fft
    use mod_atparam, only: ix

    implicit none

    private
    public wsave

    real :: wsave(2*ix+15)
end module
