!> author: Sam Hatfield
!  date: 08/07/2019
!  For setting aliases to commonly used types.
module types
    use iso_fortran_env, only: real32, real64

    implicit none

    ! Set single and double-precision variables and the actual default precision used by the model
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    integer, parameter :: p = dp
end module
