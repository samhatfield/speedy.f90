subroutine stloop(istep)
    !   subroutine stloop (istep)
    !
    !   Purpose: Perform a series of time steps calling 
    !            post-processing/output routines at selected steps
    !   Input/output : istep = time step index
    !   Updated common block : lflag2
      
    use mod_lflags, only: lradsw, lrandf
    use mod_tsteps

    implicit none

    integer, intent(inout) :: istep
    integer :: iitest = 0, j

    do j = 1, nsteps
      if (iitest == 1) print*, 'stloop: calling step ', istep

      ! Set logical flags
      lradsw = (mod(istep,nstrad) == 1)
      lrandf = ((istep <= nstrdf) .or. (nstrdf < 0))

      ! Perform one leapfrog time step
      call step(2, 2, delt2, alph, rob, wil)   

      ! Do diagnostic, post-processing and I/O tasks 
      call diagns(2, istep)

      if (mod(istep, nstppr) == 0) call tminc

      if (nstout > 0 .and. mod(istep, nstout) == 0) call tmout(1)

      istep = istep + 1
    enddo
end
