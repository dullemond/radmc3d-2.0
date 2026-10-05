module diffusion_module
use rtglobal_module
!
! This module is meant to smooth out the temperature profile in regions
! that are so optically thick that their photon statistics is bad. This
! works only well (at least for now) if these regions do not have internal
! radiation sources. 
!
! This method was originally built into RADMC (the 2-D version of this
! code), but there it was based on LU-decomposition of a dense matrix.
! This meant that in that case the diffusion was limited to a certain
! maximum number of cells. That sometimes caused annoying limitations.
! Here we use the Bi-Conjugate Gradient method for solving sparse 
! linear systems. We use no preconditioner at this stage, so the 
! convergence is slow, but since it is anyway only a minor part of the
! total CPU time, it does not matter.
!
!----------------------------------------------------------------------------

contains

!----------------------------------------------------------------------------
!            SMOOTH DUST TEMPERATURES WITH DIFFUSION ALGORITHM
!
! Note that the boundary condition is a fixed temperature, namely the
! temperature of species 1 and size 1.
!
! Note that:
!
!                 4 pi          /sigma    \
!    Flux = - ----------- Nabla |----- T^4 |
!             rho kappa_R       \ pi      /
!
! But in this routine, since Nabla.F=0, we ignore the 4*sigma factor.
!
!----------------------------------------------------------------------------
subroutine smooth_by_diffusion
  implicit none
  character*200 :: current_dir, radmc3d_exe
  character(len=:), allocatable :: dummy, radmc3d_dir
  character*400 :: cmd_line

  call get_environment_variable("PWD",current_dir)
  call get_command_argument(0, radmc3d_exe)

  dummy = trim(radmc3d_exe)
  radmc3d_dir = dummy(1:len(dummy)-11)

  write(cmd_line, "(a, a, a, a, F10.2, i2, i4)") "python3 ", &
                  radmc3d_dir, "python/diffusion.py ", trim(current_dir), &
                  rt_mcparams%nphotdiff, rt_mcparams%nphotdiff_type, setthreads
  !call execute_command_line("conda init zsh")
  call execute_command_line(trim(cmd_line))

end subroutine
end module diffusion_module





