module diffusion_module
use constants_module
use rtglobal_module
use amr_module
use bicg_module
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

integer :: matrix_ne,diffeq_nc
double precision, allocatable :: matrix_value(:)
double precision, allocatable :: diffeq_rhs(:)
integer, allocatable :: matrix_icol(:)
integer, allocatable :: matrix_irow(:)
integer, allocatable :: matrix_idiag(:)

contains

!----------------------------------------------------------------------------
!            SMOOTH DUST TEMPERATURES WITH DIFFUSION ALGORITHM
!
! Note that the boundary condition is a fixed temperature, namely the
! temperature of species 1 and size 1.
!
! Note that:
!
!                 4 pi          /sigma     \
!    Flux = - ----------- Nabla |----- T^4 |
!             rho kappa_R       \ pi      /
!
! But in this routine, since Nabla.F=0, we ignore the 4*sigma factor.
!
###! NOTE: For multi-temperature situations we take only the first dust
###!       species, and set the other 
!----------------------------------------------------------------------------
subroutine smooth_by_diffusion(errtol,ierror)
  implicit none
  integer :: ierror
  doubleprecision :: errtol
  integer :: idiff,ndiff0
  doubleprecision :: tguess,dummy,error,thghost
  doubleprecision :: dxc,dxm,dxp,dcp,dcm,solmax
  integer :: iter,ispec,icell,index,nemax,iel,irow
  logical :: gotit
  integer :: niter
  parameter(niter=20)    ! No more than 20 iteration SHOULD be necessary
  !
  ! Check
  !
  if(igrid_type.ge.100) then
     write(stdo,*) 'ERROR: Diffusion smoothing of the dust temperature is'
     write(stdo,*) '       currently only built in for regular and AMR-type'
     write(stdo,*) '       grids.'
     stop
  endif
  !
  ! Reset error
  !
  ierror = 0
  !
  ! Check validity of all temperatures
  !
  do ispec=1,dust_nr_species 
     do icell=1,nrcells
        index = cellindex(icell)
        if(index.gt.0) then
           if(number_invalid(dusttemp(ispec,index)).ne.0) then
              write(stdo,*) 'ERROR: Invalid temperature found'
              write(stdo,*) '    before diffusion algorithm.'
              write(stdo,*) dusttemp(ispec,index)
              stop
           endif
        endif
     enddo
  enddo
  !
  ! Default initial guess of the temperature for the Rosseland opacity
  !
  tguess = 100.d0
  !
  ! Count the cells which are to be included
  !
############# NOTE: MUST SWITCH ON mc_iphotcount(:) WHEN INTENDING TO DO DIFF #######
  !
  diffeq_nc = 0
  do icell=1,nrcells
     index = cellindex(icell)
     if(mc_iphotcount(index).lt.nphotdiff) then
        !
        ! (Originally, in the RADMC version, we tested here if this cell
        ! is already included, but here we don't do this because it does
        ! not seem to be necessary in this context).
        !
        ! Add the cell itself to the diffusion
        !
        diffeq_nc = diffeq_nc + 1
     endif
  enddo
  !
  ! Estimate the number of matrix elements we will need
  !
  if(amr_dim.eq.1) then
     nemax = diffeq_nc*3+10
  elseif(amr_dim.eq.2) then
     nemax = diffeq_nc*5+10
  elseif(amr_dim.eq.3) then
     nemax = diffeq_nc*7+10
  else
     write(stdo,*) 'ERROR in diffusion: Unknown dimension of problem...'
     stop
  endif
  !
  ! Now allocate the arrays for the Bi-CG method
  !
  if(allocated(matrix_value)) deallocate(matrix_value,matrix_icol,matrix_irow)
  allocate(matrix_value(nemax),matrix_icol(nemax),matrix_irow(nemax))
  if(allocated(diffeq_rhs)) deallocate(diffeq_rhs,matrix_idiag)
  allocate(diffeq_rhs(diffeq_nc),matrix_idiag(diffeq_nc))
  !
  ! And now fill these arrays
  !
  iel = 0
  irow = 0
  do icell=1,nrcells
     index = cellindex(icell)
     if(mc_iphotcount(index).lt.nphotdiff) then
        !
        ! Increase counter of row
        !
        irow = irow + 1
        !
        ! Add matrix element for 
        !
        iel = iel + 1
     endif
  enddo

                      irdiff(ndiff) = ir
                      itdiff(ndiff) = it
c
c                     Mark this cell that the PDE has to be solved here
c
                      ipde(ndiff)   = 1
c
cc                      dmatrix(ndiff,ndiff) = 1.d0
cc                      drhs(ndiff) = tdust(1,1,it,ir)**4
c
c                     Memorize current ndiff
c
                      ndiff0 = ndiff
                  endif
c
c                 Add the 4 neighboring cells to the diffusion
c                 if they have not yet been added before
c
c                 ...r-left
c
                  gotit = .false.
                  do idiff=1,ndiff
                      if((irdiff(idiff).eq.ir-1).and.
     %                   (itdiff(idiff).eq.it)) then
                          gotit=.true.
                          idiff_rleft(ndiff0) = idiff
                          idiff_rright(idiff) = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir-1
                      itdiff(ndiff) = it
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it,ir-1)**4
cc#####################################
c                          write(*,*) '** ',idiff,ir,it,drhs(ndiff)
cc#####################################
                      idiff_rleft(ndiff0) = idiff
                      idiff_rright(idiff) = ndiff0
                  endif
c
c                 ...r-right
c
                  gotit = .false.
                  do idiff=1,ndiff
                      if((irdiff(idiff).eq.ir+1).and.
     %                   (itdiff(idiff).eq.it)) then
                          gotit=.true.
                          idiff_rright(ndiff0) = idiff
                          idiff_rleft(idiff)   = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir+1
                      itdiff(ndiff) = it
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it,ir+1)**4
cc#####################################
c                          write(*,*) '$$ ',idiff,ir,it,drhs(ndiff)
cc#####################################
                      idiff_rright(ndiff0) = idiff
                      idiff_rleft(idiff)   = ndiff0
                  endif
c
c                 ...t-left
c
                  gotit = .false.
                  do idiff=1,ndiff
                      if((itdiff(idiff).eq.it-1).and.
     %                   (irdiff(idiff).eq.ir)) then
                          gotit=.true.
                          idiff_tleft(ndiff0) = idiff
                          idiff_tright(idiff) = ndiff0
                      endif
                  enddo
                  if(.not.gotit) then
                      ndiff = ndiff + 1
                      idiff = ndiff
                      if(ndiff.gt.FRSIZE_NDIFF) then
                          write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                          write(*,*) '  CRASHED BECAUSE OF TOO FEW '
                          write(*,*) '  diffusion cell places.'
                          write(*,*) ndiff,FRSIZE_NDIFF
                          stop
                      endif
                      irdiff(ndiff) = ir
                      itdiff(ndiff) = it-1
                      dmatrix(ndiff,ndiff) = 1.d0
                      drhs(ndiff) = tdust(1,1,it-1,ir)**4
                      idiff_tleft(ndiff0) = idiff
                      idiff_tright(idiff) = ndiff0
                  endif
c
c                 ...t-right (special: equator)
c
                  if(it.lt.nt) then 
                      gotit = .false.
                      do idiff=1,ndiff
                          if((itdiff(idiff).eq.it+1).and.
     %                       (irdiff(idiff).eq.ir)) then 
                              gotit=.true.
                              idiff_tright(ndiff0) = idiff
                              idiff_tleft(idiff)   = ndiff0
                          endif
                      enddo
                      if(.not.gotit) then
                          ndiff = ndiff + 1
                          idiff = ndiff
                          if(ndiff.gt.FRSIZE_NDIFF) then
                              write(*,*) 'ERROR: DIFFUSION ALGORITHM'
                              write(*,*) '  CRASHED BECAUSE OF TOO FEW'
                              write(*,*) '  diffusion cell places.'
                              write(*,*) ndiff,FRSIZE_NDIFF
                              stop
                          endif
                          irdiff(ndiff) = ir
                          itdiff(ndiff) = it+1
                          dmatrix(ndiff,ndiff) = 1.d0
                          drhs(ndiff) = tdust(1,1,it+1,ir)**4
cc#####################################
c                          write(*,*) '== ',idiff,ir,it,drhs(ndiff)
cc#####################################
                          idiff_tright(ndiff0) = idiff
                          idiff_tleft(idiff)   = ndiff0
                      endif
c
c                     Now mark the ndiff0 cell as one in which the PDE
c                     has to be solved
c
                      ipde(ndiff0) = 1
                  endif
              endif
          enddo
      enddo
c
c     If no cells need help, then return
c
      if(ndiff.eq.0) then
          return
      endif
c
c     If cells need help, then message
c
      write(*,*) 'Now fixing low-statistics ',
     %     'cells near midplane with a diffusion method'
c
c     Make initial guess of the inverse Rosseland-mean alpha
c
      do idiff=1,ndiff
          ir = irdiff(idiff)
          it = itdiff(idiff)
          diffconst(idiff) = 1.d0 / rossmeanalp(tguess,nf,freq,
     %                       alpha_a(1,it,ir),alpha_s(1,it,ir))
      enddo
c
c     Now start the loop (iteration over Rosseland mean)
c
      do iter=1,niter
c
c         Message
c
          write(*,*) 'Diffusion iteration ',iter
c
c         Implement the diffusion equation 
c
          do idiff=1,ndiff
c
c             Only include the PDE if it has been marked as such
c
              if(ipde(idiff).eq.1) then
c
c                 Distinguish between normal and equatorial cells
c
                  if(itdiff(idiff).lt.nt) then
c
c                     Normal cell
c                 
c                     In R-direction
c
                      dxc = 0.5 * ( rc(irdiff(idiff)+1) -
     %                              rc(irdiff(idiff)-1) )
                      dxp = rc(irdiff(idiff)+1) - rc(irdiff(idiff))
                      dxm = rc(irdiff(idiff)) - rc(irdiff(idiff)-1)
                      dcp = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_rright(idiff)) )
                      dcm = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_rleft(idiff)) )
                      dmatrix(idiff,idiff_rleft(idiff)) = 
     %                           (1./dxc)*(dcm/dxm)
                      dmatrix(idiff,idiff_rright(idiff)) = 
     %                           (1./dxc)*(dcp/dxp)
                      dmatrix(idiff,idiff) = 
     %                           - (1./dxc)*( (dcm/dxm) + (dcp/dxp) )
                      drhs(idiff)   = 0.d0
c
c                     In Theta-direction
c
                      dxc = 0.5 * rc(irdiff(idiff)) * 
     %                            ( tc(itdiff(idiff)+1) -
     %                              tc(itdiff(idiff)-1) )
                      dxp = rc(irdiff(idiff)) * 
     %                        ( tc(itdiff(idiff)+1) - 
     %                          tc(itdiff(idiff)) )
                      dxm = rc(irdiff(idiff)) * 
     %                        ( tc(itdiff(idiff)) - 
     %                          tc(itdiff(idiff)-1) )
                      dcp = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_tright(idiff)) )
                      dcm = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_tleft(idiff)) )
                      dmatrix(idiff,idiff_tleft(idiff)) = 
     %                           (1./dxc)*(dcm/dxm)
                      dmatrix(idiff,idiff_tright(idiff)) = 
     %                           (1./dxc)*(dcp/dxp)
                      dmatrix(idiff,idiff) = dmatrix(idiff,idiff) 
     %                           - (1./dxc)*( (dcm/dxm) + (dcp/dxp) )
                      drhs(idiff)   = 0.d0
                  else
c                 
c                     Equatorial cell
c
c                     In R-direction
c
                      dxc = 0.5 * ( rc(irdiff(idiff)+1) -
     %                              rc(irdiff(idiff)-1) )
                      dxp = rc(irdiff(idiff)+1) - rc(irdiff(idiff))
                      dxm = rc(irdiff(idiff)) - rc(irdiff(idiff)-1)
                      dcp = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_rright(idiff)) )
                      dcm = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_rleft(idiff)) )
                      dmatrix(idiff,idiff_rleft(idiff)) = 
     %                           (1./dxc)*(dcm/dxm)
                      dmatrix(idiff,idiff_rright(idiff)) = 
     %                           (1./dxc)*(dcp/dxp)
                      dmatrix(idiff,idiff) = 
     %                           - (1./dxc)*( (dcm/dxm) + (dcp/dxp) )
                      drhs(idiff)   = 0.d0
c    
c                     In Theta-direction
c    
                      thghost = 0.5*pi + abs(0.5*pi-tc(itdiff(idiff)))
                      dxc = 0.5 * rc(irdiff(idiff)) * 
     %                            ( thghost -
     %                              tc(itdiff(idiff)-1) )
                      dxp = 1d99
                      dxm = rc(irdiff(idiff)) * 
     %                        ( tc(itdiff(idiff)) - 
     %                          tc(itdiff(idiff)-1) )
                      dcp = 1d99
                      dcm = 0.5 * ( diffconst(idiff) +
     %                              diffconst(idiff_tleft(idiff)) )
                      dmatrix(idiff,idiff_tleft(idiff)) = 
     %                           (1./dxc)*(dcm/dxm)
                      dmatrix(idiff,idiff) = dmatrix(idiff,idiff)
     %                           - (1./dxc)*( (dcm/dxm) )
                      drhs(idiff)   = 0.d0
                  endif
              endif
          enddo
c
c         Find the solution 
c
          dummy = 1.d0
          do idiff=1,ndiff 
              do it=1,ndiff
                  dmat(idiff,it) = dmatrix(idiff,it)
              enddo
              dsolold(idiff) = dsol(idiff)
              dsol(idiff)    = drhs(idiff)
          enddo
#####          call ludcmp(dmat,ndiff,FRSIZE_NDIFF,indx,dummy,success)
#####          call lubksb(dmat,ndiff,FRSIZE_NDIFF,indx,dsol)
c
c         Check if we have convergence
c
          error = 0.d0
          do idiff=1,ndiff
              dummy = abs((dsol(idiff)/(dsolold(idiff)+1d-30))-1.d0)
              if(dummy.gt.error) error = dummy
          enddo
          write(*,*) 'Nr of cells = ',ndiff
c          do idiff=1,ndiff
c              ir = irdiff(idiff)
c              it = itdiff(idiff)
cc#############################
c              write(*,*) ir,it,(drhs(idiff))**0.25,(dsol(idiff))**0.25
cc#############################
c          enddo
          write(*,*) 'The error   = ',error
c
c         If we reached convergence, then exit
c
          if(error.lt.errtol) goto 100
c
c         Put the new guess of temperature back into the grid
c
          solmax=0.d0
          do idiff=1,ndiff
              if(dsol(idiff).gt.solmax) solmax=dsol(idiff)
          enddo          
          do idiff=1,ndiff
              ir = irdiff(idiff)
              it = itdiff(idiff)
              if(dsol(idiff).le.0.d0) then
                  write(*,*) 'In diffusion.F:'
                  write(*,*) 'Zero or negative temperature...'
                  write(*,*) idiff,ir,it,dsol(idiff)
                  write(*,*) 'Correcting...'
                  dsol(idiff) = 1d-2*solmax
                  ierror = 1
              endif
              tdust(1,1,it,ir) = (dsol(idiff))**0.25
          enddo
c
c         Make new guess of the inverse Rosseland-mean alpha
c
          do idiff=1,ndiff
              ir = irdiff(idiff)
              it = itdiff(idiff)
              diffconst(idiff) = 1.d0 / rossmeanalp(tdust(1,1,it,ir),
     %                       nf,freq,
     %                       alpha_a(1,it,ir),alpha_s(1,it,ir))
          enddo
          
      enddo
c
c     No convergence
c
      write(*,*) 'ERROR: No convergence in diffusion'
      ierror=2
      return
c      stop
 100  continue
      write(*,*) 'Convergence in diffusion reached in ',iter,
     %  ' iterations.'
c
c     Now plug the new dust temperatures into the arrays 
c
      do idiff=1,ndiff
          ir = irdiff(idiff)
          it = itdiff(idiff)
          dummy = (dsol(idiff))**0.25
          do ispec=1,dust_nr_species
              do isize=1,dust_nr_size(ispec)
                  tdust(isize,ispec,it,ir) = dummy
              enddo
          enddo
      enddo
c
c     Check validity of all temperatures
c
#     ifdef CHECK_NUMBERS
      do ispec=1,dust_nr_species 
          do ir=1,nr
              do  it=1,nt
                  if(number_invalid(tdust(1,ispec,it,ir)).ne.0) then
                      write(*,*) 'ERROR: Invalid temperature found'
                      write(*,*) '    after diffusion.F'
                      write(*,*) ir,it,ispec,tdust(1,ispec,it,ir)
                      stop 5827
                  endif
                  if(tdust(1,ispec,it,ir).le.0.d0) then
                      write(*,*) 'ERROR: Zero temperature found'
                      write(*,*) '   after diffusion'
                      write(*,*) ir,it,ispec,tdust(1,ispec,it,ir)
                      stop 5828
                  endif
              enddo
          enddo
      enddo
#     endif
c
      end




c     --------------------------------------------------------------
c                         ROSSELAND MEAN OPACITY
c     --------------------------------------------------------------
      function rossmeanalp(tdust,nf,freq,alpha_a,alpha_s)
      implicit none
c
      integer nf
      doubleprecision tdust,rossmeanalp
      doubleprecision freq(FRSIZE_FREQ)
      doubleprecision alpha_a(FRSIZE_FREQ)
      doubleprecision alpha_s(FRSIZE_FREQ)
c
      doubleprecision bplanckdt,integrate
      doubleprecision dumarr1(FRSIZE_FREQ)
      doubleprecision dumarr2(FRSIZE_FREQ)
      doubleprecision dum1,dum2
c
      integer inu
c
      do inu=1,nf
          dumarr1(inu) = bplanckdt(tdust,freq(inu))
          dumarr2(inu) = dumarr1(inu)/
     %               ( alpha_a(inu) + alpha_s(inu) )
      enddo
      dum1 = integrate(nf,freq,dumarr1)
      dum2 = integrate(nf,freq,dumarr2)
      rossmeanalp = dum1/dum2
      return
c     
      end


c     --------------------------------------------------------------
c           THE TEMPERATURE DERIVATIVE OF PLANCK FUNCTION 
c     
c      This function computes the temperature derivative of the
c      Blackbody function 
c      
c         dB_nu(T)     2 h^2 nu^4      exp(h nu / kT)        1 
c         --------   = ---------- ------------------------  ---
c            dT          k c^2    [ exp(h nu / kT) - 1 ]^2  T^2
c     
c      ARGUMENTS:
c         nu    [Hz]            = Frequency
c         temp  [K]             = Electron temperature
c     --------------------------------------------------------------
      function bplanckdt(temp,nu)
      implicit none
      doubleprecision temp,nu
      doubleprecision theexp,bplanckdt
c
      theexp = exp(4.7989d-11*nu/temp)
      if(theexp.lt.1.d33) then
          bplanckdt = 7.07661334104d-58 * nu**4 * theexp / 
     %       ( (theexp-1.d0)**2 * temp**2 ) + 1.d-290
      else
          bplanckdt = 7.07661334104d-58 * nu**4 /
     %       ( theexp * temp**2 ) + 1.d-290
      endif
      return
      end


c     --------------------------------------------------------------
c                          INTEGRATION ROUTINE
c     --------------------------------------------------------------
      function integrate(n,x,f)
      implicit none
      integer n,i
      doubleprecision x(n),f(n),integrate,int
      int=0.d0
      do i=2,n
          int=int+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
      enddo
      if(x(n).gt.x(1)) then
          integrate = int
      else
          integrate = -int
      endif
      return
      end

