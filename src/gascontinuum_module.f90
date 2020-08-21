!=======================================================================
!                 GAS CONTINUUM OPACITY/EMISSION MODULE
!
! In this module various gas continuum opacities and emissivities are
! implemented. Currently we have:
!  
!  - General free-free extinction and emission for a thermal electron 
!    and ion distribution (Gordon & Sorochenko 2002, Kluwer book)
!
! Planned future continua are:
! 
!  - Various bound-free continua, which may be linked to a future 
!    planned ionization module. 
!  - H- opacity (T.L.John 1988, A&A 193, 189-192), consisting of:
!    - Free-free part:                       hnu + e- + H --> H + e- 
!    - Bound-free part (photo-detachment):   hnu + H-     --> H + e-
!  - H2- opacity (W.B.Somerville 1963, ApJ 139, 192), consisting of:
!    - Free-free part:                       hnu + e- + H2 --> H2 + e-
!    - (bound-free is argued to be negligible)
!  - H2+ opacity (Mihalas 1965)
!  - He opacity (Peach 1970)
!  - He- opacity (Carbon, Gingerich & Latham 1969)
!  - H2+H2 and H2+He opacity (Borysov 2002)
! =======================================================================
module gascontinuum_module
  use rtglobal_module
  use ioput_module
  use mathroutines_module



contains

!------------------------------------------------------------------------
!                        INITIALIZE STUFF
!------------------------------------------------------------------------
subroutine gascont_init(action)
  implicit none
  integer :: action
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  endif
  !
  ! Sanity check
  !
  if(.not.rt_incl_gascont) then
     write(stdo,*) 'ERROR: Cannot initialize gas continuum opacity if rt_incl_gascont.eq..false.'
     stop
  endif
  ! 
  ! If including free-free emission, then read the ion number density and 
  ! electron number if they have not yet been read in. Also make sure the
  ! gas temperature is known
  !
  if(rt_incl_gascont_freefree) then
     call read_gas_temperature(action)
     call read_electron_numberdensity(action)
     call read_ion_numberdensity(action)
  endif
  !
end subroutine gascont_init


!------------------------------------------------------------------------
!                        CLEAN STUFF UP
!------------------------------------------------------------------------
subroutine gascont_cleanup()
  if(allocated(electron_numdens)) deallocate(electron_numdens)
  if(allocated(ion_numdens)) deallocate(ion_numdens)
end subroutine gascont_cleanup


!------------------------------------------------------------------------
!             ADD ALL GAS CONTINUUM SOURCES TO SRC AND ALP         
!------------------------------------------------------------------------
subroutine gascont_addto_jnu_alpnu(index,nrfreq,inu0,inu1,freq,src,alp)
  implicit none
  integer :: nrfreq,inu0,inu1,inu,index
  double precision :: freq(nrfreq),src(nrfreq),alp(nrfreq),dumalpha
  !
  ! Check if we are indeed in a cell
  !
  if(index.lt.1) return
  !
  ! Switch between various processes
  !
  if(rt_incl_gascont_freefree) then
     !
     ! Include thermal free-free emission.
     !
     if(debug_check_all.eq.1) then
        if(.not.allocated(electron_numdens).or. &
           .not.allocated(ion_numdens).or.&
           .not.allocated(gastemp)) then
           write(stdo,*) 'ERROR: For free-free emission we need electron and ion number densities and gastemp'
           stop
        endif
     endif
     !
     ! Now implement the free-free formula:
     !
     !   alpha_nu = 0.2120 N_e N_i / nu^2.1 T^1.35
     !
     ! which is Eq. 2.96 of Gordon & Sorochenko (2002) Kluwer Academic Publishers
     !
     do inu=inu0,inu1
        dumalpha = 0.2120d0 * electron_numdens(index) *             &
                   ion_numdens(index) * freq(inu)**(-2.1d0) *       &
                   gastemp(index)**(-1.35d0)
        alp(inu) = alp(inu) + dumalpha
        src(inu) = src(inu) + dumalpha*bplanck(gastemp(index),freq(inu))
     enddo 
     !
  endif
  !
end subroutine gascont_addto_jnu_alpnu

end module gascontinuum_module
