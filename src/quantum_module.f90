module quantum_module
use rtglobal_module

doubleprecision,allocatable :: emisquant(:,:)
doubleprecision,allocatable :: emisquant_loccum(:,:)
doubleprecision,allocatable :: emisquant_loctot(:)
doubleprecision,allocatable :: emisquant_cum(:)
doubleprecision,allocatable :: miquant(:,:)
doubleprecision :: emisquanttot

contains

!---------------------------------------------------------------------
!                 Initialize the quantum mode
!---------------------------------------------------------------------
subroutine quantum_init()
implicit none
integer :: ierr
write(stdo,*) 'Initializing quantum module'
if(nrcellsmax.le.0) then
   write(stdo,*) 'INTERNAL ERROR IN QUANTUM_MODULE: MONTE CARLO NOT INITIALIZED'
   stop
endif
allocate(emisquant(1:freq_nr,1:nrcellsmax),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate emisquant.'
   stop
endif
allocate(emisquant_loccum(1:freq_nr+1,1:nrcellsmax),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate emisquant_loccum.'
   stop
endif
allocate(emisquant_loctot(1:nrcellsmax),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate emisquant_loctot.'
   stop
endif
allocate(emisquant_cum(nrcellsmax+1),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate emisquant_cum.'
   stop
endif
allocate(miquant(1:freq_nr,1:nrcellsmax),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate miquant.'
   stop
endif
end subroutine quantum_init


!---------------------------------------------------------------------
!                   Cleanup the quantum mode
!---------------------------------------------------------------------
subroutine quantum_cleanup()
implicit none
if(allocated(emisquant)) deallocate(emisquant)
if(allocated(emisquant_loccum)) deallocate(emisquant_loccum)
if(allocated(emisquant_loctot)) deallocate(emisquant_loctot)
if(allocated(emisquant_cum)) deallocate(emisquant_cum)
if(allocated(miquant)) deallocate(miquant)
end subroutine quantum_cleanup

end module quantum_module
