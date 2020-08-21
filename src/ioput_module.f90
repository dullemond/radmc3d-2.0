module ioput_module
use amr_module
use amrray_module
use stars_module
use dust_module
use rtglobal_module
use constants_module

!character*160 :: gridfilename
!integer :: gridfilenamelen

contains



!-------------------------------------------------------------------------
!                        READ DUST DENSITY
!-------------------------------------------------------------------------
subroutine read_dust_density(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: index,ierr,irec,ispec,i,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(dustdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Now read in the density, either from a special file or from the standard
  ! dust_density.inp
  !
  inquire(file='dust_density.inp',exist=fex1)
  inquire(file='dust_density.uinp',exist=fex2)
  inquire(file='dust_density.binp',exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file dust_density.*inp'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'ERROR: Could not find any dust_density.*inp...'
     stop
  endif
  write(stdo,*) 'Reading dust densities...'
  call flush(stdo)
  if(fex1) then
     !
     ! Formatted input
     !
     style = 1
     open(unit=1,file='dust_density.inp',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_density.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
     read(1,*) kk
  elseif(fex2) then
     !
     ! Unformatted input, old fortran style (with records)
     !
     style = 2
     open(unit=1,file='dust_density.uinp',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_density.uinp is invalid/unknown.'
        write(stdo,*) 'Format number = ',iiformat
        write(stdo,*) 'Record length = ',reclen
        stop
     endif
     reclenn = reclen
     read(1) nn,kk
  else
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     style = 3
     open(unit=1,file='dust_density.binp',status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_density.binp is invalid/unknown.'
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn,kk
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: dust_density.inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  if(kk.ne.dust_nr_species) then
     write(stdo,*) 'ERROR: In the dust_density.uinp file the number'
     write(stdo,*) '  of dust species is different from that '
     write(stdo,*) '  specified in dustopac.inp'
     write(stdo,*) kk,dust_nr_species
     stop
  endif
  !
  ! Allocate the arrays
  !
  if(allocated(dustdens)) deallocate(dustdens)
  if(allocated(dust_massdust)) deallocate(dust_massdust)
  allocate(dustdens(1:dust_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust density array'
     stop
  endif
  allocate(dust_massdust(1:dust_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust mass array'
     stop
  endif
  !
  ! Now read the dust density
  !
  do ispec=1,dust_nr_species
     call read_scalarfield(1,style,precis,nrcellsinp,               &
          dust_nr_species,1,ispec,1,1d-99,reclenn,scalar1=dustdens)
  enddo
  !
  ! Close the file
  !
  close(1)
  !
  ! Compute mass in each dust species and print
  !
  call compute_dust_mass()
  if(dust_nr_species.gt.1) then
     do ispec=1,dust_nr_species
        write(stdo,*) 'Dust mass in species ',ispec,' = ',     &
             dust_massdust(ispec)/1.99d33,' Msun'
     enddo
  endif
  write(stdo,*) 'Dust mass total = ',                      &
       dust_massdusttot/1.99d33,' Msun'
  !
end subroutine read_dust_density



!-------------------------------------------------------------------------
!                       READ DUST TEMPERATURE 
!-------------------------------------------------------------------------
subroutine read_dust_temperature(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: i,irec,ispec,index,ierr,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(dusttemp)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! First check that the user has not accidently created a 
  ! .inp or .uinp file. This is a common confusion and may lead to
  ! a lot of bug hunting for nothing. The reason why this file is
  ! a .dat and not .inp file is because normally this is a file that
  ! is created by RADMC-3D by the Monte Carlo simulation, and later 
  ! re-read by RADMC-3D for the imaging/spectra. It is therefore a mix
  ! of .inp and .out file, and we thus call it a .dat (or .udat) file.
  !
  inquire(file='dust_temperature.inp',exist=fex1)
  inquire(file='dust_temperature.uinp',exist=fex2)
  inquire(file='dust_temperature.binp',exist=fex3)
  if(fex1.or.fex2.or.fex3) then
     write(stdo,*) 'ERROR: Found dust_temperature.inp or dust_temperature.uinp or dust_temperature.binp'
     write(stdo,*) '       file. But dust_temperature must be a .dat or .udat or .bdat file.'
     stop
  endif
  !
  ! Open temperature file
  !
  inquire(file='dust_temperature.dat',exist=fex1)
  inquire(file='dust_temperature.udat',exist=fex2)
  inquire(file='dust_temperature.bdat',exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file dust_temperature.*dat'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'ERROR: Could not find any dust_temperature.*dat file...'
     stop
  endif
  write(stdo,*) 'Reading dust temperatures...'
  call flush(stdo)
  if(fex1) then
     !
     ! Formatted input
     !
     style = 1
     open(unit=1,file='dust_temperature.dat',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_temperature.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
     read(1,*) kk
  elseif(fex2) then
     !
     ! Unformatted input, old fortran style (with records)
     !
     style = 2
     open(unit=1,file='dust_temperature.udat',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_temperature.udat is invalid/unknown.'
        stop
     endif
     reclenn=reclen
     read(1) nn,kk
  else
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     style = 3
     open(unit=1,file='dust_temperature.bdat',status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of dust_temperature.bdat is invalid/unknown.'
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn,kk
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: dust_temperature.inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  if(kk.ne.dust_nr_species) then
     write(stdo,*) 'ERROR: In the dust_temperature.inp file the number'
     write(stdo,*) '  of dust species is different from that '
     write(stdo,*) '  specified in dustopac.inp'
     write(stdo,*) kk,dust_nr_species
     stop
  endif
  !
  ! Allocate the array
  !
  if(allocated(dusttemp)) deallocate(dusttemp)
  allocate(dusttemp(1:dust_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate dust temperature array'
     stop
  endif
  !
  ! Now read the dust temperature
  !
  do ispec=1,dust_nr_species
     call read_scalarfield(1,style,precis,nrcellsinp,                &
          dust_nr_species,1,ispec,1,1d-99,reclenn,scalar1=dusttemp)
  enddo
  !
  ! Close the file
  !
  close(1)
  !
end subroutine read_dust_temperature


!-------------------------------------------------------------------------
!                       READ GAS DENSITY
!-------------------------------------------------------------------------
subroutine read_gas_density(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: i,index,irec,ierr,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(gasdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Open temperature file
  !
  inquire(file='gas_density.inp',exist=fex1)
  inquire(file='gas_density.uinp',exist=fex2)
  inquire(file='gas_density.binp',exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file gas_density.*dat'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'ERROR: Could not find any gas_density.*dat file...'
     stop
  endif
  write(stdo,*) 'Reading gas density...'
  call flush(stdo)
  !
  ! Now read this temperature
  !
  if(fex1) then
     !
     ! Formatted gas temperature file
     !
     style = 1
     open(unit=1,file='gas_density.inp',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of gas_density.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
  elseif(fex2) then
     !
     ! Unformatted input, old fortran style (with records)
     !
     style = 2
     open(unit=1,file='gas_density.uinp',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of gas_density.uinp is invalid/unknown.'
        stop
     endif
     reclenn=reclen
     read(1) nn
  else
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     style = 3
     open(unit=1,file='gas_density.binp',status='old')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of gas_density.binp is invalid/unknown.'
        stop
     endif
     read(1) nn
     precis = nn
     read(1,*) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: gas_density.inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Allocate the array
  !
  if(allocated(gasdens)) deallocate(gasdens)
  allocate(gasdens(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate gas density array'
     stop
  endif
  !
  ! Now read the gas density
  !
  call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn,&
                        scalar0=gasdens)
  !
  ! Close the file
  !
  close(1)
  !
end subroutine read_gas_density


!-------------------------------------------------------------------------
!                       READ GAS TEMPERATURE 
!-------------------------------------------------------------------------
subroutine read_gas_temperature(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: icell,i,ierr,irec,index,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(gastemp)) then
        !
        ! The gas temperature is already determined. It might be
        ! that also the gastmax has been calculated from
        ! this, but we don't know for sure. So let's calculate
        ! it here to make sure gastmax is up-to-date. It
        ! can be important if the gas temperature is created 
        ! using the userdef module.
        !
        ! Thanks, Rainer Rolffs, for pointing this out!
        ! 02.01.2011
        !
        gastmax = 0.d0
        do icell=1,nrcells
           index = cellindex(icell)
           if(gastemp(index).gt.gastmax) then
              gastmax = gastemp(index)
           endif
        enddo
        return
     endif
  endif
  !
  ! Default
  !
  precis = 8
  ! 
  ! If the user requests, then use one of the dust temperatures as gas
  ! temperature
  !
  if(tgas_eq_tdust.gt.0) then
     !
     ! Use dust temperature
     !
     ! First make sure the dust data are read
     ! if they are not already read yet
     !
     call read_dustdata(1)
     call read_dust_density(1)
     call read_dust_temperature(1)
     !
     ! Check something
     !
     if(tgas_eq_tdust.gt.dust_nr_species) then
        write(stdo,*) 'ERROR: tgas_eq_tdust is larger than dust_nr_species'
        stop
     endif
     !
     ! Allocate array
     !
     if(allocated(gastemp)) deallocate(gastemp)
     allocate(gastemp(1:nrcells),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate gastemp array'
        stop
     endif
     write(stdo,*) 'Using dust temperature as gas temperature...'
     call flush(stdo)
     !
     ! Now go and map the temperature 
     !
     do icell=1,nrcells
        index = cellindex(icell)
        if((index.le.0).or.(index.gt.nrcellsmax)) then
           write(stdo,*) 'ERROR: Internal error while reading gas temperature'
           stop
        endif
        gastemp(index) = dusttemp(tgas_eq_tdust,index)
     enddo
  else
     !
     ! Open temperature file
     !
     inquire(file='gas_temperature.inp',exist=fex1)
     inquire(file='gas_temperature.uinp',exist=fex2)
     inquire(file='gas_temperature.binp',exist=fex3)
     idum=0
     if(fex1) idum=idum+1
     if(fex2) idum=idum+1
     if(fex3) idum=idum+1
     if(idum.gt.1) then
        write(stdo,*) 'ERROR: Found more than one file gas_temperature.*inp'
        stop
     endif
     if(idum.eq.0) then
        write(stdo,*) 'ERROR: Could not find any gas_temperature.*inp file...'
        stop
     endif
     write(stdo,*) 'Reading gas temperature...'
     call flush(stdo)
     !
     ! Now read this temperature
     !
     if(fex1) then
        !
        ! Formatted gas temperature file
        !
        style = 1
        open(unit=1,file='gas_temperature.inp',status='old')
        read(1,*) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of gas_temperature.inp is invalid/unknown.'
           stop
        endif
        read(1,*) nn
     elseif(fex2) then
        !
        ! F77-style unformatted input of gas temperature
        !
        style = 2
        open(unit=1,file='gas_temperature.uinp',status='old',form='unformatted')
        read(1) iiformat,reclen
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of gas_temperature.uinp is invalid/unknown.'
           stop
        endif
        reclenn=reclen
        read(1) nn
     else
        !
        ! Binary input: C-compliant unformatted streaming data
        !
        style = 3
        open(unit=1,file='gas_temperature.binp',status='old',access='stream')
        read(1) iiformat
        if(iiformat.ne.1) then
           write(stdo,*) 'ERROR: Format number of gas_temperature.binp is invalid/unknown.'
           stop
        endif
        read(1) nn
        precis = nn
        read(1) nn
     endif
     !
     ! Do some checks
     !
     if(nn.ne.nrcellsinp) then
        write(stdo,*) 'ERROR: gas_temperature.*inp does not have same number'
        write(stdo,*) '       of cells as the grid.'
        write(stdo,*) nn,nrcellsinp
        stop
     endif
     !
     ! Allocate array
     !
     if(allocated(gastemp)) deallocate(gastemp)
     allocate(gastemp(1:nrcells),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR: Could not allocate gastemp array'
        stop
     endif
     !
     ! Now read the gas temperature
     !
     call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn, &
                           scalar0=gastemp)
     !
     ! Close the file
     !
     close(1)
     !
  endif
  !
  ! Compute the maximum temperature
  !
  gastmax = 0.d0
  do icell=1,nrcells
     index = cellindex(icell)
     if(gastemp(index).gt.gastmax) then
        gastmax = gastemp(index)
     endif
  enddo
  !
end subroutine read_gas_temperature


!-------------------------------------------------------------------------
!                   READ ELECTRON NUMBER DENSITY
!-------------------------------------------------------------------------
subroutine read_electron_numberdensity(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: i,index,irec,ierr,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(electron_numdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Open temperature file
  !
  inquire(file='electron_numdens.inp',exist=fex1)
  inquire(file='electron_numdens.uinp',exist=fex2)
  inquire(file='electron_numdens.binp',exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file electron_numdens.*inp'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'ERROR: Could not find any electron_numdens.*inp file...'
     stop
  endif
  write(stdo,*) 'Reading electron number density...'
  call flush(stdo)
  !
  ! Now read the electron number density
  !
  if(fex1) then
     !
     ! Formatted electron number density temperature file
     !
     style = 1
     open(unit=1,file='electron_numdens.inp',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of electron_numdens.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
  elseif(fex2) then
     !
     ! Unformatted input of electron numbder density temperature
     !
     style = 2
     open(unit=1,file='electron_numdens.uinp',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of electron_numdens.uinp is invalid/unknown.'
        stop
     endif
     reclenn=reclen
     read(1) nn
  else
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     style = 3
     open(unit=1,file='electron_numdens.binp',status='old',form='unformatted')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of electron_numdens.binp is invalid/unknown.'
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: electron_numdens.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Allocate array
  !
  if(allocated(electron_numdens)) deallocate(electron_numdens)
  allocate(electron_numdens(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate electron_numdens array'
     stop
  endif
  !
  ! Now read the electron numberdens
  !
  call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn, &
                        scalar0=electron_numdens)
  !
  ! Close the file
  !
  close(1)
  !
end subroutine read_electron_numberdensity


!-------------------------------------------------------------------------
!                   READ ION NUMBER DENSITY
!-------------------------------------------------------------------------
subroutine read_ion_numberdensity(action)
  implicit none
  integer :: action
  logical :: fex1,fex2,fex3
  double precision, allocatable :: data(:)
  double precision :: dummy
  integer(kind=8) :: iiformat,reclen,nn,kk
  integer :: i,index,irec,ierr,idum,style,precis,reclenn
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(ion_numdens)) return
  endif
  !
  ! Default
  !
  precis = 8
  !
  ! Open temperature file
  !
  inquire(file='ion_numdens.inp',exist=fex1)
  inquire(file='ion_numdens.uinp',exist=fex2)
  inquire(file='ion_numdens.uinp',exist=fex3)
  idum=0
  if(fex1) idum=idum+1
  if(fex2) idum=idum+1
  if(fex3) idum=idum+1
  if(idum.gt.1) then
     write(stdo,*) 'ERROR: Found more than one file ion_numdens.*inp'
     stop
  endif
  if(idum.eq.0) then
     write(stdo,*) 'ERROR: Could not find any ion_numdens.*inp file...'
     stop
  endif
  write(stdo,*) 'Reading ion number density...'
  call flush(stdo)
  !
  ! Now read this temperature
  !
  if(fex1) then
     !
     ! Formatted ion number density temperature file
     !
     style = 1
     open(unit=1,file='ion_numdens.inp',status='old')
     read(1,*) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of ion_numdens.inp is invalid/unknown.'
        stop
     endif
     read(1,*) nn
  elseif(fex2) then
     !
     ! Unformatted input of ion numbder density temperature
     !
     style = 2
     open(unit=1,file='ion_numdens.uinp',status='old',form='unformatted')
     read(1) iiformat,reclen
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of ion_numdens.uinp is invalid/unknown.'
        stop
     endif
     reclenn=reclen
     read(1) nn
  else
     !
     ! Binary input: C-compliant unformatted streaming data
     !
     style = 3
     open(unit=1,file='ion_numdens.binp',status='old',access='stream')
     read(1) iiformat
     if(iiformat.ne.1) then
        write(stdo,*) 'ERROR: Format number of ion_numdens.uinp is invalid/unknown.'
        stop
     endif
     read(1) nn
     precis = nn
     read(1) nn
  endif
  !
  ! Do some checks
  !
  if(nn.ne.nrcellsinp) then
     write(stdo,*) 'ERROR: ion_numdens.*inp does not have same number'
     write(stdo,*) '       of cells as the grid.'
     write(stdo,*) nn,nrcellsinp
     stop
  endif
  !
  ! Allocate array
  !
  if(allocated(ion_numdens)) deallocate(ion_numdens)
  allocate(ion_numdens(1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate ion_numdens array'
     stop
  endif
  !
  ! Now read the ion number density
  !
  call read_scalarfield(1,style,precis,nrcellsinp,1,1,1,1,1d-99,reclenn, &
                        scalar0=ion_numdens)
  !
  ! Close the file
  !
  close(1)
  !
end subroutine read_ion_numberdensity


!-------------------------------------------------------------------------
!                  WRITE RESULTS OF BJORKMAN & WOOD
!-------------------------------------------------------------------------
subroutine write_results_mc_bjorkmanwood(ifile)
  implicit none
  integer :: ifile
  !
  call write_dust_temperature()
  !
end subroutine write_results_mc_bjorkmanwood


!-------------------------------------------------------------------------
!                    WRITE GRID FILE (UNIVERSAL)
!-------------------------------------------------------------------------
subroutine write_grid_file()
  implicit none
  character*160 :: filename
  if(igrid_type.lt.100) then
     if(rto_style.eq.1) then
        filename = 'amr_grid.inp'
     elseif(rto_style.eq.2) then
        filename = 'amr_grid.uinp'
     elseif(rto_style.eq.3) then
        filename = 'amr_grid.binp'
     else
        write(stdo,*) 'ERROR: Do not know I/O style ',rto_style
        stop
     endif
     write(stdo,*) 'Writing grid file to ',trim(filename)
     call amr_write_grid(filename,rto_style,rtio_gridinfo)
  else
     write(stdo,*) 'ERROR: In write_grid_file(): unstructured grids not yet supported'
     stop
  endif
end subroutine write_grid_file


!-------------------------------------------------------------------------
!                    READ GRID FILE (UNIVERSAL)
!-------------------------------------------------------------------------
subroutine read_grid_file()
  implicit none
  character*160 :: filename
  logical :: fex_amrgrid,fex_delaunaygrid,fex_voronoigrid
  logical :: fex1,fex2,fex3
  integer :: amr,delaunay,voronoi,icell,idir,style,idum
  double precision :: dd
  !
  ! Check out which grid file is present
  !
  inquire(file='amr_grid.inp',exist=fex1)
  inquire(file='amr_grid.uinp',exist=fex2)
  inquire(file='amr_grid.binp',exist=fex3)
  fex_amrgrid=fex1.or.fex2.or.fex3
  inquire(file='triang_vertex_grid.inp',exist=fex1)
  inquire(file='triang_vertex_grid.uinp',exist=fex2)
  inquire(file='triang_vertex_grid.binp',exist=fex3)
  fex_delaunaygrid=fex1.or.fex2.or.fex3
  inquire(file='voronoi_grid.inp',exist=fex1)
  inquire(file='voronoi_grid.uinp',exist=fex2)
  inquire(file='voronoi_grid.binp',exist=fex3)
  fex_voronoigrid=fex1.or.fex2.or.fex3
  if(fex_amrgrid) then
     amr=2 
  else 
     amr=1
  endif
  if(fex_delaunaygrid) then
     delaunay=2 
  else 
     delaunay=1
  endif
  if(fex_voronoigrid) then
     voronoi=2 
  else 
     voronoi=1
  endif
  if(amr*delaunay*voronoi.gt.2) then
     write(stdo,*) 'ERROR: More than one (amr/delaunay/voronoi)_grid.inp file found.'
     write(stdo,*) '       You must choose which kind of grid you want to use.'
     stop
  elseif(amr*delaunay*voronoi.eq.1) then
     write(stdo,*) 'ERROR: No file (amr/delaunay/voronoi)_grid.inp found.'
     write(stdo,*) '       I need this file to set up a grid.'
     stop
  endif
  !
  ! For each grid type call the appropriate routine. This automatically
  ! initializes the grid and makes the linear list of cells that is 
  ! necessary for storing data. For the AMR grid this is the list that is
  ! computed by amr_compute_list_all(). The order in which the leafs (=cells)
  ! are listed in this list is also the order in which their cells are 
  ! read/written.
  !
  ! NOTE: In the future, for extremely large simulations that do not fit 
  !       onto a single node we may be forced to do this differently. We
  !       would then (in the AMR grid) have partial trees. Still, the 
  !       order of the leafs would still the same as the order of the 
  !       list. 
  !
  if(fex_amrgrid) then
     inquire(file='amr_grid.inp',exist=fex1)
     inquire(file='amr_grid.uinp',exist=fex2)
     inquire(file='amr_grid.binp',exist=fex3)
     idum=0
     if(fex1) idum=idum+1
     if(fex2) idum=idum+1
     if(fex3) idum=idum+1
     if(idum.gt.1) then
        write(stdo,*) 'ERROR: Found more than one file amr_grid.*inp'
        stop
     endif
     if(idum.eq.0) then
        write(stdo,*) 'ERROR: Could not find any amr_grid.*inp...'
        stop
     endif
     if(fex1) then
        filename  = 'amr_grid.inp'
        style = 1
     elseif(fex2) then
        filename  = 'amr_grid.uinp'
        style = 2
     else
        filename  = 'amr_grid.binp'
        style = 3
     endif
     call amr_read_grid(filename,style,checkspher=.true.)
     igrid_type = amr_style
     if(igrid_type.ge.100) stop 4091
     if(igrid_type.lt.0) stop 4092
     igrid_coord = amr_coordsystem
  else
     write(stdo,*) 'ERROR: In read_grid_file(): unstructured grids not yet supported'
     stop
  endif
  !
  ! Remember that we read this grid from a file, i.e. we did NOT 
  ! generate the grid internally in the code. This means that when 
  ! writing results of e.g. the Monte Carlo simulation, we do not
  ! need to write the grid again.
  !
  grid_was_read_from_file = .true.
  !
end subroutine read_grid_file


!-------------------------------------------------------------------------
!                        POST PROCESS GRID
!
! NOTE: You must call this routine AFTER you read the star data, because
!       this routine computes some stuff related to stars (if they are
!       present in your model.
!-------------------------------------------------------------------------
subroutine postprocess_grid_for_use(renew)
  implicit none
  double precision :: dd
  integer :: i,ierr,idum,idir,icell,ilev,ilevel
  integer :: ix,iy,iz,istar
  logical, optional :: renew
  logical :: dorenew
  character*80 :: str1,str2
  double precision :: dummy
  !
  ! Check optional argument
  ! 
  if(present(renew)) then 
     dorenew = renew 
  else
     dorenew = .false.
  endif
  !
  ! If dorenew then unallocate arrays
  !
  if(dorenew) then
     if(allocated(cellindex)) deallocate(cellindex)
     if(allocated(cellvolume)) deallocate(cellvolume)
     nrcells    = 0
     nrcellsmax = 0
  endif
  !
  ! Switch between grid forms
  !
  if((igrid_type.ge.0).and.(igrid_type.lt.100)) then
     !
     ! For backward compatibility, redefine igrid_type to be 
     ! the same as the amr_style
     !
     igrid_type = amr_style
     !
     ! Check that the coordinates are always increasing
     !
     do ix=2,amr_grid_nx+1
        if(amr_grid_xi(ix,1).le.amr_grid_xi(ix-1,1)) then
           write(stdo,*) 'ERROR: coordinates must always increase!'
           stop
        endif
     enddo
     do iy=2,amr_grid_ny+1
        if(amr_grid_xi(iy,2).le.amr_grid_xi(iy-1,2)) then
           write(stdo,*) 'ERROR: coordinates must always increase!'
           stop
        endif
     enddo
     do iz=2,amr_grid_nz+1
        if(amr_grid_xi(iz,3).le.amr_grid_xi(iz-1,3)) then
           write(stdo,*) 'ERROR: coordinates must always increase!'
           stop
        endif
     enddo
     !
     ! Import the coordinate system number
     !
     igrid_coord = amr_coordsystem
     !
     ! Check the grid
     !
     if(igrid_coord.lt.100) then
        !
        ! Cartesian grid
        !
        ! Just for safety, check that the user does not use a coordinate
        ! system (in case of Cartesian coordinates) in which the (0,0,0)
        ! point is very very far outside of the domain. That could give
        ! trouble in the imaging routines of the camera_module, because
        ! some things are not calculated with translational invariance.
        !
        if((abs(amr_grid_xi(amr_grid_nx+1,1)-amr_grid_xi(1,1)).lt.                &
             1d-4*max(abs(amr_grid_xi(amr_grid_nx+1,1)),abs(amr_grid_xi(1,1)))).or. &
             (abs(amr_grid_xi(amr_grid_ny+1,2)-amr_grid_xi(1,2)).lt.                &
             1d-4*max(abs(amr_grid_xi(amr_grid_ny+1,2)),abs(amr_grid_xi(1,2)))).or. &
             (abs(amr_grid_xi(amr_grid_nz+1,3)-amr_grid_xi(1,3)).lt.                &
             1d-4*max(abs(amr_grid_xi(amr_grid_nz+1,3)),abs(amr_grid_xi(1,3))))) then
           write(stdo,*) 'ERROR: Please put computational domain closer'
           write(stdo,*) '       to the (0,0,0) point of the coordinates'
           write(stdo,*) '       to avoid numerical accuracy problems.'
           stop
        endif
        !
        ! The normal Cartesian grid mode is 3-D, but as of version 0.31
        ! there are also special Cartesian modes: 1-D plane-parallel
        ! and 2-D pencil-parallel. Both cases involve "infinitely extended
        ! dimensions". We catch these here. 
        !
        if((igrid_coord.eq.10).or.((amr_grid_xi(1,2).le.-1d89).and. &
             (amr_grid_xi(amr_grid_ny+1,2).gt.1d89).and.            &
             (amr_ydim.eq.0).and.(amr_grid_ny.eq.1).and.            &
             (amr_grid_xi(1,1).le.-1d89).and.                       &
             (amr_grid_xi(amr_grid_nx+1,1).gt.1d89).and.            &
             (amr_xdim.eq.0).and.(amr_grid_nx.eq.1))) then
           !
           ! We have the 1-D plane parallel mode, in z-direction
           !
           if(amr_grid_nx.ne.1) then
              write(stdo,*) 'ERROR: 1-D cartesian plane-parallel mode: Must have amr_grid_nx=1.'
              stop
           endif
           if(amr_grid_ny.ne.1) then
              write(stdo,*) 'ERROR: 1-D cartesian plane-parallel mode: Must have amr_grid_ny=1.'
              stop
           endif
           igrid_coord       = 10
           amr_grid_xi(1,1)  = -1d90
           amr_grid_xi(2,1)  =  1d90
           amr_grid_xi(1,2)  = -1d90
           amr_grid_xi(2,2)  =  1d90
           amr_finegrid_xi(1:2,1,0) = amr_grid_xi(1:2,1)
           amr_finegrid_xi(1:2,2,0) = amr_grid_xi(1:2,2)
           amr_finegrid_xc(1,1,0)   = 0.5d0 * ( amr_finegrid_xi(1,1,0) + &
                                                amr_finegrid_xi(2,1,0) )
           amr_finegrid_xc(1,2,0)   = 0.5d0 * ( amr_finegrid_xi(1,2,0) + &
                                                amr_finegrid_xi(2,2,0) )
           do ilevel=1,amr_levelmax
              amr_finegrid_xi(1:2,1,ilevel) = amr_finegrid_xi(1:2,1,ilevel-1)
              amr_finegrid_xc(1,1,ilevel)   = amr_finegrid_xc(1,1,ilevel-1)
              amr_finegrid_xi(1:2,2,ilevel) = amr_finegrid_xi(1:2,2,ilevel-1)
              amr_finegrid_xc(1,2,ilevel)   = amr_finegrid_xc(1,2,ilevel-1)
           enddo
           amr_xdim          = 0
           amr_ydim          = 0
           amr_xyzdim(1)     = 0
           amr_xyzdim(2)     = 0
           amr_cyclic_xyz(1) = .false.
           amr_cyclic_xyz(2) = .false.
        elseif((igrid_coord.eq.20).or.((amr_grid_xi(1,1).le.-1d89).and. &
             (amr_grid_xi(amr_grid_nz+1,1).gt.1d89).and.                &
             (amr_xdim.eq.0).and.(amr_grid_nx.eq.1))) then
           !
           ! We have the 2-D pencil parallel mode, in y-z direction
           !
           if(amr_grid_nx.ne.1) then
              write(stdo,*) 'ERROR: 2-D cartesian pencil-parallel mode: Must have amr_grid_nx=1.'
              stop
           endif
           igrid_coord       = 20
           amr_grid_xi(1,1)  = -1d90
           amr_grid_xi(2,1)  =  1d90
           amr_finegrid_xi(1:2,1,0) = amr_grid_xi(1:2,1)
           amr_finegrid_xc(1,1,0)   = 0.5d0 * ( amr_finegrid_xi(1,1,0) + &
                                                amr_finegrid_xi(2,1,0) )
           do ilevel=1,amr_levelmax
              amr_finegrid_xi(1:2,1,ilevel) = amr_finegrid_xi(1:2,1,ilevel-1)
              amr_finegrid_xc(1,1,ilevel)   = amr_finegrid_xc(1,1,ilevel-1)
           enddo
           amr_xdim          = 0
           amr_xyzdim(1)     = 0
           amr_cyclic_xyz(1) = .false.
           grid_contsph_x = 1d99  ! Should not be necessary; in 2-D only a circle
           grid_contsph_y = 0.5d0 * ( amr_grid_xi(amr_grid_ny+1,2) + amr_grid_xi(1,2) )
           grid_contsph_z = 0.5d0 * ( amr_grid_xi(amr_grid_nz+1,3) + amr_grid_xi(1,3) )
           grid_contsph_r = 0.5001d0 * sqrt( ( amr_grid_xi(amr_grid_ny+1,2) - amr_grid_xi(1,2) )**2 +  &
                                             ( amr_grid_xi(amr_grid_nz+1,3) - amr_grid_xi(1,3) )**2 )
        else
           !
           ! We have the normal 3-D cartesian mode
           !
           ! If using Cartesian AMR gridding, then we must insist that the
           ! base cells are cubic and all of the same size. This is because
           ! the camera_module.f90 imaging routines can only find the proper
           ! recursive refinement of the image pixels if this is so. This
           ! recursive refinement is important for ensuring that all the flux
           ! is picked up properly. 
           !
           dd  = amr_grid_xi(2,1) - amr_grid_xi(1,1)
           if(dd.le.0.d0) then
              write(stdo,*) 'ERROR while reading grid: base grid is bad.'
              stop
           endif
           do i=1,amr_grid_nx
              if(abs((amr_grid_xi(i+1,1)-amr_grid_xi(i,1))/dd-1.d0).gt.1d-4) then
                 write(stdo,*) 'ERROR while reading grid: base grid is not regularly cubically spaced.'
                 stop
              endif
           enddo
           do i=1,amr_grid_ny
              if(abs((amr_grid_xi(i+1,2)-amr_grid_xi(i,2))/dd-1.d0).gt.1d-4) then
                 write(stdo,*) 'ERROR while reading grid: base grid is not regularly cubically spaced.'
                 stop
              endif
           enddo
           do i=1,amr_grid_nz
              if(abs((amr_grid_xi(i+1,3)-amr_grid_xi(i,3))/dd-1.d0).gt.1d-4) then
                 write(stdo,*) 'ERROR while reading grid: base grid is not regularly cubically spaced.'
                 stop
              endif
           enddo
           !
           ! Compute grid center and the radius of the containing sphere.
           ! These variables are useful for various Monte Carlo things, 
           ! including the external luminosity, and improved photon statistics
           ! for externally located stars.
           !
           grid_contsph_x = 0.5d0 * ( amr_grid_xi(amr_grid_nx+1,1) + amr_grid_xi(1,1) )
           grid_contsph_y = 0.5d0 * ( amr_grid_xi(amr_grid_ny+1,2) + amr_grid_xi(1,2) )
           grid_contsph_z = 0.5d0 * ( amr_grid_xi(amr_grid_nz+1,3) + amr_grid_xi(1,3) )
           grid_contsph_r = 0.5001d0 * sqrt( ( amr_grid_xi(amr_grid_nx+1,1) - amr_grid_xi(1,1) )**2 +  &
                                             ( amr_grid_xi(amr_grid_ny+1,2) - amr_grid_xi(1,2) )**2 +  &
                                             ( amr_grid_xi(amr_grid_nz+1,3) - amr_grid_xi(1,3) )**2 )
        endif
        !   
     elseif(igrid_coord.lt.200) then
        !
        ! Spherical coordinates
        !
        ! Check that the grid is always monotonically increasing
        ! 
        do i=1,amr_grid_nx
           if(amr_grid_xi(i+1,1).le.amr_grid_xi(i,1)) then
              write(stdo,*) 'ERROR while reading grid: base R grid is not monotonically increasing.'
              stop
           endif
        enddo
        do i=1,amr_grid_ny
           if(amr_grid_xi(i+1,2).le.amr_grid_xi(i,2)) then
              write(stdo,*) 'ERROR while reading grid: base Theta grid is not monotonically increasing.'
              stop
           endif
        enddo
        do i=1,amr_grid_nz
           if(amr_grid_xi(i+1,3).le.amr_grid_xi(i,3)) then
              write(stdo,*) 'ERROR while reading grid: base Phi grid is not monotonically increasing.'
              stop
           endif
        enddo
        !
        ! Compute grid center and the radius of the containing sphere.
        ! These variables are useful for various Monte Carlo things, 
        ! including the external luminosity, and improved photon statistics
        ! for externally located stars.
        !
        grid_contsph_x = 0.d0
        grid_contsph_y = 0.d0
        grid_contsph_z = 0.d0
        grid_contsph_r = 1.001 * amr_grid_xi(amr_grid_nx+1,1)
        !
     else
        stop 609
     endif
     !
     ! Because the radiative transfer code should be as independent as 
     ! possible of the details of the AMR grid (so that other types of
     ! grids can easily be included too) we make our own independent
     ! indexing array here and have some other copies of variables.
     !
     ! NOTE (bugfix 28.07.09): Here we must allocate the arrays with
     !      nrcellsmax size instead of nrcells. The reason is that 
     !      the user MAY want to refine the grid lateron.
     !
     nrcells     = amr_nrleafs
     !
     ! Now set the nrcellsinp value, which is the counter for how
     ! many cells will be read. This is equal to nrcells for regular
     ! and AMR grids of the oct-tree type. But for the layer-style
     ! AMR there is redundancy, to keep the data format easy.
     !
     if(amr_style.lt.10) then
        nrcellsinp = nrcells
     elseif(amr_style.lt.20) then
        nrcellsinp = amr_nrbranches
     else
        write(stdo,*) 'ERROR: AMR grid style ',amr_style,' does not exist.'
        stop
     endif
     !
     ! Allocate the arrays only if they have not yet been allocated
     !
     if(.not.allocated(cellindex)) then 
        !
        ! Not yet allocated, so allocate
        !
        !nrcellsmax  = amr_nrleafs       ! bugfix 28.07.09
        nrcellsmax  = amr_nrleafs_max
        ! Bug? 16.11.09: Why do I allocate with nrcells and not with nrcellsmax?
        allocate(cellindex(1:nrcells),STAT=ierr)
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cellindex array'
           stop
        endif
        if(amr_tree_present) then
           !
           ! The AMR tree is available, so get the index from the 
           ! amr_theleaf_index() array
           !
           do icell=1,nrcells
              cellindex(icell) = amr_theleaf_index(icell)
           enddo
        else
           !
           ! Regular grid, so cellindex = icell
           !
           do icell=1,nrcells
              cellindex(icell) = icell
           enddo
        endif
        ! Bug? 16.11.09: Why do I allocate with nrcells and not with nrcellsmax?
        allocate(cellvolume(1:nrcells),STAT=ierr) 
        if(ierr.ne.0) then
           write(stdo,*) 'ERROR: Could not allocate cellvolume array'
           stop
        endif
     else
        !
        ! Already allocated, so just do some checks
        !
        if(nrcellsmax.ne.amr_nrleafs_max) then
           write(stdo,*) 'ERROR: cellindex already allocated, but with '
           write(stdo,*) '       different array size...'
           stop
        endif
     endif
     !
     ! Compute the cell volumes
     !
     call compute_cellvolumes_amr()
     !
     ! Calculate the ray_nsmax
     ! The factor of 2 is necessary for spherical or cylindrical coordinates.
     !
     ray_nsmax = 0
     do idir=1,3
        idum = 2 * amr_grid_nx * 2**amr_levelmax + 1
        ray_nsmax = max(ray_nsmax,idum)
     enddo
     !
     ! Do a message
     !
     ! NOTE: Since we changed the levels to 0:amr_levelmax (from originally
     !       1:amr_levelmax), we must also start from 0 here, contrary to
     !       earlier versions.
     !
     write(stdo,*) 'Grid information (current status):'
     write(str1,*) amr_nrbranches
     write(str2,*) amr_nrleafs
     write(stdo,*) '  We have ',trim(adjustl(str1)),' branches, of which ', &
          trim(adjustl(str2)),' are actual grid cells.'
     write(str1,321) 100*amr_nrbranches/(1.d0*amr_nrbranches_max)
     write(str2,321) 100*amr_nrleafs/(1.d0*amr_nrleafs_max)
321 format(F7.3)
     write(stdo,*) '  ---> ',trim(adjustl(str1)),'% mem use for branches, and ', &
          trim(adjustl(str2)),'% mem use for actual cells.'
     ilev = 0
     if(amr_tree_present) then
        do icell=1,nrcells
           if(ilev.lt.amr_theleafs(icell)%link%level) then
              ilev = amr_theleafs(icell)%link%level
           endif
        enddo
     endif
     if(ilev.eq.0) then
        if(amr_tree_present) then
           write(stdo,*) '  No grid refinement is active (at present; userdef may change this).'
        else
           write(stdo,*) '  No grid refinement is active. The AMR tree is not allocated (this saves memory).'
        endif
     else
        write(str1,*) ilev
        write(stdo,*) '  Top refinement level of the grid = ',trim(adjustl(str1))
     endif
     !
  elseif((igrid_type.ge.100).and.(igrid_type.lt.200)) then
     write(stdo,*) 'ERROR: In postprocess_grid_for_use(): unstructured grids not yet supported'
     stop
  else
     write(stdo,*) 'ERROR: In postprocess_grid_for_use(): igrid_type ',igrid_type,' not known.'
     stop
  endif
  !
  ! Compute the farthest distance of a star from the edge of the 
  ! containing sphere
  !
  star_maxexterndist = 0.d0
  if((nstars.gt.0).and.(grid_contsph_r.gt.0.d0)) then
     do istar=1,nstars
        dummy = sqrt((star_pos(1,istar)-grid_contsph_x)**2+&
                     (star_pos(2,istar)-grid_contsph_y)**2+&
                     (star_pos(3,istar)-grid_contsph_z)**2)
        if(dummy.gt.grid_contsph_r) then
           if(dummy-grid_contsph_r+star_r(istar).gt.star_maxexterndist) then
              star_maxexterndist = dummy-grid_contsph_r+star_r(istar)
           endif
        endif
     enddo
  endif
  !
end subroutine postprocess_grid_for_use



!-------------------------------------------------------------------------
!           COMPUTE ALL CELL VOLUMES AND INDICES OF AMR GRID
!-------------------------------------------------------------------------
subroutine compute_cellvolumes_amr()
  implicit none
  !
  integer :: ileaf,leafindex,iddr
  type(amr_branch), pointer :: leaf
  doubleprecision :: axi(1:2,1:3)
  integer :: ix,iy,iz
  !
  ! Check if we can use the leaf list of the AMR module
  !
  if(amr_tree_present) then
     if(.not.allocated(amr_theleafs)) then
        write(stdo,*) 'ERROR while computing cell volume: The leaf list is not yet made.'
        stop
     endif
  endif
  !
  ! Now loop over all leafs
  !
  if(igrid_coord.lt.100) then  
     !
     ! Cartesian coordinates
     !
     do ileaf=1,amr_nrleafs
        !
        ! Get cell wall information
        !
        if(amr_tree_present) then
           leafindex = amr_theleaf_index(ileaf)
           leaf      => amr_theleafs(ileaf)%link
           do iddr=1,3
              axi(1,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr),iddr,leaf%level)
              axi(2,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr)+1,iddr,leaf%level)
           enddo
        else
           leafindex=ileaf
           call amr_regular_get_ixyz(ileaf,ix,iy,iz)
           axi(1,1) = amr_finegrid_xi(ix,1,0)
           axi(2,1) = amr_finegrid_xi(ix+1,1,0)
           axi(1,2) = amr_finegrid_xi(iy,2,0)
           axi(2,2) = amr_finegrid_xi(iy+1,2,0)
           axi(1,3) = amr_finegrid_xi(iz,3,0)
           axi(2,3) = amr_finegrid_xi(iz+1,3,0)
        endif
        !
        ! Now compute cell volume, but take special care of
        ! the 1-D plane-parallel and 2-D pencil-parallel cases
        !
        if(igrid_coord.eq.10) then
           !
           ! 1-D Plane parallel coordinates
           !
           cellvolume(leafindex) = ( axi(2,3) - axi(1,3) )
           !
        elseif(igrid_coord.eq.20) then
           !
           ! 2-D Pencil parallel coordinates
           !
           cellvolume(leafindex) = ( axi(2,2) - axi(1,2) ) *   &
                                   ( axi(2,3) - axi(1,3) )
        else
           !
           ! Normal 3-D cartesian coordinates
           !
           cellvolume(leafindex) = ( axi(2,1) - axi(1,1) ) *   &
                                   ( axi(2,2) - axi(1,2) ) *   &
                                   ( axi(2,3) - axi(1,3) )
        endif
     enddo
  elseif(igrid_coord.lt.200) then  
     !
     ! Spherical coordinates
     !
     do ileaf=1,amr_nrleafs
        if(amr_tree_present) then
           leafindex = amr_theleaf_index(ileaf)
           leaf      => amr_theleafs(ileaf)%link
           do iddr=1,3
              axi(1,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr),iddr,leaf%level)
              axi(2,iddr) = amr_finegrid_xi(leaf%ixyzf(iddr)+1,iddr,leaf%level)
           enddo
        else
           leafindex=ileaf
           call amr_regular_get_ixyz(ileaf,ix,iy,iz)
           axi(1,1) = amr_finegrid_xi(ix,1,0)
           axi(2,1) = amr_finegrid_xi(ix+1,1,0)
           axi(1,2) = amr_finegrid_xi(iy,2,0)
           axi(2,2) = amr_finegrid_xi(iy+1,2,0)
           axi(1,3) = amr_finegrid_xi(iz,3,0)
           axi(2,3) = amr_finegrid_xi(iz+1,3,0)
        endif
        cellvolume(leafindex) = (1.d0/3.d0)*                  &
             ( axi(2,1)**3 - axi(1,1)**3 ) *                  &
             ( cos(axi(1,2)) - cos(axi(2,2)) ) *              &
             ( axi(2,3) - axi(1,3) )
     enddo
  else
     stop 823
  endif
  !
end subroutine compute_cellvolumes_amr


!-------------------------------------------------------------------------
!                      COMPUTE DUST MASS
!-------------------------------------------------------------------------
subroutine compute_dust_mass()
  use stars_module
  use amr_module
  implicit none
  integer :: icell,ispec,index
  !
  ! Check if we can use the cellindex array
  !
  if(.not.allocated(cellindex)) then
     write(stdo,*) 'ERROR while computing dust mass: The cellindex list is not yet made.'
     stop
  endif
  !
  ! Reset the dust masses
  !
  do ispec=1,dust_nr_species
     dust_massdust(ispec) = 0.d0
  enddo
  !
  ! Now do a loop over the leafs and dust species
  !
  do icell=1,nrcells
     index = cellindex(icell)
     do ispec=1,dust_nr_species
        dust_massdust(ispec) = dust_massdust(ispec) +                 &
             cellvolume(index)*dustdens(ispec,index)
     enddo
  enddo
  !
  ! If mirror symmetry is used, then multiply the masses by 2
  !
  if(amrray_mirror_equator) then
     dust_massdust(:) = dust_massdust(:) * 2
  endif
  !
  ! Now compute total dust mass
  ! 
  dust_massdusttot = 0.d0
  do ispec=1,dust_nr_species
     dust_massdusttot = dust_massdusttot + dust_massdust(ispec)
  enddo
  !
end subroutine compute_dust_mass


!-------------------------------------------------------------------------
!                        WRITE DUST DENSITY
! Note: this file has .inp instead of .dat, because typically the dust
! density is a quantity that is inserted into RADMC-3D rather than 
! written out. The write_dust_density() routine is, however, useful for
! models made internally with the userdef_module.f90 for debugging or
! model analysis.
!-------------------------------------------------------------------------
subroutine write_dust_density()
  implicit none
  integer(kind=8) :: iiformat,nn,kk
  integer :: index,ierr,ispec,i,precis
  character*160 :: filename
  !
  ! Set the output style and precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Open file and write header
  !
  if(rto_style.eq.1) then
     !
     ! Formatted output
     !
     filename = 'dust_density.inp'
     open(unit=1,file=filename)
     write(1,*) 1
     write(1,*) nrcellsinp
     write(1,*) dust_nr_species
  elseif(rto_style.eq.2) then
     !
     ! F77-Unformatted output
     !
     filename = 'dust_density.uinp'
     open(unit=1,file=filename,form='unformatted')
     nn=1
     kk=rto_reclen
     write(1) nn,kk
     nn=nrcellsinp
     kk=dust_nr_species
     write(1) nn,kk
  elseif(rto_style.eq.3) then
     !
     ! C-compliant binary output
     !
     filename = 'dust_density.binp'
     open(unit=1,file=filename,status='replace',access='stream')
     nn=1
     kk=precis
     write(1) nn,kk
     nn=nrcellsinp
     kk=dust_nr_species
     write(1) nn,kk
  else
     write(stdo,*) 'ERROR: I/O Style ',rto_style,' unknown'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Writing dust densities to ',trim(filename)
  call flush(stdo)
  !
  ! Write the dust density
  !
  do ispec=1,dust_nr_species
     call write_scalarfield(1,rto_style,precis,nrcellsinp,        &
                            dust_nr_species,1,ispec,1,rto_reclen, &
                            scalar1=dustdens)
  enddo
  !
  ! Close file
  !
  close(1)
  !
end subroutine write_dust_density

!-------------------------------------------------------------------------
!                       WRITE DUST TEMPERATURE
!-------------------------------------------------------------------------
subroutine write_dust_temperature()
  implicit none
  integer :: ierror,ispec,index,inu,irec,i,ierr,precis
  integer(kind=8) :: nn,kk
  logical :: fex
  double precision, allocatable :: data(:)
  !
  ! Determine the precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Now write the dust temperature
  !
  if(igrid_type.lt.100) then
     !
     ! Regular (AMR) grid
     ! 
     ! Just make sure that the cell list is complete
     !
     if(amr_tree_present) then
        call amr_compute_list_all()
     endif
     !
     ! Do a stupidity check
     !
     if(nrcells.ne.amr_nrleafs) stop 3209
     !
     ! Open file and write the temperature to it
     !
     if(rto_style.eq.1) then 
        !
        ! Write the temperatures in ascii form
        !
        open(unit=1,file='dust_temperature.dat')
        write(1,*) 1                                   ! Format number
        write(1,*) nrcellsinp
        write(1,*) dust_nr_species
     elseif(rto_style.eq.2) then
        !
        ! Write the temperatures in f77-style unformatted form,
        ! using a record length given by rto_reclen
        !
        open(unit=1,file='dust_temperature.udat',form='unformatted')
        nn = 1
        kk = rto_reclen
        write(1) nn,kk               ! Format number and record length
        nn = nrcellsinp
        kk = dust_nr_species
        write(1) nn,kk
     elseif(rto_style.eq.3) then
        !
        ! C-compatible binary style output
        !
        open(unit=1,file='dust_temperature.bdat',status='replace',access='stream')
        nn = 1
        kk = precis
        write(1) nn,kk                ! Format number and precision
        nn = nrcellsinp
        kk = dust_nr_species
        write(1) nn,kk
     else
        write(stdo,*) 'ERROR: Do not know I/O style ',rto_style
        stop
     endif
     !
     ! Write the data
     !
     do ispec=1,dust_nr_species
        call write_scalarfield(1,rto_style,precis,nrcellsinp,dust_nr_species,1, &
                               ispec,1,rto_reclen,scalar1=dusttemp)
     enddo
     !
     ! Close file
     !
     close(1)
     !
  else
     !
     ! Other grids not yet implemented
     !
     write(stdo,*) 'ERROR: Only regular and AMR grids implemented'
     stop
  endif
  !
  ! If the grid is internally made, then we must make sure that
  ! the grid has been written to file, otherwise the output file
  ! created here makes no sense.
  !
  if((.not.grid_was_read_from_file).and.(.not.grid_was_written_to_file)) then
     call write_grid_file()
     grid_was_written_to_file = .true.     ! Avoid multiple writings
  endif
  !
end subroutine write_dust_temperature


!-------------------------------------------------------------------------
!                    WRITE VELOCITY FIELD
! Note: this file has .inp instead of .dat, because typically the velocity
! field is a quantity that is inserted into RADMC-3D rather than 
! written out. The write_velocityfield() routine is, however, 
! useful for models made internally with the userdef_module.f90 for 
! debugging or model analysis.
!-------------------------------------------------------------------------
subroutine write_velocityfield()
  implicit none
  integer(kind=8) :: iiformat,nn,kk
  integer :: index,ierr,i,precis
  character*160 :: filename
  !
  ! Set the output style and precision
  !
  if(rto_single) then
     precis = 4
  else
     precis = 8
  endif
  !
  ! Open file and write header
  !
  if(rto_style.eq.1) then
     !
     ! Formatted output
     !
     filename = 'gas_velocity.inp'
     open(unit=1,file=filename)
     write(1,*) 1
     write(1,*) nrcellsinp
  elseif(rto_style.eq.2) then
     !
     ! F77-Unformatted output
     !
     filename = 'gas_velocity.uinp'
     open(unit=1,file=filename,form='unformatted')
     nn=1
     kk=rto_reclen
     write(1) nn,kk
     nn=nrcellsinp
     write(1) nn
  elseif(rto_style.eq.3) then
     !
     ! C-compliant binary output
     !
     filename = 'gas_velocity.binp'
     open(unit=1,file=filename,status='replace',access='stream')
     nn=1
     kk=precis
     write(1) nn,kk
     nn=nrcellsinp
     write(1) nn
  else
     write(stdo,*) 'ERROR: I/O Style ',rto_style,' unknown'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Writing gas velocity field ',trim(filename)
  call flush(stdo)
  !
  ! Write the molecule number density
  !
  call write_vectorfield(1,rto_style,precis,3,3,nrcellsinp,1,1, &
                         vector0=gasvelocity)
  !
  ! Close file
  !
  close(1)
  !
end subroutine write_velocityfield



!--------------------------------------------------------------------------
!                            Main subbox routine 
!--------------------------------------------------------------------------
subroutine do_subbox_extract_and_write(nv,nc,subboxfilename,           &
     nx,ny,nz,x0,x1,y0,y1,z0,z1,phi1,theta,phi2,func,identities)
  implicit none
  integer :: nc,nx,ny,nz,nv,ix,iy,iz,iv
  integer, optional :: identities(nv)
  double precision :: x0,x1,y0,y1,z0,z1,phi1,theta,phi2
  double precision :: func(nv,nc)
  double precision :: funcslice(nx,ny,nz,nv)
  character*160 :: subboxfilename
  !
  ! Get the sub box
  !
  if(igrid_type.lt.100) then
     if(igrid_coord.lt.100) then
        call amr_subbox_3d(nv,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,           &
                           phi1,theta,phi2,func,funcslice)
     elseif(igrid_coord.lt.200) then
        call amr_subbox_3d(nv,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,           &
                  phi1,theta,phi2,func,funcslice,spherical=.true.)
     else
        stop 5092
     endif
  else
     write(stdo,*) 'ERROR: No other grid types than regular are allowed yet.'
     stop
  endif
  !
  ! Now write to a file
  !
  open(unit=1,file=trim(subboxfilename))
  write(1,*) 2                         ! Format number
  write(1,*) nx,ny,nz,nv               ! Dimensions
  write(1,*) x0,x1,y0,y1,z0,z1         ! Size and position of box
  write(1,*) phi1,theta,phi2           ! Orientation of the box
  write(1,*) 
  if(present(identities)) then
     write(1,*) (identities(iv),iv=1,nv)
  else
     write(1,*) (iv,iv=1,nv)
  endif
  write(1,*) 
  do iv=1,nv
     do iz=1,nz
        do iy=1,ny
           do ix=1,nx
              write(1,*) funcslice(ix,iy,iz,iv)
           enddo
        enddo
     enddo
  enddo
  close(1)
end subroutine do_subbox_extract_and_write


!--------------------------------------------------------------------------
!                            Main subbox routine 
!--------------------------------------------------------------------------
subroutine do_subbox_levelpop_extract_and_write(nvm1,nv1,nv2,nc,ispec,     &
     subset,subboxfilename,nx,ny,nz,x0,x1,y0,y1,z0,z1,phi1,theta,phi2,func)
  implicit none
  integer :: nc,nx,ny,nz,nvm1,nv1,nv2,ix,iy,iz,iv1,ispec
  integer :: subset(nv1)
  double precision :: x0,x1,y0,y1,z0,z1,phi1,theta,phi2
  double precision :: func(nvm1,nv2,nc)
  double precision :: funcslice(nx,ny,nz,nvm1,nv2)
  character*160 :: subboxfilename
  !
  ! Get the sub box
  !
  if(igrid_type.lt.100) then
     if(igrid_coord.lt.100) then
        call amr_subbox2_3d(nvm1,nv2,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,           &
                           phi1,theta,phi2,func,funcslice)
     elseif(igrid_coord.lt.200) then
        call amr_subbox2_3d(nvm1,nv2,nc,nx,ny,nz,x0,x1,y0,y1,z0,z1,           &
                  phi1,theta,phi2,func,funcslice,spherical=.true.)
     else
        stop 5092
     endif
  else
     write(stdo,*) 'ERROR: No other grid types than regular are allowed yet.'
     stop
  endif
  !
  ! Now write to a file
  !
  open(unit=1,file=trim(subboxfilename))
  write(1,*) 2                         ! Format number
  write(1,*) nx,ny,nz,nv1              ! Dimensions
  write(1,*) x0,x1,y0,y1,z0,z1         ! Size and position of box
  write(1,*) phi1,theta,phi2           ! Orientation of the box
  write(1,*) 
  write(1,*) (subset(iv1),iv1=1,nv1)
  do iv1=1,nv1
     write(1,*) 
     do iz=1,nz
        do iy=1,ny
           do ix=1,nx
              write(1,*) funcslice(ix,iy,iz,iv1,ispec)
           enddo
        enddo
     enddo
  enddo
  close(1)
end subroutine do_subbox_levelpop_extract_and_write


!--------------------------------------------------------------------------
!                           READ SAMPLE POINTS
!--------------------------------------------------------------------------
subroutine read_sampling_points_3d(npt,xpt,ypt,zpt)
  implicit none
  integer :: npt,iformat,ierr,ipt
  double precision, pointer :: xpt(:),ypt(:),zpt(:)
  double precision :: x,y,z
  logical :: fex
  !
  ! Check if we must first deallocate
  !
  if(associated(xpt)) deallocate(xpt)
  if(associated(ypt)) deallocate(ypt)
  if(associated(zpt)) deallocate(zpt)
  !
  ! Now check if a sampling point file exists
  !
  inquire(file='sample_points.inp',exist=fex)
  if(fex) then
     open(unit=1,file='sample_points.inp')
     read(1,*) iformat
     if(iformat.ne.1) then 
        write(stdo,*) 'ERROR: Format number of sample_points.inp not known'
        stop
     endif
     read(1,*) npt
     if(npt.le.0) stop
     allocate(xpt(1:npt),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in IOput Module: could not allocate xpt'
        stop 
     endif
     allocate(ypt(1:npt),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in IOput Module: could not allocate ypt'
        stop 
     endif
     allocate(zpt(1:npt),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in IOput Module: could not allocate zpt'
        stop 
     endif
     do ipt=1,npt
        read(1,*) x,y,z
        xpt(ipt) = x
        ypt(ipt) = y
        zpt(ipt) = z
     enddo
     close(1)
  endif
end subroutine read_sampling_points_3d

!--------------------------------------------------------------------------
!                            Main sample routine 
!--------------------------------------------------------------------------
subroutine do_sample_extract_and_write(nv,nc,samplefilename,           &
     npt,xpt,ypt,zpt,func)
  implicit none
  integer :: nc,npt,nv,ipt,iv
  double precision :: func(nv,nc)
  double precision :: values(npt,nv)
  double precision, pointer :: xpt(:),ypt(:),zpt(:)
  character*160 :: samplefilename
  !
  ! Check if the sampling point grid exists
  !
  if((.not.associated(xpt)).or.(.not.associated(ypt))             &
       .or.(.not.associated(zpt)).or.(npt.le.0)) then
     write(stdo,*) 'ERROR in sampling: Cannot sample if no sampling points are given.'
     open(unit=1,file=trim(samplefilename))
     write(1,*) 2                         ! Format number
     write(1,*) 0,0                       ! Dimensions
     close(1)
     return
  endif
  !
  ! Get the sampled values
  !
  if(igrid_type.lt.100) then
     if(igrid_coord.lt.100) then
        call amr_sample_3d(nv,nc,npt,xpt,ypt,zpt,func,values)
     elseif(igrid_coord.lt.200) then
        call amr_sample_3d(nv,nc,npt,xpt,ypt,zpt,func,values,spherical=.true.)
     else
        stop 5192
     endif
  else
     write(stdo,*) 'ERROR: No other grid types than regular are allowed yet.'
     stop
  endif
  !
  ! Now write to a file
  !
  open(unit=1,file=trim(samplefilename))
  write(1,*) 2                         ! Format number
  write(1,*) npt,nv                    ! Dimensions
  write(1,*) 
  write(1,*) (iv,iv=1,nv)
  do iv=1,nv
     write(1,*) 
     do ipt=1,npt
        write(1,*) values(ipt,iv)
     enddo
  enddo
  close(1)
end subroutine do_sample_extract_and_write

!--------------------------------------------------------------------------
!                            Main sample routine 
!--------------------------------------------------------------------------
subroutine do_sample_levelpop_extract_and_write(nvm1,nv1,nv2,nc,ispec, &
     subset,samplefilename,npt,xpt,ypt,zpt,func)
  implicit none
  integer :: nc,npt,nv1,nv2,ipt,iv1,ispec,nvm1
  integer :: subset(nv1)
  double precision :: func(nvm1,nv2,nc)
  double precision :: values(npt,nvm1,nv2)
  double precision, pointer :: xpt(:),ypt(:),zpt(:)
  character*160 :: samplefilename
  !
  ! Check if the sampling point grid exists
  !
  if((.not.associated(xpt)).or.(.not.associated(ypt))             &
       .or.(.not.associated(zpt)).or.(npt.le.0)) then
     write(stdo,*) 'ERROR in sampling: Cannot sample if no sampling points are given.'
     open(unit=1,file=trim(samplefilename))
     write(1,*) 2                         ! Format number
     write(1,*) 0,0                       ! Dimensions
     close(1)
     return
  endif
  !
  ! Get the sampled values
  !
  if(igrid_type.lt.100) then
     if(igrid_coord.lt.100) then
        call amr_sample2_3d(nvm1,nv2,nc,npt,xpt,ypt,zpt,func,values)
     elseif(igrid_coord.lt.200) then
        call amr_sample2_3d(nvm1,nv2,nc,npt,xpt,ypt,zpt,func,values,spherical=.true.)
     else
        stop 5192
     endif
  else
     write(stdo,*) 'ERROR: No other grid types than regular are allowed yet.'
     stop
  endif
  !
  ! Now write to a file
  !
  open(unit=1,file=trim(samplefilename))
  write(1,*) 2                          ! Format number
  write(1,*) npt,nv1                    ! Dimensions
  write(1,*) 
  write(1,*) (subset(iv1),iv1=1,nv1)
  do iv1=1,nv1
     write(1,*) 
     do ipt=1,npt
        write(1,*) values(ipt,iv1,ispec)
     enddo
  enddo
  close(1)
end subroutine do_sample_levelpop_extract_and_write


!--------------------------------------------------------------------------
!                  DESTROY ANY REMAINING SAMPLING POINTS
!--------------------------------------------------------------------------
subroutine destroy_sampling_points(xpt,ypt,zpt)
  implicit none
  double precision, pointer :: xpt(:),ypt(:),zpt(:)
  if(associated(xpt)) deallocate(xpt)
  if(associated(ypt)) deallocate(ypt)
  if(associated(zpt)) deallocate(zpt)
end subroutine destroy_sampling_points


end module ioput_module
