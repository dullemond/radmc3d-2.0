module sources_module
  use rtglobal_module
  use amrray_module
  use dust_module
  use lines_module
  use gascontinuum_module
  use stars_module
  use ioput_module
  use userdef_module
  use mathroutines_module
  use montecarlo_module
  use constants_module
  !
  ! Basic variables
  !
  integer :: sources_nrfreq,sources_align_munr
  logical :: sources_secondorder,sources_catch_doppler_line
  logical :: sources_interpol_jnu=.true.
  logical :: sources_localobserver=.false.  
  double precision :: sources_observer_position(1:3)
  !
  ! Some basic arrays
  !
  double precision, allocatable :: sources_frequencies(:)
  double precision, allocatable :: sources_dustdens(:),sources_dusttemp(:)
  double precision, allocatable :: sources_dustkappa_a(:,:),sources_dustkappa_s(:,:)
  double precision, allocatable :: sources_alpha_a(:),sources_alpha_s(:)
  !
  ! Arrays for the aligned grains mode
  !
  double precision, allocatable :: sources_align_mu(:)
  double precision, allocatable :: sources_align_orth(:,:,:),sources_align_para(:,:,:)
  !
  !    Arrays for templates of stellar source spectra
  !
  double precision, allocatable :: sources_stellarsrc_templates(:,:)
  !
  !    Arrays for the second order integration with corner vertices
  !
  !    NOTE: snu either means S_nu (source function) or j_nu (emissivity)
  !          depending on the setting of sources_interpol_jnu.
  !
  double precision, allocatable :: sources_vertex_snu(:,:)
  double precision, allocatable :: sources_vertex_anu(:)
  double precision, allocatable :: sources_cell_snu(:,:)
  double precision, allocatable :: sources_cell_anu(:)
  !
  !    Arrays for doppler catching algorithm
  !
  double precision, allocatable :: sources_vertex_doppler(:)
  double precision, allocatable :: sources_vertex_turb(:)
  double precision, allocatable :: sources_vertex_temp(:)
  double precision, allocatable :: sources_vertex_line_nup(:,:)
  double precision, allocatable :: sources_vertex_line_ndown(:,:)
  double precision, allocatable :: sources_cell_doppler(:)
!  double precision, allocatable :: sources_cell_turb(:)
!  double precision, allocatable :: sources_cell_temp(:)
  double precision, allocatable :: sources_cell_line_nup(:,:)
  double precision, allocatable :: sources_cell_line_ndown(:,:)
  integer :: sources_vertex_lines_nractivetot=0
  !
  double precision :: sources_local_doppler_curr,sources_local_turb_curr,sources_local_temp_curr
  double precision, allocatable :: sources_local_line_nup_curr(:)
  double precision, allocatable :: sources_local_line_ndown_curr(:)
  !
  double precision :: sources_local_doppler_prev,sources_local_turb_prev,sources_local_temp_prev
  double precision, allocatable :: sources_local_line_nup_prev(:)
  double precision, allocatable :: sources_local_line_ndown_prev(:)
  !
  double precision :: sources_local_doppler_end,sources_local_turb_end,sources_local_temp_end
  double precision, allocatable :: sources_local_line_nup_end(:)
  double precision, allocatable :: sources_local_line_ndown_end(:)

contains

!-------------------------------------------------------------------
!                   INITIALIZE THE SOURCES MODULE
!-------------------------------------------------------------------
subroutine sources_init(nrfreq,frequencies,secondorder,doppcatch)
  implicit none
  integer :: nrfreq
  logical :: secondorder,doppcatch,havealignopac
  integer :: ierr,inu,ispec,itempl,iinu,imu
  double precision :: temp,eps,orth,para
  double precision :: frequencies(nrfreq)
  !
  ! Memorize nr of frequencies and order
  !
  sources_nrfreq             = nrfreq
  sources_secondorder        = secondorder
  sources_catch_doppler_line = doppcatch
  !
  ! Allocate the sources_frequencies array
  !
  if(allocated(sources_frequencies)) deallocate(sources_frequencies)
  allocate(sources_frequencies(1:sources_nrfreq),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in sources module: Could not allocate sources_frequencies() array'
     stop
  endif
  sources_frequencies(1:nrfreq) = frequencies(1:nrfreq)
  !
  ! Allocate the arrays for the dust 
  !
  if(rt_incl_dust) then
     allocate(sources_dustdens(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_dustdens() array'
        stop
     endif
     allocate(sources_dusttemp(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_dusttemp() array'
        stop
     endif
     allocate(sources_dustkappa_a(1:nrfreq,1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_dustkappa_a()'
        stop 
     endif
     allocate(sources_dustkappa_s(1:nrfreq,1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_dustkappa_s()'
        stop 
     endif
     allocate(sources_alpha_a(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_alpha_a()'
        stop 
     endif
     allocate(sources_alpha_s(1:dust_nr_species),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate sources_alpha_s()'
        stop 
     endif
     !
     ! Get the dust opacities 
     !
     temp    = 100.d0     ! Arbitrary temperature...
     do inu=1,nrfreq
        do ispec=1,dust_nr_species
           sources_dustkappa_a(inu,ispec) =    &
                find_dust_kappa_interpol(sources_frequencies(inu),ispec,temp,1,0,0)
           sources_dustkappa_s(inu,ispec) =    &
                find_dust_kappa_interpol(sources_frequencies(inu),ispec,temp,0,1,0)
        enddo
     enddo
     !
     ! If scattering_mode is 0, then put kappa_s to zero
     ! 
     if(scattering_mode.eq.0) then
!        write(stdo,*) 'In Sources initialization: found that scattering_mode==0, so putting kappa_scat=0.'
        do inu=1,sources_nrfreq
           do ispec=1,dust_nr_species
              sources_dustkappa_s(inu,ispec) = 0.d0
           enddo
        enddo
     endif
     !
     ! Grain alignment mode
     !
     if(alignment_mode.ne.0) then
        !
        ! Do a few checks for grain alignment mode
        !
        if(rt_incl_lines) then
           write(stdo,*) 'ERROR: If alignment_mode is switched on, then gas lines are incompatible...'
           stop
        endif
        if(.not.allocated(grainalign_dir)) then
           write(stdo,*) 'ERROR: If alignment_mode.ne.0 then the alignment vector ', &
                'field must be read. Appears to be a bug in the code.'
           stop
        endif
        if(align_munr.eq.0) then
           write(stdo,*) 'ERROR: If alignment_mode.ne.0 then align_munr ', &
                'must be larger than 0. Appears to be a bug in the code.'
           stop
        endif
        if(.not.allocated(align_mui_grid)) then
           write(stdo,*) 'ERROR: If alignment_mode.ne.0 then a global align_mui_grid ', &
                'must be present. Appears to be a bug in the code.'
           stop
        endif
        havealignopac = .false.
        do ispec=1,dust_nr_species
           if(associated(dust_kappa_arrays(ispec)%alignmu)) then
              havealignopac = .true.
           endif
        enddo
        if(.not.havealignopac) then
           write(stdo,*) 'ERROR: If alignment_mode.ne.0 then at least one of the ', &
                'dust opacities must be accompanied by a file dustkapalignfact_*.inp ', &
                'which sets the strength and shape of the alignment effect.'
           stop
        endif
        !
        ! Allocate the global alignfact arrays
        !
        sources_align_munr = align_munr
        if(allocated(sources_align_mu)) deallocate(sources_align_mu)
        if(allocated(sources_align_orth)) deallocate(sources_align_orth)
        if(allocated(sources_align_para)) deallocate(sources_align_para)
        allocate(sources_align_mu(sources_align_munr))
        allocate(sources_align_orth(sources_align_munr,sources_nrfreq,dust_nr_species))
        allocate(sources_align_para(sources_align_munr,sources_nrfreq,dust_nr_species))
        !
        ! Now fill these arrays
        !
        sources_align_mu(:) = align_mui_grid(:)
        do ispec=1,dust_nr_species
           do inu=1,nrfreq
              do imu=1,sources_align_munr
                 call find_dust_alignfact_interpol(sources_frequencies(inu), &
                      sources_align_mu(imu),ispec,1,0,.true.,orth,para)
                 sources_align_orth(imu,inu,ispec) = orth
                 sources_align_para(imu,inu,ispec) = para
              enddo
           enddo
        enddo
     endif
  endif
  !
  ! If we wish to do second order integration of the formal transfer equation,
  ! then we must set up the vertex grid (corner-based grid), or at least
  ! check that it is set up.
  !
  if(sources_secondorder) then
     !
     ! Set up the corner vertex grid
     !
     if(.not.allocated(amr_vertex_cells)) then
        write(stdo,*) 'Setting up the vertices (cell corners) for second order transfer...'
        call amr_set_up_vertices()
        !#########################
        !#########################
        ! MAYBE THERE IS A BUG HERE: WE MUST MAKE A SPECIAL TREATMENT OF THE MIDPLANE IF MIRROR SYMMETRY IS PRESENT.
        ! THIS CAN ONLY EASILY BE FIXED INSIDE AMR_SET_UP_VERTICES VIA AN OPTIONAL ARGUMENT
        ! HOWEVER, THIS IS ONLY A MINOR PROBLEM, SINCE SRC/ALP REMAINS CORRECT; ONLY THE OPACITY OF THE CELL IS SMALLER
        !#########################
        !#########################
     endif
     !
     ! Allocate the j_nu and alpha_nu arrays at the vertices. Note that
     ! it would take presumably too much memory to store the entire
     ! jnu(inu,ivert) array. So we store just one frequency at a time.
     ! This means also that this works only in the one-wavelength-at-
     ! a-time raytracing mode.
     !
     if(.not.allocated(sources_vertex_snu))                    &
          allocate(sources_vertex_snu(1:amr_nr_vertices_max,1:4))
     if(.not.allocated(sources_vertex_anu))                    &
          allocate(sources_vertex_anu(1:amr_nr_vertices_max))
     if(.not.allocated(sources_cell_snu))                      &
          allocate(sources_cell_snu(1:amr_nr_vertices_max,1:4))
     if(.not.allocated(sources_cell_anu))                      &
          allocate(sources_cell_anu(1:amr_nr_vertices_max))
     !
     ! If we want to do doppler-catching of lines, then we must allocate
     ! further arrays
     !
     if(sources_catch_doppler_line) then
        if(.not.allocated(sources_vertex_doppler))                    &
             allocate(sources_vertex_doppler(1:amr_nr_vertices_max))
        if(.not.allocated(sources_cell_doppler))                    &
             allocate(sources_cell_doppler(1:amr_nrleafs_max))
        if(.not.allocated(sources_vertex_turb))                    &
             allocate(sources_vertex_turb(1:amr_nr_vertices_max))
!        if(.not.allocated(sources_cell_turb))                    &
!             allocate(sources_cell_turb(1:amr_nrleafs_max))
        if(.not.allocated(sources_vertex_temp))                    &
             allocate(sources_vertex_temp(1:amr_nr_vertices_max))
!        if(.not.allocated(sources_cell_temp))                    &
!             allocate(sources_cell_temp(1:amr_nrleafs_max))
     endif
     !
  endif
  !
  ! Allocate arrays for the stellar source templates, and
  ! map the templates onto these frequencies
  !
  if(stellarsrc_nrtemplates.gt.0) then
     allocate(sources_stellarsrc_templates(1:sources_nrfreq,1:stellarsrc_nrtemplates),STAT=ierr)
     if(ierr.ne.0) then
        write(stdo,*) 'ERROR in sources module: Could not allocate stellarsrc templates.'
        stop
     endif
     do itempl=1,stellarsrc_nrtemplates
        do inu=1,sources_nrfreq
           call hunt(freq_nu,freq_nr,sources_frequencies(inu),iinu)
           if((iinu.lt.1).or.(iinu.ge.freq_nr)) then
              sources_stellarsrc_templates(inu,itempl) = 0.d0
           else
              eps = ( sources_frequencies(inu) - freq_nu(iinu) ) / &
                    ( freq_nu(iinu+1) - freq_nu(iinu) )
              if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 674
              sources_stellarsrc_templates(inu,itempl) = &
                   (1.d0-eps)*stellarsrc_templates(iinu,itempl) + &
                          eps*stellarsrc_templates(iinu+1,itempl)
           endif
        enddo
     enddo
  endif
end subroutine sources_init

!-------------------------------------------------------------------
!                     PARTIAL CLEANUP SOURCES
!-------------------------------------------------------------------
subroutine sources_partial_cleanup()
  implicit none
  if(allocated(sources_dustdens)) deallocate(sources_dustdens)
  if(allocated(sources_dusttemp)) deallocate(sources_dusttemp)
  if(allocated(sources_dustkappa_a)) deallocate(sources_dustkappa_a)
  if(allocated(sources_dustkappa_s)) deallocate(sources_dustkappa_s)
  if(allocated(sources_alpha_a)) deallocate(sources_alpha_a)
  if(allocated(sources_alpha_s)) deallocate(sources_alpha_s)
  if(allocated(sources_align_mu)) deallocate(sources_align_mu)
  if(allocated(sources_align_orth)) deallocate(sources_align_orth)
  if(allocated(sources_align_para)) deallocate(sources_align_para)
end subroutine sources_partial_cleanup

!-------------------------------------------------------------------
!                         CLEANUP SOURCES
!-------------------------------------------------------------------
subroutine sources_cleanup()
  implicit none
  !
  ! Call the partial cleanup
  !
  call sources_partial_cleanup()
  !
  ! Deallocate the frequency array
  !
  if(allocated(sources_frequencies)) deallocate(sources_frequencies)
  sources_nrfreq = 0
  !
  ! ...Any vertex grid arrays for the sources module
  !
  if(allocated(sources_vertex_snu)) deallocate(sources_vertex_snu)
  if(allocated(sources_vertex_anu)) deallocate(sources_vertex_anu)
  if(allocated(sources_cell_snu)) deallocate(sources_cell_snu)
  if(allocated(sources_cell_anu)) deallocate(sources_cell_anu)
  !
  ! ...Any doppler catching stuff
  !
  if(allocated(sources_vertex_doppler)) deallocate(sources_vertex_doppler)
  if(allocated(sources_cell_doppler)) deallocate(sources_cell_doppler)
  if(allocated(sources_vertex_turb)) deallocate(sources_vertex_turb)
!  if(allocated(sources_cell_turb)) deallocate(sources_cell_turb)
  if(allocated(sources_vertex_temp)) deallocate(sources_vertex_temp)
!  if(allocated(sources_cell_temp)) deallocate(sources_cell_temp)
  if(allocated(sources_vertex_line_nup)) deallocate(sources_vertex_line_nup)
  if(allocated(sources_vertex_line_ndown)) deallocate(sources_vertex_line_ndown)
  if(allocated(sources_cell_line_nup)) deallocate(sources_cell_line_nup)
  if(allocated(sources_cell_line_ndown)) deallocate(sources_cell_line_ndown)
  if(allocated(sources_local_line_nup_curr)) deallocate(sources_local_line_nup_curr)
  if(allocated(sources_local_line_ndown_curr)) deallocate(sources_local_line_ndown_curr)
  if(allocated(sources_local_line_nup_prev)) deallocate(sources_local_line_nup_prev)
  if(allocated(sources_local_line_ndown_prev)) deallocate(sources_local_line_ndown_prev)
  if(allocated(sources_local_line_nup_end)) deallocate(sources_local_line_nup_end)
  if(allocated(sources_local_line_ndown_end)) deallocate(sources_local_line_ndown_end)
  !
  ! Stellar source templates
  !
  if(allocated(sources_stellarsrc_templates)) deallocate(sources_stellarsrc_templates)
end subroutine sources_cleanup

!-------------------------------------------------------------------
!           COMPUTE AND STORE POPULATIONS IN GLOBAL ARRAY
!-------------------------------------------------------------------
subroutine lines_compute_and_store_local_populations(action)
  implicit none
  integer :: action,index,ierr,ilevsub,ilevel,ispec
  integer, allocatable :: ilevmax_15(:)
  type(amr_branch), pointer :: a
  double precision, allocatable :: pop(:,:)
  double precision :: velgrad,xav,yav,zav,turb
  integer :: ixx,iyy,izz,idum
  character*80 :: sispec,slevel
  !
  ! Action=0 means do nothing, action=1 means compute if not yet computed,
  ! action=2 means re-compute.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(lines_levelpop)) return
  endif
  !
  ! Check
  !
  if(lines_nrlevels_subset_max.le.0) then
!     write(stdo,*) 'ERROR: lines_nrlevels_subset_max is zero...'
!     write(stdo,*) '       Cannot allocate level population array...'
!     stop
     return
  endif
  !
  ! Allocate local population array
  !
  allocate(pop(1:lines_maxnrlevels,1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate local pop() array'
     stop
  endif
  allocate(ilevmax_15(1:lines_nr_species),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate local ilevmax_15() array'
     stop
  endif
  !
  ! Remove any already existing array
  !
  if(allocated(lines_levelpop)) deallocate(lines_levelpop)
  !
  ! Create the level populations array
  !
  allocate(lines_levelpop(1:lines_nrlevels_subset_max,1:lines_nr_species,1:nrcells),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate the lines_levelpop() array'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Computing level populations, and storing them in global array...'
  call flush(stdo)
  !
  ! Checks
  !
  if(lines_mode.eq.1) then 
     if(.not.allocated(lines_partition_temperature)) then
        write(stdo,*) 'ERROR: wanting to compute level populations '
        write(stdo,*) '       according to LTE prescription, but'
        write(stdo,*) '       partition function array not allocated.'
        stop 
     endif
  endif
  if(.not.allocated(gastemp)) then
     write(stdo,*) 'ERROR: Wanting to compute level populations '
     write(stdo,*) '       according to LTE prescription, but'
     write(stdo,*) '       gastemp array not allocated.'
     stop
  endif
  if(.not.allocated(gas_chemspec_numberdens)) then
     write(stdo,*) 'ERROR: Wanting to compute level populations '
     write(stdo,*) '       according to LTE prescription, but'
     write(stdo,*) '       gas_chemspec_numberdens array not allocated.'
     stop
  endif
  if(lines_mode.ge.3) then
     if(.not.allocated(collpartner_numberdens)) then
        write(stdo,*) 'ERROR: Wanting to compute non-LTE level populations but'
        write(stdo,*) '       collpartner_numberdens(:,:) array not allocated.'
        stop
     endif
  endif
  !
  ! For being able to report later how many iterations for LVG or
  ! escape probability on average are used, reset a counter
  !
  lines_lvgesc_nriter   = 0.d0
  lines_lvgesc_nrsolves = 0.d0
  !
  ! Reset the maser warning flag
  !
  lines_maser_warning = .false.
  !
  ! For being able to check if we waste CPU time on always-empty
  ! levels, we reset the ilevmax_15
  !
  ilevmax_15(:) = 1
  !
  ! Now do the big loop over all cells
  !
  call amr_resetcount()
  call amr_nextcell(index)
  do while(index.ge.0)
     if(index.gt.0) then
        !
        ! Compute the local level populations
        !
        if(lines_mode.eq.1) then
           !
           ! LTE populations
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 call lines_compute_ltepop(ispec,             &
                      lines_nrlevels(ispec),gastemp(index),   &
                      gas_chemspec_numberdens(ispec,index),   &
                      pop(:,ispec))
              endif
           enddo
        elseif(lines_mode.eq.2) then
           !
           ! User-defined populations
           !
           ! Find the mid-cell position
           !
           if(amr_tree_present) then
              !
              ! If the AMR tree is available, then use that 
              !
              a => amr_index_to_leaf(index)%link
              xav = amr_finegrid_xc(a%ixyzf(1),1,a%level)
              yav = amr_finegrid_xc(a%ixyzf(2),2,a%level)
              zav = amr_finegrid_xc(a%ixyzf(3),3,a%level)
           else
              !
              ! For a regular grid, use ixx,iyy,izz
              !
              call amr_regular_get_vertex_ixyz(index,ixx,iyy,izz)
              xav = amr_finegrid_xc(ixx,1,0)
              yav = amr_finegrid_xc(iyy,2,0)
              zav = amr_finegrid_xc(izz,3,0)
           endif
           !
           ! User-defined population calculation
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 call userdef_compute_levelpop(ispec,            &
                      lines_nrlevels(ispec),index,               &
                      xav,yav,zav,                               &
                      gas_chemspec_numberdens(ispec,index),      &
                      pop(:,ispec))
              endif
           enddo
        elseif(lines_mode.eq.3) then
           !
           ! LVG (Sobolev) populations
           !
           ! First determine the velocity gradient
           !
           call lines_compute_velgradient(index,velgrad)
           !
           ! Now compute, with this velocity gradient, the
           ! LVG populations.
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 if(lines_collisions_npartners(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no collision rates are loaded for'
                    write(stdo,*) '        this molecule, ispec=',ispec
                    stop
                 endif
                 if(lines_collisions_npartners_used(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no number densities of collision partners have been'
                    write(stdo,*) '        specified for this molecule, ispec=',ispec
                    stop
                 endif
                 !
                 ! Decide whether or not to include Escape Probability in
                 ! the LVG computation (i.e. a max length scale to assure
                 ! that for low velocity gradient the optical depth does
                 ! not go to infinity).
                 !
                 if(allocated(lines_escprob_lengthscale)) then
                    if(allocated(lines_microturb)) then
                       turb = lines_microturb(index)
                    else
                       turb = 0.d0
                    endif
                    call lines_compute_lvgpop(ispec,lines_nrlevels(ispec),         &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(index),                    &
                      gas_chemspec_numberdens(ispec,index),                        &
                      collpartner_numberdens(:,index),velgrad,                     &
                      pop(:,ispec),lengthscale=lines_escprob_lengthscale(index),   &
                      turb=turb)
                 else
                    call lines_compute_lvgpop(ispec,lines_nrlevels(ispec),         &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(index),                    &
                      gas_chemspec_numberdens(ispec,index),                        &
                      collpartner_numberdens(:,index),velgrad,                     &
                      pop(:,ispec))
                 endif
              endif
           enddo
        elseif(lines_mode.eq.4) then
           !
           ! Optically thin non-LTE populations
           !
           ! Now compute the populations
           ! 
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 if(lines_collisions_npartners(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no collision rates are loaded for'
                    write(stdo,*) '        this molecule, ispec=',ispec
                    stop
                 endif
                 if(lines_collisions_npartners_used(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no number densities of collision partners have been'
                    write(stdo,*) '        specified for this molecule, ispec=',ispec
                    stop
                 endif
                 call lines_compute_optthinpop(ispec,lines_nrlevels(ispec),        &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(index),                    &
                      gas_chemspec_numberdens(ispec,index),                        &
                      collpartner_numberdens(:,index),                             &
                      pop(:,ispec))
              endif
           enddo
        else
           write(stdo,*) 'ERROR: line_mode not known'
           stop
        endif
        !
        ! Now store these populations, but only those of the level
        ! subset. This is meant to save memory. By default, however, the
        ! subset is the full set.
        !
        do ispec=1,lines_nr_species
           if(lines_species_fullmolec(ispec)) then
              do ilevsub=1,lines_nrlevels_subset(ispec)
                 ilevel = lines_levels_subset(ilevsub,ispec)
                 lines_levelpop(ilevsub,ispec,index) = pop(ilevel,ispec)
              enddo
           endif
        enddo
        !
        ! Find the max ilevel that is relatively populated more than 1d-15
        !
        do ispec=1,lines_nr_species
           if(lines_species_fullmolec(ispec)) then
              ilevel = lines_nrlevels(ispec)
              do ilevel=lines_nrlevels(ispec),1,-1
                 if(pop(ilevel,ispec).gt.1d-15*gas_chemspec_numberdens(ispec,index)) exit
              enddo
              if(ilevel.gt.ilevmax_15(ispec)) then
                 ilevmax_15(ispec) = ilevel
              endif
           endif
        enddo
        !
     endif
     call amr_nextcell(index)
  enddo
  !
  ! Message about the nr of iterations
  !
  if(lines_lvgesc_nrsolves.gt.0.d0) then
     lines_lvgesc_nriter = lines_lvgesc_nriter / lines_lvgesc_nrsolves
     write(stdo,*) '---> Average nr of iterations in the LVG/escape-probability solver: ',lines_lvgesc_nriter
  endif
  !
  ! If masers are detected, then make a message
  !
  if(lines_maser_warning) then
     write(stdo,*) 'WARNING: Masering (negative opacity) was detected in the line transfer.'
     write(stdo,*) '         In those cases the opacity was set to 0 in the non-LTE iteration.'
  endif
  !
  ! Message about nr of levels that are non-zero populated
  ! (only for non-LTE stuff, because only then is it useful
  ! to know, because only then one can save substantial
  ! amount of computing time by trimming the molecular
  ! input files
  !
  if(lines_mode.gt.1) then
     do ispec=1,lines_nr_species
        if(lines_species_fullmolec(ispec)) then
           if(ilevmax_15(ispec).lt.lines_nrlevels(ispec)) then
              call integer_to_string(ispec,sispec)
              idum = ilevmax_15(ispec)
              call integer_to_string(idum,slevel)
              write(stdo,*) '---> For molecule/atom '//trim(sispec)//' only levels up to ilevel='//trim(slevel)// &
                   ' are populated above 1E-15 (fractional population).' 
              idum = lines_nrlevels(ispec)-ilevmax_15(ispec)
              call integer_to_string(idum,slevel)
              write(stdo,*) '     This means that the top '//trim(slevel)//' levels are as good as empty everywhere.'
              write(stdo,*) '     This is not a problem, but if the code is slow, maybe it is because'
              write(stdo,*) '     it wastes time doing matrix LU-decomposition for more levels than necessary.'
              write(stdo,*) '     By trimming these top levels from the molecular/atomic data file you could'
              write(stdo,*) '     speed up the code. But beware: if you later have a hotter setup and still'
              write(stdo,*) '     use the trimmed molecular data file, you might get wrong results. So this'
              write(stdo,*) '     speed-up advice is at your own risk!'
           endif
        endif
     enddo
  endif
  !
  ! Deallocate
  !
  deallocate(pop)
  deallocate(ilevmax_15)
  !
end subroutine lines_compute_and_store_local_populations

!-------------------------------------------------------------------------
!                      GET SRC AND ALP LOCALLY 
!-------------------------------------------------------------------------
subroutine sources_get_src_alp(inu0,inu1,nf,src,alp,inclstokes)
  implicit none
  integer :: ispec,inu,itemplate,inu0,inu1,nf,ilevel,ilevsub
  double precision :: dum,dummy,src(nf,1:4),alp(nf)
  double precision :: xav,yav,zav
  double precision :: rlen,cosp,sinp,cost,sint,vx,vy,vz
  double precision :: velgrad,turb
  double precision :: levent,phievent
  integer :: iphievent
  logical :: inclstokes
  !
  ! Find the local dust densities
  !
  if(rt_incl_dust) then
     if(ray_index.ge.1) then
        sources_dustdens(:) = dustdens(:,ray_index)
        sources_dusttemp(:) = dusttemp(:,ray_index)
     else
        sources_dustdens(:) = 0.d0
        sources_dusttemp(:) = 1.d0
     endif
  endif
  !
  ! Find the local line-related quantities
  !
  if(rt_incl_lines) then
     if(ray_index.ge.1) then
        if(debug_check_all.eq.1) then
           if(.not.allocated(lines_ray_levpop)) stop 3023
           if(.not.allocated(gas_chemspec_numberdens)) stop 3024
           if((lines_mode.gt.0).and..not.allocated(lines_levelpop)) stop 3025
        endif
        !
        ! Compute the local level populations
        !
        if(lines_mode.gt.0) then
           !
           ! Get the level populations from the big stored array,
           ! meaning that the (should) have been calculated 
           ! beforehand.
           !
           ! NOTE: Only the levels from the subset that is stored
           !       are copied. The populations other levels remain
           !       undefined. 
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then           
                 !
                 !##########################################################
                 ! NOTE: For testing purposes (can later be removed)
                 !       we will put these unused level populations
                 !       to 1D90, to enforce an error if they are
                 !       accidently used.
                 !
                 lines_ray_levpop(:,ispec,1) = 1d90
                 !##########################################################
                 !
                 ! Now copy the levels from the big array to the local
                 ! array
                 !
                 do ilevsub=1,lines_nrlevels_subset(ispec)
                    ilevel = lines_levels_subset(ilevsub,ispec)
                    lines_ray_levpop(ilevel,ispec,1) =            &
                         lines_levelpop(ilevsub,ispec,ray_index)
                 enddo
              endif
           enddo
           !
        elseif(lines_mode.eq.-1) then
           !
           ! On-the-fly LTE population calculation
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 call lines_compute_ltepop_subset(lines_nrlevels(ispec),    &
                      active_nrlevels(ispec),                               &
                      active_levels(:,ispec),                               &
                      ispec,gastemp(ray_index),                             &
                      gas_chemspec_numberdens(ispec,ray_index),             &
                      lines_ray_levpop(:,ispec,1))
              endif
           enddo
        elseif(lines_mode.eq.-2) then
           !
           ! On-the-fly user-defined populations
           ! 
           ! First find the "average" position along this ray segment
           !
           xav = 0.5d0 * ( ray_cart_x + ray_prev_x )
           yav = 0.5d0 * ( ray_cart_y + ray_prev_y )
           zav = 0.5d0 * ( ray_cart_z + ray_prev_z )
           !
           ! On-the-fly user-defined population calculation
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 call userdef_compute_levelpop(ispec,                &
                      lines_nrlevels(ispec),ray_index,               &
                      xav,yav,zav,                                   &
                      gas_chemspec_numberdens(ispec,ray_index),      &
                      lines_ray_levpop(:,ispec,1))
              endif
           enddo
        elseif(lines_mode.eq.-3) then
           !
           ! On-the-fly LVG (Sobolev) populations
           !
           ! First determine the velocity gradient
           !
           call lines_compute_velgradient(ray_index,velgrad)
           !
           ! Now compute, with this velocity gradient, the
           ! LVG populations.
           !
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 if(lines_collisions_npartners(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no collision rates are loaded for'
                    write(stdo,*) '        this molecule, ispec=',ispec
                    stop
                 endif
                 if(lines_collisions_npartners_used(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no number densities of collision partners have been'
                    write(stdo,*) '        specified for this molecule, ispec=',ispec
                    stop
                 endif
                 !
                 ! Decide whether or not to include Escape Probability in
                 ! the LVG computation (i.e. a max length scale to assure
                 ! that for low velocity gradient the optical depth does
                 ! not go to infinity).
                 !
                 if(allocated(lines_escprob_lengthscale)) then
                    if(allocated(lines_microturb)) then
                       turb = lines_microturb(ray_index)
                    else
                       turb = 0.d0
                    endif
                    call lines_compute_lvgpop(ispec,lines_nrlevels(ispec),         &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(ray_index),                &
                      gas_chemspec_numberdens(ispec,ray_index),                    &
                      collpartner_numberdens(:,ray_index),                         &
                      velgrad,lines_ray_levpop(:,ispec,1),                         &
                      lengthscale=lines_escprob_lengthscale(ray_index),            &
                      turb=turb)
                 else
                    call lines_compute_lvgpop(ispec,lines_nrlevels(ispec),         &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(ray_index),                &
                      gas_chemspec_numberdens(ispec,ray_index),                    &
                      collpartner_numberdens(:,ray_index),                         &
                      velgrad,lines_ray_levpop(:,ispec,1))
                 endif
              endif
           enddo
        elseif(lines_mode.eq.-4) then
           !
           ! On-the-fly optically thin non-LTE populations
           !
           ! Now compute the populations
           ! 
           do ispec=1,lines_nr_species
              if(lines_species_fullmolec(ispec)) then
                 if(lines_collisions_npartners(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no collision rates are loaded for'
                    write(stdo,*) '        this molecule, ispec=',ispec
                    stop
                 endif
                 if(lines_collisions_npartners_used(ispec).eq.0) then
                    write(stdo,*) 'ERROR in line transfer: Cannot do non-LTE transfer because'
                    write(stdo,*) '        no number densities of collision partners have been'
                    write(stdo,*) '        specified for this molecule, ispec=',ispec
                    stop
                 endif
                 call lines_compute_optthinpop(ispec,lines_nrlevels(ispec),        &
                      lines_nrlines(ispec),lines_collisions_npartners_used(ispec), &
                      lines_collisions_npartmax,gastemp(ray_index),                &
                      gas_chemspec_numberdens(ispec,ray_index),                    &
                      collpartner_numberdens(:,ray_index),                         &
                      lines_ray_levpop(:,ispec,1))
              endif
           enddo
        elseif(lines_mode.eq.-10) then
           !
           ! Entirely general userdef routine for computation of level populations
           ! NOTE: This mode is only available as on-the-fly calculation.
           !
           call userdef_general_compute_levelpop(ray_index, lines_ray_levpop(:,:,1))
        else
           write(stdo,*) 'ERROR: line_mode not known'
           stop
        endif
        !
        ! Copy the temperature and turbulence to local variable
        !
        lines_ray_temp(1)       = gastemp(ray_index)
        if(allocated(lines_microturb)) then
           lines_ray_turb(1)    = lines_microturb(ray_index)
        else
           lines_ray_turb(1)    = 1.d-30
        endif
        !
        ! If linelist molecules are present, then we must also store the
        ! number densities of the molecules
        !
        if(allocated(lines_linelist_eup)) then
           do ispec=1,lines_nr_species
              lines_ray_nrdens(ispec,1) = gas_chemspec_numberdens(ispec,ray_index)
           enddo
        endif
        !
        ! Compute the Lorentz delta parameter in case we are using the Voigt profile
        ! (Added by Thomas Peters 2011)
        !
        if(lines_profile.eq.1) call userdef_compute_lorentz_delta(ray_index)
        !
        ! Compute the doppler shift
        !
        ! ********** CHECK THE SIGN! ***********
        !
        if(igrid_coord.lt.100) then
           !
           ! Cartesian coordinates:
           !
           dummy = gasvelocity(1,ray_index) * ray_cart_dirx +   &
                   gasvelocity(2,ray_index) * ray_cart_diry +   &
                   gasvelocity(3,ray_index) * ray_cart_dirz
        elseif(igrid_coord.lt.200) then
           !
           ! Spherical coordinates: Here component 1 is the v_r component,
           ! component 2 is the v_theta component (in the same unit as
           ! v_r, i.e. cm/s) and component 3 is the v_phi component (also
           ! in the same unit: cm/s).
           !
           xav   = 0.5d0 * ( ray_cart_x + ray_prev_x )
           yav   = 0.5d0 * ( ray_cart_y + ray_prev_y )
           zav   = 0.5d0 * ( ray_cart_z + ray_prev_z )
           dummy = xav**2+yav**2
           rlen  = sqrt(dummy+zav**2)+1d-90
           dummy = sqrt(dummy)+1d-90
           cosp  = xav/dummy
           sinp  = yav/dummy
           cost  = zav/rlen
           sint  = dummy/rlen
           dummy = sint * gasvelocity(1,ray_index) +   &
                   cost * gasvelocity(2,ray_index) 
           vx    = cosp * dummy - sinp * gasvelocity(3,ray_index)
           vy    = sinp * dummy + cosp * gasvelocity(3,ray_index)
           vz    = cost * gasvelocity(1,ray_index) -   &
                   sint * gasvelocity(2,ray_index)
           dummy = vx * ray_cart_dirx +   &
                   vy * ray_cart_diry +   &
                   vz * ray_cart_dirz
        else
           stop 501
        endif
        lines_ray_doppler(1) = dummy / cc
        !
     else
        !
        ! Not in the grid, so put everything to zero
        !
        lines_ray_levpop(:,:,1) = 0.d0
        lines_ray_temp(1)       = 1d-2
        lines_ray_turb(1)       = 1d-30
        lines_ray_doppler(1)    = 1d-30
     endif
  endif
  !
  ! Compute the source terms and extinction coefficients
  !
  ! ...Reset the alp and src
  !
  alp(inu0:inu1) = 0.d0
  if(inclstokes) then
     src(inu0:inu1,1:4) = 0.d0
  else
     src(inu0:inu1,1) = 0.d0
  endif
  !
  ! Now only add sources if ray_index.gt.0
  !
  if(ray_index.gt.0) then
     !
     ! ...Add the dust emission and extinction
     !
     ! IMPORTANT NOTE: If we include grain alignment, then the dust
     !    opacities must be dealt with in a more complex manner:
     !    the opacity then becomes a matrix which, if the coordinate
     !    basis in the image plane is appropriately rotated, 
     !    becomes diagonal. Nevertheless the diagonal components
     !    of this Mueller matrix, i.e. the opacities for the
     !    different Stokes components are different for the Stokes
     !    components. This cannot be done here, because in this
     !    subroutine we assume that the opacity is isotropic. 
     !    So if grain alignment is included, we skip the dust
     !    opacity and emissivity here; assuming they are dealt
     !    with properly in the other parts of the code. 
     !
     if(rt_incl_dust) then
        do inu=inu0,inu1
           if(alignment_mode.eq.0) then
              !
              ! Find the the dust continuum extinction coefficients
              !
              do ispec=1,dust_nr_species
                 sources_alpha_a(ispec) = sources_dustdens(ispec) * sources_dustkappa_a(inu,ispec)
                 sources_alpha_s(ispec) = sources_dustdens(ispec) * sources_dustkappa_s(inu,ispec)
                 alp(inu) = alp(inu) + sources_alpha_a(ispec) + sources_alpha_s(ispec)
              enddo
              ! 
              ! The source function. First the thermal part.
              !
              do ispec=1,dust_nr_species
                 src(inu,1) = src(inu,1) + sources_alpha_a(ispec) *                       &
                      bplanck(sources_dusttemp(ispec),sources_frequencies(inu))
              enddo
           else
              !
              ! For alignment mode: only do scattering part of angle-averaged opacity
              !
              do ispec=1,dust_nr_species
                 sources_alpha_s(ispec) = sources_dustdens(ispec) * sources_dustkappa_s(inu,ispec)
                 alp(inu) = alp(inu) + sources_alpha_s(ispec)
              enddo
           endif
           !
           ! Then the scattering part; this part also should be done for aligned grains
           !
           if(scattering_mode.ge.1) then
              if(.not.dust_2daniso) then
                 !
                 ! Normal 3-D case
                 !
                 if(inclstokes) then
                    if(camera_mcscat_monochromatic) then
                       src(inu,1:4) = src(inu,1:4) + mcscat_scatsrc_iquv(1,ray_index,1:4,mcscat_current_dir)
                    else
                       src(inu,1:4) = src(inu,1:4) + mcscat_scatsrc_iquv(inu,ray_index,1:4,mcscat_current_dir)
                    endif
                 else
                    if(camera_mcscat_monochromatic) then
                       src(inu,1) = src(inu,1) + mcscat_scatsrc_iquv(1,ray_index,1,mcscat_current_dir)
                    else
                       src(inu,1) = src(inu,1) + mcscat_scatsrc_iquv(inu,ray_index,1,mcscat_current_dir)
                    endif
                 endif
              else
                 !
                 ! Special 2-D axisymmetric mode for spherical coordinates
                 !
                 ! First find the "average" position along this ray segment
                 !
                 xav = 0.5d0 * ( ray_cart_x + ray_prev_x )
                 yav = 0.5d0 * ( ray_cart_y + ray_prev_y )
                 zav = 0.5d0 * ( ray_cart_z + ray_prev_z )
                 !
                 ! Find the phi angle of position
                 !
                 levent   = sqrt(xav**2+yav**2)
                 if(abs(xav).gt.1d-8*levent) then
                    phievent = atan(yav/xav)
                    if(phievent.lt.0.d0) then
                       phievent = phievent + pi
                    endif
                    if(yav.lt.0.d0) then
                       phievent = phievent + pi
                    endif
                 else
                    if(yav.gt.0.d0) then
                       phievent = pihalf
                    else
                       phievent = pi+pihalf
                    endif
                 endif
                 !
                 ! Find the index of this phi-cell 
                 !
                 iphievent = floor(dust_2daniso_nphi*phievent/twopi) + 1
                 if(iphievent.eq.dust_2daniso_nphi+1) then
                    iphievent = dust_2daniso_nphi
                 endif
                 !
                 ! Now add the scattering source function
                 !
                 if(inclstokes) then
                    if(camera_mcscat_monochromatic) then
                       src(inu,1:4) = src(inu,1:4) + 0.5d0 * mcscat_scatsrc_iquv(1,ray_index,1:4,iphievent) + &
                                                     0.5d0 * mcscat_scatsrc_iquv(1,ray_index,1:4,iphievent+1)
                    else
                       src(inu,1:4) = src(inu,1:4) + 0.5d0 * mcscat_scatsrc_iquv(inu,ray_index,1:4,iphievent) + &
                                                     0.5d0 * mcscat_scatsrc_iquv(inu,ray_index,1:4,iphievent+1)
                    endif
                 else
                    if(camera_mcscat_monochromatic) then
                       src(inu,1) = src(inu,1) + 0.5d0 * mcscat_scatsrc_iquv(1,ray_index,1,iphievent) + &
                                                 0.5d0 * mcscat_scatsrc_iquv(1,ray_index,1,iphievent+1)
                    else
                       src(inu,1) = src(inu,1) + 0.5d0 * mcscat_scatsrc_iquv(inu,ray_index,1,iphievent) + &
                                                 0.5d0 * mcscat_scatsrc_iquv(inu,ray_index,1,iphievent+1)
                    endif
                 endif
              endif
           endif
        enddo
     endif
     !
     ! ...Then the continuous stellar source part
     !    Note: the fourpi factor is not needed here, because the stellarsrc_templates
     !          are already in "per steradian" form.
     !
     if(incl_stellarsrc.ne.0) then
        do inu=inu0,inu1
           dum = 0.d0
           do itemplate=1,stellarsrc_nrtemplates
              stellarsrc_cum(itemplate) = dum
              dum = dum + stellarsrc_dens(itemplate,ray_index) *      &
                          sources_stellarsrc_templates(inu,itemplate)
           enddo
           src(inu,1) = src(inu,1) + dum 
        enddo
     endif
     !
     ! ...Then the PAH part 
     !
     !########## STILL TO DO ###########
     ! 
     ! ...Add gas continuum sources, if applicable
     !
     if(rt_incl_gascont) then
        call gascont_addto_jnu_alpnu(ray_index,sources_nrfreq,inu0,inu1,   &
             sources_frequencies(:),src(:,1),alp(:))
     endif
     !
     ! ...Add the line sources and extinction, unless we use
     !    the doppler catching algorith in which case we must
     !    do the final adding of the line emission and extinction
     !    on-the-fly during the ray-tracing (and thus possibly
     !    during the sub-cell steps along the ray).
     !
     if(rt_incl_lines.and.(.not.sources_catch_doppler_line)) then
        call lines_serial_addto_jnu_alpnu(sources_nrfreq,inu0,inu1,   &
             sources_frequencies(:),src(:,1),alp(:))
     endif
     !
     ! ...Add the user-defined emissivity and absorptivity
     !
     if(rt_incl_userdef_srcalp) then
        call userdef_srcalp(ray_index,sources_nrfreq,inu0,inu1,     &
                            sources_frequencies(:),src(:,:),alp(:))
     endif
     !
     ! Now divide by the total opacity to obtain the source function
     !
     if(inclstokes) then
        src(inu0:inu1,1) = src(inu0:inu1,1) / ( alp(inu0:inu1) + 1d-199 )
        src(inu0:inu1,2) = src(inu0:inu1,2) / ( alp(inu0:inu1) + 1d-199 )
        src(inu0:inu1,3) = src(inu0:inu1,3) / ( alp(inu0:inu1) + 1d-199 )
        src(inu0:inu1,4) = src(inu0:inu1,4) / ( alp(inu0:inu1) + 1d-199 )
     else
        src(inu0:inu1,1)   = src(inu0:inu1,1) / ( alp(inu0:inu1) + 1d-199 )
     endif
  endif
  !
end subroutine sources_get_src_alp


!==========================================================================
!        SUBROUTINES FOR SECOND ORDER INTEGRATION ON VERTEX GRID
!==========================================================================

!--------------------------------------------------------------------------
!            COMPUTE SOURCE FUNCTION AND ALPHA AT VERTICES
!
! This subroutine is used for the second-order integration of the transfer
! equation. For this we need to know the sources at the corner points.
!
! NOTE: If sources_interpol_jnu is set, then snu = jnu rather than the
!       source function jnu/anu. 
!
! NOTE: Before calling this subroutine, please set the ray_cart_dirx,y,z
!       to the correct values: the direction in which the rays toward the
!       observer point. If you use local observer mode, then you must
!       let the sources_module know by setting sources_localobserver=.true.
!       and setting sources_observer_position(:) to the position of the
!       observer.
!--------------------------------------------------------------------------
subroutine sources_compute_snualphanu_at_vertices(inu,inclstokes)
  implicit none
  integer :: inu,icell,ivt,ix,iy,iz,cnt,icnt
  integer :: ixx,iyy,izz
  double precision :: nu,snu(1:4),anu,dpl,trb,tmp
  double precision :: r,theta,phi
  double precision :: dirx,diry,dirz,xbk,ybk,zbk
  double precision, allocatable :: src(:,:),alp(:)
  double precision, allocatable :: nup(:),ndown(:)
  type(amr_branch), pointer :: b
  logical :: inclstokes
  !
  ! Check if frequency index is in range
  !
  if(inu.gt.sources_nrfreq) then
     write(stdo,*) 'ERROR when computing snu and alphanu at vertices:'
     write(stdo,*) '      inu out of range.'
     stop
  endif
  !
  ! Allocate temporary arrays
  !
  allocate(src(1:sources_nrfreq,1:4))
  allocate(alp(1:sources_nrfreq))
  if(sources_catch_doppler_line) then
     allocate(nup(1:sources_vertex_lines_nractivetot))
     allocate(ndown(1:sources_vertex_lines_nractivetot))
  endif
  !
  ! Get frequency
  !
  nu = sources_frequencies(inu)
  !
  ! Fill the cell center values of j_nu and alpha_nu
  !
  do icell=1,nrcells
     !
     ! Find the position of the cell
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is present, so we get our information from there
        !
        ! Get the current cell index
        !
        ray_index = cellindex(icell)
        !
        ! Set the coordinates
        !
        b => amr_index_to_leaf(ray_index)%link
        if(.not.associated(b)) stop 3130
        ray_cart_x = amr_finegrid_xc(b%ixyzf(1),1,b%level)
        ray_cart_y = amr_finegrid_xc(b%ixyzf(2),2,b%level)
        ray_cart_z = amr_finegrid_xc(b%ixyzf(3),3,b%level)
     else
        !
        ! We have a regular grid, so no AMR tree available
        !
        ray_index = icell
        call amr_regular_get_ixyz(ray_index,ixx,iyy,izz)
        ray_cart_x = amr_finegrid_xc(ixx,1,0)
        ray_cart_y = amr_finegrid_xc(iyy,2,0)
        ray_cart_z = amr_finegrid_xc(izz,3,0)
     endif
     !
     ! If spherical coordinates, then convert to cartesian
     ! BUGFIX: 06.03.2015-2
     !
     if((igrid_coord.ge.100).and.(igrid_coord.lt.200)) then
        r          = ray_cart_x
        theta      = ray_cart_y
        phi        = ray_cart_z
        ray_cart_x = r*sin(theta)*cos(phi)
        ray_cart_y = r*sin(theta)*sin(phi)
        ray_cart_z = r*cos(theta)
     endif
     !
     ! Copy also to prev
     !
     ray_prev_x = ray_cart_x
     ray_prev_y = ray_cart_y
     ray_prev_z = ray_cart_z
     !
     ! If local observer, set the direction vector. For an observer
     ! at infinity this direction will be the same for all cells, and
     ! will thus be set beforehand, at the start of this routine.
     !
     if(sources_localobserver) then
        write(stdo,*) 'Local observer for now not possible in second order integration'
        stop
     endif
     !
     ! Compute the j_nu and alpha_nu, as well as (if requested)
     ! the line level populations and other line-related stuff.
     ! Note that if 'sources_catch_doppler_line' is switched on,
     ! we only calculate what is needed to compute the line 
     ! j_nu and alpha_nu; we do not add these to the src and alp
     ! yet. If this option is not set, then we DO add them to
     ! the src and alp.
     !
     call sources_get_src_alp(inu,inu,sources_nrfreq,src,alp,inclstokes)
     !
     ! Fill the anu and jnu
     !
     sources_cell_anu(ray_index) = alp(inu)
     if(inclstokes) then
        if(sources_interpol_jnu) then
           sources_cell_snu(ray_index,1:4) = src(inu,1:4)*alp(inu)
        else
           sources_cell_snu(ray_index,1:4) = src(inu,1:4)
        endif
     else
        if(sources_interpol_jnu) then
           sources_cell_snu(ray_index,1) = src(inu,1)*alp(inu)
        else
           sources_cell_snu(ray_index,1) = src(inu,1)
        endif
     endif
     !
     ! If sources_catch_doppler_line option is set, then fill the line stuff
     !
     if(sources_catch_doppler_line) then
        !
        ! The line-of-sight Doppler shift at this point
        !
        sources_cell_doppler(ray_index) = lines_ray_doppler(1)
!        !
!        ! The local microturbulence at this point
!        !
!        sources_cell_turb(ray_index) = lines_ray_turb(1)
!        !
!        ! The local temperature at this point
!        ! NOTE: THIS IS NOT REALLY NECESSARY. DOUBLE STORAGE...
!        ! 
!        sources_cell_temp(ray_index) = gastemp(ray_index)
!
!       Get the nup and ndown        
!
        call lines_get_active_nup_ndown(sources_vertex_lines_nractivetot, &
             nup,ndown)
        do icnt=1,sources_vertex_lines_nractivetot
           sources_cell_line_nup(icnt,ray_index)   = nup(icnt)
           sources_cell_line_ndown(icnt,ray_index) = ndown(icnt)
        enddo
     endif
  enddo
  !
  ! Now average these onto the corner cells
  !
  do ivt=1,amr_nr_vertices
     snu(1:4)=0.d0
     anu=0.d0
     cnt=0
     if(sources_catch_doppler_line) then
        dpl        = 0.d0
        trb        = 0.d0
        tmp        = 0.d0
        nup(:)     = 0.d0
        ndown(:)   = 0.d0
     endif
     if(.not.amr_tree_present) then
        call amr_regular_get_vertex_ixyz(ivt,ixx,iyy,izz)
     endif
     do iz=1,1+amr_zdim
        do iy=1,1+amr_ydim
           do ix=1,1+amr_xdim
              !
              ! Find the cell connecting to the vertex in this
              ! octant
              !
              if(amr_tree_present) then
                 !
                 ! If there is an AMR tree, then it is easy
                 !
                 if(associated(amr_vertex_cells(ix,iy,iz,ivt)%link)) then
                    ray_index = amr_vertex_cells(ix,iy,iz,ivt)%link%leafindex
                 else
                    ray_index = 0
                 endif
              else
                 !
                 ! If the grid is regular, then we must do some more work
                 !
                 if(amr_xdim.eq.1) then
                    amrray_ix_curr = ixx+ix-2
                 else
                    amrray_ix_curr = 1
                 endif
                 if(amr_ydim.eq.1) then
                    amrray_iy_curr = iyy+iy-2
                 else
                    amrray_iy_curr = 1
                 endif
                 if(amr_zdim.eq.1) then
                    amrray_iz_curr = izz+iz-2
                 else
                    amrray_iz_curr = 1
                 endif
                 !
                 ! Special treatment if phi is cyclic
                 !
                 if(amr_cyclic_xyz(3)) then
                    if(amrray_iz_curr.lt.1) amrray_iz_curr=amr_grid_nz
                    if(amrray_iz_curr.gt.amr_grid_nz) amrray_iz_curr=1
                 endif
                 !
                 ! Special treatment of theta if near the midplane
                 !
                 if(amrray_iy_curr.gt.amr_grid_ny) then
                    if(amrray_mirror_equator) then
                       amrray_iy_curr = amr_grid_ny
                    endif
                 endif
                 !
                 ! For the rest simply check if we are on the grid or not
                 !
                 if((amrray_ix_curr.ge.1).and.(amrray_ix_curr.le.amr_grid_nx).and. &
                    (amrray_iy_curr.ge.1).and.(amrray_iy_curr.le.amr_grid_ny).and. &
                    (amrray_iz_curr.ge.1).and.(amrray_iz_curr.le.amr_grid_nz)) then
                    ray_index = amrray_ix_curr+((amrray_iy_curr-1)+(amrray_iz_curr-1)*amr_grid_ny)*amr_grid_nx
                 else
                    ray_index =  0
                 endif
              endif
              if(ray_index.gt.0) then
                 if(inclstokes) then
                    snu(1:4) = snu(1:4) + sources_cell_snu(ray_index,1:4)
                 else
                    snu(1) = snu(1) + sources_cell_snu(ray_index,1)
                 endif
                 anu = anu + sources_cell_anu(ray_index)
                 if(sources_catch_doppler_line) then
                    dpl = dpl + sources_cell_doppler(ray_index)
                    if(allocated(lines_microturb)) trb = trb + lines_microturb(ray_index)
                    tmp = tmp + gastemp(ray_index)
                    do icnt=1,sources_vertex_lines_nractivetot
                       nup(icnt)   = nup(icnt)   + sources_cell_line_nup(icnt,ray_index)
                       ndown(icnt) = ndown(icnt) + sources_cell_line_ndown(icnt,ray_index)
                    enddo
                 endif
                 cnt = cnt + 1
              endif
           enddo
        enddo
     enddo
     if(cnt.ge.1) then
        snu(1:4) = snu(1:4) / cnt
        anu = anu / cnt
        if(sources_catch_doppler_line) then
           dpl     = dpl / cnt
           trb     = trb / cnt
           tmp     = tmp / cnt
           nup   = nup / cnt
           ndown = ndown / cnt
        endif
     endif
     if(inclstokes) then
        sources_vertex_snu(ivt,1:4) = snu(1:4)
     else
        sources_vertex_snu(ivt,1) = snu(1)
     endif
     sources_vertex_anu(ivt) = anu
     if(sources_catch_doppler_line) then
        sources_vertex_doppler(ivt) = dpl
        sources_vertex_turb(ivt)    = trb
        sources_vertex_temp(ivt)    = tmp
        do icnt=1,sources_vertex_lines_nractivetot
           sources_vertex_line_nup(icnt,ivt)   = nup(icnt)
           sources_vertex_line_ndown(icnt,ivt) = ndown(icnt)
        enddo
     endif
  enddo
  !
  deallocate(src,alp)
  if(allocated(nup)) deallocate(nup,ndown)
  !
end subroutine sources_compute_snualphanu_at_vertices


!-------------------------------------------------------------------------
!             DO LINEAR INTERPOLATION OF ALPHA AND SOURCE
!
! This subroutine does the bilinear interpolation of the s_nu and alpha_nu
! (snu and anu) at the cell interfaces, based on the values at the cell
! corners (vertices).
!
! NOTE: If sources_interpol_jnu is set, then snu = jnu rather than the
!       source function jnu/anu. 
! 
! NOTE: This routine also uses information from the latest call to
!       amrray_find_next_location_***() to determine on which cell 
!       interface we are currently.
!-------------------------------------------------------------------------
subroutine sources_find_srcalp_interpol(x,y,z,snu,anu,inclstokes)
  implicit none
  double precision :: x,y,z,snu(1:4),anu
  logical :: inclstokes
  integer :: index,icross,ivt00,ivt01,ivt10,ivt11
  integer :: ivt000,ivt001,ivt010,ivt011
  integer :: ivt100,ivt101,ivt110,ivt111
  type(amr_branch), pointer :: a
  double precision :: epsx,epsx1,epsy,epsy1,epsz,epsz1
  double precision :: x0,x1,y0,y1,z0,z1
  integer :: ixx,iyy,izz,nnx,nny
  !
  ! Small computation
  !
  nnx = amr_grid_nx + amr_xdim
  nny = amr_grid_ny + amr_ydim
  !
  ! Check which cell to use
  !
  ! New 2017.03.18: Allow inside-cell point 
  ! 
  if(amrray_icross.gt.0) then
     index = ray_index
     icross = amrray_icross
  elseif(amrray_icross.lt.0) then
     index = ray_indexnext
     icross = -amrray_icross
  else
     if(ray_index.le.0) then
        snu(1:4) = 0.d0
        anu = 1d-99
        return
     else
        index = ray_index
        icross = 0
     endif
  endif
  !
  ! Check
  !
  if(index.le.0) then
     write(stdo,*) 'ERROR in srcalp interpol: index le 0'
     stop
  endif
  if(amr_tree_present) then
     if(.not.associated(amr_index_to_leaf(index)%link)) then
        write(stdo,*) 'ERROR in srcalp interpol: link to leaf not associated'
        stop
     endif
  endif
  !
  ! Find the pointer to the cell
  !
  if(amr_tree_present) then
     !
     ! We have an AMR tree, so get the information from the tree
     !
     a => amr_index_to_leaf(index)%link
  else
     call amr_regular_get_ixyz(index,ixx,iyy,izz)
  endif
  !
  ! Check out which interface
  !
  select case(icross)
  case(0)
     !
     ! Somewhere inside the cell. We have to do a trillinear interpolation
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(a%ixyzf(1),1,a%level)
           x1 = amr_finegrid_xi(a%ixyzf(1)+1,1,a%level)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(a%ixyzf(2),2,a%level)
           y1 = amr_finegrid_xi(a%ixyzf(2)+1,2,a%level)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(a%ixyzf(3),3,a%level)
           z1 = amr_finegrid_xi(a%ixyzf(3)+1,3,a%level)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt000 = amr_cell_corners(1,1,1,index)
        ivt010 = amr_cell_corners(1,1+amr_ydim,1,index)
        ivt001 = amr_cell_corners(1,1,1+amr_zdim,index)
        ivt011 = amr_cell_corners(1,1+amr_ydim,1+amr_zdim,index)
        ivt100 = amr_cell_corners(1+amr_xdim,1,1,index)
        ivt110 = amr_cell_corners(1+amr_xdim,1+amr_ydim,1,index)
        ivt101 = amr_cell_corners(1+amr_xdim,1,1+amr_zdim,index)
        ivt111 = amr_cell_corners(1+amr_xdim,1+amr_ydim,1+amr_zdim,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(ixx,1,0)
           x1 = amr_finegrid_xi(ixx+1,1,0)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(iyy,2,0)
           y1 = amr_finegrid_xi(iyy+1,2,0)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(izz,3,0)
           z1 = amr_finegrid_xi(izz+1,3,0)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt000 = ixx+((iyy-1)+(izz-1)*nny)*nnx
        ivt010 = ixx+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt001 = ixx+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt011 = ixx+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
        ivt100 = ixx+amr_xdim+((iyy-1)+(izz-1)*nny)*nnx
        ivt110 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt101 = ixx+amr_xdim+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt111 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsx1*(epsz1*(epsy1*sources_vertex_snu(ivt000,1:4)+epsy*sources_vertex_snu(ivt010,1:4)) + &
                          epsz*(epsy1*sources_vertex_snu(ivt001,1:4)+epsy*sources_vertex_snu(ivt011,1:4))) + &
                   epsx*(epsz1*(epsy1*sources_vertex_snu(ivt100,1:4)+epsy*sources_vertex_snu(ivt110,1:4)) +  &
                          epsz*(epsy1*sources_vertex_snu(ivt101,1:4)+epsy*sources_vertex_snu(ivt111,1:4)))
     else
        snu(1) = epsx1*(epsz1*(epsy1*sources_vertex_snu(ivt000,1)+epsy*sources_vertex_snu(ivt010,1)) + &
                        epsz*(epsy1*sources_vertex_snu(ivt001,1)+epsy*sources_vertex_snu(ivt011,1))) + &
                 epsx*(epsz1*(epsy1*sources_vertex_snu(ivt100,1)+epsy*sources_vertex_snu(ivt110,1)) + &
                        epsz*(epsy1*sources_vertex_snu(ivt101,1)+epsy*sources_vertex_snu(ivt111,1))) 

     endif
     anu   = epsx1*(epsz1*(epsy1*sources_vertex_anu(ivt000)+epsy*sources_vertex_anu(ivt010)) + &
                    epsz*(epsy1*sources_vertex_anu(ivt001)+epsy*sources_vertex_anu(ivt011))) + &
             epsx*(epsz1*(epsy1*sources_vertex_anu(ivt100)+epsy*sources_vertex_anu(ivt110)) + &
                   epsz*(epsy1*sources_vertex_anu(ivt101)+epsy*sources_vertex_anu(ivt111)))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
         epsx1*(epsz1*(epsy1*sources_vertex_doppler(ivt000)+epsy*sources_vertex_doppler(ivt010)) + &
                epsz *(epsy1*sources_vertex_doppler(ivt001)+epsy*sources_vertex_doppler(ivt011)))+ &
          epsx*(epsz1*(epsy1*sources_vertex_doppler(ivt100)+epsy*sources_vertex_doppler(ivt110)) + &
                epsz *(epsy1*sources_vertex_doppler(ivt101)+epsy*sources_vertex_doppler(ivt111)))
        sources_local_turb_curr = &
         epsx1*(epsz1*(epsy1*sources_vertex_turb(ivt000)+epsy*sources_vertex_turb(ivt010)) + &
                epsz *(epsy1*sources_vertex_turb(ivt001)+epsy*sources_vertex_turb(ivt011)))+ &
          epsx*(epsz1*(epsy1*sources_vertex_turb(ivt100)+epsy*sources_vertex_turb(ivt110)) + &
                epsz *(epsy1*sources_vertex_turb(ivt101)+epsy*sources_vertex_turb(ivt111)))
        sources_local_temp_curr = &
         epsx1*(epsz1*(epsy1*sources_vertex_temp(ivt000)+epsy*sources_vertex_temp(ivt010)) + &
                epsz *(epsy1*sources_vertex_temp(ivt001)+epsy*sources_vertex_temp(ivt011)))+ &
          epsx*(epsz1*(epsy1*sources_vertex_temp(ivt100)+epsy*sources_vertex_temp(ivt110)) + &
                epsz *(epsy1*sources_vertex_temp(ivt101)+epsy*sources_vertex_temp(ivt111)))
        sources_local_line_nup_curr(:) = &
         epsx1*(epsz1*(epsy1*sources_vertex_line_nup(:,ivt000)+epsy*sources_vertex_line_nup(:,ivt010)) + &
                epsz *(epsy1*sources_vertex_line_nup(:,ivt001)+epsy*sources_vertex_line_nup(:,ivt011)))+ &
          epsx*(epsz1*(epsy1*sources_vertex_line_nup(:,ivt100)+epsy*sources_vertex_line_nup(:,ivt110)) + &
                epsz *(epsy1*sources_vertex_line_nup(:,ivt101)+epsy*sources_vertex_line_nup(:,ivt111)))
        sources_local_line_ndown_curr(:) = &
         epsx1*(epsz1*(epsy1*sources_vertex_line_ndown(:,ivt000)+epsy*sources_vertex_line_ndown(:,ivt010)) + &
                epsz *(epsy1*sources_vertex_line_ndown(:,ivt001)+epsy*sources_vertex_line_ndown(:,ivt011)))+ &
          epsx*(epsz1*(epsy1*sources_vertex_line_ndown(:,ivt100)+epsy*sources_vertex_line_ndown(:,ivt110)) + &
                epsz *(epsy1*sources_vertex_line_ndown(:,ivt101)+epsy*sources_vertex_line_ndown(:,ivt111)))
     endif
  case(1)
     !
     ! Left x-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(a%ixyzf(2),2,a%level)
           y1 = amr_finegrid_xi(a%ixyzf(2)+1,2,a%level)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(a%ixyzf(3),3,a%level)
           z1 = amr_finegrid_xi(a%ixyzf(3)+1,3,a%level)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = amr_cell_corners(1,1,1,index)
        ivt10 = amr_cell_corners(1,1+amr_ydim,1,index)
        ivt01 = amr_cell_corners(1,1,1+amr_zdim,index)
        ivt11 = amr_cell_corners(1,1+amr_ydim,1+amr_zdim,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(iyy,2,0)
           y1 = amr_finegrid_xi(iyy+1,2,0)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(izz,3,0)
           z1 = amr_finegrid_xi(izz+1,3,0)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = ixx+((iyy-1)+(izz-1)*nny)*nnx
        ivt10 = ixx+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt01 = ixx+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt11 = ixx+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsz1*(epsy1*sources_vertex_snu(ivt00,1:4)+epsy*sources_vertex_snu(ivt10,1:4)) + &
                    epsz*(epsy1*sources_vertex_snu(ivt01,1:4)+epsy*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsz1*(epsy1*sources_vertex_snu(ivt00,1)+epsy*sources_vertex_snu(ivt10,1)) + &
                  epsz*(epsy1*sources_vertex_snu(ivt01,1)+epsy*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsz1*(epsy1*sources_vertex_anu(ivt00)+epsy*sources_vertex_anu(ivt10)) + &
              epsz*(epsy1*sources_vertex_anu(ivt01)+epsy*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsz1*(epsy1*sources_vertex_doppler(ivt00)+epsy*sources_vertex_doppler(ivt10)) + &
             epsz *(epsy1*sources_vertex_doppler(ivt01)+epsy*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsz1*(epsy1*sources_vertex_turb(ivt00)+epsy*sources_vertex_turb(ivt10)) + &
             epsz *(epsy1*sources_vertex_turb(ivt01)+epsy*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsz1*(epsy1*sources_vertex_temp(ivt00)+epsy*sources_vertex_temp(ivt10)) + &
             epsz *(epsy1*sources_vertex_temp(ivt01)+epsy*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsz1*(epsy1*sources_vertex_line_nup(:,ivt00)+epsy*sources_vertex_line_nup(:,ivt10)) + &
             epsz *(epsy1*sources_vertex_line_nup(:,ivt01)+epsy*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsz1*(epsy1*sources_vertex_line_ndown(:,ivt00)+epsy*sources_vertex_line_ndown(:,ivt10)) + &
             epsz *(epsy1*sources_vertex_line_ndown(:,ivt01)+epsy*sources_vertex_line_ndown(:,ivt11))
     endif
  case(2)
     !
     ! Right x-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(a%ixyzf(2),2,a%level)
           y1 = amr_finegrid_xi(a%ixyzf(2)+1,2,a%level)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(a%ixyzf(3),3,a%level)
           z1 = amr_finegrid_xi(a%ixyzf(3)+1,3,a%level)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = amr_cell_corners(2,1,1,index)
        ivt10 = amr_cell_corners(2,1+amr_ydim,1,index)
        ivt01 = amr_cell_corners(2,1,1+amr_zdim,index)
        ivt11 = amr_cell_corners(2,1+amr_ydim,1+amr_zdim,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(iyy,2,0)
           y1 = amr_finegrid_xi(iyy+1,2,0)
           epsy  = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2091
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(izz,3,0)
           z1 = amr_finegrid_xi(izz+1,3,0)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = ixx+amr_xdim+((iyy-1)+(izz-1)*nny)*nnx
        ivt10 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt01 = ixx+amr_xdim+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt11 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsz1*(epsy1*sources_vertex_snu(ivt00,1:4)+epsy*sources_vertex_snu(ivt10,1:4)) + &
                    epsz*(epsy1*sources_vertex_snu(ivt01,1:4)+epsy*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsz1*(epsy1*sources_vertex_snu(ivt00,1)+epsy*sources_vertex_snu(ivt10,1)) + &
                  epsz*(epsy1*sources_vertex_snu(ivt01,1)+epsy*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsz1*(epsy1*sources_vertex_anu(ivt00)+epsy*sources_vertex_anu(ivt10)) + &
              epsz*(epsy1*sources_vertex_anu(ivt01)+epsy*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsz1*(epsy1*sources_vertex_doppler(ivt00)+epsy*sources_vertex_doppler(ivt10)) + &
             epsz *(epsy1*sources_vertex_doppler(ivt01)+epsy*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsz1*(epsy1*sources_vertex_turb(ivt00)+epsy*sources_vertex_turb(ivt10)) + &
             epsz *(epsy1*sources_vertex_turb(ivt01)+epsy*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsz1*(epsy1*sources_vertex_temp(ivt00)+epsy*sources_vertex_temp(ivt10)) + &
             epsz *(epsy1*sources_vertex_temp(ivt01)+epsy*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsz1*(epsy1*sources_vertex_line_nup(:,ivt00)+epsy*sources_vertex_line_nup(:,ivt10)) + &
             epsz *(epsy1*sources_vertex_line_nup(:,ivt01)+epsy*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsz1*(epsy1*sources_vertex_line_ndown(:,ivt00)+epsy*sources_vertex_line_ndown(:,ivt10)) + &
             epsz *(epsy1*sources_vertex_line_ndown(:,ivt01)+epsy*sources_vertex_line_ndown(:,ivt11))
     endif
  case(3)
     !
     ! Left y-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(a%ixyzf(1),1,a%level)
           x1 = amr_finegrid_xi(a%ixyzf(1)+1,1,a%level)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(a%ixyzf(3),3,a%level)
           z1 = amr_finegrid_xi(a%ixyzf(3)+1,3,a%level)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = amr_cell_corners(1,1,1,index)
        ivt10 = amr_cell_corners(1+amr_xdim,1,1,index)
        ivt01 = amr_cell_corners(1,1,1+amr_zdim,index)
        ivt11 = amr_cell_corners(1+amr_xdim,1,1+amr_zdim,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(ixx,1,0)
           x1 = amr_finegrid_xi(ixx+1,1,0)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(izz,3,0)
           z1 = amr_finegrid_xi(izz+1,3,0)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = ixx+((iyy-1)+(izz-1)*nny)*nnx
        ivt10 = ixx+amr_xdim+((iyy-1)+(izz-1)*nny)*nnx
        ivt01 = ixx+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt11 = ixx+amr_xdim+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsz1*(epsx1*sources_vertex_snu(ivt00,1:4)+epsx*sources_vertex_snu(ivt10,1:4)) + &
                    epsz*(epsx1*sources_vertex_snu(ivt01,1:4)+epsx*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsz1*(epsx1*sources_vertex_snu(ivt00,1)+epsx*sources_vertex_snu(ivt10,1)) + &
                  epsz*(epsx1*sources_vertex_snu(ivt01,1)+epsx*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsz1*(epsx1*sources_vertex_anu(ivt00)+epsx*sources_vertex_anu(ivt10)) + &
              epsz*(epsx1*sources_vertex_anu(ivt01)+epsx*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsz1*(epsx1*sources_vertex_doppler(ivt00)+epsx*sources_vertex_doppler(ivt10)) + &
             epsz *(epsx1*sources_vertex_doppler(ivt01)+epsx*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsz1*(epsx1*sources_vertex_turb(ivt00)+epsx*sources_vertex_turb(ivt10)) + &
             epsz *(epsx1*sources_vertex_turb(ivt01)+epsx*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsz1*(epsx1*sources_vertex_temp(ivt00)+epsx*sources_vertex_temp(ivt10)) + &
             epsz *(epsx1*sources_vertex_temp(ivt01)+epsx*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsz1*(epsx1*sources_vertex_line_nup(:,ivt00)+epsx*sources_vertex_line_nup(:,ivt10)) + &
             epsz *(epsx1*sources_vertex_line_nup(:,ivt01)+epsx*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsz1*(epsx1*sources_vertex_line_ndown(:,ivt00)+epsx*sources_vertex_line_ndown(:,ivt10)) + &
             epsz *(epsx1*sources_vertex_line_ndown(:,ivt01)+epsx*sources_vertex_line_ndown(:,ivt11))
     endif
  case(4)
     !
     ! Right y-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(a%ixyzf(1),1,a%level)
           x1 = amr_finegrid_xi(a%ixyzf(1)+1,1,a%level)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(a%ixyzf(3),3,a%level)
           z1 = amr_finegrid_xi(a%ixyzf(3)+1,3,a%level)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = amr_cell_corners(1,2,1,index)
        ivt10 = amr_cell_corners(1+amr_xdim,2,1,index)
        ivt01 = amr_cell_corners(1,2,1+amr_zdim,index)
        ivt11 = amr_cell_corners(1+amr_xdim,2,1+amr_zdim,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(ixx,1,0)
           x1 = amr_finegrid_xi(ixx+1,1,0)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_zdim.eq.1) then
           z0 = amr_finegrid_xi(izz,3,0)
           z1 = amr_finegrid_xi(izz+1,3,0)
           epsz = (z-z0)/(z1-z0)
           epsz1 = 1.d0-epsz
           if((epsz.lt.-1d-6).or.(epsz.gt.1.000001d0)) stop 2092
        else
           epsz  = 0.d0
           epsz1 = 1.d0
        endif
        ivt00 = ixx+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt10 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt01 = ixx+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
        ivt11 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsz1*(epsx1*sources_vertex_snu(ivt00,1:4)+epsx*sources_vertex_snu(ivt10,1:4)) + &
                    epsz*(epsx1*sources_vertex_snu(ivt01,1:4)+epsx*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsz1*(epsx1*sources_vertex_snu(ivt00,1)+epsx*sources_vertex_snu(ivt10,1)) + &
                  epsz*(epsx1*sources_vertex_snu(ivt01,1)+epsx*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsz1*(epsx1*sources_vertex_anu(ivt00)+epsx*sources_vertex_anu(ivt10)) + &
              epsz*(epsx1*sources_vertex_anu(ivt01)+epsx*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsz1*(epsx1*sources_vertex_doppler(ivt00)+epsx*sources_vertex_doppler(ivt10)) + &
             epsz *(epsx1*sources_vertex_doppler(ivt01)+epsx*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsz1*(epsx1*sources_vertex_turb(ivt00)+epsx*sources_vertex_turb(ivt10)) + &
             epsz *(epsx1*sources_vertex_turb(ivt01)+epsx*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsz1*(epsx1*sources_vertex_temp(ivt00)+epsx*sources_vertex_temp(ivt10)) + &
             epsz *(epsx1*sources_vertex_temp(ivt01)+epsx*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsz1*(epsx1*sources_vertex_line_nup(:,ivt00)+epsx*sources_vertex_line_nup(:,ivt10)) + &
             epsz *(epsx1*sources_vertex_line_nup(:,ivt01)+epsx*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsz1*(epsx1*sources_vertex_line_ndown(:,ivt00)+epsx*sources_vertex_line_ndown(:,ivt10)) + &
             epsz *(epsx1*sources_vertex_line_ndown(:,ivt01)+epsx*sources_vertex_line_ndown(:,ivt11))
     endif
  case(5)
     !
     ! Left z-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(a%ixyzf(1),1,a%level)
           x1 = amr_finegrid_xi(a%ixyzf(1)+1,1,a%level)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(a%ixyzf(2),2,a%level)
           y1 = amr_finegrid_xi(a%ixyzf(2)+1,2,a%level)
           epsy = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2092
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        ivt00 = amr_cell_corners(1,1,1,index)
        ivt10 = amr_cell_corners(1+amr_xdim,1,1,index)
        ivt01 = amr_cell_corners(1,1+amr_ydim,1,index)
        ivt11 = amr_cell_corners(1+amr_xdim,1+amr_ydim,1,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(ixx,1,0)
           x1 = amr_finegrid_xi(ixx+1,1,0)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(iyy,2,0)
           y1 = amr_finegrid_xi(iyy+1,2,0)
           epsy = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2092
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        ivt00 = ixx+((iyy-1)+(izz-1)*nny)*nnx
        ivt10 = ixx+amr_xdim+((iyy-1)+(izz-1)*nny)*nnx
        ivt01 = ixx+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
        ivt11 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsy1*(epsx1*sources_vertex_snu(ivt00,1:4)+epsx*sources_vertex_snu(ivt10,1:4)) + &
                    epsy*(epsx1*sources_vertex_snu(ivt01,1:4)+epsx*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsy1*(epsx1*sources_vertex_snu(ivt00,1)+epsx*sources_vertex_snu(ivt10,1)) + &
                  epsy*(epsx1*sources_vertex_snu(ivt01,1)+epsx*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsy1*(epsx1*sources_vertex_anu(ivt00)+epsx*sources_vertex_anu(ivt10)) + &
              epsy*(epsx1*sources_vertex_anu(ivt01)+epsx*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsy1*(epsx1*sources_vertex_doppler(ivt00)+epsx*sources_vertex_doppler(ivt10)) + &
             epsy *(epsx1*sources_vertex_doppler(ivt01)+epsx*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsy1*(epsx1*sources_vertex_turb(ivt00)+epsx*sources_vertex_turb(ivt10)) + &
             epsy *(epsx1*sources_vertex_turb(ivt01)+epsx*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsy1*(epsx1*sources_vertex_temp(ivt00)+epsx*sources_vertex_temp(ivt10)) + &
             epsy *(epsx1*sources_vertex_temp(ivt01)+epsx*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsy1*(epsx1*sources_vertex_line_nup(:,ivt00)+epsx*sources_vertex_line_nup(:,ivt10)) + &
             epsy *(epsx1*sources_vertex_line_nup(:,ivt01)+epsx*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsy1*(epsx1*sources_vertex_line_ndown(:,ivt00)+epsx*sources_vertex_line_ndown(:,ivt10)) + &
             epsy *(epsx1*sources_vertex_line_ndown(:,ivt01)+epsx*sources_vertex_line_ndown(:,ivt11))
     endif
  case(6)
     !
     ! Right z-wall...
     !
     ! Find the vertex indices and the linear interpolation epsilons
     !
     if(amr_tree_present) then
        !
        ! The AMR tree is available, so get the information from there
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(a%ixyzf(1),1,a%level)
           x1 = amr_finegrid_xi(a%ixyzf(1)+1,1,a%level)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(a%ixyzf(2),2,a%level)
           y1 = amr_finegrid_xi(a%ixyzf(2)+1,2,a%level)
           epsy = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2092
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        ivt00 = amr_cell_corners(1,1,2,index)
        ivt10 = amr_cell_corners(1+amr_xdim,1,2,index)
        ivt01 = amr_cell_corners(1,1+amr_ydim,2,index)
        ivt11 = amr_cell_corners(1+amr_xdim,1+amr_ydim,2,index)
     else
        !
        ! Regular grid, so construct the information from ixx,iyy,izz
        !
        if(amr_xdim.eq.1) then
           x0 = amr_finegrid_xi(ixx,1,0)
           x1 = amr_finegrid_xi(ixx+1,1,0)
           epsx  = (x-x0)/(x1-x0)
           epsx1 = 1.d0-epsx
           if((epsx.lt.-1d-6).or.(epsx.gt.1.000001d0)) stop 2091
        else
           epsx  = 0.d0
           epsx1 = 1.d0
        endif
        if(amr_ydim.eq.1) then
           y0 = amr_finegrid_xi(iyy,2,0)
           y1 = amr_finegrid_xi(iyy+1,2,0)
           epsy = (y-y0)/(y1-y0)
           epsy1 = 1.d0-epsy
           if((epsy.lt.-1d-6).or.(epsy.gt.1.000001d0)) stop 2092
        else
           epsy  = 0.d0
           epsy1 = 1.d0
        endif
        ivt00 = ixx+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt10 = ixx+amr_xdim+((iyy-1)+(izz-1+amr_zdim)*nny)*nnx
        ivt01 = ixx+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
        ivt11 = ixx+amr_xdim+((iyy-1+amr_ydim)+(izz-1+amr_zdim)*nny)*nnx
     endif
     !
     ! Do the interpolations
     !
     if(inclstokes) then
        snu(1:4) = epsy1*(epsx1*sources_vertex_snu(ivt00,1:4)+epsx*sources_vertex_snu(ivt10,1:4)) + &
                    epsy*(epsx1*sources_vertex_snu(ivt01,1:4)+epsx*sources_vertex_snu(ivt11,1:4))
     else
        snu(1) = epsy1*(epsx1*sources_vertex_snu(ivt00,1)+epsx*sources_vertex_snu(ivt10,1)) + &
                  epsy*(epsx1*sources_vertex_snu(ivt01,1)+epsx*sources_vertex_snu(ivt11,1))
     endif
     anu   = epsy1*(epsx1*sources_vertex_anu(ivt00)+epsx*sources_vertex_anu(ivt10)) + &
              epsy*(epsx1*sources_vertex_anu(ivt01)+epsx*sources_vertex_anu(ivt11))
     if(sources_catch_doppler_line) then
        sources_local_doppler_curr = &
             epsy1*(epsx1*sources_vertex_doppler(ivt00)+epsx*sources_vertex_doppler(ivt10)) + &
             epsy *(epsx1*sources_vertex_doppler(ivt01)+epsx*sources_vertex_doppler(ivt11))
        sources_local_turb_curr = &
             epsy1*(epsx1*sources_vertex_turb(ivt00)+epsx*sources_vertex_turb(ivt10)) + &
             epsy *(epsx1*sources_vertex_turb(ivt01)+epsx*sources_vertex_turb(ivt11))
        sources_local_temp_curr = &
             epsy1*(epsx1*sources_vertex_temp(ivt00)+epsx*sources_vertex_temp(ivt10)) + &
             epsy *(epsx1*sources_vertex_temp(ivt01)+epsx*sources_vertex_temp(ivt11))
        sources_local_line_nup_curr(:) = &
             epsy1*(epsx1*sources_vertex_line_nup(:,ivt00)+epsx*sources_vertex_line_nup(:,ivt10)) + &
             epsy *(epsx1*sources_vertex_line_nup(:,ivt01)+epsx*sources_vertex_line_nup(:,ivt11))
        sources_local_line_ndown_curr(:) = &
             epsy1*(epsx1*sources_vertex_line_ndown(:,ivt00)+epsx*sources_vertex_line_ndown(:,ivt10)) + &
             epsy *(epsx1*sources_vertex_line_ndown(:,ivt01)+epsx*sources_vertex_line_ndown(:,ivt11))
     endif
  end select
end subroutine sources_find_srcalp_interpol




end module sources_module
