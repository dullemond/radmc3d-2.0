@wavelength_vac_air.pro
;-----------------------------------------------------------
; Produce a Leiden LAMDA style atomic data file from the
; Chianti database. It will automatically produce a file
; with filename "molecule_XXXX.inp" where XXXX is the 
; appropriate name of the atom. The word "molecule" is
; because RADMC-3D does not distinguish between atoms and
; molecules. And since the LAMDA database format was 
; originally meant for molecules, the name "molecule"
; was chosen.
;
; NOTE: Before using this subroutine you MUST type
;       @chianti_startup in IDL (possibly you may need
;       also to adjust the chianti_startup.pro file to
;       set the paths right).
;
; ARGUMENTS:
;  iz        Atomic number (e.g. oxygen has iz=8)
;  ion       Ionization degree (e.g. O++ = OIII ---> ion=3)
;  nlevels   Max number of levels to include
;  temp      Array of temperatures for which the collision
;            rates are to be computed and tabulated.
;
; KEYWORDS:
;  nprec     Precision of the collision rate tables (decimals
;            behind the comma). Default = 1.
;  full      If set, more information will be included in
;            the output file.
;
;-----------------------------------------------------------
pro convert_chianti_to_leiden,iz,ion,nlevels,temp,$
                              nprec=nprec,full=full
  COMMON elvlc,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
  COMMON wgfa, wvl,gf,a_value
  cc  = 2.9979245800000d10      ; Light speed             [cm/s]
  kk  = 1.3807d-16              ; Bolzmann's constant     [erg/K]
  hh  = 6.6262d-27              ; Planck's constant       [erg.s]
  ;;
  ;; Get atomic weight and name for this atom
  ;; 
  openr,1,'atom_names'
  str=' '
  for i=1,iz do begin
     readf,1,str
  endfor
  close,1
  molweight=strmid(str,0,8)*1.d0
  print,'Atomic weight = ',molweight
  name=strlowcase(strcompress(strmid(str,34,2),/remove_all))
  namechianti=name+'_'+strcompress(string(ion),/remove_all)
  if ion gt 1 then begin
     name = name+'+'
     if ion gt 2 then begin
        name = name+strcompress(string(ion-1),/remove_all)
     endif
  endif
  ;;
  ;; Get information about the levels
  ;;
  show_pops,iz,ion,popstr
  ;;
  ;; Now the information of the levels is in the common
  ;; blocks (presumably there is a more elegant way, but
  ;; this works).
  ;;
  if nlevels gt n_elements(ecm) then nlevels=n_elements(ecm)
  ;;
  ;; Count the lines 
  ;;
  nlines = 0
  for ilevel=0,nlevels-1 do begin
     for k=ilevel+1,nlevels-1 do begin
        if a_value[ilevel,k] gt 0.d0 then begin
           nlines = nlines + 1
        endif
     endfor
  endfor
  ;;
  ;; Get the collision rates in cm^3/s for the e- collision
  ;; partner. 
  ;;
  ntemp = n_elements(temp)
  colrate = dblarr(nlevels,nlevels,ntemp)
  for itemp=0,ntemp-1 do begin
     q=rate_coeff(namechianti,temp[itemp])
     colrate[*,*,itemp] = q[0:nlevels-1,0:nlevels-1]
  endfor
  ;;
  ;; If "full" is set, let us compute the critical density
  ;; of each level at 1E4 K, where we focus only on downward
  ;; collisional rates.
  ;;
  if keyword_set(full) then begin
     clr=rate_coeff(namechianti,1d4)
     ncrit=dblarr(nlevels)
     for ilevel=0,nlevels-1 do begin
        rad = 0.d0
        col = 0.d0
        for k=0,ilevel-1 do begin
           rad = rad + a_value[k,ilevel]
           col = col + clr[ilevel,k]
        endfor
        if col gt 0.d0 then begin
           ncrit[ilevel] = rad/col
        endif
     endfor
  endif
  ;;
  ;; Now write the leiden database
  ;;
  fileout='molecule_'+name+'.inp'
  openw,1,fileout
  printf,1,'! Atom'
  printf,1,name
  printf,1,'! Atomic weight'
  printf,1,molweight
  printf,1,'! Number of levels'
  printf,1,nlevels
  if keyword_set(full) then begin
     printf,1,'! Level + Energy[cm^-1] + weight + J + Term + N_crit [cm^-3] @ 1d4K'
     for ilevel=0,nlevels-1 do begin
        printf,1,ilevel+1,ecm[ilevel],2*jj[ilevel]+1,jj[ilevel],term[ilevel],ncrit[ilevel],$
               format='(I3,1X,E13.6,1X,E9.2,1X,E9.2,5X,A15,1X,E13.6)'
     endfor
  endif else begin
     printf,1,'! Level + Energy[cm^-1] + weight + J + Term'
     for ilevel=0,nlevels-1 do begin
        printf,1,ilevel+1,ecm[ilevel],2*jj[ilevel]+1,jj[ilevel],'   ',term[ilevel]
     endfor
  endelse
  printf,1,'! Nr of radiative transitions'
  printf,1,nlines
  if keyword_set(full) then begin
     printf,1,'! Trans + up + low + EinsteinA[s^-1] + Freq[GHz] + E_u[K] + Lam [micron] + ObsLam [micron] (Freq=0 means two-phot continuum)'
  endif else begin
     printf,1,'! Trans + up + low + EinsteinA[s^-1] + Freq[GHz] + E_u[K] (Freq=0 means two-phot continuum)'
  endelse
  iline=1
  for ilevel=0,nlevels-1 do begin
     for k=ilevel+1,nlevels-1 do begin
        if a_value[ilevel,k] gt 0.d0 then begin
           if wvl[ilevel,k] gt 0.d0 then begin
              nughz=cc/(wvl[ilevel,k]*1d-8)/1d9
           endif else begin
              nughz=0.d0   ;; Two-photon continuum
           endelse
           if keyword_set(full) then begin
              lamvac=wvl[ilevel,k]*1d-4
              if lamvac ge 0.2d0 and lamvac le 1.69d0 then begin
                 convert_lamvac_to_lamair,lamvac,lamair
                 slamair = strcompress(string(lamair,format='(E13.6)'),/remove_all)
              endif else begin
                 slamair = ''
              endelse
              printf,1,iline,k+1,ilevel+1,a_value[ilevel,k],nughz,hh*cc*ecm[k]/kk,wvl[ilevel,k]*1d-4,slamair,$
                     format='(I4,1X,I3,1X,I3,1X,E13.6,1X,E13.6,1X,E13.6,1X,E13.6,1X,A13)'
           endif else begin
              printf,1,iline,k+1,ilevel+1,a_value[ilevel,k],nughz,hh*cc*ecm[k]/kk
           endelse
           iline = iline+1
        endif
     endfor
  endfor
  printf,1,'! Nr of collision partners'
  printf,1,1
  printf,1,'! Collisions with'
  printf,1,'e-'
  printf,1,'! Nr of collisional transitions'
  printf,1,nlevels*(nlevels-1)/2
  printf,1,'! Nr of temperatures'
  printf,1,ntemp
  if not keyword_set(nprec) then nprec=1
  precis='E'+trim(string(7+nprec))+'.'+trim(string(nprec))
  printf,1,'! Temperatures'
  printf,1,temp,format='(15X,'+trim(string(ntemp))+'('+precis+'))'
  printf,1,'! Trans + up + low + Collrates[cm^3 s^-1]'
  itr=1
  for iup=0,nlevels-1 do begin
     for ilo=0,iup-1 do begin
        printf,1,itr,iup+1,ilo+1,colrate[iup,ilo,*],$
               format='(I6,1X,I3,1X,I3,1X,'+trim(string(ntemp))+'('+precis+'))'
        itr = itr + 1
     endfor
  endfor
  printf,1,'!'
  printf,1,'! Data from the CHIANTI database: http://www.chiantidatabase.org/'
  printf,1,'! Automatically converted into Leiden LAMDA format using the IDL routine '
  printf,1,'! convert_chianti_to_leiden.pro written by C.P. Dullemond, which uses IDL'
  printf,1,'! library routines from the CHIANTI distribution.'
  printf,1,'! Note: This file may contain only a subset of all possible levels and lines.'
  printf,1,'! Note: For finer sampling of collision rates in temperature, please redo conversion.'
  close,1
end

