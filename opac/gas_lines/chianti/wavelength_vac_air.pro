;-----------------------------------------------------------------------------
; Reference: B. Edlen (1966) Metrologia, Vol 2., pp. 71
;-----------------------------------------------------------------------------
pro convert_lamvac_to_lamair,lamvac_mic,lamair_mic
  ;;
  ;; Check range of validity
  ;;
  if min(lamvac_mic) lt 0.2d0 or max(lamvac_mic) gt 1.69d0 then begin
     print,'Unfortunately the vac to air conversion of lambda is only possible between 0.2 and 1.69 microns.'
     print,'Putting all values outside of this domain to zero.'
  endif
  ;;
  ;; Calculate refractive index of air
  ;;
  sigma  = 1.d0 / lamvac_mic
  nmin1s = 1d-8*(8342.13d0+2406030d0*(130d0-sigma^2)^(-1)+15997d0*(38.9d0-sigma^2)^(-1))
  ;;
  ;; Convert wavelength
  ;;
  lamair_mic = lamvac_mic / (1.d0+nmin1s)
  ;;
  ;; All wavelengths outside domain are put to zero
  ;;
  ii = where(lamvac_mic lt 0.2d0 or lamvac_mic gt 1.69d0)
  if ii[0] ge 0 then begin
     lamair_mic[ii] = 0.d0
  endif
  ;;
end
