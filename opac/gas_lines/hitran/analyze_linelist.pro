;-----------------------------------------------------------------------
;    Routine for calculating opacity and emissivity of LTE molecule
;
; INPUT:
;   linelistfile     Which linelist to plot
;   temp             Temperature of the gas in Kelvin
;   akms             Microturbulence in km/s
; 
; RETURNS:
;   lambdamic        Wavelengths of lines in microns
;   a0               alpha_nu0 extinction at line center per molecule
;                    in units of cm^2.
;   j0               j_nu0 emissivity at line center per molecule
;                    in units of erg sc^-1 Hz^-1 ster^-1
;-----------------------------------------------------------------------
@bplanck.pro
@read_linelist.pro
pro analyze_linelist,linelistfile,temp,akms,lambdamic,a0,j0
@natconst.pro
const1=hh/(4*pi)
const2=cc*cc/(2*hh)
;
; Read the molecular data file
;
q    = read_linelist(linelistfile)
nlin = n_elements(q.aud)
lambdamic = q.lambda
;
; Line width
;
a    = sqrt((akms*1d5)^2 + 2*kk*temp/q.molweight)
;
; Interpolate the partition sum
;
if temp lt min(q.temp) or temp gt max(q.temp) then stop
psum = interpol(q.psum,q.temp,temp)
print,'Partition sum = ',psum
;
; Compute the line center extinction coefficient per molecule and
; the line center emissivity per molecule
;
a0   = dblarr(nlin)
j0   = dblarr(nlin)
i    = 0LL
while i lt nlin do begin
   ;;
   ;; Compute the Bdu and Bud
   ;;
   gratio = q.g_up[i] / q.g_lo[i]
   bud    = q.aud[i] * const2 / q.nu[i]^3
   bdu    = bud * gratio
   ;;
   ;; Get the populations
   ;;
   ndown  = q.g_lo[i]*exp(-q.e_lo[i]/(kk*temp)) / psum
   nup    = q.g_up[i]*exp(-q.e_up[i]/(kk*temp)) / psum
   ;;
   ;; Determine the value of the line profile function
   ;; at line center
   ;;
   phi0  = cc / (a*q.nu[i]*sqrt(pi))
   ;;
   ;; The j_nu is now:
   ;;
   ;;                 h nu_0
   ;;          j_nu = ------ n_up A_ud phi(nu)
   ;;                  4 pi
   ;;
   j0[i] = const1 * q.nu[i] * phi0 * q.aud[i] * nup 
   ;;
   ;; The alpha_nu is now:
   ;;
   ;;                 h nu_0
   ;;      alpha_nu = ------ ( n_down B_du - n_up B_ud ) phi(nu)
   ;;                  4 pi
   ;;
   if(ndown gt 1d-40 and nup gt 1d-40) then begin
      a0[i] = const1 * q.nu[i] * phi0 * ( bdu * ndown - bud * nup )
   endif else begin
      a0[i] = 0.d0
   endelse
   ;;
   i = i+1LL
endwhile
;
end
