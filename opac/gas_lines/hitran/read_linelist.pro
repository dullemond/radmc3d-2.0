;-------------------------------------------------------------
;              READ LINE LIST (RADMC-3D STYLE)
;
; ARGUMENT:
;   filename           File to read
;
; RETURNS:
;   Structure containing, for each of the lines:
;     id               Counter for the lines
;     nu               Frequency of line in Hz
;     lambda           Wavelength of line in micron
;     aud              Einstein A_ud coefficient in sec^-1
;     e_lo             Energy of lower level in erg
;     e_up             Energy of upper level in erg
;     g_lo             Statistical weight of lower level
;     g_up             Statistical weight of upper level
;     n_air            (if present: HITRAN n_air parameter)
;     delta            (if present: HITRAN delta parameter)
;     gam_air          (if present: HITRAN gamma_air parameter)
;     gam_self         (if present: HITRAN gamma_self parameter)
;   and (if available in the linelist file) the partition sum:
;     temp             Temperature grid points in Kelvin
;     psum             Partition sum at those temperatures
;-------------------------------------------------------------
function read_linelist,filename
  @natconst
  openr,1,filename
  str=''
  readf,1,str
  readf,1,str
  iformat=0
  readf,1,iformat
  if iformat ne 1 then stop
  readf,1,str
  molecule=''
  readf,1,molecule
  readf,1,str
  readf,1,str
  molweight=0.d0
  readf,1,molweight
  readf,1,str
  ipart=0
  readf,1,ipart
  readf,1,str
  ivoigt=0
  readf,1,ivoigt
  if ipart ne 0 then begin
     readf,1,str
     ntemp=0
     readf,1,ntemp
     readf,1,str
     data=dblarr(2,ntemp)
     readf,1,data
     temp=transpose(data[0,*])
     psum=transpose(data[1,*])
  endif else begin
     temp=0
     psum=0
  endelse
  readf,1,str
  nlines=0LL
  readf,1,nlines
  id=intarr(nlines)*0LL
  lambda=dblarr(nlines)
  aud=dblarr(nlines)
  e_lo=dblarr(nlines)
  e_up=dblarr(nlines)
  g_lo=dblarr(nlines)
  g_up=dblarr(nlines)
  if ivoigt ne 0 then begin
     n_air = dblarr(nlines)
     delta = dblarr(nlines)
     gam_air = dblarr(nlines)
     gam_self = dblarr(nlines)
     idum1 = 0LL
     dat1  = dblarr(10)
  endif else begin
     n_air = 0
     delta = 0
     gam_air = 0
     gam_self = 0
     idum1 = 0LL
     dat1  = dblarr(6)
  endelse
  readf,1,str
  i=0LL
  while i lt nlines do begin
     readf,1,idum1,dat1
     id[i]      = idum1
     lambda[i]  = dat1[0]        ;; In micron
     aud[i]     = dat1[1]        ;; In sec^-1
     e_lo[i]    = dat1[2]*cc*hh  ;; In erg
     e_up[i]    = dat1[3]*cc*hh  ;; In erg
     g_lo[i]    = dat1[4]
     g_up[i]    = dat1[5]
     if ivoigt ne 0 then begin
        n_air[i] = dat1[6]
        delta[i] = dat1[7]
        gam_air[i] = dat1[8]
        gam_self[i]= dat1[9]
     endif
     i=i+1LL
  endwhile
  close,1
  nu=1d4*cc/lambda
  return,{molecule:molecule,molweight:molweight,$
          id:id,nu:nu,lambda:lambda,aud:aud,e_lo:e_lo,e_up:e_up,g_lo:g_lo,g_up:g_up,$
          n_air:n_air,delta:delta,gam_air:gam_air,gam_self:gam_self,temp:temp,psum:psum}
end
