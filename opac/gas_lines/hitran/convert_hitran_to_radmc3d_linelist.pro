@read_hitran
;-----------------------------------------------------------
;    CONVERT HITRAN TO RADMC-3D STANDARD LINE LIST STYLE
;
; ARGUMENTS:
;  filein            Input file from the HITRAN database
;  imol              Index of the molecule:
;                    (01=H2O,02=CO2,03=O3, etc, see 
;                    page 36 of HAWKSmanual.pdf)
;  iso               Isotopomer number (see molparam.txt,
;                    though that file is a bit tricky to
;                    fully understand; good luck ;-))
;  nameout           Name of the output molecule
;
; KEYWORDS:
;  ivoigt            If 0, then no Voigt profile info is
;                    included, if 1, then the HITRAN 
;                    style n_air, delta, gam_air, gam_self
;                    are included.
;  ipart             If 1, then a table of the partition
;                    sum as a function of temperature is
;                    included in the file.
;-----------------------------------------------------------
pro convert_hitran_to_radmc3d_linelist,filein,imol,iso,nameout,$
     ivoigt=ivoigt,ipart=ipart
@natconst
if n_elements(ivoigt) eq 0 then ivoigt=0
if n_elements(ipart) eq 0 then ipart=0
;
; Read HITRAN file
;
print,'Reading RATRAN file '+strcompress(filein,/remove_all)+'...'
a=read_hitran(filein)
;
; Names of molecules the HITRAN database
;
molecules = ["H2O","CO2","O3","N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl","HBr","HI","ClO","OCS","H2CO","HOCl","N2","HCN","CH3Cl","H2O2","C2H2","C2H6","PH3","COF2","SF6","H2S","HCOOH","HO2","O","ClONO2","NO+","HOBr","C2H4","CH3OH","CH3Br","CH3CN","CF4"]
print,'Base molecule name (irrespective of isotopomer) = ',molecules[imol-1]
;
; Number of isotopologues for each of these molecules, as listed in 
; the molparam.txt file
;
niso      = [6    ,10   ,5   ,5    ,6   ,4    ,3   ,3   ,2    ,1    ,2    ,1     ,3   ,1   ,2    ,2    ,1   ,2    ,5    ,3     ,2     ,1   ,3    ,2      ,1     ,2     ,2     ,1    ,1     ,1    ,3    ,1      ,1    , 1 ,2       ,1    ,2     ,2     ,1      ,2      ,1      ,1    ]
;
; Get the molecular weight of this isotopomer
;
print,'Reading molparam.txt...'
openr,1,'molparam.txt'
str=""
readf,1,str
for i=0,imol-1 do begin
   readf,1,str
   for k=0,niso[i]-1 do begin
      if i eq imol-1 and k eq iso-1 then begin
         d=dblarr(5)
         readf,1,d
         molweight=d[4]
      endif else begin
         readf,1,str
      endelse
   endfor
endfor
close,1
print,'Molecular weight = ',molweight
;
; Count number of lines for this molecule and this isotopomer
;
nlinmax = n_elements(a.nu)
nlin = 0LL
i=0LL
while i lt nlinmax do begin
   if a.imol[i] eq imol and a.iiso[i] eq iso then nlin=nlin+1
   i=i+1
endwhile
;
; If requested, read also the partition sum data
;
if ipart gt 0 then begin
   print,'Reading partition function data from parsum.dat...'
   format='(A17,'+strcompress(string(ipart),/remove_all)+'A27)'
   str=strarr(ipart+1)
   openr,1,'parsum.dat'
   readf,1,str,format=format
   print,'Taking parition sum for column '+strcompress(str[ipart],/remove_all)+' of parsum.dat'
   nt=3000-70+1
   temp=dblarr(nt)
   psum=dblarr(nt)
   data=dblarr(ipart+1)
   for i=0,nt-1 do begin
      readf,1,data
      temp[i]=data[0]
      psum[i]=data[ipart]
   endfor
   close,1
endif
;
; Now write RADMC-3D standard line list style for this 
; molecule and this isotopomer
;
fileout='linelist_'+nameout+'.inp'
print,'Writing output file '+fileout
openw,1,fileout
printf,1,'! RADMC-3D Standard line list'
printf,1,'! Format number:'
printf,1,'1'
printf,1,'! Molecule name:'
printf,1,nameout
printf,1,'! Reference: From the HITRAN Database (see below for more info)'
printf,1,'! Molecular weight (in atomic units)'
smolweight=strcompress(string(molweight),/remove_all)
printf,1,'  ',smolweight
printf,1,'! Include table of partition sum? (0=no, 1=yes)'
if ipart eq 0 then begin
   printf,1,'  0'
endif else begin
   printf,1,'  1'
endelse
printf,1,'! Include information on Voigt profile? (0=no, 1=yes)'
sivoigt=strcompress(string(ivoigt),/remove_all)
printf,1,'  ',sivoigt
if ipart ne 0 then begin
   printf,1,'! Nr of temperature points for the partition sum'
   printf,1,nt
   printf,1,'!  Temp [K]      PartSum'
   format='(E13.6,1X,E13.6)'
   for i=0,nt-1 do begin
      printf,1,temp[i],psum[i],format=format
   endfor
endif
printf,1,'! Nr of lines'
snlin=strcompress(string(nlin),/remove_all)
printf,1,'  ',snlin
if ivoigt eq 0 then begin
   printf,1,'! ID    Lambda [mic]  Aud [sec^-1]  E_lo [cm^-1]  E_up [cm^-1]  g_lo  g_up  '
   format='(I6,1X,E13.6,1X,E13.6,1X,E13.6,1X,E13.6,1X,F5.0,1X,F5.0)'
endif else begin
   printf,1,'! ID    Lambda [mic]  Aud [sec^-1]  E_lo [cm^-1]  E_up [cm^-1]  g_lo  g_up   n_air delta    gair   gself'
   format='(I6,1X,E13.6,1X,E13.6,1X,E13.6,1X,E13.6,1X,F5.0,1X,F5.0,2X,F5.2,1X,F9.6,1X,F6.4,1X,F6.4)'
endelse
ilin=1LL
i=0LL
while i lt nlinmax do begin
   if a.imol[i] eq imol and a.iiso[i] eq iso then begin
      if ivoigt eq 0 then begin
         printf,1,ilin,1d4*cc/a.nu[i],a.aud[i],a.edown[i],a.eup[i],a.gdown[i],a.gup[i],format=format
      endif else begin
         printf,1,ilin,1d4*cc/a.nu[i],a.aud[i],a.edown[i],a.eup[i],a.gdown[i],a.gup[i],$
                a.n_air[i],a.delta[i],a.gam_air[i],a.gam_self[i],format=format
      endelse
      ilin = ilin + 1
   endif
   i=i+1
endwhile
printf,1,'!--------------------------------------------------------'
printf,1,'! More information on the source of the data:'
printf,1,'! The HITRAN database is maintained by L. Rothman.'
printf,1,'! The web site is: http://www.cfa.harvard.edu/HITRAN/'
printf,1,'! Note that there is a special JQRST issue on HITRAN:'
printf,1,'! Journal of Quantitative Spectroscopy and Radiative '
printf,1,'! Transfer (JQSRT) volume 110, numbers 9-10, June/July 2009'
printf,1,'! These data were converted to RADMC-3D style by C.P. Dullemond'
close,1
end
