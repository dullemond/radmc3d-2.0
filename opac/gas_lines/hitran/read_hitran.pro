function read_hitran,file,nnl=nnl
@natconst
;file='01_hit08.par'
openr,1,file
if keyword_set(nnl) then begin
   nnl = nl
endif else begin
   nl=0LL
   readf,1,nl
endelse
nu=dblarr(nl)
aud=dblarr(nl)
eup=dblarr(nl)
edown=dblarr(nl)
gup=dblarr(nl)
gdown=dblarr(nl)
imol=intarr(nl)
iiso=intarr(nl)
n_air=dblarr(nl)
delta=dblarr(nl)
gam_air=dblarr(nl)
gam_self=dblarr(nl)
;
format='(I2,I1,F12.6,2E10.3,2F5.4,F10.4,F4.2,F8.6,2A15,2A15,6I1,6I2,A1,2F7.1)'
idum1=[0,0]
dum1=dblarr(8)
str1=strarr(4)
idum2=intarr(12)
str2=''
dum2=dblarr(2)
i=0LL
while i lt nl do begin
   readf,1,idum1,dum1,str1,idum2,str2,dum2,format=format
   imol[i]=idum1[0]               ;; ID of molecule
   iiso[i]=idum1[1]               ;; ID of isotopologue of molecule
   nu[i]=cc*dum1[0]               ;; Frequency in Hertz
   aud[i]=dum1[2]                 ;; A_ud in sec^-1
;   edown[i]=cc*hh*dum1[5]         ;; E_down in erg
;   eup[i]=cc*hh*(dum1[5]+dum1[0]) ;; E_up in erg
   edown[i]=dum1[5]               ;; E_down in cm^-1
   eup[i]=dum1[5]+dum1[0]         ;; E_up in cm^-1
   gup[i]=dum2[0]                 ;; g_up (statistical weight)
   gdown[i]=dum2[1]               ;; g_down (statistical weight)
   n_air[i]=dum1[6]               ;; n parameter of HITRAN for Voigt profile
   delta[i]=dum1[7]               ;; delta parameter of HITRAN for Voigt profile
   gam_air[i]=dum1[3]             ;; gamma_air of HITRAN for Voigt profile
   gam_self[i]=dum1[4]            ;; gamma_self of HITRAN for Voigt profile
   i=i+1LL
endwhile
close,1
return,{nu:nu,aud:aud,edown:edown,eup:eup,gdown:gdown,gup:gup,imol:imol,iiso:iiso,$
       n_air:n_air,delta:delta,gam_air:gam_air,gam_self:gam_self}
end
