@readradmc.pro

print,'Which opacity name (only the base name, e.g. pyrmg70)?'
name=''
read,name
o=readopac(spec=name)
print,'Which wavelength index?'
read,inu
print,'   ----> lambda = ',2.9979d14/o.freq[inu],' micron'
mu = cos(o.angle*!dpi/180.d0)
phasefunc = o.zmat[0,*,inu]
;
; Normalize the phase function to 1
;
sum = 0.d0
for i=1,n_elements(mu)-1 do begin
   sum = sum + 0.25d0 * (phasefunc[i]+phasefunc[i-1]) * abs(mu[i]-mu[i-1])
endfor
phasefunc = phasefunc/sum
;
; Compute the g, just to compare
;
sum = 0.d0
for i=1,n_elements(mu)-1 do begin
   sum = sum + 0.25d0 * (phasefunc[i]*mu[i]+phasefunc[i-1]*mu[i-1]) * abs(mu[i]-mu[i-1])
endfor
print,'g computed = ',sum,' while g from file = ',o.g[inu]
;
; Make the Henyey-Greenstein phase function for comparison
;
g = o.g[inu]
phasefunc_henyey = (1.d0-g^2)/(1.d0+g^2-2*g*mu)^1.5d0

plot,o.angle,phasefunc,/yl,xtitle='scattering angle',ytitle='Phasefunction'
oplot,o.angle,phasefunc_henyey,line=2

end
