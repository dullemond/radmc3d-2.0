@readradmc.pro

print,'Which opacity name (only the base name, e.g. pyrmg70)?'
name=''
read,name
if name eq '' then name='pyrmg70'
o=readopac(spec=name,/scatmat)
print,'Which wavelength index?'
read,inu
print,'   ----> lambda = ',2.9979d14/o.freq[inu],' micron'
mu = cos(o.angle*!dpi/180.d0)
print,'Which scattering angle index?'
read,iang
print,'   ----> theta = ',o.angle[iang],' degrees'

z = dblarr(4,4)
z[0,0] = o.zmat[0,iang,inu]
z[0,1] = o.zmat[1,iang,inu]
z[1,0] = o.zmat[1,iang,inu]
z[1,1] = o.zmat[2,iang,inu]
z[2,2] = o.zmat[3,iang,inu]
z[2,3] = o.zmat[4,iang,inu]
z[3,2] = -o.zmat[4,iang,inu]
z[3,3] = o.zmat[5,iang,inu]

print,'Z = '
print,transpose(z)    ;; Transpose because in IDL the rightmost index is rows

;; Phase matrix:
P = 4*!dpi*z/o.kappa_scat[inu]   ;; Convention: isotropic P[0,0] = 1.0
print,'P = '
print,transpose(P)    ;; Transpose because in IDL the rightmost index is rows


print,'Give the 4 Stokes param (I,Q,U,V) for incident radiation'
stokes_in = dblarr(4)
read,stokes_in

stokes_out = dblarr(4)
for i=0,3 do for k=0,3 do stokes_out[i] = stokes_out[i] + p[i,k]*stokes_in[k]

print,'P#incident = '
print,transpose(stokes_out)

end
