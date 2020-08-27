;
; Some natural constants
;
AU  = 1.49598d13     ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]
;
; Grid parameters
;
nx       = 32L
ny       = 32L
nz       = 32L
sizex    = 10*AU
sizey    = 10*AU
sizez    = 10*AU
;
; Model parameters
;
radius   = 5*AU
rho0     = 1d-16
;
; Make the coordinates
;
xi       = -sizex + 2*sizex*dindgen(nx+1)/(1.d0*nx)
yi       = -sizey + 2*sizey*dindgen(ny+1)/(1.d0*ny)
zi       = -sizez + 2*sizez*dindgen(nz+1)/(1.d0*nz)
xc       = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
yc       = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
zc       = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
;
; Make the dust density model
;
xx       = rebin(xc,nx,ny,nz)
yy       = transpose(rebin(yc,ny,nx,nz),[1,0,2])
zz       = transpose(rebin(zc,nz,ny,nx),[2,1,0])
rr       = sqrt(xx^2+yy^2+zz^2)
rhog     = rho0 * exp(-(rr^2/radius^2)/2.d0)
;
; Write the wavelength_micron.inp file
;
lambda1 = 0.1d0
lambda2 = 7.0d0
lambda3 = 25.d0
lambda4 = 1.0d4
n12     = 20
n23     = 100
n34     = 30
lam12   = lambda1 * (lambda2/lambda1)^(dindgen(n12)/(1.d0*n12))
lam23   = lambda2 * (lambda3/lambda2)^(dindgen(n23)/(1.d0*n23))
lam34   = lambda3 * (lambda4/lambda3)^(dindgen(n34)/(1.d0*(n34-1.d0)))
lambda  = [lam12,lam23,lam34]
nlam    = n_elements(lambda)
;
; Write the wavelength file
;
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
close,1
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,0                      ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,nx,ny,nz               ; Size of grid
for i=0,nx do printf,1,xi[i]    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i]    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'gas_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhog[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'incl_userdef_srcalp = 1'
close,1
;
; Set up the default file for the transfer.inp 
;
openw,1,'transfer.inp.default'
printf,1,'transfer_density_plaw_r = 1.0'
printf,1,'transfer_density_plaw_g = 2.0'
printf,1,'transfer_density_plaw_b = 3.0'
close,1
;
end
