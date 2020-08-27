@readopac.pro
;
; In this example we convert the opacity from:
;
;   http://hera.ph1.uni-koeln.de/~ossk/Jena/tables/thin5
;
; (see also http://hera.ph1.uni-koeln.de/~ossk/Jena/tables.html)
;
; into RADMC-3D format. So please first download that file.
;
file='thin5'

data=[0.,0.]
openr,1,file
nlines=0
while not eof(1) do begin
   readf,1,data
   nlines=nlines+1
endwhile
close,1
data=dblarr(2,nlines)
openr,1,file
readf,1,data
close,1

;
; Write dustkappa file, but add one extra point at small
; wavelengths, to extrapolate the opacity into the optical.
; NOTE: This is very ad-hoc!!! But I see no better way.
;
openw,1,'dustkappa_osshenn.inp'
printf,1,1
printf,1,nlines+1
printf,1,0.1d0,data[1,0]
for i=0,nlines-1 do printf,1,data[0,i],data[1,i]
close,1

;
; Now read again in the "standard" RADMC-3D way
;
o=readopac(spec='osshenn')
plotopac,o

end

