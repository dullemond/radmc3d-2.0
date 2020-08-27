@convert_chianti_to_leiden.pro


print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
print,'Please first type:'
print,'@chianti_startup'
print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

;
; Temperatures at which to tabulate the collision rates with e-
;
ntemp = 5*2+1
tmin  = 1d1
tmax  = 1d6
temp  = tmin * (tmax/tmin)^(dindgen(ntemp)/(ntemp-1.d0))

;
; Now make the atomic data files (called "molecule_***.inp")
;
;convert_chianti_to_leiden,1,1,100,temp,nprec=2,/full      ;; Full H-atom
;convert_chianti_to_leiden,2,1,100,temp,nprec=2,/full      ;; Full He-atom
;convert_chianti_to_leiden,2,2,100,temp,nprec=2,/full      ;; Full He+
convert_chianti_to_leiden,8,3,6,temp,nprec=2,/full        ;; O++, only first 6 levels

end
