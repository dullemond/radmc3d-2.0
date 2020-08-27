@analyze_linelist
;@bplanck
@natconst

temp=300.
akms=2.d0

;analyze_linelist,'linelist_h2o.inp',temp,akms,lambdamic,a0,j0
analyze_linelist,'linelist_co2.inp',temp,akms,lambdamic,a0,j0
nlin = n_elements(lambdamic)

plot,lambdamic,a0,/xl,/yl,ps=6,/nodata
i=0LL
while i lt nlin do begin
   oplot,[lambdamic[i],lambdamic[i]],[1d-50,a0[i]]
   i = i+1LL
endwhile

;plot,lambdamic,j0,/xl,/yl,ps=6,/nodata
;i=0LL
;while i lt nlin do begin
;   oplot,[lambdamic[i],lambdamic[i]],[1d-50,j0[i]]
;   i = i+1LL
;endwhile

;plot,ener/ev>1d-5,npop,/xl,/yl

;plot,ener/(kk*temp)>1d-2,npop>1d-10,/xl,/yl

end
