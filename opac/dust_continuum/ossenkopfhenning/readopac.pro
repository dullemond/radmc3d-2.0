;------------------------------------------------------;
; READ THE DUST OPACITY FILES                          ;
;------------------------------------------------------;
function readopac,spec=spec,used=used
close,1
if not keyword_set(spec) then spec=1
filename = 'dustkappa_'+strcompress(string(spec),/remove_all)+'.inp'
if keyword_set(used) then filename=filename+'.used'
print,"Reading ",filename
openr,1,filename
icnt=0
str=' '
flag=0
while flag eq 0 do begin
   readf,1,str
   i=0
   while strmid(str,i,1) eq ' ' do i=i+1
   if((strmid(str,i,1) eq '#') or (strmid(str,i,1) eq ';') or (strmid(str,i,1) eq '!')) then begin
      icnt=icnt+1
   endif else begin
      flag = 1
   endelse
endwhile
close,1
openr,1,filename
for i=1,icnt do readf,1,str
iformat=0
readf,1,iformat
nf=0
readf,1,nf
case iformat of
   1: ncol=2
   2: ncol=3
   3: ncol=4
endcase
data=dblarr(ncol,nf)
readf,1,data
close,1
lambda    = dblarr(nf)
cabs      = dblarr(nf)
csca      = dblarr(nf)
heny      = dblarr(nf)
lambda[*] = data[0,*]
cabs[*]   = data[1,*]
if iformat ge 2 then csca[*]   = data[2,*]
if iformat ge 3 then heny[*]   = data[3,*]
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
freq      = 1d4*cc/lambda
;
return,{nf:nf,ns:1,freq:freq,wave:lambda,cabs:cabs,csca:csca,nrt:1,trange:0.d0}
end




;------------------------------------------------------;
; PLOT THE DUST OPACITIES                              ;
;------------------------------------------------------;
pro plotopac,a,abs=abs,scat=scat,is=is,oplot=oplot,$
             ylin=ylin,irt=irt,mult=mult,xlin=xlin,$
             _extra=_extra
ipl=0
ylog=1
xlog=1
if not keyword_set(irt) then irt=0
if irt ge a.nrt then begin
   print,'irt too large'
   return
endif
xtitle='!4k!X [!4l!Xm]'
ytitle='!4j!X!D!4k!X!N [cm!U2!Ng!U-1!N]'
if keyword_set(xlin) then xlog=0
if keyword_set(ylin) then ylog=0
if not keyword_set(oplot) then oplot=0
if not keyword_set(mult) then mult=1.d0
if keyword_set(abs) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
	    if oplot eq 0 then begin
                plot,a.wave,mult*a.cabs(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.cabs(*,is,irt),_extra=_extra
            endelse
        endif else begin
            dum = mult*a.cabs(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.cabs(*,is,irt)
            endfor
	    if oplot eq 0 then begin
               plot,a.wave,dum,ylog=ylog,xlog=xlog,_extra=_extra
	    endif else begin
               oplot,a.wave,dum,_extra=_extra
            endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.cabs(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.cabs(*,0,irt),_extra=_extra
	endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.cabs(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if keyword_set(scat) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
            if oplot eq 0 then begin
                plot,a.wave,mult*a.csca(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.csca(*,is,irt),_extra=_extra
	    endelse
        endif else begin
            dum = mult*a.csca(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.csca(*,is,irt)
            endfor
            if oplot eq 0 then begin
                plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,dum,_extra=_extra
	    endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.csca(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.csca(*,0,irt),_extra=_extra
        endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.csca(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if ipl eq 0 then begin
    dum = mult*a.csca(*,0,irt) + mult*a.cabs(*,0,irt)
    for is=1,a.ns-1 do begin
        dum = dum + mult*a.csca(*,is,irt) + mult*a.cabs(*,is,irt)
    endfor
    if oplot eq 0 then begin
       plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
    endif else begin
       oplot,a.wave,dum,_extra=_extra
    endelse
endif
end
