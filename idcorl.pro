pro idcorl,ri,ti,idsh,filter,quality,fit,inc=inc,f_pnt=f_pnt
;  Measures the relative shift of a template spectrum with
;  respect to an another spectrum; output shift = outsh, with associated
;  chi sq value, minval. Must Supply a guess of the shift, insh good to .1 pxl.
;  If quality = 'bad', then minimum not found!
;  This code runs idcorlsh.pro.
	ln=n_elements(ri)
	test=0
	quality='good'
        index=where(filter gt 0.,num_good)
   if num_good ge 6 then begin
	a=ri(index)
	a=a/mean(a)
	if n_elements(inc) eq 0 then inc=0.02
	iodsh=idsh
	fti=fft(ti,1)
	sig=findgen(ln)/ln - .5
	sig=shift(sig,ln/2) * 2. * !pi
	idcorlsh,5,inc,a,fti,sig,index,iodsh,test
        if test eq 1 then idcorlsh,7,inc,a,fti,sig,index,iodsh,test
	if test eq 2 then begin
		   iodsh=idsh  &  idcorlsh,9,inc,a,fti,sig,index,iodsh,test
        endif
;	if test eq 3 then begin
;	   oldinc=inc
;          inc=inc*5. & iodsh=idsh & idcorlsh,9,inc,a,fti,sig,index,iodsh,test
;	   if test ne 4 then begin
;              inc=oldinc  &   idcorlsh,9,inc,a,fti,sig,index,iodsh,test
;          endif ;test ne 4
;        endif   ;test eq 3
	if test ge 3 then begin
	   if n_elements(f_pnt) eq 0 then f_pnt=0
	   lo=strtrim(index(0)+f_pnt,2)
	   hi=strtrim(index(num_good-1)+f_pnt,2)
	   quality='bad'
           talk='did not find minimum in idcorl, from pixel '
           talk=talk+lo+' to '+hi+',  inc = '+strtrim(inc,2)
           print,talk
        endif ;test eq 4

;calculating chi-sq fit
        idsh=iodsh
	shfour,ti,idsh,sti
	sti=sti(index)   &   sti=sti/mean(sti)
	d=sti-a
	fit=(total(d*d))/num_good
	fit=fit^0.5
   endif else begin
       print,'not enough good elements in idcorl'
       quality='bad'
       fit=9
   endelse
return
end ;idcorl
