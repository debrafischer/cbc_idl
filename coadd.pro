pro coadd,im1,im2,filt_pix,order,sh_pix,imadd,sh_im2,inc=inc,vd=vd
;im1 and im2 should be two different spectra of the same object.  This
;routine measures the pixel shift between the spectra and returns imadd
;which is the addition of spectra im1 and the shifted spectra of im2.
;sh_pix is the pixel displacement between im1 and im2.  (PB 5/25/88).
;Revised to accept 'benitz' style spectra on 8/8/88 (PB)
;bo=22 for iodine-cepheid \ bo=0 for everything else
;eo=22 for iodine spectra \ eo=24 for stellar spectra
;eo=52 for cepheid spectra \ eo=47 for iodine-cepheid spectra
	ord=[order]   &  n_ord=n_elements(ord)-1
	bo=ord(0)     &  eo=ord(n_ord)
        sh_pix=fltarr(n_elements(im1(0,*)))
	imdum=im2
	if n_elements(inc) ne 1 then inc=0.02 else inc=abs(inc)
;	chicorl,im1(*,bo),im2(*,bo),pnt,30
;print,'First guess shift (order '+strtrim(bo,2)+'): '+strtrim(pnt,2)
        for ord_ct=0,n_ord do begin
	   n=ord(ord_ct)
           sp1=im1(*,n)  &   sp2=im2(*,n)
	   filter=filt_pix(*,n)
	   chicorl,sp1,sp2,pnt,30
           idcorl,sp1,sp2,pnt,filter,check,inc=inc
	   if check eq 'bad' then begin
print,'In the bad-lands, order '+strtrim(n,2)+'   shift: '+strtrim(pnt,2)
		   if n gt (bo+3) then pnt=mean(sh_pix(n-4:n-1))
		   if n le (bo+3) and n gt bo then pnt=mean(sh_pix(bo:n-1))
		   if n eq bo then pnt=0
print,'Out of the bad-lands, order '+strtrim(n,2)+'   shift: '+strtrim(pnt,2)
	   endif
	   shfour,sp2,pnt,sh_i2
	   sh_pix(n)=pnt
	   imdum(*,n)=sh_i2
           talk='order '+strtrim(n,2)+' thru '+strtrim(eo,2)
	   print,talk+'     shift = '+strtrim(pnt,2)
	   if n_elements(vd) gt 0 then begin
	      ind=where(vd.order eq n,nind)
	      if nind gt 0 then vd(ind).pixel=vd(ind).pixel-pnt
           endif
        endfor
imadd=imdum+im1
sh_im2=imdum
return
end ;coadd
