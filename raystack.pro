pro  raystack,strr,filter,ordr,snr_arr=snr_arr, auto=auto,mdtck=mdtck, px0=px0, px1=px1

   ordr  = [ordr]  &  n_ord=n_elements(ordr)-1
   col   = n_elements(strr[*,0,0])-1
   nobs  = n_elements(strr[0,0,*])
   mdcts = fltarr(nobs)
   nstr  = fltarr(col+1,nobs) 
   nflt  = nstr
   cotb  = [221,141,191,111,81,41]
   plord = 3
   if col gt 1000 then plord = 5

   for n=0,n_ord do begin 
;      ;for q=0,nobs-1 do mdcts[q]=median(strr[*,ordr[n],q])
;      mdobs=first_el(where(mdcts eq median(mdcts)))
      mdobs = where(snr_arr^2 eq median(snr_arr^2), nmdobs)
      sp1=reform(strr[*,ordr[n],mdobs]) ;median observation
      sp1 = sp1[px0:px1]
      for q=0,nobs-1 do begin             ;normalize observations to median observation
         sp2=reform(strr[*,ordr[n],q])
         sp2 = sp2[px0:px1]
         fltr=reform(filter[*,ordr[n]])
	 divdif=sp1/sp2   ; bug???  28Jul99  RPB
         if n_elements(where(fltr gt 0)) gt 100 then for m=0,9 do begin
;toss the 10 greatest and 10 least pixels
;	   divdif=sp1/sp2  ; bug???  28Jul99  RPB
           dum=maxloc(divdif,/first)
           fltr[dum]=-1.
	   divdif[dum]=median(divdif)
           dum=minloc(divdif,/first)
           fltr[dum]=-1.
	   divdif(dum)=median(divdif)
         endfor ;m
         ind=where(fltr gt 0)
	 cof = pl_fit(ind,divdif[ind],plord)
         if n_elements(cof) eq 1 then begin   ; rare BOMB condition, 28Jul99 PB
print,'Order '+strtrim(ordr[n],2)+': Using simply ratio to compare spectra'
            nflt[*,q]=fltarr(col+1)*0.+median(divdif[ind])
         endif else nflt[*,q] = poly_fat(findgen(col+1),cof)
	 nstr[*,q]=reform(strr[*,ordr[n],q])*reform(nflt[*,q])
      endfor ;q                   
      rayslave,nstr,ordr(n),auto=auto,mdtck=mdtck
      for q=0,nobs-1 do strr[*,ordr[n],q]=nstr[*,q]/nflt[*,q]
   endfor ;n=0,n_ord do begin 

return
end
