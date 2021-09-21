pro  rayzap,sp1,sp2,filter,ordr,mindiff,median=median, px0=px0, px1=px1

   if n_params() lt 5 then mindiff=900
   if keyword_set(median) and mindiff gt 5 then mindiff=5
   dumtalk='n'
   ordr=[ordr]  &  n_ord=n_elements(ordr)-1
   col=n_elements(sp1(*,0))-1
   for n=0,n_ord do begin 
      s1=reform(sp1[*,ordr[n]]) & s1 = s1[px0:px1]
      s2=reform(sp2[*,ordr[n]]) & s2 = s2[px0:px1] 
      fltr=reform(filter[*,ordr[n]])
      if n_elements(where(fltr gt 0)) gt 100 then begin
         for m=0,9 do begin       ;toss the 10 greatest and 10 least pixels
           ind=where(fltr gt 0)
           dum1=maxloc(s1[ind],/first)
           dum2=maxloc(s2[ind],/first)
           fltr(ind[dum1])=-1.
           fltr(ind[dum2])=-1.
           dum1=minloc(s1[ind],/first)
           dum2=minloc(s2[ind],/first)
           fltr(ind[dum1])=-1.
           fltr(ind[dum2])=-1.
         endfor ;m
         ind=where(fltr gt 0)
         dum1=mean(s1[ind])
         dum2=mean(s2[ind])
	 plord=3
	 if n_elements(s1) gt 1000 then plord=5
         xfltr = where(fltr gt 0)    ;good pixels
         if dum2 lt dum1 then begin
	    spud = s1/s2
;	    smcof=pl_fit(indgen(col+1),spud,plord)
	    smcof=pl_fit(xfltr,spud[xfltr],plord)
	    smfit=poly(indgen(col+1),smcof)
	    s2=s2*smfit
	 endif else begin
	    spud = s2/s1
;	    smcof=pl_fit(indgen(col+1),spud,plord)
	    smcof=pl_fit(xfltr,spud(xfltr),plord)
	    smfit=poly(indgen(col+1),smcof)
	    s1=s1*smfit
	 endelse
         dffr=abs(s1-s2)
	 mdffr=median(dffr)
;cosmic rays are assumed to cause at least mindiff dn difference
         if keyword_set(median) then $
             ind=where(dffr gt mindiff*mdffr and filter[*,ordr[n]] gt 0,n_ind) $
	     else ind=where(dffr gt mindiff and filter[*,ordr[n]] gt 0,n_ind)
;tack on one extra pixel on each side of bad pixels 
         if n_ind gt 0 then begin
	   ind = [ind-1,ind,ind+1]
	   ind = ind(rem_dup(ind))
	   ind = ind(sort(ind))
	   ind = ind(where(ind ge 0 and ind le col))
	   n_ind = n_elements(ind)
           for count=1,n_ind do begin
              bdpx=ind(count-1) 
              fpix=max([0,bdpx-20])
              lpix=min([col,bdpx+20])
              xx=indgen(lpix-fpix+1)+fpix
              plot,xx,s1(fpix:lpix),/xstyle, $
                 yra=[min(s1(fpix:lpix))-1000.,max(s1(fpix:lpix))+1000.], $
                 title='Order # '+strtrim(ordr(n),2)
              oplot,xx,s2(fpix:lpix),co=121
              oplot,[bdpx],[max([s1(bdpx),s2(bdpx)])],psym=1,co=121,symsize=1.5,thick=1.5
              oplot,[bdpx],[min([s1(bdpx),s2(bdpx)])],psym=4,co=221,symsize=1.5,thick=1.5
              talk='Possible cosmic ray hit at order '+strtrim(ordr(n),2)
              talk=talk+', pixel '+strtrim(bdpx,2)
              print,talk
              if s1(bdpx) gt s2(bdpx) then begin
                 read,'Do you wish to fix spectrum #1 (red line)  y/n/u? ',dumtalk
                 if dumtalk eq 'y' then s1(bdpx)=s2(bdpx)
                 if dumtalk eq 'u' then s2(bdpx)=s1(bdpx)
              endif else begin
                 read,'Do you wish to fix spectrum #2 (green line)  y/n/u? ',dumtalk
                 if dumtalk eq 'y' then s2(bdpx)=s1(bdpx)
                 if dumtalk eq 'u' then s1(bdpx)=s2(bdpx)
              endelse
           endfor
         endif
;         if dum2 lt dum1 then s2=s2*(dum2/dum1) else s1=s1*(dum1/dum2)
         if dum2 lt dum1 then s2=s2/smfit else s1=s1/smfit
         sp1(px0:px1,ordr(n))=s1
         sp2(px0:px1,ordr(n))=s2
      endif
   endfor ;n   
return
end
