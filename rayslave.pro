pro  rayslave,nstr,ordr,auto=auto,mdtck=mdtck

   dumtalk='n'
   col   = n_elements(nstr[*,0])-1
   nobs  = n_elements(nstr[0,*])
   divdif= fltarr(col+1)
;   cotb  = [221,141,191,111,81,41,61,101,121,161,201]  ;colors, assumes color_table 13
;   autobad = 1600.          ;differences of 1600 are automatically bad
;   autobad = 5000.          ;PB Kludge Dec 12 1996 for Keck
   autobad = 1200.          ;PB Kludge Dec 12 1996 for Keck
   if max(nstr) lt 5.0 then autobad = 0.2  ; normalized spectra

   if n_elements(mdtck) ne 1 and keyword_set(auto) then mdtck=5.0
   if n_elements(mdtck) ne 1 and not keyword_set(auto) then mdtck=4.5
   if mdtck gt 5 then mdtck = 5.0
   if max(nstr) lt 5.0 then mdtck = 7.5   ; normalized spectrum

      for q=0,col do begin
	 bd=maxloc(abs(nstr[q,*]-median(nstr[q,*])),/first)
	 gd=indgen(nobs)
	 gd=gd(where(gd ne bd))
	 divdif[q]=abs(nstr[q,bd]-mean(nstr[q,gd])) 
      endfor  ;q
      mdstdev=median(divdif)
      ind=where(divdif gt autobad,n_ind)  ;automatically bad pixels
      if n_ind gt 0 then for count=0,n_ind-1 do begin
          bdpx=ind[count] 
          bdobs=maxloc(abs(nstr[bdpx,*]-median(nstr[bdpx,*])),/first)
	  gdobs=indgen(nobs)
	  gdobs=gdobs(where(gdobs ne bdobs))
          fpix=max([0,bdpx-20])
          lpix=min([col,bdpx+20])
          xx=indgen(lpix-fpix+1)+fpix
          if not keyword_set(auto) then begin
              ypad=1000
              if max(nstr[fpix:lpix,0]) lt 2000 then ypad=200
              if max(nstr[fpix:lpix,0]) lt 2 then ypad=0.1  ; normalized spectrum
	      plot,xx,nstr[fpix:lpix,0],/xstyle, $
                   title='Order # '+strtrim(ordr,2), $
                 yra=[max([min(nstr[fpix:lpix,0])-ypad,0]), $
                      max(nstr[fpix:lpix,0])+ypad]
              for q=1,nobs-1 do oplot,xx,nstr[fpix:lpix,q],co=q*10
              talk='Bad Pixel at order '+strtrim(ordr,2)
              talk=talk+', pixel '+strtrim(bdpx,2)
              print,talk
	      oplot,[bdpx],[nstr[bdpx,bdobs]],psym=1,symsize=1.5,thick=1.5,co=151
              oplot,[bdpx],[mean(nstr[bdpx,gdobs])],psym=4,thick=2.5,symsize=1.5,co=211
	      wait,0.1
          endif
          nstr[bdpx,bdobs]=mean(nstr[bdpx,gdobs])
      endfor   ;count

      for q=0,col do begin
	 bd=maxloc(abs(nstr[q,*]-median(nstr[q,*])),/first)
	 gd=indgen(nobs)
	 gd=gd(where(gd ne bd))
	 divdif[q]=abs(nstr[q,bd]-mean(nstr[q,gd])) 
      endfor  ;q
      mdstdev=median(divdif)
      ind=where(divdif gt (mdtck*mdstdev),n_ind)
;also add any negative pixels  5August99  PB
      if n_ind gt 0 then begin 
         for q=0,col do if min(nstr[q,*]) le 0 then ind=[ind,q]
         ind = where(histogram(ind) gt 0) + min(ind)
         ind = ind(where(ind gt 0,n_ind))
      endif
;end add any negative pixels  5August99  PB
      if n_ind gt 0 then for count=0,n_ind-1 do begin
          bdpx=ind[count] 
          bdobs=maxloc(abs(nstr[bdpx,*]-median(nstr[bdpx,*])),/first)
; the following line forces negative pixels to be the bad observation
          if min(nstr[bdpx,*]) lt 0 then bdobs = minloc(nstr[bdpx,*],/first)
	  gdobs=indgen(nobs)
	  gdobs=gdobs(where(gdobs ne bdobs))
          fpix=max([0,bdpx-20])
          lpix=min([col,bdpx+20])
          xx=indgen(lpix-fpix+1)+fpix
	  if not keyword_set(auto) then begin
             ypad=1000
             if max(nstr[fpix:lpix,0]) lt 2000 then ypad=200
             if max(nstr[fpix:lpix,0]) lt 2. then ypad=0.1
             plot,xx,nstr[fpix:lpix,0],/xstyle, $
                 yra=[max([min(nstr[fpix:lpix,0])-ypad,0]), $
                      max(nstr[fpix:lpix,0])+ypad], $
                 title='Order # '+strtrim(ordr,2)
             for q=1,nobs-1 do oplot,xx,nstr[fpix:lpix,q],co=q*5
             talk='Possible cosmic ray hit at order '+strtrim(ordr,2)
             talk=talk+', pixel '+strtrim(bdpx,2)
             print,talk
	     oplot,[bdpx],[nstr(bdpx,bdobs)],psym=1,symsize=1.5,thick=1.5,co=151
             oplot,[bdpx],[mean(nstr(bdpx,gdobs))],psym=4,thick=2.5,symsize=1.5,co=211
             read,'Do you wish to fix the marked (+) pixel (y/n)? ',dumtalk
             if dumtalk eq 'y' or dumtalk eq '1' $
                then nstr[bdpx,bdobs]=mean(nstr[bdpx,gdobs])
          endif else nstr[bdpx,bdobs]=mean(nstr[bdpx,gdobs])   ;auto bad pixel
      endfor   ;count

return
end
