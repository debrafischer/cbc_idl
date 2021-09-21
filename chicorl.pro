pro   chicorl, star,template,pnt,np,shtem,plot=plot

if n_params() lt 4 then begin
	print,'CHICORL, star,template,pnt,range,tempsh'
return 
endif

;Like Cross-correlation routine, but instead this finds
;the shift to nearest pixel such that chi-square is minimized.
;PNT is the shift of TEMPLATE relative to STAR.
;NP is +/- shift range, SHTEM is the shifted version of template
;Precursor to MINCORL

	cs = fltarr(2*np+1)
	clip = np
	ln = n_elements(star)
	amean = mean(star(clip+np:ln-(clip+1)-np))
	bmean = mean(template(clip+np:ln-(clip+1)-np))
	al = star/amean
	b = template(clip:ln-clip-1)/bmean
;So b is shorter than al by clip*2 pixels.
;This is where the initial chi-squared function gets computed
	for n=0,2*np do begin
		a = al(clip-np+n:ln-(clip+1)-np+n)   ;shift star
		d = (b-a)
		cs(n) = total(d*d)		     ;cs is chi sq.
	endfor
	q1 = where(cs eq min(cs))
	pnt = q1(0)-np
	shfour, template,pnt,shtem		     ;call shfour

if keyword_set(plot) then begin
    loadct,39
    plot,template,/xsty,/ysty,/nodata,col=0
    oplot,template,col=255
    oplot,shtem,col=240
stop
endif

end

