pro pr_RV, star=star, tag=tag, dace=dace

;if ~keyword_set(star) then ff = file_search('vst*_'+tag+'.dat',count=nstars)
if keyword_set(star) then ff=file_search('vel/vst'+star+'_'+tag+'.dat', count=nstars)


for i=0, nstars-1 do begin
   restore,ff[i]
   s1 = strpos(ff[i],'_')
   starnm = strmid(ff[i],4,s1-4)
   fname = starnm+'_'+tag+'.txt'
   if keyword_set(dace) then fname = 'HD'+star+'_EXPRES_'+tag+'.rdb'
stop
   nobs=n_elements(cf3)
   openw,1,fname
   printf,1,"rjd","vrad","svrad","fwhm","sig_fwhm","bis_span","sig_bis_span","s_mw","sig_s","ha","sig_ha",$
          format='(a3,TR17,a4,TR5,a5, TR6, a4,TR7, a8,TR3,a8,TR2,a12,TR3,a4,TR3,a5,TR6,a2,TR6,a6)'
   printf,1,"---","----","-----","----","--------","--------","------------","----","-----","--","------",$
          format='(a3,TR17,a4,TR5,a5, TR6, a4,TR7, a8,TR3,a8,TR2,a12,TR3,a4,TR3,a5,TR6,a2,TR6,a6)'
;          format='(a3,TR5,a4,TR5,a5)'
   for j=0, nobs-1 do printf,1, cf3[j].jd+40000.d, cf3[j].mnvel,cf3[j].errvel,cf3[j].fwhm_ccf, "0.003",$
                             cf3[j].bis,"0.001",cf3[j].sval,"0.001",cf3[j].halpha_em,"0.001",$
                             f='(d13.6,TR3, f9.3, TR3, f6.3,TR3,f10.3,TR3,f8.3,TR3,f9.4,TR6,f6.3,TR3,f6.3,TR3,f5.3,TR3, f7.4,TR3,f7.4)'
                                ; rjd        vrad      svrad    fwhm   sfwhm     bis     sbis      sval      ssval     ha      sha
   close,1
endfor 


end
