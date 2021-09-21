pro ck_precision, starnm=starnm, tag=tag, fname=fname, verbose=verbose, all=all, pre=pre 

;ff = fname
if keyword_set(all) then ff=file_search('vel/vst*'+tag+'.dat', count=nstars)

if keyword_set(starnm) then begin
   fstar = 'vst'+starnm+'_'+tag+'.dat'
   ff = file_search('vel/'+fstar, count=nstars) 
   ff = ff[0]
   nstars=1
endif 
starnm = strarr(nstars)

;starnm='3651' 
;stop

;ndate=40
pre = {starnm:' ', nights:0, nobs:0, rms_nightly:0.0, expt:0.0, am:0.0, std_am:0.0, $
       rms_all:0.0, mean_err:0.0, mean_rmsnightly:0.0}
pre = replicate(pre, nstars) 

;if keyword_set(verbose) then $
;   print, '     Starnm', '     Date ','     Expt', '      AM', ' $
;   std(AM)', '  RMS_n' ;,$
;   '        RMS_all', '      <Err>'
for i=0, nstars-1 do begin
   restore,ff[i]   ; restore the vst
   indx = strpos(ff[i],'_')
   starnm[i] = strmid(ff[i],7, indx-7)
;   x5 = where(cf3.epoch eq 5, nx5)
;   cf = cf3[x5]
   cf = cf3
   date = strmid(cf.obnm,0,6)
   ii = uniq(date) 
   kdate = date[ii]
   ndate = n_elements(kdate)
   pre[i].starnm = starnm[i]
   pre[i].nights = ndate
   pre[i].nobs = n_elements(cf)
   pre[i].rms_all = stddev(cf.mnvel)
   pre[i].mean_err = mean(cf.errvel) 
;   print, starnm[i], ndate
   for j=0, ndate-1 do begin                               ; j loops thru dates 
      kk=where(strmid(cf.obnm, 0, 6) eq kdate[j],nkk)      ; kk loops through obs on date
;       if keyword_set(verbose) then print, $
;          starnm, '  ',kdate[j], mean(cf[kk].expt), mean(cf[kk].am),$
;          stddev(cf[kk].am), stddev(cf[kk].mnvel),$
;          f='(a10,a3,a8,f9.2, f9.2, f9.2, f9.3)'
;                                           stddev(cf.mnvel), mean(cf.errvel)      
       if nkk ge 2 then begin
          rms_nightly=stddev(cf[kk].mnvel) 
          expt = mean(cf[kk].expt)
          am = mean(cf[kk].am)
          std_am = stddev(cf[kk].am)
       endif
    endfor
   pre[i].rms_nightly = rms_nightly
   pre[i].expt = expt
   pre[i].am = am
   pre[i].std_am = std_am
;       if keyword_set(verbose) then $
;          for l = 0, nkk-1 do print, starnm, '  ',cf[kk[l]].obnm, stddev(cf[kk].mnvel), $
;                                    stddev(cf.mnvel), mean(cf.errvel)
;         if nkk ge 2 then begin
;            pre[i].rms_nightly[j]=stddev(cf[kk].mnvel) 
;         endif
;      endfor
;   endif
endfor  ; i


 
print, ' '
;openw, 1, 'starlist'
print, 'Star', ' Nobs ', ' Nnights ',' <err>', ' RMS_all ', ' <RMS_nightly> ', $
       f='(a8, a10, a7, a7, a8, a14)'
for i=0, nstars-1 do begin
;   x=where(pre[i].rms_nightly ne 0.0, nx)
;   if nx gt 1 then pre[i].mean_rmsnightly = mean(pre[i].rms_nightly[x])
   print, i, pre[i].starnm, pre[i].nobs, pre[i].nights, pre[i].mean_err, $
          pre[i].rms_all, pre[i].rms_nightly, $
          f='(i2,a10, i6, i6, f9.2, f7.2, f10.2)'
endfor
;close,1   
print,' ****************************************************************'
print, ' ' 


;stop
end

