pro dr_cbc_star, obj_nm=obj_nm, obsnm=obsnm, date=date, templ_tag=templ_tag, $
                 templ_nm=templ_nm, vdtag=vdtag, demo=demo, div_telluric=div_telluric,$
                 ddir=ddir, excalibur=excalibur, mincts=mincts

; dr_cbc_star, obj_nm='101501', templ_tag='am10', vdtag='am10a', ddir='nn' 

path = '/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'
if strlowcase(obj_nm) eq 'sun' then path = ddir


  ff = file_search(path+obj_nm+'_??????.????.fits', count=nobs)
  if keyword_set(obsnm) then $
     ff = file_search(path+obj_nm+'_'+obsnm+'.fits', count=nobs)
  if keyword_set(date) then $
     ff = file_search(path+obj_nm+'_'+date+'*.fits', count=nobs)

  if ~keyword_set(templ_nm) then templ_nm = obj_nm+'_templ_'+templ_tag+'.dat'  ; morphed template 

  if ~keyword_set(demo) then demo=0
  if ~keyword_set(div_telluric) then div_telluric=0

  for i = 0, nobs-1 do begin 
     x1=strpos(ff[i],'_', /reverse_search)
     obsnm = strmid(ff[i], x1+1, 11)
;     print, templ_nm, '   ',obsnm
     if strlowcase(obj_nm) eq 'sun' then begin
        snr=ck_snr(dir=path, obnm = obsnm)
        if snr lt 300 then stop,'reject and remove this observation by hand?'
        print, obsnm, ' ', snr
     endif
     xfound = file_search(path+'vd'+vdtag+'_'+obj_nm+'.'+obsnm+'.dat', count=nfound) 

     if nfound eq 1 then print, 'Found: ',xfound
     if nfound eq 0 then begin

     ; SETUP SPECTRA INTO VD CHUNKS
        if ~keyword_set(vdtag) then vdtag=templ_tag 
       cbc_chunk_setup, obj_nm, obsnm, templ_nm, templ_tag=templ_tag, $
                         vdtag=vdtag, ddir=ddir, excalibur=excalibur, div_telluric=div_telluric
                                ; RUN L-M FITTING OF MODEL B/T TEMPLATE AND OBS 
                                ; FLAGGING TELLURIC LINES 
       cbc_main, obj_nm, obsnm, templ_nm, demo=demo, vdtag=vdtag, ddir=ddir, excalibur=excalibur
       print, 'Ran: '+'vd'+vdtag+'_'+obj_nm+'.'+obsnm+'.dat'
    endif 
  end                           ; loop through nobs

  cbc_vank, obj_nm, vdtag, ddir=ddir, mincts=mincts


end   ; pro
