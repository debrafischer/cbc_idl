pro dr_all, obj_nm=obj_nm, tag_in=tag_in, ddir=ddir, start=start, $
            thresh=thresh, ngau=ngau, night=night, all=all, demo=demo, div_telluric=div_telluric, $
            excalibur=excalibur, mincts=mincts

; TAG_IN 'amz' for exalibur with poly-fitting
; DDIR directory name for FITS and VD files default 'excalibur'
; NIGHT '210104' finds and runs all observations for that night
; THRESH threshold for significant residuals for solar morph code 
; NGAU number of Gaussians for solar morph code
; DEMO plot template construction
; DIV_TELLURIC divide out telluric model in FITS files of program observations
; EXCALIBUR use LZhao's excalibur wavelengths with polynomial
;           fitting of PCA amplitudes (nearest neighbor weighting)
; MINCTS mininum count cut used in cbc_vank for including RVs

; dr_all, tag_in='amz', night='210219', /excalibur, dir='v3'
; dr_all, tag_in='amz', obj_nm='10700', /excalibur, dir='v3'

if ~keyword_set(tag_in) then tag_in='amz'
if ~keyword_set(ddir) then ddir='/Volumes/G/expres/extracted/v3/'
if ~keyword_set(thresh) then thresh=0.01 
if ~keyword_set(ngau) then ngau=130
if keyword_set(demo) then demo=1 else demo=0 

if ~keyword_set(excalibur) then begin
   ansb='y'
   read, 'Default is DM-9. Use excalibur? (y/n) ',ansb
   if ansb eq 'y' then begin
      excalibur = 1
      print, 'Using excalibur wavelengths for RVs' 
      tag_in = 'amz' 
   endif
   if ansb eq 'n' then begin
      excalibur = 0
      print, 'Using DM-9 wavelengths for RVs' 
      tag_in = 'amd' 
   endif
endif

tag_out = tag_in+strcompress(string(ngau),/remove_all)

restore,'cat_drive.dat'    ; cat structure
;xw=where(strmid(cat.comment,0,4) eq 'Wait', nxw)
;if nxw gt 0 then print, 'Skipping stars with comment = Wait;',cat[xw].starnm 
xx=where(cat.comment ne 'Wait',nstars)
ncat = cat[xx]

if keyword_set(all) and ~keyword_set(start) then begin
   starlist = ncat.starnm
   nstars = n_elements(starlist)
endif

if keyword_set(all) and keyword_set(start) then begin
   xnbr = n_elements(ncat)
   xbeg = where(ncat.starnm eq start)
   ncat = ncat[xbeg:xnbr-1]
   starlist = ncat.starnm
   nstars = n_elements(starlist)
endif

if keyword_set(obj_nm) then begin
   xx=where(ncat.starnm eq obj_nm, nstars)
   ncat = ncat[xx]
   starlist = ncat.starnm
   nstars = n_elements(ncat)
endif

if keyword_set(night) then begin  ; find stars observed on given night
   ss = file_search('/Volumes/G/expres/extracted/'+ddir+'/*/*'+night+'*.fits', count = ncount)
   namestarlist = strarr(ncount) 
   for j=0,ncount-1 do begin
      s1 = strpos(ss[j], '/',/reverse_search) 
      s2 = strpos(ss[j], '_',/reverse_search) 
      namestarlist[j] = strmid(ss[j], s1+1, (s2-s1)-1)
   endfor
   ii=uniq(namestarlist) 
   namestarlist = namestarlist[ii] 
   nstars = n_elements(namestarlist) 
   for k=0,nstars-1 do begin
      xx=where(namestarlist[k] eq cat.starnm, nxx) 
      if nxx eq 0 then namestarlist[k] = '0' 
   endfor 
   xkp = where(namestarlist ne '0',nstars)
   starlist = namestarlist[xkp]
   print, starlist
   dum = ncat[0]                     ; seed the structure
   for j = 0, nstars-1 do begin
      xx=where(ncat.starnm eq starlist[j],nxx)
      if nxx eq 0 then print, 'no match for starlist', starlist[j] 
      if nxx ne 0 then dum = [dum, ncat[xx]]
   endfor
   ncat = dum[1:n_elements(dum)-1]   ; remove the seed
   nstars = n_elements(ncat)
   starlist = ncat.starnm
endif 

;print, nstars 
;stop

for i=0, nstars-1 do begin 
   obj_nm = ncat[i].starnm 
   print, obj_nm 
   ffound1 = file_search('/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'+obj_nm+'_templ_'+tag_in+'.dat', count=c1)
   ffound2 = file_search('/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'+obj_nm+'_templ_'+tag_out+'.dat', count=c2)
   if c1 + c1 eq 1 then stop,'what?'

; NO MORPHED TEMPLATE - MAKE ONE
    if c1 eq 0 and c2 eq 0 then begin
      xt = where(ncat[i].templ_obnm ne 'n', nxt)
      templ_obnm = ncat[i].templ_obnm[xt]
      if ncat[i].morph_rv ne 0.0 then templ_rv = ncat[i].morph_rv
        ;if obj_nm eq 'TOI-1518' then templ_obnm = '200803.'+strcompress(string(indgen(35)+1093),/remove_all)
      print, 'template obs: ',templ_obnm
      if nxt eq 2 then t_obj_nm = templ_obnm[0]
      if nxt eq 3 then t_obj_nm = templ_obnm[1]
      if nxt eq 4 then t_obj_nm = templ_obnm[1]
      if nxt eq 5 then t_obj_nm = templ_obnm[2]
      date = strmid(t_obj_nm, 0, 6)
      print, 'Running: ',obj_nm
      coadd_templ_nm = obj_nm+'_coadd_'+date+'.dat'
      rayclean, templ_obnm, obstack, star, starnm=obj_nm, outfname=coadd_templ_nm, ddir=ddir,$
                excalibur=excalibur, /auto
      cbc_star_templ, obj_nm, t_obj_nm, tag=tag_in, coadd_obnm=coadd_obnm, $
                      coadd_templ_nm=coadd_templ_nm, ddir=ddir, excalibur=excalibur
      cbc_templ_morph, obj_nm, t_obj_nm, tag_in=tag_in, tag_out=tag_out, ngau=ngau, templ_rv=templ_rv, $
                       thresh=thresh, /coadd, ddir=ddir, excalibur=excalibur,fix_rotbro=fix_rotbro, $
                       init_rotbro=init_rotbro, fix_rv=fix_rv, demo=demo
   endif 

   ; ANALYZE PROGRAM OBSERVATIONS
   if ~keyword_set(mincts) then mincts=10000.
   dr_cbc_star, obj_nm=obj_nm, templ_tag=tag_out, vdtag=tag_in, $
                ddir=ddir, excalibur=excalibur, demo=demo, div_telluric=div_telluric, $
                mincts=mincts

endfor


end 
