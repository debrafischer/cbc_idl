;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
; Pro CBC_CHUNK_SETUP  
;   
; PROCEDURES CALLED: CBC_INIT, RDNSO
;     
; PURPOSE:   
;   Setup chunk structure for observations
;    
; INPUT: 
;
;  
; OPTIONAL INPUT: 
;
; OUTPUT:
;  VD: IDL structure holding chunks of the observations
; 
; PROCEDURE:  
;   1. collect observations
;   2. call cbc_init to setup cbcenv structure with global parameters
;      all the EXPRES parameters are hard-wired there for consistency
;      in all the cbc programs
;   3. read template structure with padded chunks
;   4. format with the same pixels and wavelength range as EXPRES order
;      divide into 160-pixel chunks with 40-pixel padding on each side 
;   5. write observations 
; 
; CALLING SEQUENCE
;   cbc_101501_setup, '190518.1122', '101501_templ_am10.dat', templ_tag=templ_tag
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

pro cbc_chunk_setup, obj_nm, obsnm, templ_nm, templ_tag=templ_tag, $
                     vdtag=vdtag, div_telluric=div_telluric, ddir=ddir, excalibur=excalibur

; COLLECT OBSERVATIONS 
  osamp = 1                       ; osamp of spectra
  date = strmid(obsnm,0,6) 

; CALL CBC_INIT TO CREATE CBCENV
; WITH GLOBAL AND EXPRES-SPECIFIC PARAMETERS
; WORK IN PIXEL SPACE FOR CONSTANT ARRAY SIZES
; SETUP YOUR DIRECTORIES AND RENAME THE PATHS BELOW
;     ORD 42 [5154, 5218]  ORD 75 [7156, 7245]  ; epoch eq 5
;     ORD 35 [4865, 4926]  ORD 78 [7418, 7510]  ; epoch le 4
;     
  obnm = obj_nm+'_'+obsnm          ; obsname to read
  if ~keyword_set(templ_nm) then stop
; get cbcenv
  cbcenv = cbc_init(obsnm, templ_nm, obj_nm=obj_nm, osamp=osamp, ddir=ddir, excalibur=excalibur)  
  c_light = cbcenv.c_light        ; speed of light
  templ_pad = cbcenv.templ_pad    ; pixels for padding template
  npix_chunk = cbcenv.n_pix       ; 160  pixels on each vd chunk   
  npix_ord = cbcenv.npix_ord      ; 6600 per extracted order 
;  edge_pad = cbcenv.edge_pad      ; number of pixels to skip at blaze edges  
  ord_init = cbcenv.st_ord        ; initial order for CBC (LFC) 
  ord_final = cbcenv.fin_ord      ; final order for CBC
  nord = ord_final - ord_init + 1 ; number of orders for CBC
  nchunk_ord = cbcenv.nch_ord     ; number of chunks per order 
  nchunks = cbcenv.n_chunks       ; 40 chunks / order * nord (2160)
  epoch = cbcenv.epoch            ; instrument status epoch 
  px0_init = cbcenv.st_pix   

; RESTORE TEMPLATE FILE
  restore, cbcenv.templ_dir+templ_nm ; templ: restore template for CBC

; RESTORE OBSERVATIONS AND LOAD INTO VD STRUCTURE
  fname = cbcenv.obs_dir+obnm+'.fits'
  spec = mrdfits(fname, cbcenv.file_ext, header, /silent)
  jd = strcompress(sxpar(header, 'BARYMJD'),/remove_all)

  vd = {obnm:'',                                         $
        blaze:dblarr(npix_chunk),                        $
        contin:dblarr(npix_chunk),                       $  ; continuum
        cts: 1000L,                                      $  ; max((sobs/unc)^2)
        cont_offset:0.0,                                 $  ; vertical scaling factor
        epoch:epoch,                                     $  ; instrument epoch
        gdpix:intarr(npix_chunk),                        $
        jd:0.000d,                                       $  ; julian date of observation
;        keep:0,                                          $  ; retain chunk in vank
        npix_chunk:npix_chunk,                           $
        ord:ord_init,                                    $  ; relative order
        perror:dblarr(2),                                $  ; from mpfit
        pixt:cbcenv.st_pix,                              $  ; first pixel for template
        pixob:cbcenv.st_pix,                             $  ; first pixel for observation (drifts wrt template) 
        red_chi:0.0,                                     $  
        rms:0.0,                                         $  ; stddev(abs(smod-sobs))
        smod:dblarr(npix_chunk),                         $  ; model spectrum 
        sobs:dblarr(npix_chunk),                         $  ; observed spectrum
        stel:dblarr(npix_chunk),                         $  ; telluric spectrum 
        stmp:dblarr((2.*templ_pad)+npix_chunk),          $  ; padded template spectrum
        stmp_morph:dblarr((2.*templ_pad)+npix_chunk),    $  ; padded morphed template spectrum
        templ_wt:0.0d,                                   $  ; template weight (line slopes) 
        templ_nm: templ_nm,                              $  ; file name of template used in CBC model   
        unc:dblarr(npix_chunk),                          $  ; RP uncertainties
        wav_lab:dblarr(npix_chunk),                      $  ; lab wavelengths for tellurics
        wav:dblarr(npix_chunk),                          $  ; bary-corrected wavelengths
        wtel:dblarr(npix_chunk),                         $  ; 
        wtmp:dblarr((2.*templ_pad)+npix_chunk),          $  ; padded template wavelengths
        vel:0.0d,                                        $  ; fitted RV
        weight:0.0d}
 
 vd = replicate(vd, nchunks)

; FIND PIXOB, WOBS AND SOBS
  for i = 0, nord-1 do begin
     ord = ord_init + i         ; order index
     sp = double(spec[ord].spectrum)                         ; stellar spectrum - NOT NORMALIZED
     contin = double(spec[ord].continuum)                    ; continuum: spectrum / continuum = normalized spec
     stel = spec[ord].tellurics
     wtel = spec[ord].wavelength
     blz1 = spec[ord].blaze
     blz1_trim = blz1[px0_init:npix_ord-1]
     contf, blz1_trim, c, sbin=50, nord=6
     blz1[px0_init:npix_ord-1] = c                           ; smoothed blaze function (remove QE variations)
     if keyword_set(excalibur) then $
        wtel = spec[ord].excalibur                           ; excalibur lab wavelength
     unc = double(spec[ord].uncertainty)                     ; RP uncertainty - sp/unc to get SNR
     wv = double(spec[ord].bary_wavelength)                  ;;     barycentric wavelengths
     if keyword_set(excalibur) then $
        wv = double(spec[ord].bary_excalibur)
     if max(wv) - min(wv) eq 0 then stop,'need BC wavelengths' 
     wv_lb = double(spec[ord].wavelength)                    ; lab frame wavelengths
     if keyword_set(excalibur) then $                        ; if excalibur, replace with exalibur lab wav
        wv_lb = double(spec[ord].excalibur)

;     if ord eq 41 and keyword_set(excalibur) then begin
;        nchunk_ord = 37
;        px0_init = 1470
;     endif
     if ord gt 41 and keyword_set(excalibur) then begin
        nchunk_ord = cbcenv.nch_ord
        px0_init = cbcenv.st_pix
     endif

     for j = 0, nchunk_ord-1 do begin                        ; loop thru chunks in each order
        chnk = j + (i*nchunk_ord)                            ; chunk index
;        if ord_init eq 41 and keyword_set(excalibur) and i ge 1 then chnk = j + ((i-1)*nchunk_ord) + 37
        vd[chnk].obnm = obsnm                                ; obnm (e.g. 190704.1003)
        vd[chnk].jd = double(jd)-40000.+0.5                  ; BARYMJD from ext 1 FITS header => RJD
        vd[chnk].ord = ord                                   ; relative order
        vd[chnk].wtmp = templ[chnk].wav                      ; templ wavelength
        if strmid(templ_tag,1,1) eq 'm' then $
           vd[chnk].stmp_morph = templ[chnk].morph_spec      ; morphed templ spec
        vd[chnk].stmp = templ[chnk].spec                     ; not morphed 
        pix0 = px0_init + (j*npix_chunk)                     ; pix0 of unshifted, unpadded templ chunk
        pix1 = pix0 + npix_chunk-1                           ; last pix of unshifted, unpadded templ chunk
        disp = (wv[pix1] - wv[pix0])/npix_chunk              ; template dispersion: dlambda / dpix

        ; WAVELENGTH OF SHIFTED OBSERVATION (SPECTRA DRIFT ON DETECTOR) 
;        d_lambda1 = mean([templ[chnk].wav[templ_pad:templ_pad+npix_chunk] - wv[pix0:pix1]]) 
        d_lambda = templ[chnk].wav[templ_pad] - wv[pix0]
        dpix = round(d_lambda / disp)                        ; pixel offset for observed chunk
        vd[chnk].pixt = templ[chnk].pixt                     ; unpadded first pixel in chunk template (aka pix0)
        vd[chnk].pixob = vd[chnk].pixt + dpix                ; pixob is location of pix0 in obs
        p0_obs = pix0 + dpix                                 ; pixel location of observed spectrum
        p1_obs = p0_obs + npix_chunk-1                       ; pixel location of observed spectrum
        vd[chnk].wav_lab = double(wv_lb[p0_obs:p1_obs])      ; lab wavelength solution of obs
        vd[chnk].wav = double(wv[p0_obs:p1_obs])             ; bc-corrected wavelength of obs
        vd[chnk].sobs = sp[p0_obs:p1_obs]/contin[p0_obs:p1_obs]  ; normalized observed spectrum           
        vd[chnk].blaze = blz1[p0_obs:p1_obs]                ; blaze function
        vd[chnk].unc = unc[p0_obs:p1_obs]                    ; RP uncertainty for chunk
        vd[chnk].contin = contin[p0_obs:p1_obs]
;print, pix0, dpix, chnk, p0_obs, vd[chnk].wtmp[templ_pad],
; vd[chnk].wav[0], 
; disp, (vd[chnk].wav[129]-vd[chnk].wav[0])/npix_chunk
;print, d_lambda1, d_lambda
        ; LOCATION OF TELLURIC CHUNK
;;        vd[chnk].wtel = double(wtel[pix0:pix1])          ; tel wavelengths for tell chunk
;;        vd[chnk].stel = stel[pix0:pix1]                  ; telluric spectrum for tell chunk
        vd[chnk].wtel = double(wtel[p0_obs:p1_obs])          ; tel wavelengths for tell chunk
        vd[chnk].stel = stel[p0_obs:p1_obs]                  ; telluric spectrum for tell chunk
        if keyword_set(div_telluric) then begin
;           if min(vd[chnk].stel) gt 0.5 then $
           vd[chnk].sobs = vd[chnk].sobs/vd[chnk].stel
        endif

        vd[chnk].templ_wt = templ[chnk].weight           ; chunk weight for template
        vd[chnk].weight = templ[chnk].weight             ; initial weight for chunk 

      ; WHAT ARE THE CONTINUUM COUNTS?
        allcounts = (sp[p0_obs:p1_obs] / vd[chnk].unc)^2 ; photon counts
        ii = sort(allcounts)  
        allcounts=allcounts[ii] 
        revall=reverse(allcounts)
        counts = median(revall[0:10])                    ; median of top ten points in spectrum
        vd[chnk].cts = counts
        vd[chnk].gdpix[*] = 1                            ; all pixels initially marked as good
     endfor                                              ; step through chunks in each order
  endfor                                                 ; step through each order

  if ~keyword_set(vdtag) then $
     fname = 'vd'+templ_tag+'_'+obj_nm+'.'+obsnm+'.dat' ;'vd_NSO.'+obnm+'.dat'
  if keyword_set(vdtag) then $
     fname = 'vd'+vdtag+'_'+obj_nm+'.'+obsnm+'.dat' ;'vd_NSO.'+obnm+'.dat'

  save, vd, f=cbcenv.obs_dir+fname

end                             ; pro
