;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
; PRO CBC_STAR_TEMPL  
;   
; PROCEDURES CALLED: 
;   CBC_INIT
;   MRDFITS
;     
; PURPOSE:   
;   Generate and begin to fill a template structure for CBC or morphing
;    
; INPUT: 
;  OBJ_NM: star name (required
;  OBSNM: 190531.1107 (dummy for coadded, e.g.: 190210.1153)
;  TAG: tag for template structure ("amd60") - "m" = will be morphed
;
; OPTIONAL INPUT
;  COADD_TEMPL_NM: 
;  OSAMP: pixel sampling of the template 
;  
; OUTPUT:
;  TEMPL: IDL structure holding chunks of the (co-added) template
; 
; PROCEDURE:  
;   1. call cbc_init to setup cbcenv structure with global parameters
;      all the EXPRES parameters are hard-wired there for consistency
;      in all the cbc programs
;   2. read NSO atlas to fill orders
;   3. format with same orders, pixels and wavelength range as EXPRES order
;   4. divide into 160-pixel chunks with 40-pixel padding on each side 
;   5. derive a weight for the chunk based on spectral line slopes 
; 
; CALLING SEQUENCE
;   cbc_star_templ, '101501','190531.1107', tag='amf' 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

pro cbc_star_templ, obj_nm, obsnm, tag=tag, osamp=osamp, demo=demo, ddir=ddir, $
                    excalibur=excalibur, coadd_templ_nm=coadd_templ_nm, coadd_obnm=coadd_obnm

osamp = 1

; CALL CBC_INIT TO CREATE CBCENV
; WITH GLOBAL AND EXPRES-SPECIFIC PARAMETERS
; WORK IN PIXEL SPACE FOR CONSTANT ARRAY SIZES
; FOR EPOCHS 3, 4, 5 USE PIXELS 350 TO 7150
;     ST_ORD 42 [5154, 5218]  ORD 75 [7156, 7245]  ; epoch eq 5
;     ST_ORD 35 [4865, 4926]  ORD 78 [7418, 7510]  ; epoch ge 3
;     ST_ORD 29 [4643, 4700]  ORD 78 [7418, 7510]  ; epoch le 4
  
  if n_params () lt 2 then begin
     print, ' IDL> cbc_star_templ, obj_nm, obsnm, tag=tag'
     return
  endif

  templ_nm = obj_nm+'_templ_'+tag+'.dat' 
  cbcenv = cbc_init(obsnm, templ_nm, obj_nm=obj_nm, osamp=osamp, ddir=ddir, excalibur=excalibur)
  c_light = cbcenv.c_light
  date = strmid(obsnm, 0, 6)       ; e.g. 190531
  edge_pad = cbcenv.edge_pad       ; toss 100 pixels from each side of order  
  epoch = cbcenv.epoch             ; 1 - 5
  nchunk_ord = cbcenv.nch_ord      ; number of chunks per order
  nchunks = cbcenv.n_chunks        ; 40 chunks / order * nord (+ 28? DF dec25, 2020)
  ord_init = cbcenv.st_ord         ; initial order for CCF (LFC)
  ord_final = cbcenv.fin_ord       ; final order for CCF
  nord = ord_final - ord_init + 1  ; number of orders for CCF
  npix_chunk = cbcenv.n_pix        ; 140*osamp pixels per chunk 
  px0_init = cbcenv.st_pix         ; trim NaNs from spectrum
  templ_dir = cbcenv.templ_dir     ; directory for template obs
  templ_pad = cbcenv.templ_pad     ; 40*osamp padding on each side of chunk

  tpix_chunk = npix_chunk + (2.*templ_pad) 

; TEMPL IS IDL STRUCTURE FOR EACH CHUNK 
; SAME STRUCTURE IF MORPHED
  if ~keyword_set(coadd_obnm) then coadd_obnm=' '
  templ = {obnm:'', coadd:'n',             $
           coadd_spec:dblarr(tpix_chunk),  $
           coadd_obnm:coadd_obnm,          $
           contin:dblarr(tpix_chunk),      $
           morph_spec:dblarr(tpix_chunk),  $ 
           morph_resid:dblarr(tpix_chunk), $ 
           morph_rv:0.0,                   $
           morph_vsini:0.0,                $
           morph_lor:0.0,                  $
           morph_ngau:0,                   $
           morph_tau:0.0,                  $
           ord:ord_init,                   $
           pixt:cbcenv.st_pix,             $
           pixob:cbcenv.st_pix,            $
           redchi:0.0d,                    $
           re_blaze:dblarr(tpix_chunk),    $
           snr:dblarr(tpix_chunk),         $
           spec:dblarr(tpix_chunk),        $
           stel:dblarr(tpix_chunk),        $
           tell:dblarr(tpix_chunk),        $
           unc:dblarr(tpix_chunk),         $
           wav:dblarr(tpix_chunk),         $
           weight:1.0,                     $
           wv_lab:dblarr(tpix_chunk)}

  templ = replicate(templ, nchunks)

;; READ TEMPLATE OBSERVATION AND FILL CHUNKS
  fname = cbcenv.obs_dir+obj_nm+'_'+obsnm+'.fits'       ; file name of observed spectrum
  spec = mrdfits(fname, cbcenv.file_ext)                ; observed spectrum
  if keyword_set(coadd_templ_nm) then $
     restore,templ_dir+coadd_templ_nm                   ; star[pix, ord]

; TEL_FLUX: MOLECFIT TELLURIC SPECTRUM (SC)
  for i = 0, nord-1 do begin
     ord = ord_init + i 
     sp = spec[ord].spectrum/spec[ord].continuum ; normalized spectrum
     if keyword_set(coadd_templ_nm) then $
        cosp = reform(star[*,ord])
     ws = spec[ord].bary_wavelength                     ;;     barycentric wavelength
     w_lab = spec[ord].wavelength                       ; wavelength in rest frame of lab
     if keyword_set(excalibur) then begin
        ws = spec[ord].bary_excalibur                   ; barycentric wavelength
        w_lab = spec[ord].excalibur                     ; wavelength in rest frame of lab
     endif

     tell = spec[ord].tellurics
     if mean(tell eq 0.0) then stop, 'No telluric model'
     uncert = spec[ord].uncertainty 
     continuum = spec[ord].continuum 
     snr = 1./spec[ord].uncertainty     ; RP noise-weighted extraction SNR (norm spec)

;     if ord eq 41 and keyword_set(excalibur) then begin
;        nchunk_ord = 37
;        px0_init = 1470 
;     endif 
;     if ord gt 41 and keyword_set(excalibur) then begin
;        nchunk_ord = cbcenv.nch_ord
;        px0_init = cbcenv.st_pix
;     endif

     for j = 0, nchunk_ord-1 do begin                   ; fill each chunk
        chnk = j + (i*nchunk_ord)                       ; chunk index
;        if ord_init eq 41 and keyword_set(excalibur) and i ge 1 then chnk = j + ((i-1)*nchunk_ord) + 37
;print, i, j, ord, chnk
        pix0 = px0_init + (j * npix_chunk)              ; pix0 of chunk
        pix1 = pix0 + npix_chunk - 1                    ; last pix of chunk
        templ[chnk].pixt = pix0                         ; pix0 of unpadded chunk
        templ[chnk].pixob = pix0                        ; pix0 of unpadded chunk
        ppix0 = pix0 - templ_pad                        ; pix0 - padding
        ppix1 = pix1 + templ_pad                        ; last pix plus padding
        templ[chnk].obnm = obj_nm+'_'+obsnm             ; template obnm
        templ[chnk].tell = tell[ppix0:ppix1]            ; wavelengths for telluric spectrum
        templ[chnk].wav = ws[ppix0:ppix1]               ; wavlength array padded chunk
        templ[chnk].spec = sp[ppix0:ppix1]              ; spectrum array padded chunk
        xnan = where( finite(templ[chnk].spec, /NAN), nxnan )   ; find any lingering NaNs
        if nxnan gt 0 then begin
           print, 'number of NaNs: ',nxnan              ;  ...report them,
           templ[chnk].spec[xnan] = 0.01                ;  ...and set them to 0.01
        endif
        templ[chnk].wv_lab = w_lab[ppix0:ppix1]         ; wavlength array padded chunk
        if keyword_set(coadd_templ_nm) then begin
           templ[chnk].coadd_spec = cosp[ppix0:ppix1]   ; coadded template spectrum
           templ[chnk].coadd = 'y'
        endif
        templ[chnk].ord = ord                           ; relative order 
        templ[chnk].snr = snr[ppix0:ppix1]              ; sqrt counts - SNR
        templ[chnk].unc = uncert[ppix0:ppix1]
        templ[chnk].contin = continuum[ppix0:ppix1]

;plot, templ[chnk].wav, templ[chnk].spec, /xsty 
;oplot, templ[chnk].wav, templ[chnk].tell, col=222
;oplot, templ[chnk].wav, templ[chnk].spec / templ[chnk].tell, col=90
;wait, 2 

; DERIVE TEMPLATE CHUNK WEIGHT FROM SLOPE
        sp_ck = templ[chnk].spec                                ; normalized spectrum
        xneg = where(sp_ck lt 0.0, nxneg)                       ; enforce positivity
        if nxneg gt 0 then sp_ck[xneg]=0.01                     ; positivity
        lam = templ[chnk].wav                                   ; wavelength array of chunk 
        disp_arr = lam[1:npix_chunk-1] - lam[0:npix_chunk-2]    ; dispersion
        di_dpix = sp_ck[1:npix_chunk-1] - sp_ck[0:npix_chunk-2] ; dI/dpix
        di_dv = di_dpix * lam/(c_light*disp_arr)                ; dI/dvel
        templ[chnk].weight = total((di_dv / sp_ck)^2)           ; weight for chunk
        if stddev(templ[chnk].tell) gt 0.02 then templ[chnk].weight=0  ; too many tellurics in chunk
        if keyword_set(demo) then begin
           str_ckwt = 'Order: '+strcompress(string(templ[chnk].ord),/remove_all)+$
                      ' Pixel: '+strcompress(string(templ[chnk].pixt),/remove_all)
           plot, templ[chnk].wav, templ[chnk].spec, /xsty, /ysty, yra=[0, 1.05], $
                 xtitl='!6 Wavelength', titl=str_ckwt
           if keyword_set(coadd_templ_nm) then oplot, templ[chnk].wav, templ[chnk].coadd_spec, col=222
           oplot, templ[chnk].wav, templ[chnk].tell, col=80
           wait, 0.5
        endif
     endfor                     ; step through chunks in order
  endfor                        ; step through orders

  templ.weight = templ.weight/total(templ.weight) ; sum of weights = 1
  save, templ, f=cbcenv.templ_dir+templ_nm        ; save templ structure in templ directory

end ; pro
