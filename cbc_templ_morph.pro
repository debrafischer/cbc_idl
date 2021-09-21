;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRO CBC_TEMPL_MORPH
;                 
; fischer 20 Jul 2019 
; morph the NSO spectrum to templ observation
; 
; PROCEDURES CALLED: 
;       CBC_INIT
;       RDNSO
;       MRDFITS
;       CONTF
;       MPFITFUN
;       TEMPL_FIT
;     
; PURPOSE:
;   Generate a NSO-derived template for CBC
;    
; INPUT:
;  CBCENV: an IDL structure, like a common block passed b/t programs
;  OBJ_NM: starname (required)
;  OBSNM: placeholder - one obsnm for cbc_init
;  TAG_IN: tag of cbc_spec_templ structure, e.g., 'ame'
;  TAG_OUT: add on ngau, e.g., 'ame60' 
;  OSAMP: pixel sampling of the template
;  COADD: if this set, then templ[ichnk].spec = coadd_spec from
;                                                 cbc_star_templ 
;  DEMO: optional plots
;      
; OUTPUT:
;  TEMPL: IDL structure holding chunks of the template
;   
; PROCEDURE:
;   1. before running this program, run cbc_star_templ.pro to generate 
;      the template structure 
;   2. load in the NSO spectrum
;   3. call cbc_init to setup cbcenv structure with global parameters
;      all the EXPRES parameters are hard-wired there for consistency
;      in all the cbc programs
;   4. find the templ chunk in the observation 
;      (will be shifted from pixt to pixob) 
;   5. carry out L-M fitting for Doppler shift, normalization, 
;      Lorentzian scaling, gaussian heights to reshape the NSO 
;      spectrum to match the observed spectrum in templ_nm_in
;                                                                                          
; CALLING SEQUENCE
;   cbc_templ_morph, '101501', '190531.1107',tag_in='a', tag_out='am'
;   osamp=1, ngau=10
;                                                                                          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
pro cbc_templ_morph, obj_nm, obsnm, tag_in=tag_in, tag_out=tag_out, ddir=ddir, $
                     coadd=coadd, osamp=osamp, ngau=ngau, demo=demo, $
                     thresh=thresh, div_telluric=div_telluric, templ_rv=templ_rv, $
                     excalibur=excalibur, fix_rotbro=fix_rotbro, $
                     init_rotbro=init_rotbro, fix_rv=fix_rv, hdcopy=hdcopy

  if ~keyword_set(osamp) then osamp=1

  strngau = strcompress(string(ngau),/remove_all) 
  templ_nm_in = obj_nm+'_templ_'+tag_in+'.dat' 
  templ_nm_out = obj_nm+'_templ_'+tag_out+'.dat' 
  cbcenv = cbc_init(obsnm, templ_nm_in, obj_nm=obj_nm, ddir=ddir, excalibur=excalibur)
  apod = cbcenv.apod               ; ignore these edge pixels when morphing templ
  c_light = cbcenv.c_light         ; speed of you-know-who
  templ_pad = cbcenv.templ_pad     ; 150*osamp padding on each side of chunk
  npix_chunk = cbcenv.n_pix        ; 140*osamp pixels per chunk 
;  edge_pad = cbcenv.edge_pad       ; toss 100 pixels from each side of order  
  gauwid = cbcenv.gauwid           ; 
  ord_init = cbcenv.st_ord         ; initial order for CBC (LFC)
  ord_final = cbcenv.fin_ord       ; final order for CBC 
  nord = ord_final - ord_init + 1  ; number of orders for CBC 
  nchunk_ord = cbcenv.nch_ord      ; number of chunks per order
  nchunks = cbcenv.n_chunks        ; 40 chunks / order * nord  
  tpix_chunk = npix_chunk + (2.*templ_pad)  ; number of pixels in padded templ with osamp=1

;; READ TEMPLATE OBSERVATION
  restore, cbcenv.templ_dir+templ_nm_in    ; templ structure
  if keyword_set(coadd) then templ.spec = templ.coadd_spec

; SETUP THE PARINFO STRUCTURE FOR FREE PARAMETERS
  names = ['RV', 'tau', 'vsini', 'gamma']
  nnames = n_elements(names)
  npar = n_elements(names)+ngau 
  parinfo = {value: 0.0d,       $ ; double precision
             fixed: 0,          $
             limited: [0,0],    $ ; use with caution
             limits: fltarr(2), $
             parname: '?',      $
             step: 0.01d,       $
             relstep: 0.00,     $
             mpside: 2}         ; 0, 1, -1, 2
  parinfo=replicate(parinfo, npar)

  for i=0, n_elements(names)-1 do parinfo[i].parname = names[i]

  ; PARINFO[0] = RV
  parinfo[0].value = 0.d       ; doppler shift m/s
  parinfo[0].step = 2.0d
;parinfo[0].limited=[1,1]        ; for 17156 191101
;parinfo[0].limits=[-24250,-5000] ; for 17156 191101
  if keyword_set(fix_rv) then begin
     parinfo[0].fixed = 1
     parinfo[0].value = fix_rv
  endif

  ; PARINFO[1] = TAU
  parinfo[1].value = 1.0       ; tau factor (1.0 = no change)
  parinfo[1].step = 0.005

  ; PARINFO[2] = VSINI
  if keyword_set(init_rotbro) then parinfo[2].value = init_rotbro else parinfo[2].value = 1.0
  parinfo[2].step = 0.5 ;0.01       ; vsini
  parinfo[2].limited = [1,1]   ; no negative vsini
  parinfo[2].limits = [0.1, 10.0] ; min and max vsini
  if keyword_set(fix_rotbro) then begin
     parinfo[2].fixed = 1
     parinfo[2].value = fix_rotbro
  endif

  ; PARINFO[3] = GAMMA
  parinfo[3].value = 0.1       ; gamma for Lorentzian
  parinfo[3].fixed = 1

  ; PARINFO[4:*] = GAU HEIGHTS
  for jj = 0, (ngau)-1 do begin   ; for first pass, zero out heights of gaussians
     nm = 'gau'+strcompress(string(jj),/remove_all)
     parinfo[jj+nnames].parname = nm
     parinfo[jj+nnames].value = 0.0
     parinfo[jj+nnames].fixed = 1.0        ; initially, gaussians not a free parameter
  endfor

  if ~keyword_set(templ_rv) then begin
  ; MANUALLY GET THE APPROX DOPPLER SHIFT OF STAR W.R.T. NSO
  ; INITIAL GRAB OF NSO SPECTRUM NEAR NA-D LINES

     sp = mrdfits(cbcenv.obs_dir+cbcenv.obj_nm+'_'+cbcenv.obsnm+'.fits',1)
     w0 = 5892.0
     w1 = 5903.0

; EXPRES SPECTRA ARE VACUUM, BUT THE NSO IS AIR 
     vactoair, w0, w0air        ; convert EXPRES wavelengths to air for NSO query
     vactoair, w1, w1air        ; ditto
     rdnso, wavsun, ssun, w0air-2.3, w1air+2.3   ; padded NSO chunk
     airtovac, wavsun                            ; NSO vacuum wavelengths         
     if keyword_set(excalibur) then begin
        wsnip0 = sp[56].bary_excalibur
        ssnip0 = sp[56].spectrum/sp[56].continuum
        kp_pix = where(wsnip0 gt w0 and wsnip0 lt w1) 
        wsnip = wsnip0[kp_pix]
        ssnip = ssnip0[kp_pix] 
     endif
     if ~keyword_set(excalibur) then begin
        wsnip0 = sp[56].bary_wavelength
        ssnip0 = sp[56].spectrum/sp[56].continuum
        kp_pix = where(wsnip0 gt w0 and wsnip0 lt w1)
        wsnip = wsnip0[kp_pix]
        ssnip = ssnip0[kp_pix] 
     endif
     plot, wsnip, ssnip, /xsty     ; vacuum wavelengths 
     oplot, wavsun, ssun, col=222  ;vacuum wavelength NSO
     print, 'Click on any line in the observation (black)' 
     cursor, wobs, y
     print, 'Click on the same line in the NSO (red)' 
     wait, 1
     cursor, wsun, y 
  ; INITAL GUESS FOR DOPP SHIFT TO APPLY TO NSO 
     parinfo[0].value = c_light * (wobs - wsun) / wsun      
     print, parinfo[0].value 
  endif
  if keyword_set(templ_rv) then begin
     parinfo[0].value=templ_rv
  endif

  ; LIMIT RV SHIFT TO W/I 1500 M/S OF INITIAL VALUE
  parinfo[0].limited=[1,1]
  parinfo[0].limits=[parinfo[0].value-3000., parinfo[0].value+3000.]


; REFINED GRAB OF NSO SPECTRUM AND 
; MORPH the NSO WITH L-M LEAST SQUARES FIT
  for ichnk = 0, nchunks-1 do begin
     w0 = templ[ichnk].wav[0]      &   vactoair, w0, w0air
     w1 = max(templ[ichnk].wav)    &   vactoair, w1, w1air
     rdnso, wavsun, ssun, w0air-2.3, w1air+2.3             ; fat padded NSO chunk, shifted wavelengths 
     airtovac, wavsun                                      ; and back to vacuum wavelengths

   ; ALIGN CONTINUUM OF TEMPL OBS AND NSO
     contf, templ[ichnk].spec, ct, sbin=10, nord=1
     templ[ichnk].spec = templ[ichnk].spec / ct
     contf, ssun, cs, sbin=10, nord=1 
     ssun = ssun / cs

   ; DIVIDE OUT THE TELLURIC SPECTRUM
     if keyword_set(div_telluric) then begin
        if min(templ[ichnk].tell) ge 0.5 then $
           templ[ichnk].spec = templ[ichnk].spec / templ[ichnk].tell
     endif

   ; SETUP THE EXTRA ARRAY AND INPUT FOR MPFITFUN 
   ; WITH NSO ON OSAMP=1 SCALE
     wtmpl = templ[ichnk].wav             
     stmpl = templ[ichnk].spec
     disp_arr = wtmpl[1:tpix_chunk-1] - wtmpl[0:tpix_chunk-2]    ; templ obs dispersion
     dispfine = mean(disp_arr)/osamp                             ; mean dispersion
     npixf = (max(wavsun) - wavsun[0])/dispfine                  ; divide NSO into dispfine pixels
     wavfine = wavsun[0] + dindgen(npixf-1)*dispfine             ; osamp pix size of obs imposed on NSO
     sunfine = spline(wavsun, ssun, wavfine, /double)            ; fine pix NSO spectrum
;stop
; check it out
;plot, (wavfine + parinfo[0].value * wsun / c_light), sunfine, /xsty
;oplot, wtmpl, stmpl, col=222
;wait,1
;stop

     ; PIXEL WEIGHTS
     noise_flat=0.0055
     tt = (templ[ichnk].unc / templ[ichnk].contin)^2            ; normalized 
     wt = (1. / tt) / (1. + tt * noise_flat^2)                  ; add noise from flat
     xneg=where(wt lt 0.,nxneg) & if nxneg gt 0 then wt[xneg]=0.d
     xhi=where(wt gt 5.*median(wt),nxhi) & if nxhi gt 0 then wt[xhi]=0.d

     ; TELLURIC PIXEL WEIGHT - TEST Jun19, 2021 DAF
;test     if min(templ[ichnk].tell) lt 0.5 then templ[ichnk].weight = 0.0
;     xz = where(templ[ichnk].tell lt 0.0, nxz)
;     if nxz gt 0 then wt[xz] = 0.0 
;     xp = where(templ[ichnk].tell gt 0.0, nxp) 
;     if nxp gt 0 then wt = wt * templ[ichnk].tell^2
;     wt = wt / total(wt) 

;     ; TEST BLOCK FOR DOWN-WEIGHTING PIXELS AFFECTED BY TELLURICS
;     ; 30 Jul 2020 - DAF
;     down_weight = templ[ichnk].tell 
;     xtst = where(down_weight lt 0.9, nxtst)    ; CLeet used -1.0 for saturated tellurics
;     if nxtst gt 0 then begin
;        down_weight[xtst]=0.0
;        wt[xtst] = 0.0
;     endif

     parinfo1 = parinfo[0:nnames-1]              ; parinfo[0:3] - really just Dopp shift
     functargs = {stempl: stmpl,         $
                  wfinenso: wavfine,     $
                  sfinenso: sunfine,     $
                  cbcenv: cbcenv,        $
                  nnames: nnames,        $
                  osamp: osamp,          $
                  c_light:c_light}
     x = wtmpl 
     y = stmpl 

   ; L-M FIT FOR DOPPLER SHIFT ONLY
;     if ~keyword_set(fix_rv) then begin
        newpar = mpfitfun('templ_fit', x, y, parinfo=parinfo1,      $ ; just the Dopp shift w/parinfo1
                          functargs=functargs, maxiter=200, /nan,   $
                          errmsg=errmsg, /iterstop,                 $
                          bestnorm=bestnorm, perror=perror,         $
                          yfit=syn_fit, weight=wt, /quiet)
     
        nfree = 1
        dof = n_elements(syn_fit) - nfree 
        redchisq = total( wt * (syn_fit - templ[ichnk].spec)^2)
        redchisq = redchisq / dof 
        diff = stmpl - syn_fit

; commented out Nov 25, 2020 Fischer 
; no syn_fit
;        if keyword_set(hdcopy) then ps_open,'morph_fit',/color
;        if keyword_set(demo) then begin
;           !p.charsize = 1.8
;           !x.omargin = [8,2]
;           !y.omargin = [2,2]
;           !x.margin = [3, 1]
;           !y.margin = [2, 1]
;           plot, wtmpl, stmpl, /xsty, /ysty, yr=[-0.2,1.2],xtitl='!6 Wavelength', thick=2
;           oplot, wtmpl, syn_fit, col=222, linesty=2, thick=2
;           oplot, wtmpl, diff
;        endif
;     print, newpar
;     endif


; FULL BLOWN FIT, INCLUDING GAUSSIANS
   ; NOW CALCULATE DIFFERENCE B/T MODEL AND TEMPL OBS
   ; LOCATE NGAU GAUSSIANS AT THE BUMPS IN THE DIFF SPECTRUM 
     if ~keyword_set(thresh) then thresh = 0.01                ; threshold for significant residuals to be fit
     node_pix = fltarr(ngau)      ; pixel locations of ngau Gaussians
     gau_ht = fltarr(ngau)        ; pixel heights of ngau Gaussians
     wid = gauwid                 ; fixed to the LSF width
     parinfo.fixed = 0            ; gaussians no longer fixed
     if keyword_set(fix_rv) then begin    ; par[0] is RV
        parinfo[0].fixed = 1
        parinfo[0].value = fix_rv 
     endif
     if ~keyword_set(fix_rv) then parinfo[0].value = newpar[0]

     parinfo[1].value = 1.0                ; par[1] is tau   

     if keyword_set(fix_rotbro) then begin ; par[2] is broadening
        parinfo[2].fixed = 1
        parinfo[2].value = fix_rotbro 
     endif
     if ~keyword_set(fix_rotbro) then begin
        if keyword_set(init_rotbro) then parinfo[2].value = init_rotbro else parinfo[2].value = 0.5
        parinfo[2].limited = [1,1] ; vsini must be positive and < 10.0 km/s
        parinfo[2].limits = [0.0, 5.0]
     endif

     parinfo[3].fixed = 1                   ; par[3] is gamma for lorentzian
     parinfo[3].value = 0.1                 ; gamma for Lorentzian

     j = -1                       ; counter
     resid_ht = 2.0*thresh        ; start big for while loop below
     len = n_elements(stmpl)      ; template chunk size
     xarr = indgen(len)           ; index array
     orig_diff = diff             ; save it 
     diff = smooth(diff,3)        ; smooth by 3 pixels to reduce fitting noise

   ; FIND THE NODE LOCATIONS AND INITIAL HEIGHTS FOR GAUSSIANS
   ; WHILE LOOP FINDS THE LARGEST RESIDUALS (GREATER THAN THRESH) 
   ; AND DROPS NODES THERE UNTIL NGAU NODES HAVE BEEN USED UP. 
   ; DOF = 480 (pixels) - NGAU 
     while abs(resid_ht) gt thresh and j lt ngau-1 do begin
        j = j+1                   ; counting
        diff[0:apod] = 0          ; zero diff inner edge pixels
        diff[len-(apod-1):len-1] = 0     ; zero diff for outer edge pixels
        resid_ht = max(abs(diff), n_ind) ; find greatest differences in loop
        resid_ht = diff[n_ind]    ; value of diff at peak or dip
        node_pix[j] = n_ind       ; pixel index in diff for Gaussian node
        gau_ht[j] = resid_ht      ; pos or neg

      ; SUBARRAY IN XARR FOR BUILDING GAUSSIAN
        ind = where(xarr ge n_ind - 15.0 * wid and xarr le n_ind + 15.0 * wid)
        gau = resid_ht * exp(-0.5 * ((xarr[ind] - n_ind)/wid)^2)
;        oplot, ind, gau, ps=-8, co=100 
        diff[ind] = diff[ind] - gau
     end ; while 

     ii = where(node_pix ne 0, nxgau) ; found the significant bumps in diff

   ; WHERE RESID_HT GT THRESHOLD, 
   ; SET UP NODES AND GAUSSIAN HEIGHT PARS FOR L-M FIT 
;if nxgau lt 1 then stop
     if nxgau ge 1 then begin
        node_pix = node_pix[ii] 
        gauwav = wtmpl[node_pix]
      ; SHIFT TO NSO REST WAVELENGTHS
        gauwav = gauwav * (1.d - parinfo[0].value/c_light)  ; w' = w * (1 + v/c)
        gau_ht = gau_ht[ii]
        noise_flat=0.002
        tt = (templ[ichnk].unc / templ[ichnk].contin)^2
        wt = (1. / tt) / (1. + tt * noise_flat^2)
        xneg=where(wt lt 0.,nxneg) & if nxneg gt 0 then wt[xneg]=0.d
        xhi=where(wt gt 5.*median(wt),nxhi) & if nxhi gt 0 then wt[xhi]=0.d

     ; TELLURIC PIXEL WEIGHT - TEST Juneteenth, 2021 DAF
;test     if min(templ[ichnk].tell) lt 0.5 then templ[ichnk].weight = 0.0
;     xz = where(templ[ichnk].tell lt 0.0, nxz)
;     if nxz gt 0 then wt[xz] = 0.0 
;     xp = where(templ[ichnk].tell gt 0.0, nxp) 
;     if nxp gt 0 then wt = wt * templ[ichnk].tell^2
;     wt = wt / total(wt) 

;     ; TEST BLOCK FOR DOWN-WEIGHTING PIXELS AFFECTED BY TELLURICS
;     ; 30 Jul 2020 - DAF
;     down_weight = templ[ichnk].tell 
;     xtst = where(down_weight lt 0.9, nxtst)    ; CLeet used -1.0 for saturated tellurics
;     if nxtst gt 0 then begin
;        down_weight[xtst]=0.0
;        wt[xtst] = 0.0
;     endif

        parinfo[nnames:nnames+nxgau-1].value = gau_ht

        functargs = {stempl: stmpl,         $  ; observed spectrum
                     wfinenso: wavfine,     $  ; NSO wavelengths
                     sfinenso: sunfine,     $  ; NSO spectrum
                     cbcenv: cbcenv,        $  ; parameter structure
                     osamp: osamp,          $  ; oversample?
                     nnames: nnames,        $  ; number of non-gaussian free pars
                     gauwid: wid,           $  ; width of gaussians
                     node_pix: node_pix,    $  ; pixels of gaussian nodes
                     gauwav: gauwav,        $  ; wavelengths of gaussian nodes
                     c_light:c_light}
        x = wtmpl 
        y = stmpl 
;help,/st, parinfo
;stop
        newpar = mpfitfun('templ_fit', x, y, parinfo=parinfo,       $ ; pars include gaussians
                          functargs=functargs, maxiter=200, /nan,   $
                          errmsg=errmsg, /iterstop,weight=wt,       $
                          bestnorm=bestnorm, perror=perror,         $
                          yfit=syn_fit, /quiet)   

        templ[ichnk].morph_spec = syn_fit
        templ[ichnk].morph_resid = syn_fit - stmpl
        templ[ichnk].morph_rv = newpar[0]
        templ[ichnk].morph_tau = newpar[1]
        templ[ichnk].morph_vsini = newpar[2]
        templ[ichnk].morph_ngau = nxgau

        nfree = nnames + n_elements(node_pix)
        dof = n_elements(syn_fit) - nfree - 4        ; 2 pixels fixed at each edge of chunk
        redchisq = total( wt * (syn_fit - templ[ichnk].spec)^2)
        redchisq = redchisq / dof 
        templ[ichnk].redchi=redchisq 

        print, 'chunk: ', ichnk, ' Ngau: ', nxgau, ' Mean(resid): ', mean(templ[ichnk].morph_resid)

   ; PLOTS, if keyword_set(demo) 
        if keyword_set(demo) then begin
           !p.charsize = 1.6
           !x.omargin = [8,2]
           !y.omargin = [2,2]
           !x.margin = [3, 1]
           !y.margin = [2, 1]
           plot, wtmpl, stmpl, /xsty, /ysty, yr=[-0.1,1.1]
           oplot, wtmpl, syn_fit, col=222, thick=2
           diff = stmpl - syn_fit
           oplot, wtmpl, diff, col=90, thick=2
           ;wait, 0.5
        endif                   ; demo
     endif                      ; nxgau fits

     if keyword_set(hdcopy) then ps_close

  endfor                        ; step through chunks in order

fout = cbcenv.templ_dir+templ_nm_out 
;stop 

  save, templ, f=fout ; cbcenv.templ_dir+templ_nm_out        ; save templ structure in templ directory

end ; pro


