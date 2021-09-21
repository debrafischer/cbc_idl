;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION CBC_INIT
;
; Initial coding: D. Fischer July 4, 2019
;
; CALLED FROM:
;   cbc_main, cbc_nso_templ, cbc_obs_templ 
;
; PURPOSE: 
;   One-stop editing for changing root paths, spectrographs
;   Establishes the Doppler code paths, global parameters, 
;       default EXPRES parameters.
;
; HARDWIRED DIRECTORY PATHS
;   ROOT_DIR: directory where code is located
;   OBS_DIR: directory where observations are located 
;   TEMPL_DIR: directory where templates located
;   OUT_DIR: directory for vd files 
;
; INPUTS:
;   OBSNM: observation ID (e.g., '190704.1000') 
;   TEMPL_NM: template file name 
;   OBJ_NM: name of target (e.g., 'NSO' or '10700') 
;   
; OPTIONAL INPUTS (DEFAULTS ARE HARDWIRED IN THIS CODE)

; USE PIXELS 100 TO 6500
;     ORD 42 [5154, 5218]  ORD 75 [7156, 7245]  ; epoch eq 5
;     ORD 35 [4865, 4926]  ORD 78 [7418, 7510]  ; epoch le 4
;
;   NPIX: number of pixels per chunk (not including padding) 
;   N_CHUNKS: number of chunks used for CBC in the spectrum
;   ST_PIX: initial pixel in order
;   ST_ORD: starting order in spectrum for CBC
;   N_ORD: number of orders used for CBC
;   NCH_ORD: number of chunks per order
;   OSAMP: oversampling of spectrum (template?)
;   PFILT_FILE: file with masks for bad pixels 
;
; OUTPUTS:
;   CBCENV: structure with relevant paths and parameters
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION CBC_INIT, obsnm, templ_nm, obj_nm = obj_nm, $
                   npix=npix,n_chunks = n_chunks,    $
                   st_pix = st_pix, st_ord = st_ord, $
                   n_ord = n_ord, nch_ord = nch_ord, excalibur=excalibur, $
                   osamp = osamp, ddir=ddir  ;, $
;                   pfilt = pfilt 

  root_path='/Users/debrafischer/research/cbc/'                     ; the code lives here
  if obj_nm eq 'sun' then begin
     obs_dir = '/Users/debrafischer/research/cbc/data/'             ; observations live here
     templ_dir = '/Users/debrafischer/research/cbc/template/'       ; template files here 
  endif else begin
     obs_dir = '/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'
     templ_dir = '/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'
  endelse
  out_dir = obs_dir

  
  date=strmid(obsnm,0,6)
  if date lt 180405 then epoch = 0
  if date ge 180426 and date le 180611 then epoch = 1
  if date ge 180621 and date le 180715 then epoch = 2
  if date ge 181005 and date le 190130 then epoch = 3
  if date ge 190201 and date le 190710 then epoch = 4
  if date gt 190710 then epoch = 5

if ~keyword_set(ddir) then ddir='excalibur'

; SETUP YOUR DIRECTORIES AND RENAME THE PATHS BELOW
;     LFC orders 
;     ORD 40 [5067, 5134]
;     ORD 41 [5109, 5177]  ORD 74 [7069, 7163]
;     ORD 42 [5154, 5218]  ORD 75 [7156, 7245]  ; epoch eq 5
;     ORD 35 [4865, 4926]  ORD 78 [7418, 7510]  ; epoch le 4
;     ORD 29 to ORD 78 (match JO CCF)           ; using ThAr wavelengths

; NEW FORMAT: 6850 pixels / order for rel orders 29 - 78 
     
; ALL EPOCHS FIT IN THIS TRIMMED PART OF THE SPECTRUM
  if ~keyword_set(osamp) then osamp = 1
  edge_pad = 250*osamp                               ; skip pixels at blaze edges to allow for BC shifts 

      ; DM-9 wavelengths (e.g., ame)
      st_ord = 37               
      fin_ord = 75
      px0_init = 350 
      npix_ord = 6800
      if keyword_set(excalibur) then begin
         st_ord = 42            ; for excalibur wavelengths
         fin_ord = 74
         px0_init = 630
         npix_ord = 6380 ; amt 
;         npix_ord = 6100
; test 25 Dec 2020 DF - add one and a half more orders - crashes (NaN wavelengths)
;         st_ord = 40            ; test this with 28 chunks in order 40 
;         st_ord = 41            ; for excalibur wavelengths
;         px0_init = 770
;         npix_ord = 6380
;         fin_ord = 74
      endif 


;  pfilt_file = 'pfilt.dat'                                          ; bad pixel file 
  st_pix = (px0_init+edge_pad)*osamp                 ; will be updated
  nord = fin_ord - st_ord + 1                        ; number of orders with LFC lines
  n_pix = 140*osamp                                  ; number pixels / chunk
  templ_pad = 150*osamp                              ; extra pixels to pad template chunks
  apod = 20                                          ; npix to ignore at edges of morphed chunks
  nch_ord = fix((npix_ord-2*edge_pad)*osamp / n_pix) ; 42 / 45 chunks / order (excalibur vs DM-9)
  n_chunks = nord * nch_ord                          ; total chunks / spec for epoch 4
;  if st_ord eq 40 and keyword_set(excalibur) then $
;     n_chunks = ((nord-1) * nch_ord) + 28   ; half of order 40, px0_init=2730
;  if st_ord eq 41 and keyword_set(excalibur) then $
;     n_chunks = ((nord-1) * nch_ord) + 37   ; half of order 41, px0_init=1470

;print, n_chunks 

  ; Declaration of the 'cbcenv' structure. 
   cbcenv = {cbcenv,                          $
             obsnm: obsnm,                    $ ; observation ID, input to dop_init
             templ_nm: templ_nm,              $ ; name of template structure file
             obj_nm: obj_nm,                  $ ; starname
             apod: apod,                      $ ; number of pix to ignore at edges of templ morph
             c_light: 2.99792458d8,           $ ; speed of light
             edge_pad: edge_pad,              $ ; npixels to skip at edges of order 
             epoch: epoch,                    $ ; instrument status epoch
             file_ext:1,                      $ ; file extension for mrdfits
             fin_ord: fin_ord,                $ ; final order in obs
             gauwid: 3.1,                     $ ; width for gaussians to morph templ
             n_chunks: n_chunks,              $ ; number of chunks in vd structure 
             nch_ord: nch_ord,                $ ; number of chunks per order
             n_ord: nord,                     $ ; number of orders
             n_pix: n_pix,                    $ ; number of pixels per chunk
             npix_ord: npix_ord,              $ ; num pixels per order
             osamp: osamp,                    $ ; oversampling for model
             obs_dir: obs_dir,                $ ; observation dir
             out_dir: out_dir,                $ ; output dir
             st_pix: st_pix,                  $ ; first pixel
             st_ord: st_ord,                  $ ; first chunk
             templ_dir: templ_dir,            $
             templ_pad: templ_pad,            $ ; extra template pixels on edges of chunks 
             pfilt_file: ''                   $ ; pixel filter file for observations
            }

   if n_elements(obsnm) gt 0 then cbcenv.obsnm =  obsnm
   if n_elements(templ_nm) gt 0 then cbcenv.templ_nm =  templ_nm
   if n_elements(obj_nm) gt 0 then cbcenv.obj_nm =  obj_nm
   if n_elements(npix) gt 0 then cbcenv.n_pix = n_pix
   if n_elements(n_chunks) gt 0 then cbcenv.n_chunks = n_chunks
   if n_elements(st_pix) gt 0 then cbcenv.st_pix = st_pix
   if n_elements(st_ord) gt 0 then cbcenv.st_ord = st_ord
   if n_elements(n_ord) gt 0 then cbcenv.n_ord = nord
   if n_elements(nch_ord) gt 0 then cbcenv.nch_ord = nch_ord
   if n_elements(osamp) gt 0 then cbcenv.osamp = osamp
;   if n_elements(file_ext) eq 0 then cbcenv.file_ext = file_ext
   if n_elements(obs_dir) gt 0 then cbcenv.obs_dir =  obs_dir
   if n_elements(templ_dir) gt 0 then cbcenv.templ_dir = templ_dir
   if n_elements(templ_pad) gt 0 then cbcenv.templ_pad = templ_pad
   if n_elements(out_dir) gt 0 then cbcenv.out_dir = out_dir
   if n_elements(pfilt_file) gt 0 then cbcenv.pfilt_file = pfilt_file  

   return, cbcenv


end
