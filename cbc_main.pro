;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
; PRO CBC_MAIN 
;   
; PROCEDURES CALLED: CBC_INIT, CBC_MARQ, CBC_FIT
;     
; PURPOSE:   
;   Engine for chunk RV analysis    
;     
; PROCEDURE:  
;   1) set up the pathnames and input variables: 
;      CBC_INIT sets up the program-specific info 
;   2) set up initial velocity data (VD) structure (will be CBC_CHUNK_SETUP) 
;   3) starting with the template, 
;      derive best Doppler fit to observations with Levenberg-Marquardt 
;        - restore the observation 
;        - restore the template
;          - step through each order and chunk
;              - set pixel weights and filter out bad pixels 
;      - for each chunk, derive shift, continuum offset
;
; CALLING SEQUENCE:  
;   cbc_main, /demo
;   
; INPUTS:  
;   STARNAME  
;   CBCENV structure:
;       OBSNM:  (obs name) stellar spectrum to be analyzed  
;       TEMPL_NM: name of template 
;           structure has keywords:   
;             pixt, pixob, ord, wav, spec, weight  
;   TAG: for filename of output vd structure    
;     
; OPTIONAL INPUT: 
;  
; OUTPUTS  
;   VD: velocity data structure
;   
; Written by Debra Fischer, Yale, 4 July 2019
;   
; OUTSTANDING: 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

pro cbc_main, obj_nm, obsnm, templ_nm, morph_spec=morph_spec, demo=demo, $
              vdtag=vdtag, verbose=verbose, ddir=ddir, excalibur=excalibur

; fischer 29 Jun 2019 

; some basics
  loadct, 39, /silent
  !p.background = 255
  !p.color = 1

  if keyword_set(demo) then demo=1 else demo=0 
  if keyword_set(verbose) then verbose=1 else verbose=0 
;  verbose=1 

; CBCENV STRUCTURE SETS UP PATHS ETC
  cbcenv = cbc_init(obsnm, templ_nm, osamp=osamp, obj_nm=obj_nm, ddir=ddir, excalibur=excalibur) 
  nchunk = cbcenv.n_chunks
  osamp = cbcenv.osamp

; RESTORE THE OBSERVED DATA
  vd_fname=cbcenv.obs_dir+'vd'+vdtag+'_'+cbcenv.obj_nm+'.'+obsnm+'.dat'
  restore, vd_fname  ; vd structure 

; FIT EACH CHUNK
  for ch = 0, nchunk - 1 do begin
     if vd[ch].templ_wt gt 0.0d then begin
        chunk = vd[ch]
      ; MODIFY THE CHUNK STRUCTURE BY RUNNING LM FITTING 
        mod_chunk = cbc_marq(chunk, cbcenv, vdtag=vdtag, demo=demo, verbose=verbose)
        vd[ch] = mod_chunk      ; update the chunk with newpars 
        if keyword_set(verbose) then $
           print, vd[ch].ord, vd[ch].pixt, vd[ch].pixob, vd[ch].vel, vd[ch].red_chi, vd[ch].rms
     endif
  endfor   ; ch = chunk

  save, vd, f=vd_fname 
;stop
end ; pro 
