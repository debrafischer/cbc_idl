;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION CBC_FIT
;
; initial coding: Debra Fischer, Jun 2019
;
; CALLED FROM: CBC_MAIN
;
; PURPOSE:  
;   Creates the synthetic observation with free parameters 
;   (doppler shift and continuum offset)
;   driven by a Levenberg-Marquardt algorithm
;
; PROCEDURE:
;   1. unpack the variables and chunk
;   2. doppler shift the template wavelengths
;   3. oversample the template
;   4. cubic spline interpolation of template onto wavelength scale of observation
;   5. vertical scaling of spectrum 
; 
; INPUTS: 
;   X: independent variable (wavelength of observation)
;   PAR: initial values for free parameters 
;
; -EXTRA:
;   WTEMPL: wavelengths of template
;   STEMPL: padded template spectrum
;   SOBS: observed spectrum 
;   CBCENV: structure 
;   C_LIGHT: speed of light
; 
; OUTPUTS:
;   NEWPAR: best fit free parameters from mpfit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION CBC_FIT, X, PAR, _EXTRA=FUNCTARGS

COMMON MPFIT_ERROR, ERROR_CODE 

  ; UNPACK THE VARIABLES 
    wobs = x 
    extra=functargs          ; wtmpl, stmpl, sobs
    wtempl = extra.wtempl    ; template padded wavelength
    stempl = extra.stempl    ; template spectrum
    sobs = extra.sobs        ; observed chunk spectrum
    c_light = extra.c_light  ; speed of you-know-who
    osamp = 1

   ; print, par
   ; print, 'n_elements in par: ', n_elements(par)

  ; DOPPLER SHIFT THE OVERSAMPLED TEMPLATE 
    dopshift = double(par[0] / c_light)  ; trial for RV
    offset = par[1]                      ; trial for offset

    if osamp gt 1 then begin
      ; OVERSAMPLE THE TEMPLATE (MODEL) 
       w0 = wtempl[0]     
       w1 = max(wtempl)
       disp = double(w1 - w0)/(n_elements(wtempl))
       disp_fine = disp / osamp
       nfine = fix( (w1 - w0) / disp_fine)
       wtempl_fine = w0 + dindgen(nfine) * disp_fine

  ; CHECK THAT SHIFTED TEMPLATE IS STILL WIDER THAN OBSERVATION
    error_code = 0
    if wtempl_fine[0] gt wobs[0] then error_code=1
    if max(wtempl_fine) lt max(wobs) then error_code=2 
    if error_code ne 0 then stop

     ; DOUBLE PRECISION CUBIC SPLINE INTERPOLATION OF SPECTRUM ONTO FINE
     ; WAVELENGTH SCALE
       stempl_fine = spline(wtempl, stempl, wtempl_fine, /double) 

     ; SHIFT TEMPL WAVELENGTHS: lambda_new
                                ; = lambda*(1+ v/c)*sqrt(1.-(v/c)^2)

       wtempl_shft = double(wtempl_fine)*(1.0 + dopshift)*sqrt(1.0 - 0.5*dopshift^2)
;       wtempl_shft = double(wtempl_fine) +
;       double(wtempl_fine*dopshift) ;so non-relativistic

     ; SYNTH_NEW IS THE MODEL, SPLINED ONTO THE WAVELENGTH 
     ; SCALE OF THE OBSERVATION 
       rebin, wtempl_shft, stempl_fine, wobs, synth_new ; cubic spline interpolation
    endif

    if osamp eq 1 then begin
       wtempl_shft = double(wtempl) + double(wtempl*dopshift) 
       rebin, wtempl_shft, stempl, wobs, synth_new ; cubic spline interpolation
    endif

  ; MATCH CONTINUUM SHIFT FOR THE OBSERVED AND SYNTHETIC SPECTRA 
    syn_fit = par[1]*synth_new
 
    return, syn_fit

end
