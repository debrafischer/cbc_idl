;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION TEMPL_FIT
;
; initial coding: Debra Fischer, 21 Jul 2019
;
; CALLED FROM: CBC_TEMPL_MORPH
;
; PURPOSE:  
;   Creates the synthetic NSO template with free parameters 
;   (doppler shift, tau, gamma)
;   driven by a Levenberg-Marquardt algorithm
;
; PROCEDURE:
;   1. unpack the variables and spectra
;   2. doppler shift the template wavelengths
;   3. oversample the template
;   4. morph Gaussians to difference spectrum to generate model
;   5. cubic spline interpolation of model onto wavelength scale of observation
;   6. tweak vertical scaling of spectrum 
; 
; INPUTS: 
;   X: independent variable (wavelength of template observation)
;   PAR: initial values for free parameters defined in parinfo
;
; -EXTRA:
;   WTEMPL:
;   STEMPL: padded observed template spectrum
;   WAVFINE: NSO wavelengths
;   SUNFINE: normalized NSO
;   CBCENV: parameter structure 
;   C_LIGHT: speed of you know who
; 
; OUTPUTS:
;   NEWPAR: best fit free parameters from mpfit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION TEMPL_FIT, X, PAR, _EXTRA=FUNCTARGS

COMMON MPFIT_ERROR, ERROR_CODE 

  ; UNPACK THE VARIABLES 
    wtempl = x 
    extra=functargs                  ; stmpl, wnso, snso
    stempl = extra.stempl            ; template spectrum
    wavfine = extra.wfinenso         ; template padded wavelength
    sunfine = extra.sfinenso         ; observed chunk spectrum
    cbcenv = extra.cbcenv            ; parameter structure
    c_light = extra.c_light          ; speed of you-know-who
    osamp = extra.osamp              ; currently hardwired in
    nnames = extra.nnames            ; number of non-Gaussian free pars
    if n_elements(par) gt nnames then begin
       gauwid = extra.gauwid         ; pixel width of Gaussians
       node_pix = extra.node_pix     ; pix location of Gaussian nodes
       gauwav = extra.gauwav         ; wavelengths of Gaussians nodes
    endif

   ; print, par  
   ; print, 'n_elements in par: ', n_elements(par)

  ; PAR[0]: SHIFT TEMPL WAVELENGTHS: lambda_new = (1+ v/c)*lambda
    dopshift = double(par[0] / c_light)  ; trial for RV
    wfine_shft = double(wavfine) + double(wavfine*dopshift)

  ; CHECK THAT SHIFTED TEMPLATE IS STILL WIDER THAN OBSERVATION
    error_code = 0
    if wfine_shft[0] gt wtempl[0] then error_code=1
    if max(wfine_shft) lt max(wtempl) then error_code=2 
    if error_code ne 0 then stop

  ; PAR[1]: ADJUST TAU (GLOBAL)
    tau = alog(1.0d/sunfine)               ; initial ln(spectrum)
    newtau = tau * par[1]                  ; adjust global optical depth
    spec1 = exp(-newtau)                   ; sunfine with new opt depth

  ; MORPHING NSO TO DIFFERENCES 
    if n_elements(par) gt 4 then begin     ; adding in Gaussians to morph

     ; ENHANCE SPEC1
     ; PAR[NNAMES:NNAMES+NGAU-1]: ADD IN GAUSSIANS
       ngau = n_elements(node_pix)         ; initial number of Gaussians
       npar = n_elements(par)              ; all the free pars
       xarr = findgen(n_elements(wavfine)) ; 
       wid = gauwid  ;0.8 * gauwid                  ; think thin
       for j = 0, ngau-1 do begin
          gau_wv = min(abs(wavfine - gauwav[j]), wv_indx) ; closest wavelength
          wv_indx = wv_indx[0] 
          ind = where(xarr gt wv_indx - 15. * wid and xarr le wv_indx + 15. * wid)
          ht = par[j+nnames]
          gau = ht * exp(-0.5 * ((xarr[ind] - wv_indx)/wid)^2)
          spec1[ind] = spec1[ind] + gau    ; update sunfine model with Gaussian perturbations
       endfor
    endif

  ; PAR[2]: ROTATIONAL BROADENING SPEC1 => SPEC2
    lcen = mean(wavfine) 
    w = wavfine - lcen
    s = spec1
    vsini = par[2] > 0.0 
    if vsini gt 0.1 then spec1 = rotbro(w, s, lcen, vsini)      ; rotationally broadened, bumped up sunfine

;  ; PAR[3]: CONVOLVE WITH LORENTZIAN SPEC2 => SPEC3 
;    xlor_arr = findgen(100)-50             ; pixel grid for Lorentzian
;    gam = par[3] 
;    if abs(gam) lt 0.001 then begin 
;       lor = abs(par[3])/(xlor_arr^2 + par[3]^2) 
;       sp = spec2 
;       num_conv, sp, lor, spec2
;    endif

  ; SYNTH_NEW IS THE MODEL, SPLINED ONTO THE WAVELENGTH 
  ; SCALE OF THE OBSERVATION 
    rebin, wfine_shft, spec1, wtempl, synth_new ; cubic spline interpolation

;stop

;plot, wfine_shft, sunfine, /xsty
;oplot, wtempl, stempl, col=222 

  ; MATCH FLUX FOR THE OBSERVED AND SYNTHETIC SPECTRA 
    scale = stempl / synth_new
    contf,scale,c_scale, sbin=20, nord=3, frac=0.5
    syn_fit = c_scale*synth_new
 
    return, syn_fit

 end

