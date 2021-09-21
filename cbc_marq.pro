;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION CBC_MARQ
;
; CALLED FROM: CBC_MAIN
;
; CALLS TO: CBC_FIT 
; 
; PURPOSE:  
;   Drives Levenberg-Marquardt search for best fit parameters:
;   Doppler shift and continuum offset between Template and Stellar Spectrum
;    
; PROCEDURE:
;   1) weight the pixels 
;   2) set up the parinfo structure for free parameters
;   3) call LM fitting
; 
; INPUTS: 
;      chunk structure
;      cbcenv structure 
; 
; OUTPUTS:
;   MODEL: chunk.free_pars (updated RV shift, continuum offset) 
;   
; OUTSTANDING ISSUES: 
;
; Written by: Debra Fischer, Yale, Jul 2019
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION CBC_MARQ, vd, cbcenv, vdtag=vdtag, demo=demo, $
                   verbose=verbose, pfilt=pfilt


; USE RP UNCERTAINTIES INSTEAD
; PIXEL WEIGHTS FOR CALCULATING CHISQ OF THE MODEL	
          noise_flat = 0.004              ; typical s/n in flatfield is ~250-500 (DF jun30 2019) 
          ss = (vd.unc/vd.contin)^2       ; normalized uncertainty per RP
          wt = (1. / ss) / (1. + ss * noise_flat^2)
          xneg=where(wt lt 0.,nxneg) & if nxneg gt 0 then wt[xneg]=0.d
          xhi=where(wt gt 5.*median(wt),nxhi) & if nxhi gt 0 then wt[xhi]=0.d 

;test     if min(vd.stel) lt 0.5 then vd.weight = 0.0
;     ; TELLURIC PIXEL WEIGHT - TEST Jun19, 2021 DAF
;     xz = where(vd.stel lt 0.0, nxz)
;     if nxz gt 0 then wt[xz] = 0.10 
;     xp = where(vd.stel gt 0.0, nxp) 
;     if nxp gt 0 then wt = wt * vd.stel
;     wt = wt / total(wt) 


;  TELLURIC WEIGHTING DONE IN TEMPLATE CONSTRUCTION - NO NEED TO REPEAT
;     ; TELLURIC PIXEL WEIGHT
;     if min(templ[ichnk].tell) lt 0.5 then templ[ichnk].weight = 0.0

;     ; TEST BLOCK FOR DOWN-WEIGHTING PIXELS AFFECTED BY TELLURICS
;     ; 30 Jul 2020 - DAF
;     down_weight = vd.stel   ; tellurics
;     xtst = where(down_weight lt 0.0, nxtst)    ; CLeet used -1.0 for saturated tellurics
;     if nxtst gt 0 then down_weight[xtst]=0.0
;     wt = wt*down_weight^2 


;         later, do this: 
;          wt = wt*pfilt*1.0d   ; zero out the bad / telluric pixels 

    ; BEGIN L-M FITTING
    ; INFORMATION AND LIMITS FOR FREE PARAMETERS: PARINFO STRUCTURE

         if keyword_set(verbose) then begin
            print,'Parameters before Pre-fit: '
            print, vd.vel 
            print, vd.cont_offset
         endif
         
       ; initialize mod_chunk
        vd.red_chi = 100.   ; start bad, make fit improve this
    	parinfo = cbc_parinfo(vd, cbcenv)
        mod_chunk = vd
                      
; SETUP THE EXTRA ARRAY AND INPUT FOR MPFITFUN
        wtmpl = vd.wtmp
        if strmid(vdtag,1,1) eq 'm' then stmpl = vd.stmp_morph else stmpl = vd.stmp
        wobs = vd.wav
        sobs = vd.sobs

        functargs = {wtempl:wtmpl, $
                     stempl:stmpl, sobs:sobs, $
                     c_light:cbcenv.c_light}
        x = wobs
        y = sobs
        err = 1./wt

        if keyword_set(demo) then begin 
           plot, wtmpl, stmpl, /xsty
           oplot, wobs, sobs, col=222  
;           oploterr, wobs, sobs, err
        endif 

        newpar = mpfitfun('cbc_fit', x, y, parinfo=parinfo, $
                          functargs=functargs, maxiter=200, /nan, $
                          errmsg=errmsg, /iterstop, weight=wt, $
                          bestnorm=bestnorm, nfree=nfree, perror=perror, $
                          yfit=syn_fit, /quiet)

;stop
        dof = n_elements(sobs) - nfree
        chi1 = total(wt * (syn_fit - sobs)^2) 
        chi1 = chi1 / dof
        mod_chunk.red_chi = chi1
        mod_chunk.vel = newpar[0]
        mod_chunk.cont_offset = newpar[1]
        parinfo.value = newpar
        mod_chunk.smod = syn_fit 
        mod_chunk.rms = stddev(syn_fit-sobs)
        mod_chunk.perror = perror 

        if keyword_set(demo) then begin
           !p.charsize = 1.6
           !x.charsize = 1.8
           !y.charsize = 1.8
           !x.omargin = [8, 2]
           !y.omargin = [3, 1]
           !x.margin = [3, 1]
           !y.margin = [2, 1]
           ;plot, wtmpl, stmpl, /xsty, /ysty
           plot, wobs, sobs, /xsty, /ysty, yra=[0,1.01] ;col=222
           oplot, wobs, syn_fit, col=90, linesty=2
           oplot, wobs, sobs-syn_fit+0.2
           wait,0.5
;print, mod_chunk.perror
;if mod_chunk.perror[0] eq 0.0 then stop
        endif                   ; demo   

        return, mod_chunk

end

