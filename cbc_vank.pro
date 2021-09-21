pro cbc_vank, obj_nm, tag, mincts=mincts, ddir=ddir, pre=pre, skip_ord=skip_ord

; Fischer 6 July 2019 
; combine "obs" velocity data structures after chunk CBC 
; with cbc_main.pro

; vank, 'NSO', 'a' 

; INPUTS
;     obj_nm: 'NSO' '10700'  
;     tag: 'a' is the "template" tag to locate the vd files
;        a: osamp=1
;        b: osamp=2
;        c: osamp=4
;     mincts: minimum acceptable counts (sqrt(counts) = SNR) 
;     ddir: directory for fits spectra and vdfiles
;     pre: include epochs before 2019 (i.e., epochs 1, 2, 3)
;          only recommended for very large amplitude signals like 47 UMa
;
; OUTPUT
;   velocity structure: vst10700_tag.dat in /research/cbc/vel/
; 

if ~keyword_set(ddir) then ddir='excalibur'

if keyword_set(skip_ord) then skip_ord = [57, 66] ; badly contaminated with tellurics 

; GATHER THE INPUT VDS
  file_dir = '/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'  ; '/Users/debrafischer/research/cbc/files/' 
  if strlowcase(obj_nm) eq 'sun' then file_dir = ddir 
  froot = 'vd'+tag+'_'+obj_nm
  vdfiles = file_search(file_dir+froot+'*.dat', count=nobs)  ; CBC output of indiv obs
;print, vdfiles
;stop
; SPECIFY THE OUTPUT VST
  vel_dir = '/Users/debrafischer/research/cbc/vel/' 
  vstnm = vel_dir + 'vst'+strlowcase(obj_nm)+'_'+tag+'.dat'           ; output velocity structure 

; BASIC PARAMETERS 
   c_light = 2.99792458d8

; NUMBER OF CHUNKS
   restore, vdfiles[0]
   xx=where(vd.ord eq 56, nxx)             ; nxx = nchunks per order
   nchunk = (max(vd.ord) - min(vd.ord) +1) * nxx  ; total number of chunks

   ; EPOCHS
   ; 1: 180426 - 180611
   ; 2: 180621 - 180710
   ; 3: 181005 - 190130
   ; 4: 190201 - 190702
   ; 5: 190801 - 

; TOSS THE BAD BOYZ 
  if ~keyword_set(mincts) then begin
     mincts = 10000.            ; low SNR

     case obj_nm of 
        '102121': mincts=1000
;     '101501': mincts=5000
;     '85380': mincts=1000
;     '82885': mincts=1000
;     '91204': mincts=1000
;     'BD+62_1237': mincts=1000
        else: mincts=10000
     endcase 
  endif


  if obj_nm eq '17156' then mincts = 1000 
  maxcts = 122000. ; too high SNR - saturated? 
  maxchi = 5.0
  minchi = 0.2 

  if obj_nm eq 'Sun' then begin
     mincts = 40000
     maxcts = 200000
  endif

;print, nobs
;stop 

; PRINT HEADER                                                                                                          
     print,' '
     print,'--------------------------------------------------------------------------'                                 
     print,' Restored File  Median Velocity   Median Photons  Initial'                                                
     print,'                 *All Chunks*      *All Chunks*    chisq'                                                  
     print,'-------------------------------------------------------------------------'                                 

; LOOP THROUGH ALL OBSERVATIONS AND CHECK FOR BAD COUNTS OR BAD RED_CHI
     medcts_arr = fltarr(nobs)
     medchi_arr = fltarr(nobs)
     for iobs = 0, nobs-1 do begin
        restore, vdfiles[iobs]
        if keyword_set(skip_ord) then begin
           for ijk = 0, n_elements(skip_ord)-1 do begin
              xrej = where(vd.ord eq skip_ord[ijk])
              vd[xrej].templ_wt = 0.00
           endfor
        endif

        xx = where(vd.templ_wt ne 0.00000 and abs(vd.vel) le 499.9, nchunk)  ;toss templ_wt = 0.0 (tellurics)
;stop
        vd = vd[xx]
        medcts_arr[iobs] = long(median(vd.cts))  ; good chunks
        medchi_arr[iobs] = median(vd.red_chi)    ; good chunks
     ;  SUFFICIENT COUNTS AND RED_CHI LOOKS OK - KEEP the OBS

        if medcts_arr[iobs] gt mincts and medcts_arr[iobs] lt maxcts $
           ;and medchi_arr[iobs] gt minchi $
           ;and medchi_arr[iobs] lt maxchi $
        then $
           print, iobs, vd[iobs].obnm, median(vd.vel), medcts_arr[iobs], medchi_arr[iobs], f='(i6, a14, f12.2, i18, f12.2)' 

        if medcts_arr[iobs] lt mincts or medcts_arr[iobs] gt maxcts $
           ;or medchi_arr[iobs] lt minchi $
           ;or medchi_arr[iobs] gt maxchi $
        then begin
           print, iobs, vd[iobs].obnm, median(vd.vel), medcts_arr[iobs], 'tossed', f='(i6, a14, f12.2, i18, a15)'
          ; stop
        endif
     endfor ; looping through observations

; VDCUBE - STACK OF VD STRUCTURES FOR STARNAME AND TAG
   ; XKEEP: OBSERVATIONS WITH GOOD COUNTS AND GOOD CHISQ FITS 
     xkeep = where(medcts_arr gt mincts and medcts_arr lt maxcts,nxkeep); $
                   ;and medchi_arr lt maxchi $
                   ;and medchi_arr gt minchi, nxkeep) 
;print, nxkeep
;help, vdfiles
;stop
     medvel_arr = fltarr(nxkeep) 
     medcts_arr = fltarr(nxkeep)
     medchi_arr = fltarr(nxkeep) 
     vdfiles = vdfiles(xkeep)

     restore, vdfiles[0]               ; vd structure
     vdarr = replicate(vd[0], nxkeep, nchunk) ; array of good vd structures 
     diff = fltarr(nxkeep, nchunk)            ; dummy value easy to identify

     for ikp = 0, nxkeep-1 do begin  ; good chunks (counts, chisq)
        restore, vdfiles[ikp]
        xx = where(vd.templ_wt ne 0.00, nxx)
        medvel_arr[ikp] = median(vd[xx].vel) 
        medcts_arr[ikp] = long(median(vd[xx].cts))
        medchi_arr[ikp] = median(vd[xx].red_chi)
        for jcnk = 0, nchunk-1 do vdarr[ikp, jcnk] = vd[xx[jcnk]]  ; just with templ_wt ne 0.0 
     endfor 

;help, vdarr, /st
;stop

; VDCLEAN - UPDATE VDARR TO REJECT BAD CHUNKS 

; SETUP THE VST STRUCTURE (COMBINES AND WEIGHTS THE VDS FOR A SINGLE RV) 
     cf={obnm:'?',          $   ; '190704.1003'
         agitator:'?',      $   ; agitator service on?
         am:0.0,            $   ; air mass (average)
         cts:long(0),       $   ; counts in chunk 900 
         cryotemp0:0.0,    $   ;
         cryotemp1:0.0,    $   ;
         cryotemp2:0.0,    $   ;
         cryopress:0.0d,    $   ; 
         dewar: 1,          $   ; 
         epoch:0,           $   ;
         errvel:0.,         $   ;
         expt:0.0,          $   ; exposure time
         expmtr:'?',        $   ; exposure meter on?
         focpos:100L,       $  
         ha:0.0,            $
         jd:0.d,            $   ;
         mdvel:0.,          $   ; median weighted chunk velocity
         mnvel:0.,          $   ; mean weighted chunk velocity
         mdchi:0.,          $   ; 
         moondist:0.0,      $   ; degrees 
         nchunk:nchunk,     $   ; 2160
         rms:0.0,           $   ; mean RMS of residuals
         pwv:0.0,           $   ; telluric water strength
         sval:0.,           $   ; ca II H&K
         halpha_em:0.0,     $
         halpha_ew:0.0,     $
         fwhm_ccf:0.0,      $
         bis:0.0,           $
         sp1:0.,            $   ; dummy variable
         sp2:0.,            $   ; dummy variable
         wave_cal:'',       $   ; ThAr or LFC? 
         weight:0.0,        $   ;
         z:0.}                  ; v/c

     cf=replicate(cf, nxkeep)     ; cf3 structure name used by KFME etc.

   ; DON'T USE CHUNKS WITH TEMPL_WT = 0, GET STDEV FOR CHUNK SET
     mdchi = fltarr(nchunk)
     for ichk = 0, nchunk-1 do begin
        mdchi[ichk] = median( vdarr[*,ichk].red_chi )
     endfor

     err = reform(vdarr[0,*].perror[0])                       ; obs[0], all chunks
     ierr = sort(err)   &   therr = err[ierr[0.97*nchunk]]    ; error cut
     ichi = sort(mdchi) &   thchi = mdchi[ichi[0.97*nchunk]]  ; red_chi cut
     if keyword_set(maxchi) then thchi = min([maxchi, thchi])

;     igood = where(err lt therr and mdchi lt thchi, n_igood)  ;
;     n_igood replaces nchunk
     igood = where(mdchi lt thchi, n_igood)  ; n_igood replaces nchunk
     vdarr = vdarr[*,igood]

;plothist, vdarr.vel, bin=5

   ; CENTER THE CHUNK SET VELS
     for ii = 0, n_igood -1 do begin
        vset = vdarr[*,ii].vel
        vind = where(vdarr[*,ii].red_chi lt thchi, nvind)
        mean_vel = mean(vset[vind])
        ;median_vel = median(vset[vind])
        vdarr[*,ii].vel = vset - mean_vel              ; each chunk vel centered on mean
     endfor
;plothist, vdarr.vel, bin=5, col=222, /overplot

   ; START LOOKING FOR OUTLIERS - FIRST CALC DIFF B/T CHUNK RV AND MEDVEL
     diff = fltarr(nxkeep, n_igood) 
     for ob = 0, nxkeep-1 do begin
        medvel = median(vdarr[ob,*].vel)
        diff[ob, *] = vdarr[ob,*].vel - medvel 
     endfor

   ; KEEP CHUNKS WITH RED_CHI < MAXCHI AND CALC SIGMA AS STDDEV OF DIFF 
     sigma = fltarr(n_igood)    ; sigma of chunk vel - median vel
     for n = 0, n_igood - 1 do begin
        igood = where(vdarr[*,n].red_chi lt maxchi)
        sigma[n] = stddev(diff[igood,n])
     endfor 

   ; CALCULATE MAGNITUDE OF RATIO OF DIFF / SIGMA FOR EACH CHUNK - 
   ; USE THIS AS WEIGHT
     for ob = 0, nxkeep-1 do begin
        const = median(abs(diff[ob,*])/sigma)
        sigmaob = const * sigma
;print, sigma, sigmaob  
;;plothist, sigmaob 
;     sigmaob = sigma 
      ; DEFINE WEIGHT FROM SIGMA IN CHUNK RV
        for n=0, n_igood - 1 do begin
           vdarr[ob, n].weight = 1./sigmaob[n]^2
           if vdarr[ob, n].red_chi gt maxchi then vdarr[ob,n].weight=0.0d ; toss one or two bad chunks
        endfor 
     endfor

     dum = abs(diff)
     idum = sort(dum) 
     dum = dum[idum] 
     max_ind = round(0.97 * n_elements(dum))   ; keep up to 97th percentile (3-sigma)  
     lim_dum = dum[max_ind]
;     plothist, diff, xra=[-lim_dum, lim_dum], bin=5 
     x = where(abs(diff) gt lim_dum, nx)       ; zero weight for outliers
;     if nx gt 0 then stop
     if nx gt 0 then vdarr[x].weight = 0. 
 
     numobs = n_elements(vdarr[*,0])

     for ikp = 0, numobs-1 do begin
        ch_vel = reform(vdarr[ikp,*].vel)           ; chunk velocities for obs[ikp]
        diff = abs(ch_vel - median(ch_vel))         ; diff b/t chunk vel and median(vel)
        gd1=where(diff le 2.0*stddev(diff), ngd)    ; find chunks with diff < 2 sig 
;        vdarr[ikp,gd1].keep = 1
;           vd = vdarr[ikp]
;           save, vd, file=vdfiles[ikp]
        ch_vel = reform(vdarr[ikp, gd1].vel)  
        wt = reform(vdarr[ikp, gd1].weight)
        fit = reform(vdarr[ikp, gd1].red_chi)
        if strlowcase(obj_nm) ne 'sun' then $
           junk = mrdfits('/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'+obj_nm+'_'+vdarr[ikp,0].obnm+'.fits', 1, header,/silent) 
        if strlowcase(obj_nm) eq 'sun' then $
           junk = mrdfits(ddir+obj_nm+'_'+vdarr[ikp,0].obnm+'.fits', 1, header,/silent) 

        pwv = sxpar(header,'PWV', count=nfound)
        wave_cal = sxpar(header, 'WAVE_CAL', count=wfound)
        bjd = sxpar(header,'BARYMJD', count=bfound)
        sval = sxpar(header, 'S-VALUE', count=n1found) 
        if n1found eq 1 then cf[ikp].sval = sval 
        halpha_em = sxpar(header, 'HALPHA', count=n2found)
        if n2found eq 1 then cf[ikp].halpha_em = halpha_em 
        halpha_ew = sxpar(header, 'HWIDTH', count=n3found)
        if n3found eq 1 then cf[ikp].halpha_ew = halpha_ew 
        fwhm_ccf = sxpar(header, 'CCFFWHM', count=n4found)
        if n4found eq 1 then cf[ikp].fwhm_ccf = fwhm_ccf
        bis = sxpar(header, 'BIS', count=n5found)
        if n5found eq 1 then cf[ikp].bis = bis

        if strlowcase(obj_nm) ne 'sun' then $
           junk0 = mrdfits('/Volumes/G/expres/extracted/'+ddir+'/'+obj_nm+'/'+obj_nm+'_'+vdarr[ikp,0].obnm+'.fits', 0, header0,/silent) 
        if strlowcase(obj_nm) eq 'sun' then $
           junk0 = mrdfits(ddir+obj_nm+'_'+vdarr[ikp,0].obnm+'.fits', 0, header0,/silent) 
        am = sxpar(header0, 'AIRMASS', count=nf0a)
        expt = sxpar(header0, 'AEXPTIME', count=nf0b)
        expmtr = sxpar(header0, 'EXPMTR', count=nf0c)
        agitator = sxpar(header0, 'AGITATOR', count=nf0d)
        crytemp0 = sxpar(header0, 'CRYTEMP0', count=nf0e)
        crytemp1 = sxpar(header0, 'CRYTEMP1', count=nf0f)
        crytemp2 = sxpar(header0, 'CRYTEMP2', count=nf0g)
        crypress = sxpar(header0, 'CRYPRESS', count=nf0h)
        focpos = sxpar(header0, 'FOC-POS', count=nf0i)
        hourang = sxpar(header0, 'HA', count=nf0j)
        moondist = sxpar(header0, 'MOONDIST', count=nf0k)

;help, ngd, ch_vel, wt, fit 

;        vind = where(fit lt maxchi)
;        ch_vel = ch_vel[vind]
;        wt = wt[vind]
;        fit = fit[vind] 
        cf[ikp].obnm = vdarr[ikp,0].obnm              ; cf summarizes each observation
        if bfound eq 1 then cf[ikp].jd = bjd-40000.+0.5          ; baryMJD => baryRJD
        cf[ikp].mdvel = median(ch_vel)                ; median vels 
        cf[ikp].epoch = vdarr[ikp,0].epoch            ;
        if nfound eq 1 then cf[ikp].pwv = pwv
        if wfound eq 1 then cf[ikp].wave_cal = strcompress(wave_cal,/remove_all)
        cf[ikp].expt = expt
        cf[ikp].expmtr = strcompress(expmtr,/remove_all)
        cf[ikp].agitator = strcompress(agitator,/remove_all)
        cf[ikp].cryotemp0 = float(crytemp0)
        cf[ikp].cryotemp1 = float(crytemp1)
        cf[ikp].cryotemp2 = float(crytemp2)
        cf[ikp].cryopress = double(crypress)
        cf[ikp].focpos = long(focpos)
        cf[ikp].ha = ten(hourang)
        cf[ikp].moondist = float(moondist)
        cf[ikp].am = float(am)
        cf[ikp].mnvel = total(wt*ch_vel)/total(wt)    ; weighted mean
        cf[ikp].mdchi = median(vdarr[ikp, gd1].red_chi)
        cf[ikp].cts = median(vdarr[ikp, gd1].cts)
        cf[ikp].nchunk = ngd
        cf[ikp].z = median(vdarr[ikp, gd1].vel / c_light)
        cf[ikp].errvel = 1./sqrt(total(wt))           ; stddev(ch_vel/sqrt(ngd-1))
;        print, cf[ikp].obnm, '  Nchunks: ',n_elements(ch_vel)
     endfor       ; looping through good observations


; TOSS spectra w/o proper MJD in headers
 yb=where(cf.jd lt 10.,nyb)
if nyb gt 0 then begin
   print, 'Tossing the following observations: ',cf[yb].obnm
   yg = where(cf.jd gt 1.0, nyg)
   cf = cf[yg]
;stop
endif

;cf=cf[x]
xp = where(cf.ha gt 18., nxp)
if nxp gt 0 then cf[xp].ha = cf[xp].ha-24.

starname = obj_nm
velplot, cf, starname, 1./24.
cf3 = cf

; REMOVE CORRELATION WITH CRYOTEMP0

;   coef = [-3825.06, -20.3555]  ;[-3754.60, -19.9758] 
;   ny=poly(cf3.cryotemp0, coef)
;   cf3.mnvel = cf3.mnvel - ny

;   coef2 = [-44.1038, -0.624569]
;   ny2 = poly(cf3.cryotemp2, coef2)

;    if keyword_set(cryo_detrend) then begin
; 28 July 2020 - no longer think this is warrented
;       coef = [-3871.80, -20.6029] 
;       coef = [-2719.46, -14.4685]
;       coef = [-1927.76, -9.596]
;       ny2 = poly(cf3.cryotemp0, coef)
;       cf3.mnvel = cf3.mnvel - ny2
;    endif

;if ~keyword_set(excalibur) then begin
;   coef = [-3311.02, -17.6157]
;   ny=poly(cf3.cryotemp0, coef)
;   cf3.mnvel = cf3.mnvel - ny

;   coef2 = [-1741.24, -9.264]
;   nny=poly(cf3.cryotemp0, coef2) 
;   cf3.mnvel = cf3.mnvel - nny
;endif 

if keyword_set(moon) then begin
if obj_nm ne 'Sun' then begin
   x1=where(cf3.moondist lt 15.0,nx1)
   if nx1 gt 0 then for i = 0, nx1-1 do print, 'Rejecting obs closer than 15 degrees from moon: ', $
      cf3[x1[i]].obnm, cf3[x1[i]].moondist
   x11=where(cf3.moondist ge 15.0, nx11) 
   if nx11 gt 0 then cf3 = cf3[x11] 
endif
endif


; TOSS OBSERVATIONS W/O AGITATOR 
   x2a=where(strcompress(cf3.agitator,/remove_all) eq 'True',nx2a)
   x2b=where(strcompress(cf3.agitator,/remove_all) ne 'True',nx2)
   if nx2a gt 0 then cf3 = cf3[x2a] 

; TOSS OBSERVATIONS W/O EXPOSURE METER
   x3a=where(strcompress(cf3.expmtr,/remove_all) eq 'True',nx3a)
   x3b=where(strcompress(cf3.expmtr,/remove_all) ne 'True',nx3)
   if nx3a gt 0 then cf3 = cf3[x3a]

   if nx2 gt 0 then print, 'no agitator? Tossing these obs: ', cf3[x2b].obnm
   if nx3 gt 0 then print, 'no exposure meter? Tossing these obs: ', cf3[x3b].obnm

;   ffact = file_search('activity/act'+obj_nm+'.dat', count=n_cact)
;   if n_cact eq 1 then begin
;      restore, ffact
;      for i=0, n_elements(cf3)-1 do begin
;         xm1=where(cf3[i].obnm eq cact.obnm,nxm1)
;         if nxm1 gt 0 then begin
;            cf3[i].sval = cact[xm1].sval
;            cf3[i].halpha_em = cact[xm1].halpha_em
;            cf3[i].halpha_ew = cact[xm1].halpha_ew
;            cf3[i].fwhm_ccf = cact[xm1].fwhm_ccf
;            cf3[i].bis = cact[xm1].bis
;         endif
;      endfor
;   endif
   
   if obj_nm eq '16160' then pre=1
   if obj_nm eq '95128' then pre=1
   if obj_nm eq '221354' then pre=1
   if keyword_set(pre) then begin
      xeh = where(cf3.epoch lt 4, nxeh)
      if nxeh gt 0 then cf3[xeh].errvel = sqrt(cf3[xeh].errvel^2 + 9.) ; add 3 m/s errors
      save, cf3, f='vel/vst'+obj_nm+'_'+tag+'.dat'
   endif

   if ~keyword_set(pre) then begin
      xep = where(cf3.epoch ge 4,nxep)
      cf3=cf3[xep]
      save, cf3, f=vstnm
   endif

;   pr_RV, star=obj_nm, tag=tag, /dace

end ; pro 


