pro  rayclean,obsnm,obstack,star,nocoadd=nocoadd,diff=diff,auto=auto, $
              date=date, averg=averg, starnm=starnm, excalibur=excalibur, $
              epoch=epoch, outfname=outfname, ddir=ddir, solar=solar

;obsnm [input]  string array    observation names used to be "ray cleaned
;                               and coadded to form template star
;star  [output] floating array  coadded stellar template
;        
;Created PB, Nov 27, 1995

;starnm = '101501'
;obsnm=['190210.1153','190210.1154','190210.1155']
;obsnm=['190425.1090','190425.1091']
;obsnm=['190425.1090','190425.1091','190427.1068','190503.1087','190505.1069','190518.1122','190531.1107']
ordr = indgen(50)+29   ;  ;just operate on EXPRES LFC orders
nnord=n_elements(ordr)
ordrstnd = 55

px0=350 
px1 = px0+6800-1  ; epochs 3 - 5


path='/Volumes/G/expres/extracted/'+ddir+'/'+starnm+'/'
tpath=path

if keyword_set(excalibur) then begin 
   px0 = 630
   px1 = 630+6380-1  
endif

nobs=n_elements(obsnm)
fnm = path+starnm+'_'+obsnm[0]+'.fits'

ob=mrdfits(path+starnm+'_'+obsnm[0]+'.fits',1)
sp=ob.spectrum/ob.continuum  ; normalized spectrum
w = ob.bary_wavelength
;w = ob.wavelength
if keyword_set(excalibur) then w = ob.bary_excalibur

pfilt = sp*1.0   ; perfect chip
snr_arr = fltarr(nobs)
;plot, w[*,55],sp[*,55], /xsty
;stop

for i=0,nobs-1 do snr_arr[i] = ck_snr(dir=path, obnm=obsnm[i])
xx=where(snr_arr lt 150, nxx) 
if nxx gt 0 then stop, 'reject and start over? '

snr = ck_snr(dir=path, obnm=obsnm[0])
PRINT,'median counts at order '+strtrim(ordrstnd,2)+' = ',snr^2
obstack = fltarr(n_elements(sp[*,0]),n_elements(sp[0,*]),nobs)
obstack[*,*,0] = sp

for n=1,nobs-1 do begin
   ob=mrdfits(path+starnm+'_'+obsnm[n]+'.fits',1)
   sp=(ob.spectrum/ob.continuum)
   if ~keyword_set(excalibur) then w=ob.bary_wavelength
;   w=ob.wavelength
   if keyword_set(excalibur) then w = ob.bary_excalibur
   snr_test=ck_snr(dir=path, obnm=obsnm[n]) ;[mdcts,median(sp[*,ordrstnd])]
   mdcts = snr_test^2
   print,'Median counts at order '+strtrim(ordrstnd,2)+' for '+obsnm[n]+'    '+strtrim(mdcts,2)
   obstack[*,*,n]=sp
endfor

if nobs eq 2 then begin                         ; 2 observations, use RAYZAP
   print,'Cleaning with RAYZAP'
   ob=reform(obstack[*,*,0])
   dum=reform(obstack[*,*,1])
   if n_elements(diff) ne 1 then diff=3.5
   mincts=2000 ; min([2000,0.25*median(mdcts)])  ;2000 DN or 1/4 of median exposure
   rayzap,ob,dum,pfilt,ordr,mincts, px0=px0, px1=px1  
   obstack[*,*,0] = ob
   obstack[*,*,1] = dum
endif else raystack,obstack,pfilt,ordr,snr_arr=snr_arr, mdtck=mdtck,px0=px0, px1=px1, auto=auto ; >2 observations, use RAYSTACK

;coadding the "cleaned" observations
if not keyword_set(nocoadd) and nobs gt 2 then begin
   print,'Now creating a median template observation'
   npix = n_elements(obstack[*,0,0])
   ncol = n_elements(obstack[0,*,0])
   star = fltarr(npix, ncol) 
   for i = 0, npix-1 do begin
      for j=0, ncol-1 do begin
         star[i, j]=median(obstack[i,j,*])
      endfor ; j (col) 
   endfor ; i (pix)
endif

;coadding the "cleaned" observations
if not keyword_set(nocoadd) and nobs eq 2 then begin
   print,'Now creating an average template observation'
   npix = n_elements(obstack[*,0,0])
   ncol = n_elements(obstack[0,*,0])
   star = fltarr(npix, ncol) 
   for i = 0, npix-1 do begin
      for j=0, ncol-1 do begin
         star[i, j]=mean(obstack[i,j,*])
      endfor ; j (col) 
   endfor ; i (pix)
endif

save, star, f=tpath+outfname 

;stop
end
