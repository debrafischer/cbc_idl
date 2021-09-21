pro ck_moon 

restore,'cat_drive.dat' 

nstars = n_elements(cat.starnm) 
fitsdir = '/Volumes/G/expres/extracted/excalibur/'

for i=0, nstars-1 do begin
   ss=where(cat[i].templ_obnm ne 'n',nss)
   for j=0,nss-1 do begin 
      sp = mrdfits(fitsdir+'/'+cat[i].starnm+'/'+cat[i].starnm+'_'+cat[i].templ_obnm[j]+'.fits',0,hd0, /silent)
      moond = sxpar(hd0,'moondist')
      if moond lt 30 then print, cat[i].starnm+'_'+cat[i].templ_obnm[j], ' ',moond
   endfor
endfor




end 
