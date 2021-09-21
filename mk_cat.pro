pro mk_cat, addstar=addstar, tag=tag 

; add stars to cat_drive.dat 
; cat_drive.dat is used by dr_all.pro for CBC analysis.

; mk_cat, addstar=['3651', '4628']

tag1 = 'ame'
;restore,'cat_drive.dat' 
starnm = ['101177', '16160',   '37394', '101501', '161797',  '38858', '102121', '164922',  '43587', $
          '102195', '165341',   '4614', '103095',    '166',   '4628', '103287', '166620',  '50692', $
           '10476', '168009',  '52711', '105631',  '17156',  '55575',  '10700', '182488',  '62044', $
           '10780', '183123',   '6582', '109358', '183756', '110897', '185144',  '69830', '111395', $                  
          '185603',  '71148', '112060', '186408',  '72905', '114783', '186427',  '73732', '115617', $                  
          '187123',  '75732', '117043',  '18803',  '76151', '120066', '189733',  '80606', '120136', $                  
          '190404',  '81009', '122064', '190406',  '82885', '124755', '191785',  '84737', '126053', $                  
          '193664',  '85380', '127334', '195689',  '86728', '131156A','197076',  '88133', '131511', $                
          '199960',  '88230', '135599', '201033',  '89221', '136923', '201091',  '89269', '140538', $                  
          '201092',  '89744', '141004', '206374',  '91204', '142373', '209458',   '9407', '143761', $                  
          '210277',  '95128', '144579', '217014',  '95735', '145148', '218868',  '99491', '145675', $                  
          '219134',  '99492', '146233',  '22049',  '99995', '149661', '221354', '152391', $                  
           '23249', 'K2100b', '154345',  '25680', '157214',  '26965', '157347',   '32147', '158259', $                  
           '32923','TOI-1518','158614',  '33643', '158633',  '34411',  '159062', '159222', '3651']

ii=sort(starnm)  & starnm = starnm[ii]
num = n_elements(starnm) 

ddir='excalibur/'
fdir = '/Volumes/G/expres/extracted/'+ddir

onecat = {starnm:'',                          $
           templ_obnm:['n','n','n','n','n'],  $
           morph_rv: 0.0,                     $
           comment:'  '}

copycat = replicate(onecat, num) 

for i=0, num-1 do begin
   copycat[i].starnm = starnm[i]
;   copycat[i].templ_obnm = cat[i].templ_obnm 
   ff=file_search(fdir+copycat[i].starnm+'/'+copycat[i].starnm+'_templ_'+tag1+'60.dat', count=nf1) 
   if nf1 eq 1 then restore,ff 
   if nf1 eq 0 then begin
      ff = file_search(fdir+copycat[i].starnm+'/'+copycat[i].starnm+'_templ_*60.dat', count=nf2)
      if nf2 gt 0 then ff = ff[0]
   endif
   if nf1 gt 0 or nf2 gt 0 then begin
      restore, ff
      tobnm = templ[0].obnm
      ft = file_search(fdir+copycat[i].starnm+'/'+strmid(tobnm,0,strlen(tobnm)-5)+'*.fits', count=ntempl)
      if ntempl gt 0 then begin
         if ntempl gt 5 then ntempl = 5
         for j=0, ntempl-1 do begin
            i0 = strpos(ft[j],'_')
            obsnm = strmid(ft[j], i0+1, 11)
            copycat[i].templ_obnm[j]= obsnm
         endfor
      endif
      xgd = where(templ.morph_rv ne 0.0)
      copycat[i].morph_rv = median(templ[xgd].morph_rv)
;      stop   
   endif else copycat[i].comment = 'Wait'
endfor

if keyword_set(addstar) then begin
   nnum = n_elements(addstar) 
   ncat = replicate(onecat, nnum)
   for i = 0, nnum-1 do ncat[i].starnm = addstar[i]
endif

if keyword_set(addstar) then cat = [copycat, ncat] else cat = copycat

save, cat, f='cat_drive_new.dat'


end


