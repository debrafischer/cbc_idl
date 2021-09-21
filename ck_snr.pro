function ck_snr, dir=dir, obnm=obnm

;snr = ck_snr(dir=dir, obnm=obnm) 

ff = file_search(dir+'*'+obnm+'*', count=count) 

if count eq 0 then stop, 'file not found' 

;print, ff[0]
obs = mrdfits(ff[0] , 1, /silent) 
sp = (obs[55].spectrum/obs[55].uncertainty)^2 

snip = sp[3000:5000] 
ii=sort(snip)
iirev = reverse(ii) 
snip = snip[iirev] 

snr = sqrt(median(snip[0:10]) )
;print, 'SNR at center of order 55: ', snr 
if snr lt 100 then print, 'Consider rejecting: ',ff[0], ' SNR: ',snr

return, snr 

end ;pro
