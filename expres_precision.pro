pro expres_precision, hdcopy=hdcopy 


ck_precision, tag='amz', /all, pre=pre

if keyword_set(hdcopy) then ps_open, 'test',/color

plothist, pre.rms_all, bin=0.3, xra=[0,4], yra=[0,40], fcolor=0, /fill, xtitl='!6 meters per second'
plothist, pre.rms_nightly, bin=0.2, /overplot, col=90, /fline, forientation=-45, /fill, fcolor=90, thick=5, fthick=6
plothist, pre.mean_err, bin=0.1, /overplot, col=222, /fline, forientation=45, /fill, fcolor=222, thick=5, fthick=6

xyouts, 1,36,/data,'!6 Single Observation Precision [m/s]', color=222
xyouts, 1,32,/data,'!6 RMS of consecutive Observations [m/s]', color=90 
xyouts, 1,28,/data,'!6 RMS over 2.5 years [m/s]', color=0

if keyword_set(hdcopy) then ps_close 

end ; pro 
