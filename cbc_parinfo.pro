FUNCTION CBC_PARINFO, obs, cbcenv

 ; CREATE A TEMPORARY STRUCTURE FOR FREE PARAMETERS
  
c_light = cbcenv.c_light

; SETUP THE PARINFO STRUCTURE 
  names = ['rv', 'offset']
  npar = n_elements(names)

  parinfo = {value: 0.0d,       $ ; double precision 
             fixed: 0,          $
             limited: [0,0],    $ ; use with caution 
             limits: fltarr(2), $
             parname: '?',      $
             step: 0.01d,       $
             relstep: 0.00,     $
             mpside: 2}         ; 0, 1, -1, 2    
  parinfo=replicate(parinfo, npar)

  for i=0, npar-1 do parinfo[i].parname = names[i]

 ; vel
  parinfo[0].value = 0.  ; obs.vel 
  parinfo[0].step = 10.d
  parinfo[0].limited = [1, 1]
  parinfo[0].limits = [-1000., +1000.] 

 ; continuum 
  parinfo[1].value = 1.0   ; obs.cont
  parinfo[1].step = 0.001
                             
;print, parinfo.value
 
   return, parinfo

end

