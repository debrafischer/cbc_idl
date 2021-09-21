;
;Plot Radial Velocity vs. Time, based on data in VST structure.
;VST structure is generated in VEL.PRO
;
;INPUT:
;  velst  (input structure)  from VEL.PRO.  Contains Time of Obs (VELST.JD),
;                        mean and median velocities  (velst.mnvel, velst.mdvel)
;   
;OPTIONAL INPUT:
;  starname       (string)     Name of Star
;  timebin        (float)      Time interval (DAYS) to "bin" velocities
;  tspan          fltarr(2)    vector of starting and ending time (JD or yrs)
;  cai            float        Calcium IR triplet index
;
;OPTIONAL KEYWORDS:
;  yrs      (0 or 1)  /yrs  forces abscissa to be in years, not JD   
;  nocolor  (0 or 1)  /nocolor for a monochrome terminal
;  nosig    (0 or 1)   if nosig eq 0 then standard default, else don't xyout,sig
;  nobox    (0, 1, 2)  if nobox eq 0 then standard default
;                      if nobox eq 1 then don't plot the timebin averages
;                      if nobox eq 2 then don't plot the error bars
;  errcut   (real, 1.5 or 2 or 2.5) toss points with cf.errvel worse than (errcut)*median(cf.errvel)
;
pro velplot, $
velst, $
starname, $
timebin, $
dates, $
filen=filen, $
speed, $
errv, $
cai, $
nav, $
dumcf, $
yrs=yrs, $
tspan=tspan, $
nocolor=nocolor, $
color=color, $
yra=yra, $
nobox=nobox, $
median=median, $
nosig=nosig, $
subtt=subtt, $
noplot=noplot, $
nozero=nozero, $
nave=nave, $
lzr=lzr, $
errcut=errcut, $
szfac=szfac, $
linear=linear, $
zeroline=zeroline, $
table=table, $
two_obs=two_obs


IF n_params () lt 1 then begin
  print,'Syntax: '
  print,' VELPLOT,velst, [starname, timebin (days), dates, speed, errv], yrs=yrs, tspan=[start,end],nocolor=nocolor'
  return
ENDIF


if n_elements(szfac) lt 1 then szfac=0.5
dd=velst.jd

color_orig = !p.color
back_orig  = !p.background
!p.font=1
   black = 0
   darkblue = 15
   blue = 30
   purple = 50
   red = 70
   orange = 85
   lightorange = 100
   algaegreen = 122
   yellow = 130
   white = 255
if keyword_set(color) then begin   ;COLOR SETTINGS
   loadct,39
  !p.background = 255  ;white  background
  !p.color = 90      ;DEFAULT COLOR
  pcolor=!p.color
  col = [30, 50, 70, 85, 100, 120, 130] 
  nocolor=0
END ELSE begin                   ;BLACK & WHITE SETTINGS
  loadct,0
  !p.color = 0        ;black
  !p.background = 255 ;white background
  col = [0, 0, 0, 0, 0, 0, 0]   ;black 
  nocolor = 1         ;in case "nocolor" arises.
ENDELSE

if min(dd) gt 2440000. then dd=dd-2440000.              ;reduced JD
;	

   av = velst              ;Protect Structure containing results from VEL.PRO
tagnam=tag_names(av)
obdex=first_el(where(tagnam eq 'OBNM'))
if obdex eq -1 then obdex=first_el(where(tagnam eq 'OBNAM'))

   IF n_elements(nave) eq 1 then if nave eq 1 then begin   ;apply nightly corrections
      print,'Making **Primitive** Nightly Correction'
print,'here here here'
      rascii,nave,4,'~/idle/vel/nave.ascii',skip=1
      for n=0,n_elements(av)-1 do begin
         nnd=minloc(abs(reform(nave(0,*)) - dd(n)),/first) ;make correction
         if abs(dd(n)-nave(0,nnd)) lt 0.6 then av(n).mnvel=av(n).mnvel-nave(1,nnd) $
	    else print,'No nightly correction for the night of: '+strtrim(dd(n),2)
      endfor
   ENDIF

   vel = av.mnvel
   errvel = av.errvel
   cat = intarr(n_elements(errvel))             ;CAT initially set to 0
   cair = av.sp1                                ;Calcium IR triplet index
;   dewar = av.dewar 

   if keyword_set(median) then vel=av.mdvel
   if not keyword_set(nozero) then vel=vel-median(vel)
   numobs = n_elements(av.jd)
   if n_elements(starname) lt 1 then starname=' '
   if n_elements(timebin) ne 1 then timebin=2./24.
   if timebin gt (max(dd)-min(dd)) then timebin = 1.e-10
   if n_elements(nobox) ne 1 then nobox=0
   if n_elements(nosig) ne 1 then nosig=0
   if n_elements(subtt) ne 1 then subtt=' '

;TIME MATTERS (ABSCISSA: JD or years, TIMEBIN UNITS)
      xarr = av.jd
      xtit = '!6 JD - 2440000'
      timetrim = 0.
;      prepost=9675.                              ;pre/post Julian Date
       IF keyword_set(yrs) then begin                             ; YEARS?
;         prepost=1994.87                         ;pre/post year
;         xarr = 1987. + (xarr - 6795.5)/365.25   ;1987.0 = 2446795.5 JD
;   the above conversion is off by 24 hours  DF 1-24-00
;         xarr = 1987. + (xarr - 6796.5)/365.25   ;1987.0 = 2446796.5 JD
;  New Constant:juldate,[2000,1,1,0,0],jd  --> jd = 51544.5
;;;       tconst = 11544.5
;;;       xarr = 2000. + (xarr - tconst)/365.25   ;1987.0 = 2446795.5 JD
;       tconst = 15561.5       

      tconst = 15197.5          ;jd for 2010 at UT midnight: 2455197.5
       xarr = 2010. + (xarr - tconst)/365.25   ;jd for 2010 at UT midnight: 2455561.5
         xtit = '!6Time (Years)'
         timebin = timebin/365.25                   ;timebin in yrs.

         if ( max(xarr) - min(xarr) ) lt 3.0 then begin
	    xtit='!6Time  (Year-'+strtrim(fix(min(xarr)),2)+')'
;	    prepost=prepost-fix(min(xarr))
	    xarr=xarr-fix(min(xarr))
         endif  
      endif else if max(xarr) lt (min(xarr) + 300) then begin    ; <300 days?
	 dumtime = fix(xarr[0]/100.) * 100.
         timetrim = dumtime             ;time trimmed off original JD's
	 xarr = av.jd-dumtime
	 dumst = strtrim(long(2400000)+long(dumtime),2)
         xtit = '!6 Julian Date - '+dumst
      ENDIF

      t0 = long(min(xarr))
      if max(av.jd) lt (min(av.jd)+1) then begin                 ; <1 day?
	 xt = 'Julian Date - '+strtrim(long(2440000)+long(min(av.jd)),2)
;         timetrim = -13200.   ;-(min(av.jd)-5.)
	 xarr = av.jd - min(av.jd)
      endif

;PLOT CHARACTERISTICS
      !p.charsize=1.5
      !x.charsize=2
      !y.charsize=2
      !p.thick=2
      yt='!6Velocity  (m s!e-1!n)'
      starname = '!6'+string(starname)  ;make duplex Roman Font
      dt = max(xarr) - min(xarr)                     ;time range
      xr = [min(xarr)-0.05*dt, max(xarr)+0.05*dt ]   ;plot limits
      if keyword_set(tspan) then xr = tspan
      yr = [-10,10]
      if max(vel + errvel) gt yr(1) or min(vel - errvel) lt yr(0) then $
         yr = [min(vel-errvel)-10, max(vel+errvel)+10]
;      if (max(vel)-min(vel) lt 40) then yr=[-30,30]
      if max(yr) lt max(vel) then yr=yr+10   ;PB kludge, 24 Sept 1997
      if n_elements(yra) eq 2 then yr = yra
      a=findgen(32) * (!pi*2/32.)
      usersym,cos(a),sin(a),/fill      ;circle plot sym #8
;
;COLOR SETTINGS
;   IF not keyword_set(nocolor) then BEGIN     ;for Color Monitor
;     loadct,13
;     col = [100, 180, 80, 90, 260, 240, 141]
;   END ELSE begin
;;     loadct,0
;     col = [300, 300, 300, 300, 300, 300]                  ;for B+W Monitor
;   ENDELSE
;

if keyword_set(filen) then begin
set_plot,'ps'
device, filename = filen
endif
if not keyword_set(noplot) then begin
if keyword_set(nocolor) then begin
   if keyword_set(lzr) then ps_open,'idl'
   !x.margin=[7,5]
   !y.margin=[6,7]
   PLOT, xarr, vel, xtit=xtit, ytit=yt, titl=starname,charsize=1.0, $
            xr=xr, yra=yr, /xsty, /ysty, psym=8,syms=0.4*szfac, subtitle=subtt 
endif else PLOT, xarr, vel, xtit=xtit, ytit=yt, titl=starname, co=90,charsize=1.0, $
            xr=xr, yra=yr, /xsty, /ysty, psym=8,syms=0.4*szfac, subtitle=subtt 
;
;DEWAR OVERPLOT
;if not keyword_set(nocolor) then begin
;    i = where(av.dewar eq 1 or av.dewar eq 2,n) 
;      if n gt 0 then OPLOT,xarr(i),[vel(i)],psym=8,co=col(1),syms=1.3
;    i = where(av.dewar eq 6,n)
;      if n gt 0 then OPLOT,xarr(i),[vel(i)],psym=8,co=col(2),syms=1.4
;    i = where(av.dewar eq 8,n)
;      if n gt 0 then OPLOT,xarr(i),[vel(i)],psym=8,co=col(3),syms=1.3
;    i = where(av.dewar eq 39,n)
;      if n gt 0 then OPLOT,xarr(i),[vel(i)],psym=8,co=col(6),syms=1.3
;    i = where(av.dewar eq 13,n)
;      if n gt 0 then OPLOT,xarr(i),[vel(i)],psym=8,co=col(4),syms=1.3
;         IF n ge 1 then begin               ;Blink Dewar 13 points
;            for blink=0,2 do begin
;               wait,.2
;               oplot,xarr(i),[vel(i)],psym=8,co=0,syms=1.3 & wait,.2
;               oplot,xarr(i),[vel(i)],psym=8,co=col(5),syms=1.3 
;            endfor
;         ENDIF
;endif  ;not keyword_set(nocolor)

;TELESCOPE OVERPLOT
   for qq=0,n_elements(cat)-1 do begin
      dum=chip(av(qq).(obdex),gain,telescope)
      if telescope eq 'CAT' or telescope eq 'SHANE' then cat(qq)=1     ;cat is telescope flag
      if telescope eq 'KECK' then cat(qq)=2
   endfor
;   ind=where(cat eq 1,nind)
;   if nind gt 0 then if keyword_set(nocolor) then $ 
;      oplot, xarr(ind), vel(ind), psym=8,syms=1.3 else $
;      oplot, xarr(ind), vel(ind), co=211 , psym=8,syms=1.3


   sigall = stdev(vel)
       xo = max(xarr)-0.45*dt
       dy = yr(1) - yr(0)
       yo = yr(1) - 0.05*dy
   text = strmid(strtrim(sigall,2),0,4)+' ms!u-1!n'
   IF not keyword_set(nocolor) then !p.color=col(0)
   ;if szfac eq 1 then if nosig eq 0 then XYOUTS,xo,yo,'!7r!6!d ALL!n ='+ text,size=1.5
;   oplot,[xo-0.025*dt],[yo+0.007*dy],psym=8,syms=1.2
;   oplot,[xo-0.025*dt,xo-0.025*dt],[yo+0.04*dy,yo-0.03*dy]
   print,format='(a15,F6.2,a4)','Sigma (all) = ',sigall,' m s!u-1!n'
endif  ;if not keyword_set(noplot) 
 
dewar=av.dewar
;ob=strarr(n_elements(dewar))
dew=fltarr(n_elements(dewar))
;for i=0,numobs-1 do dewar(i)=0

;Time BINS
IF timebin eq 0 then timebin = 1.e-10    ;tiny time bin
If timebin gt (max(dd)-min(dd)) then timebin = 1.e-10
     ct=0
WHILE ct lt numobs do begin
      ind = where(xarr ge xarr(ct) and xarr lt xarr(ct)+timebin, num) 
      test=av(ind).dewar
;print,'test: ',test
;print,'ind: ',ind
      ntest=n_elements(test)
        if min(test) ne mean(test) then begin
print,'trouble! trying to bin observations from different telescopes.'
  	ind1=where(xarr ge xarr(ct) and xarr lt xarr(ct)+timebin and strmid(av.obnm,0,2) ne 'rk')  ;select all Lick obs
        ind2=where(xarr ge xarr(ct) and xarr lt xarr(ct)+timebin and strmid(av.obnm,0,2) eq 'rk')  ;select all Keck obs
print,'ind1: ',ind1
print,'ind2: ',ind2

        if ind2(0) ne -1 then begin
          ind=ind2					;keck observations only
          tmp=av(ind).dewar
          wt = (1./errvel(ind))^2.       ;weights based on internal errors
          wt = wt/total(wt)              ;normalized weights
          if ct eq 0 then catday = [max(cat(ind))] $    ;telescope flag
	    else catday = [catday,max(cat(ind))]
          if ct eq 0 then ob = [av(ind(0)).obnm] $
	    else ob = [ob,av(ind(0)).obnm]
          if ct eq 0 then nav = [n_elements(ind)] $
            else nav = [nav,n_elements(ind)]
          if ct eq 0 then dates = [total(xarr(ind)*wt)]      $
	    else dates = [dates,total(xarr(ind)*wt)]
          if ct eq 0 then speed = [total(vel(ind)*wt)] $
	    else speed=[speed,total(vel(ind)*wt)]
          if ct eq 0 then counts = [total(av(ind).cts)] $
	    else counts = [counts,total(av(ind).cts)]
          if ct eq 0 then cai = [total(cair(ind)*wt)] $
	    else cai=[cai,total(cair(ind)*wt)]
          dew(ct)=min(tmp)
          inverr=0. & for qq=min(ind),max(ind) do inverr=inverr+1./errvel(qq)^2.
	  errv = [errv,sqrt(1./inverr)]
          ct=ct+num
          ind=ind1					;Lick only
        endif
      endif						;end trouble 

          tmp=av(ind).dewar
          wt = (1./errvel(ind))^2.       ;weights based in internal errors
          wt = wt/total(wt)              ;normalized weights
          if ct eq 0 then catday = [max(cat(ind))] $
	    else catday = [catday,max(cat(ind))]
          if ct eq 0 then ob = [av[ind[0]].obnm] $
	    else ob = [ob,av[ind[0]].obnm]
          if ct eq 0 then nav = [n_elements(ind)] $
            else nav = [nav,n_elements(ind)]
          if ct eq 0 then dates = [total(xarr(ind)*wt)]      $
	    else dates = [dates,total(xarr(ind)*wt)]
          if ct eq 0 then speed = [total(vel(ind)*wt)] $
	    else speed=[speed,total(vel(ind)*wt)]
          if ct eq 0 then counts = [total(av(ind).cts)] $
	    else counts = [counts,total(av(ind).cts)]
          if ct eq 0 then cai = [total(cair(ind)*wt)] $
	    else cai=[cai,total(cair(ind)*wt)]
          dew(ct)=min(tmp)
          if ct eq 0 then begin
            inverr=0. & for qq=min(ind),max(ind) do inverr=inverr+1./errvel(qq)^2.
	    errv = [sqrt(1./inverr)]
          endif else begin
            inverr=0. & for qq=min(ind),max(ind) do inverr=inverr+1./errvel(qq)^2.
	    errv = [errv,sqrt(1./inverr)]
          endelse
          ct=ct+num 
    END  ;while
;print,'nav=',nav

x=where(dew gt 0.0)
dwr=fix(dew(x))
     badgate=0
     if n_elements(errcut) eq 1 then if errcut gt 0. then begin
	i0=where(dates lt prepost,ni0)
        if ni0 gt 1 then begin
        dumz=where(errv ge (errcut*median(errv(i0))) and dates lt prepost)
	if dumz(0) ge 0 then begin
	   baddate=dates(dumz)
	   badspeed=speed(dumz)
	   baderrv=errv(dumz)
	   badgate=1
        endif
        dumx=where(errv lt (errcut*median(errv(i0))) and dates lt prepost)
	dates0=dates(dumx)
  	if max(dates0) le 300. then dates0=dates0+dumtime
	ob0=ob(dumx)
	speed0=speed(dumx)
	errv0=errv(dumx)
	counts0=counts(dumx)
	cai0=cai(dumx)
        dwr0=dwr(dumx)
	endif   ;ni0 gt 1

	i1=where(dates gt prepost,ni1)
	if ni1 gt 1 then begin
          dumz=where(errv ge (errcut*median(errv(i1))) and dates gt prepost)
	  if dumz(0) ge 0 then begin
	     if n_elements(baddate) gt 0 then baddate=[baddate,dates(dumz)] $
	        else baddate=dates(dumz)
	     if n_elements(badspeed) gt 0 then badspeed=[badspeed,speed(dumz)] $
	        else badspeed=speed(dumz)
	     if n_elements(baderrv) gt 0 then baderrv=[baderrv,errv(dumz)] $
	        else baderrv=errv(dumz)
	     badgate=1
          endif
          dumx=where(errv lt (errcut*median(errv(i1))) and dates gt prepost)
	  if n_elements(dates0) gt 0 then dates=[dates0,dates(dumx)] else dates=dates(dumx)
	  if n_elements(speed0) gt 0 then speed=[speed0,speed(dumx)] else speed=speed(dumx)
	  if n_elements(errv0) gt 0 then errv=[errv0,errv(dumx)] else errv=errv(dumx)
	  if n_elements(ob0) gt 0 then ob=[ob0,ob(dumx)] else ob=ob(dumx)
	  if n_elements(counts0) gt 0 then counts=[counts0,counts(dumx)] else counts=counts(dumx)
	  if n_elements(cai0) gt 0 then cai=[cai,cai(dumx)] else cai=cai(dumx)

	  if n_elements(dwr0) gt 0 then dwr=[dwr0,dwr(dumx)] else dwr=dwr(dumx)
	endif else begin  ;if ni1 gt 1
          dates=dates0
          speed=speed0
          errv=errv0
          counts=counts0
          cai=cai0
          dwr=dwr0
        endelse
     endif
     dumcf=replicate(av(0),n_elements(dates))
     if n_elements(ob) ne n_elements(dates) then $
	dumcf.obnm = ob0 else dumcf.obnm = ob
     dumcf.jd     = dates ;+ timetrim
     if n_elements(dwr) eq n_elements(dates) then dumcf.dewar = dwr
     dumcf.mnvel  = speed
     dumcf.errvel = errv
     dumcf.cts    = counts
;     dumcf.sp1    = cai
;   endif   ;timebin
 
    sigbin = stdev(speed)

if not keyword_set(noplot) then if nobox lt 1 then begin 
       IF keyword_set(nocolor) then OPLOT, dates, speed, psym=8, symsize=1.3*szfac, thick=2 $
           ELSE OPLOT, dates, speed, psym=8, symsize=1.1*szfac, thick=2, co=90 ;col(5)
       If badgate eq 1 then begin
          IF keyword_set(nocolor) then OPLOT, baddate, badspeed, psym=7, symsize=1.3*szfac, thick=2 $
           ELSE OPLOT, baddate, badspeed, psym=7, symsize=1.3*szfac, thick=2, co=90
	   for bb=0,n_elements(baddate)-1 do begin
	      bdd=[baddate(bb),baddate(bb)]
	      bvv=[badspeed(bb)-baderrv(bb),badspeed(bb)+baderrv(bb)]
              IF keyword_set(nocolor) then OPLOT, bdd, bvv $
                 ELSE OPLOT, bdd, bvv, thick=2, co=90,symsize=0.5*szfac
           endfor
       endif  ;badgate
	 ind1=where(catday eq 1,nind1)
	 ind2=where(catday eq 1,nind2)
	 if nind1 gt 0 then if keyword_set(nocolor) then $
	   OPLOT, dates(ind1), speed(ind1), psym=8, symsize=1.3*szfac, thick=2 ELSE $
           OPLOT, dates(ind1), speed(ind1), psym=8, symsize=1.3*szfac, thick=2, co=90
	 if nind2 gt 0 then if keyword_set(nocolor) then $
	   OPLOT, dates(ind2), speed(ind2), psym=8, symsize=1.3*szfac, thick=2 ELSE $
           OPLOT, dates(ind2), speed(ind2), psym=8, symsize=1.3*szfac, thick=2, co=90
       yo = yo - 0.05*dy
       text = strmid(strtrim(sigbin,2),0,4)+' ms!u-1!n'
       IF not keyword_set(nocolor) then !p.color=90  ;col(0)
       if nosig eq 0 and keyword_set(lzr) then XYOUTS,xo,yo,cgGreek('sigma')+'!6!d BIN!n ='+ text,size=1.1
       if nosig eq 0 and ~keyword_set(lzr) then XYOUTS,xo,yo,'!9s !6!d BIN!n ='+ text,size=1.2
       print,format='(a15,F6.2,a4)','Sigma (binned) = ',sigbin,' m/s'
endif      ; if not keyword_set(noplot), nobox lt 1 

;openw,4,'6623.dat'
;for jd instead of years, remove /yrs in the call to velplot
;for i=0,n_elements(speed)-1 do printf,4,dates(i)+timetrim,speed(i),errv(i)
;close,4 
form='(F10.3, a3, F8.2, a3, F8.2,a3,a6,a4)'

if keyword_set(table) and not keyword_set(two_obs) then begin
  openw,4,'table.tex'
  n_tbl=n_elements(speed)
  chksz=n_tbl/35.  & npages=fix(n_tbl/35.)
  regpag=npages  ;number of pages with 35 entries
  extr=n_tbl-(regpag*35)  ; number of entries for last page
  if chksz gt npages then npages = npages+1
  for j=0,npages-1 do begin
     ;table header
     printf,4,'\begin{table}'
     printf,4,'\tablenum{?}'
     printf,4,'\caption{Radial Velocities for ??Star??}'
     printf,4,'\begin{tabular}{rrrc}'
     printf,4,'\tableline'
     printf,4,'\tableline'
     printf,4,'       JD  &   RV    & Uncertainties  \\'
     printf,4,'           &  (\ms)  &   (\ms)        \\'
     printf,4,'\tableline'
     ;print full pages with 35 lines
     if j lt regpag then begin
       for i=0, 34 do begin
	  ii=i+j*35
	  printf,4, format=form,dates(ii),' & ',speed(ii),' & ',errv(ii),'  \\ '
       endfor
     endif
     ;print last, partial page
     if j eq regpag then begin
       for i=0,extr-1 do begin
	  ii=i+j*35
	  printf,4, format=form,dates(ii),' & ',speed(ii),' & ',errv(ii),'  \\ '
       endfor
     endif
     ;print endlines on table
     printf,4,'\tableline'				  
     printf,4,'\end{tabular}'
     printf,4,'\end{table}'
     printf,4,'\clearpage'
     printf,4,' '
  endfor
  close,4 
endif

if keyword_set(table) and keyword_set(two_obs) then begin
  openw,4,'table.tex'
  n_tbl=n_elements(speed)
  chksz=n_tbl/35.  & npages=fix(n_tbl/35.)
  regpag=npages  ;number of pages with 35 entries
  extr=n_tbl-(regpag*35)  ; number of entries for last page
  if chksz gt npages then npages = npages+1
  for j=0,npages-1 do begin
     ;table header
     printf,4,'\begin{table}'
     printf,4,'\tablenum{??}'
     printf,4,'\caption{Radial Velocities for ??Star??}'
     printf,4,'\begin{tabular}{rrrc}'
     printf,4,'\tableline'
     printf,4,'\tableline'
     printf,4,'       JD    &      RV       & Uncertainties & Observatory \\'
     printf,4,' -2440000    &     (\ms)     &   (\ms)       &         \\'
     printf,4,'\tableline'
     ;print full pages with 35 lines
     if j lt regpag then begin
       for i=0, 34 do begin
	  ii=i+j*35
	  if velst(ii).dewar eq 102 then obs='Keck' else obs='Lick'
	  printf,4, format=form,dates(ii),' & ',speed(ii),' & ',errv(ii),' & ',obs, ' \\ '
       endfor
     endif
     ;print last, partial page
     if j eq regpag then begin
       for i=0,extr-1 do begin
	  ii=i+j*35
	  if velst(ii).dewar eq 102 then obs='Keck' else obs='Lick'
	  printf,4, format=form,dates(ii),' & ',speed(ii),' & ',errv(ii),' & ',obs, ' \\ '
       endfor
     endif
     ;print endlines on table
     printf,4,'\tableline'				  
     printf,4,'\end{tabular}'
     printf,4,'\end{table}'
     printf,4,'\clearpage'
     printf,4,' '
  endfor
  close,4 
endif

      if n_elements(errv) lt 2 then begin
	 errv = av.errvel                         ;error in the mean
	 dates = xarr 
	 speed = av.mnvel
	 if n_elements(catday) ne n_elements(speed) then catday=fix(speed*0)
      endif

 if not keyword_set(noplot) then begin
      FOR j=0,n_elements(errv)-1 do begin
        err = errv(j)                           ;error in the mean
        x = [dates(j), dates(j)]
        bar = [speed(j)-err, speed(j)+err]
        if nobox ne 2 then begin
	   if keyword_set(nocolor) then oplot,x, bar $
               else oplot,x, bar, co=90 ;col(5)
           if catday(j) eq 1 then $
	       if not keyword_set(nocolor) then oplot,x, bar,co=90
        endif
      ENDFOR
     yo = yo - 0.05*dy
     text = strmid(strtrim(median(errv),2),0,4)+' ms!u-1!n'
;     glet = "163B
;     gsig = '!9'+string(glet)
     if nosig eq 0 and keyword_set(lzr) then XYOUTS,xo,yo,cgGreek('sigma')+'!6!d INT!n ='+ text,size=1.1
     if nosig eq 0 and ~keyword_set(lzr) then XYOUTS,xo,yo,'!9s !6!d INT!n ='+ text,size=1.2
  endif
if keyword_set(linear) then begin
	coef=poly_fit(dates,speed,1)
	thv=poly(dates,coef)
	linvel=speed-thv
	oplot,dates,thv
endif

if keyword_set(zeroline) then plots,[min(dates),max(dates)],[0.,0.]

if keyword_set(lzr) then ps_close
if keyword_set(filen) then begin
device, /close
set_plot, 'x'
endif ;postscript plot because KW:filen set. 

return
end



