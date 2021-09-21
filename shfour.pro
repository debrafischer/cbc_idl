pro shfour,sp,shift,newsp
ln=n_elements(sp)
nsp=sp
fourtr=fft(nsp,1)
sig=findgen(ln)/ln-.5
sig=shift(sig,(ln/2))
sh=sig*2.*!pi*shift
shfourtr=complex(cos(sh),sin(sh))*fourtr
newsp=fft(shfourtr,-1)
newsp=float(newsp(0:ln-1))
end  ;shfour
