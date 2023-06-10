%%%% SEGMENTATION


m = mirmfcc('valse_triste_happy','Rank',2:10,'Frame',0.05,1)
mx=mirgetdata(m)
pause

sim = mirsimatrix(m)
n = mirnovelty(sim,'KernelSize',150)
nx=mirgetdata(n)
pause


p = mirpeaks(n,'Contrast',.1,'Total',Inf,'NoBegin','NoEnd')
zcx=mirgetdata(p)
pause

seg = mirsegment('valse_triste_happy',p)
mirplay(seg)


display('Strike any key to continue...');
pause
close all

[seg p m a] = mirsegment('valse_triste_happy','MFCC',2:10,...
                                'KernelSize',150,'Contrast',.1)
       
display('Strike any key to continue...');
pause
close all

                                
%%%% TEMPO
                                
fb = mirfilterbank('czardas')
%mirplay(fb)
e = mirenvelope(fb) 
ex=mirgetdata(e)
pause

de = mirenvelope(e,'Diff','Halfwave')
dex=mirgetdata(de)
pause


s = mirsum(de,'Centered') 
f = mirframe(s,3,.2);
ac = mirautocor(s,'Resonance','Enhanced') 
p = mirpeaks(ac,'Total',1) 

px=mirgetdata(p)
pause

t = mirtempo(p)
tx=mirgetdata(t)
pause

display('Strike any key to continue...');
pause
close all

[t,p] = mirtempo('czardas','Periodicity','Frame')
h = mirhisto(t)
hx=mirgetdata(h)
pause

display('Strike any key to continue...');
pause
close all

%%%% TONALITY

c = mirchromagram('vivaldi','Frame',2) 
cx=mirgetdata(c)
pause


k = mirkeystrength(c) 
p = mirpeaks(k,'Total',1) 
px=mirgetdata(p)
pause

[k,p] = mirkey('vivaldi','Frame',1)



