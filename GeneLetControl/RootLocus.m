%%
clear all
g1 = .3;g2 = .5-g1;g3=g1/g2*0.2;
p1 = 0.15; p2 = 0.1; p3 = 0.05;

pc = [1 g1+g2 g1*g2 0];
ps = conv([1 p1+p2 p1*p2],[1 p3]);
psc = conv([1  g3],[1 0.08]);
p1 = tf(psc,conv(pc,ps));

k1 = -.4; k2 =.1;

hFig=figure(2);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2*1])



rlocus(p1)
xlabel(''),ylabel(''),title('')
xlim([k1 k2])
ylim([-1 1]*0.2)