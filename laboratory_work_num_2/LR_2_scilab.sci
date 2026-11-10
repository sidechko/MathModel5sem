
W1=1.2
W2=0.6

t=1.0:0.001:3.0
A=sin(t.*W1)+cos(t.*W2)

h0=0.4
t0=1:h0:3
A0=sin(t0.*W1)+cos(t0.*W2)

p3=polyfit(t0,A0,3)
p5=polyfit(t0,A0,5)
p7=polyfit(t0,A0,7)
p9=polyfit(t0,A0,9)

xx = linspace(t(1), t($), 100)

Yap3 = horner(p3, xx)
Yap5 = horner(p5, xx)
Yap7 = horner(p7, xx)
Yap9 = horner(p9, xx)
clf();
plot(xx,Yap3)
xgrid()

