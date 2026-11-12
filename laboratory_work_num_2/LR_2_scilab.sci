clear()
W1=1.2
W2=0.6

t=1.0:0.001:3.0
A=sin(t.*W1)+cos(t.*W2)

N=6
h=(3-1)/(N-1)
t0=1:h:3
A0=sin(t0.*W1)+cos(t0.*W2)

_p3=poly(t0,"x",'coeff')
p3=polyfit(t0,A0,3)
p5=polyfit(t0,A0,5)
p7=polyfit(t0,A0,7)
p9=polyfit(t0,A0,9)

Yap3 = polyval(p3, t0)
Yap5 = polyval(p5, t0)
Yap7 = polyval(p7, t0)
Yap9 = polyval(p9, t0)

plot(t,A,t0,A0,'o',t0,Yap3,t0,Yap5,t0,Yap7,t0,Yap9)
legend('A','A0','Yap3','Yap5','Yap7','Yap9')

