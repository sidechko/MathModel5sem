clear()
clf()
t=1.0:0.001:3.0
W1=1.2*2*%pi
W2=0.6*2*%pi
A=sin(t.*W1)+cos(t.*W2)
W0=20*%pi
Famp=(A).*cos(t.*W0)
K=10
Ffm=cos((W0+K.*A).*t)

subplot(3,1,1)
plot(t,A,t,-A,t,Famp)
legend('A','-A','F AM')
xtitle('F amp')

subplot(3,1,2)
plot(t,A,t,Ffm)
legend('A','F FM')
xtitle('F fm')

subplot(3,1,3)
plot(t,A,t,-A,t,Famp,t,Ffm)
xtitle('Mixed')
