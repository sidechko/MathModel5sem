t=1.0:0.001:3.0
W1=1.2
W2=0.6
A=sin(t.*W1)+cos(t.*W2)
W0=20*%pi
Famp=(A+1).*cos(t.*W0)
K=10
Ffm=cos((W0+K.*A).*t)
plot(t,A,t,-A,t,Famp,t,Ffm)
legend('A','-A','F AM',"F FM")
