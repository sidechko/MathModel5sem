clear()
W1=1.2
W2=0.6

t=1.0:0.001:3.0
A1=sin(t.*W1)
A2=cos(t.*W2)

N=6
h=3-1/(N-1)
t0=1:h:3
A0=sin(t0.*W1)+cos(t0.*W2)
