clear()
clf()

W1=1.2*2*%pi
W2=0.6*2*%pi

t=1.0:0.001:3.0
A=sin(t.*W1)+cos(t.*W2)

N=6
h=(t($)-t(1))/(N-1)
t0=1:h:3
A0=sin(t0.*W1)+cos(t0.*W2)

tt = linspace(t(1),t($),100)

//Пункт 4
V3 = [ones(t0') t0' t0'.^2 t0'.^3]
V5 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5]
V7 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7]
V9 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7 t0'.^8 t0'.^9]

c3 = V3 \ A0'
c5 = V5 \ A0'
c7 = V7 \ A0'
c9 = V9 \ A0'

Vt3 = [ones(tt') tt' tt'.^2 tt'.^3]
Vt5 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5]
Vt7 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7]
Vt9 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7 tt'.^8 tt'.^9]

Yap3 = (Vt3 * c3)'
Yap5 = (Vt5 * c5)'
Yap7 = (Vt7 * c7)'
Yap9 = (Vt9 * c9)'

subplot(3,2,1)
plot(t, A, t0, A0, 'o', tt, Yap3, tt, Yap5, tt, Yap7, tt, Yap9)
legend('A(t)','A0(t0)','Yap3','Yap5','Yap7','Yap9')
xtitle('Пункт 4. МНК-аппроксимация (h=0.4)')

//Пункт 5
// Yint0 — линейная
Yint0 = interp1(t0, A0, tt, 'linear')

// Yint1 — кубическая сплайн-интерполяция (spline из Scilab)
Yint1 = interp1(t0, A0, tt, 'spline')

// Yint2 — nearest (ступенчатая)
Yint2 = interp1(t0, A0, tt, 'nearest')

// Yint3 — глобальный полином Лагранжа через polyfit степени N-1 = 5
pL = polyfit(t0, A0, length(t0)-1)
Yint3 = horner(pL, tt)

subplot(3,2,2)
plot(t, A, t0, A0, 'o', tt, Yint0, tt, Yint1, tt, Yint2)
legend('A(t)','A0(t0)','Yint0 (лин.)','Yint1 (сплайн)','Yint2 (nearest)','Yint3 (полином)')
xtitle('Пункт 5. Интерполяция (h=0.4)')

//Пункт 6
N=12
h=(t($)-t(1))/(N-1)
t0=1:h:3
A0=sin(t0.*W1)+cos(t0.*W2)

//h=0.2
V3 = [ones(t0') t0' t0'.^2 t0'.^3]
Vt3 = [ones(tt') tt' tt'.^2 tt'.^3]
c3 = V3 \ A0'
Yap3 = (Vt3 * c3)'

V5 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5]
Vt5 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5]
c5 = V5 \ A0'
Yap5 = (Vt5 * c5)'

V7 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7]
Vt7 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7]
c7 = V7 \ A0'
Yap7 = (Vt7 * c7)'

V9 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7 t0'.^8 t0'.^9]
Vt9 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7 tt'.^8 tt'.^9]
c9 = V9 \ A0'
Yap9 = (Vt9 * c9)'

Yint0 = interp1(t0, A0, tt, 'linear')
Yint1 = interp1(t0, A0, tt, 'spline')
Yint2 = interp1(t0, A0, tt, 'nearest')
pL = polyfit(t0, A0, length(t0)-1)
Yint3 = horner(pL, tt)

subplot(3,2,3)
plot(t, A, t0, A0, 'o', tt, Yap3, tt, Yap5, tt, Yap7, tt, Yap9, tt, Yint0, tt, Yint1, tt, Yint2)
legend('A(t)','A0(t0)','Yap3','Yap5','Yap7','Yap9','Yint0','Yint1','Yint2')
xtitle('Пункт 6.1 МНК-аппроксимация (h=0.2)')

N=24
h=(t($)-t(1))/(N-1)
t0=1:h:3
A0=sin(t0.*W1)+cos(t0.*W2)

//h=0.1
V3 = [ones(t0') t0' t0'.^2 t0'.^3]
Vt3 = [ones(tt') tt' tt'.^2 tt'.^3]
c3 = V3 \ A0'
Yap3 = (Vt3 * c3)'

V5 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5]
Vt5 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5]
c5 = V5 \ A0'
Yap5 = (Vt5 * c5)'

V7 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7]
Vt7 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7]
c7 = V7 \ A0'
Yap7 = (Vt7 * c7)'

V9 = [ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7 t0'.^8 t0'.^9]
Vt9 = [ones(tt') tt' tt'.^2 tt'.^3 tt'.^4 tt'.^5 tt'.^6 tt'.^7 tt'.^8 tt'.^9]
c9 = V9 \ A0'
Yap9 = (Vt9 * c9)'

Yint0 = interp1(t0, A0, tt, 'linear')
Yint1 = interp1(t0, A0, tt, 'spline')
Yint2 = interp1(t0, A0, tt, 'nearest')
pL = polyfit(t0, A0, length(t0)-1)
Yint3 = horner(pL, tt)

subplot(3,2,4)
plot(t, A, t0, A0, 'o', tt, Yap3, tt, Yap5, tt, Yap7, tt, Yap9, tt, Yint0, tt, Yint1, tt, Yint2)
legend('A(t)','A0(t0)','Yap3','Yap5','Yap7','Yap9','Yint0','Yint1','Yint2')
xtitle('Пункт 6.2 МНК-аппроксимация (h=0.4)')

//Пункт 7
t_ext = 0.5:0.01:3.5
A_ext = sin(t_ext.*W1) + cos(t_ext.*W2)

// МНК-полиномы (коэффициенты берём с h=0.1)
Vte = [ones(t_ext') t_ext' t_ext'.^2 t_ext'.^3]
Yap3_e = (Vte * ([ones(t0') t0' t0'.^2 t0'.^3] \ A0'))'
Vte = [ones(t_ext') t_ext' t_ext'.^2 t_ext'.^3 t_ext'.^4 t_ext'.^5]
Yap5_e = (Vte * ([ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5] \ A0'))'
Vte = [ones(t_ext') t_ext' t_ext'.^2 t_ext'.^3 t_ext'.^4 t_ext'.^5 t_ext'.^6 t_ext'.^7]
Yap7_e = (Vte * ([ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7] \ A0'))'
Vte = [ones(t_ext') t_ext' t_ext'.^2 t_ext'.^3 t_ext'.^4 t_ext'.^5 t_ext'.^6 t_ext'.^7 t_ext'.^8 t_ext'.^9]
Yap9_e = (Vte * ([ones(t0') t0' t0'.^2 t0'.^3 t0'.^4 t0'.^5 t0'.^6 t0'.^7 t0'.^8 t0'.^9] \ A0'))'

// Сплайн (вне диапазона — NaN по умолчанию)
Yint1_e = interp1(t0, A0, t_ext, 'spline')

// Глобальный полином
Yint3_e = horner(pL1, t_ext)

subplot(3,2,5)
plot(t_ext, A_ext, t0, A0, 'o',t_ext, Yap3_e, t_ext, Yap5_e,t_ext, Yap7_e, t_ext, Yap9_e, t_ext, Yint3_e)
legend('A(t)','A0(t0)','Yap3','Yap5','Yap7','Yap9','Yint3 (полином)')
xtitle('Экстраполяция на [0.5 3.5] (h=0.1)')
xgrid()
