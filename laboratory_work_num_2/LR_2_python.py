import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

W1 = 1.2 * 2 * np.pi
W2 = 0.6 * 2 * np.pi

t  = np.linspace(1.0, 3.0, 2001)
A  = np.sin(t * W1) + np.cos(t * W2)
tt = np.linspace(t[0], t[-1], 100)

def cubic_spline_notaknot(x, y, x_new):
    x     = np.asarray(x,     dtype=float)
    y     = np.asarray(y,     dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    n = x.size

    # Вырожденные случаи
    if n == 2:
        out = np.interp(x_new, x, y)
        out[(x_new < x[0]) | (x_new > x[-1])] = np.nan
        return out
    if n == 3:
        coef = np.polyfit(x, y, 2)
        out = np.polyval(coef, x_new)
        out[(x_new < x[0]) | (x_new > x[-1])] = np.nan
        return out

    h = np.diff(x)
    d = np.diff(y) / h

    # Система A·M = rhs для вторых производных M_i
    A_   = np.zeros((n, n))
    rhs  = np.zeros(n)

    # Внутренние строки: непрерывность первой производной в x_i, i=1..n-2
    for i in range(1, n - 1):
        A_[i, i-1] = h[i-1]
        A_[i, i]   = 2.0 * (h[i-1] + h[i])
        A_[i, i+1] = h[i]
        rhs[i]     = 6.0 * (d[i] - d[i-1])

    # Строка 0: not-a-knot в x_1  →  h1·M0 − (h0+h1)·M1 + h0·M2 = 0
    A_[0, 0] =  h[1]
    A_[0, 1] = -(h[0] + h[1])
    A_[0, 2] =  h[0]
    rhs[0]   =  0.0

    # Строка n−1: not-a-knot в x_{n-2}  →  −h_{n-2}·M_{n-3} + (h_{n-3}+h_{n-2})·M_{n-2} − h_{n-3}·M_{n-1} = 0
    A_[n-1, n-3] = -h[n-2]
    A_[n-1, n-2] =  h[n-3] + h[n-2]
    A_[n-1, n-1] = -h[n-3]
    rhs[n-1]     =  0.0

    M = np.linalg.solve(A_, rhs)

    # Вычисление сплайна на x_new
    out = np.full(x_new.shape, np.nan, dtype=float)
    for k, xn in enumerate(x_new):
        if xn < x[0] or xn > x[-1]:
            continue
        i = n - 2 if xn == x[-1] else int(np.searchsorted(x, xn, side='right') - 1)
        hi = h[i]
        dx = xn - x[i]
        out[k] = (M[i]   * (x[i+1] - xn)**3 / (6.0 * hi)
                  + M[i+1] * dx**3              / (6.0 * hi)
                  + (y[i]   - M[i]   * hi*hi / 6.0) * (x[i+1] - xn) / hi
                  + (y[i+1] - M[i+1] * hi*hi / 6.0) * dx             / hi)
    return out


def interp_nearest(t_nodes, y_nodes, t_new):
    t_nodes = np.asarray(t_nodes, dtype=float)
    y_nodes = np.asarray(y_nodes, dtype=float)
    t_new   = np.asarray(t_new,   dtype=float)
    idx = np.abs(t_new[:, None] - t_nodes[None, :]).argmin(axis=1)
    return y_nodes[idx]


def vander(t_arr, degree):
    return np.vander(np.asarray(t_arr), degree + 1, increasing=True)

def ls_poly_fit(t_nodes, y_nodes, degree):
    V = vander(t_nodes, degree)
    c, *_ = np.linalg.lstsq(V, y_nodes, rcond=None)
    return c

def poly_eval(coef, t_new):
    return vander(t_new, len(coef) - 1) @ coef

fig, axes = plt.subplots(3, 2, figsize=(16, 12))

#  Пункты 4 и 5: N=6, h=0.4
N  = 6
t0 = np.linspace(1.0, 3.0, N)
A0 = np.sin(t0 * W1) + np.cos(t0 * W2)

# Пункт 4 — МНК
Yap3 = poly_eval(ls_poly_fit(t0, A0, 3), tt)
Yap5 = poly_eval(ls_poly_fit(t0, A0, 5), tt)
Yap7 = poly_eval(ls_poly_fit(t0, A0, 7), tt)
Yap9 = poly_eval(ls_poly_fit(t0, A0, 9), tt)

ax = axes[0, 0]
ax.plot(t,  A,        label='A(t)')
ax.plot(t0, A0,  'o', label='A0(t0)')
ax.plot(tt, Yap3,     label='Yap3')
ax.plot(tt, Yap5,     label='Yap5')
ax.plot(tt, Yap7,     label='Yap7')
ax.plot(tt, Yap9,     label='Yap9')
ax.legend()
ax.set_title('Пункт 4. МНК-аппроксимация (h=0.4)')
ax.grid(True)

# Пункт 5 — интерполяции
Yint0 = np.interp(tt, t0, A0)
Yint1 = cubic_spline_notaknot(t0, A0, tt)
Yint2 = interp_nearest(t0, A0, tt)
pL    = np.polyfit(t0, A0, N - 1)
Yint3 = np.polyval(pL, tt)

ax = axes[0, 1]
ax.plot(t,  A,        label='A(t)')
ax.plot(t0, A0,  'o', label='A0(t0)')
ax.plot(tt, Yint0,    label='Yint0 (лин.)')
ax.plot(tt, Yint1,    label='Yint1 (сплайн)')
ax.plot(tt, Yint2,    label='Yint2 (nearest)')
ax.plot(tt, Yint3,    label='Yint3 (полином)')
ax.legend()
ax.set_title('Пункт 5. Интерполяция (h=0.4)')
ax.grid(True)

#  Пункт 6.1: N=12
N  = 12
t0 = np.linspace(1.0, 3.0, N)
A0 = np.sin(t0 * W1) + np.cos(t0 * W2)

Yap3 = poly_eval(ls_poly_fit(t0, A0, 3), tt)
Yap5 = poly_eval(ls_poly_fit(t0, A0, 5), tt)
Yap7 = poly_eval(ls_poly_fit(t0, A0, 7), tt)
Yap9 = poly_eval(ls_poly_fit(t0, A0, 9), tt)

Yint0 = np.interp(tt, t0, A0)
Yint1 = cubic_spline_notaknot(t0, A0, tt)
Yint2 = interp_nearest(t0, A0, tt)
pL    = np.polyfit(t0, A0, N - 1)
Yint3 = np.polyval(pL, tt)

ax = axes[1, 0]
ax.plot(t,  A,        label='A(t)')
ax.plot(t0, A0,  'o', label='A0(t0)')
ax.plot(tt, Yap3,     label='Yap3')
ax.plot(tt, Yap5,     label='Yap5')
ax.plot(tt, Yap7,     label='Yap7')
ax.plot(tt, Yap9,     label='Yap9')
ax.plot(tt, Yint0,    label='Yint0')
ax.plot(tt, Yint1,    label='Yint1')
ax.plot(tt, Yint2,    label='Yint2')
ax.legend(fontsize=8)
ax.set_title('Пункт 6.1 МНК-аппроксимация (h=0.2)')
ax.grid(True)

#  Пункт 6.2: N=24
N  = 24
t0 = np.linspace(1.0, 3.0, N)
A0 = np.sin(t0 * W1) + np.cos(t0 * W2)

Yap3 = poly_eval(ls_poly_fit(t0, A0, 3), tt)
Yap5 = poly_eval(ls_poly_fit(t0, A0, 5), tt)
Yap7 = poly_eval(ls_poly_fit(t0, A0, 7), tt)
Yap9 = poly_eval(ls_poly_fit(t0, A0, 9), tt)

Yint0 = np.interp(tt, t0, A0)
Yint1 = cubic_spline_notaknot(t0, A0, tt)
Yint2 = interp_nearest(t0, A0, tt)
pL    = np.polyfit(t0, A0, N - 1)
Yint3 = np.polyval(pL, tt)

# Сохраняем узлы/полином для Пункта 7
t0_24, A0_24, pL_24 = t0, A0, pL

ax = axes[1, 1]
ax.plot(t,  A,        label='A(t)')
ax.plot(t0, A0,  'o', label='A0(t0)')
ax.plot(tt, Yap3,     label='Yap3')
ax.plot(tt, Yap5,     label='Yap5')
ax.plot(tt, Yap7,     label='Yap7')
ax.plot(tt, Yap9,     label='Yap9')
ax.plot(tt, Yint0,    label='Yint0')
ax.plot(tt, Yint1,    label='Yint1')
ax.plot(tt, Yint2,    label='Yint2')
ax.legend(fontsize=8)
ax.set_title('Пункт 6.2 МНК-аппроксимация (h=0.4)')
ax.grid(True)

#  Пункт 7: экстраполяция на [0.5, 3.5]
#  Коэффициенты МНК и pL — из блока N=24
t_ext = np.linspace(0.5, 3.5, 301)      # 0.5:0.01:3.5
A_ext = np.sin(t_ext * W1) + np.cos(t_ext * W2)

Yap3_e = poly_eval(ls_poly_fit(t0_24, A0_24, 3), t_ext)
Yap5_e = poly_eval(ls_poly_fit(t0_24, A0_24, 5), t_ext)
Yap7_e = poly_eval(ls_poly_fit(t0_24, A0_24, 7), t_ext)
Yap9_e = poly_eval(ls_poly_fit(t0_24, A0_24, 9), t_ext)

Yint1_e = cubic_spline_notaknot(t0_24, A0_24, t_ext) 
Yint3_e = np.polyval(pL_24, t_ext)

ax = axes[2, 0]
ax.plot(t_ext, A_ext,        label='A(t)')
ax.plot(t0_24, A0_24,   'o', label='A0(t0)')
ax.plot(t_ext, Yap3_e,       label='Yap3')
ax.plot(t_ext, Yap5_e,       label='Yap5')
ax.plot(t_ext, Yap7_e,       label='Yap7')
ax.plot(t_ext, Yap9_e,       label='Yap9')
ax.plot(t_ext, Yint3_e,      label='Yint3 (полином)')
ax.legend(fontsize=8)
ax.set_title('Экстраполяция на [0.5 3.5] (h=0.1)')
ax.grid(True)

axes[2, 1].axis('off')

fig.tight_layout()
plt.show()