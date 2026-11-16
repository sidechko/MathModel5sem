"""
Эквивалент LR_3_scilab.sci на Python с графикой matplotlib.
Без scipy.

Пункты:
  1 — графики A1(t), A2(t)
  2 — сетка hp=0.4, узлы в серединах отрезков
  3 — гистограмма Ap2−Ap1 + график A2(t)−A1(t)
  4 — метод средних прямоугольников, сходимость
  5 — Монте-Карло
  6 — аналитический метод
  Таблица 3 — итоговые результаты

Зависимости: pip install numpy matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# ============================================================
#  Функции
# ============================================================
W1 = 1.2 * 2 * np.pi
W2 = 0.6 * 2 * np.pi

def A1(t):
    return np.sin(t * W1)

def A2(t):
    return np.cos(t * W2) + 2.0

def F_antider(t):
    return np.sin(t * W2) / W2 + 2.0 * t + np.cos(t * W1) / W1

# ============================================================
#  ПУНКТ 1. Графики A1(t) и A2(t)
# ============================================================
t = np.linspace(1.0, 3.0, 2001)          # t = 1.0:0.001:3.0

fig1 = plt.figure(figsize=(10, 6))
plt.plot(t, A1(t), 'b-', label='A1(t) = sin(W1*t)')
plt.plot(t, A2(t), 'r-', label='A2(t) = cos(W2*t) + 2')
plt.legend()
plt.title('Пункт 1. Функции A1(t) и A2(t)')
plt.grid(True)

# ============================================================
#  ПУНКТ 2. Сетка hp=0.4, узлы в серединах отрезков
# ============================================================
hp = 0.4
ap = 1.0 + hp / 2.0                      # 1.2
bp = 3.0 - hp / 2.0                      # 2.8
n_pts = int(round((bp - ap) / hp)) + 1   # 5 точек
tp = np.linspace(ap, bp, n_pts)

Ap1 = A1(tp)
Ap2 = A2(tp)

print(f'Пункт 2. Сетка (hp = {hp:.1f}):')
print(f'tp  = {tp}')
print(f'Ap1 = {Ap1}')
print(f'Ap2 = {Ap2}')

# ============================================================
#  ПУНКТ 3. Гистограмма Ap2−Ap1 + график A2(t)−A1(t)
# ============================================================
fig2 = plt.figure(figsize=(10, 6))
plt.bar(tp, Ap2 - Ap1, width=hp, edgecolor='black', color='#1f77b4', alpha=0.7)
plt.plot(t, A2(t) - A1(t), 'r-')
plt.grid(True)
plt.title('Пункт 3. Гистограмма Ap2-Ap1 и график A2(t)-A1(t)')

# ============================================================
#  ПУНКТ 4. Метод средних прямоугольников
# ============================================================
print('\n=== ПУНКТ 4. Метод средних прямоугольников ===')
print(f'{"№":<4} {"hp":<10} {"S":<14} {"Delta":<14}')

hp_vals    = []
S_vals     = []
Delta_vals = []

hp_cur   = 0.4
S_prev   = np.nan
iter_num = 1
conv4    = False

while iter_num <= 20:
    ap_cur = 1.0 + hp_cur / 2.0
    bp_cur = 3.0 - hp_cur / 2.0
    n = int(round((bp_cur - ap_cur) / hp_cur)) + 1
    tp_cur = np.linspace(ap_cur, bp_cur, n)

    Ap1_cur = A1(tp_cur)
    Ap2_cur = A2(tp_cur)
    S_cur   = hp_cur * np.sum(Ap2_cur - Ap1_cur)

    if iter_num == 1:
        Delta_cur = np.nan
        print(f'{iter_num:<4d} {hp_cur:<10.5f} {S_cur:<14.6f} {"--":<14}')
    else:
        Delta_cur = abs(S_cur - S_prev)
        print(f'{iter_num:<4d} {hp_cur:<10.5f} {S_cur:<14.6f} {Delta_cur:<14.6f}')

    hp_vals.append(hp_cur)
    S_vals.append(S_cur)
    Delta_vals.append(Delta_cur)

    if iter_num > 1 and Delta_cur < 0.001:
        conv4 = True
        print(f'--> Сходимость достигнута: Delta = {Delta_cur:.6f} < 0.001')
        break

    S_prev   = S_cur
    hp_cur   = hp_cur / 2.0
    iter_num = iter_num + 1

S_rect = S_vals[-1]

# Кривая сходимости в логарифмическом масштабе
fig3 = plt.figure(figsize=(10, 6))
idx    = np.arange(2, len(Delta_vals) + 1)      # 2:length(Delta_vals)
y_vals = np.log10(np.array(Delta_vals[1:]))     # log10(Delta_vals(2:$))
plt.plot(idx, y_vals, 'b-')
plt.plot(idx, y_vals, 'bo')
plt.xlabel('Номер сетки i')
plt.ylabel('log10(Delta)')
plt.title('Пункт 4. Кривая сходимости (лог. масштаб)')
plt.grid(True)

# ============================================================
#  ПУНКТ 5. Метод Монте-Карло
#  random-числа зависят от RNG, поэтому конкретные значения
#  не совпадают со Scilab'ом; для воспроизводимости
#  результатов между запусками Python фиксирован seed.
# ============================================================
print('\n=== ПУНКТ 5. Метод Монте-Карло ===')

rng = np.random.default_rng(42)          # seed для воспроизводимости

N0        = 100
max_iters = 12
N_vals    = []
S_means   = []
S_ranges  = []
N         = N0
conv5     = False

for it in range(1, max_iters + 1):
    S_exps = np.zeros(5)
    for j in range(5):
        xm = 2.0 * rng.random(N) + 1.0    # [1, 3]
        ym = 4.0 * rng.random(N) - 1.0    # [-1, 3]
        M = np.sum((ym >= A1(xm)) & (ym <= A2(xm)))
        S_exps[j] = 8.0 * M / N           # площадь прямоугольника = 2·4 = 8

    S_mean  = np.mean(S_exps)
    S_range = np.max(S_exps) - np.min(S_exps)

    N_vals.append(N)
    S_means.append(S_mean)
    S_ranges.append(S_range)

    print(f'N = {N:<8d}  S_mean = {S_mean:<12.6f}  Delta = {S_range:<12.6f}')

    if S_range < 0.001 and it > 1:
        conv5 = True
        print(f'--> Сходимость достигнута при N = {N}')
        break
    N = 2 * N

if not conv5:
    print(f'--> Внимание: Delta < 0.001 не достигнуто за {max_iters} итераций.')
    print('    Для чистого Монте-Карло нужно N ~ 10^8 (нецелесообразно).')

S_MC = S_means[-1]

# Визуализация: 1000 точек в прямоугольнике + линии A1, A2
fig4 = plt.figure(figsize=(10, 6))
N_vis = 1000
xm_v  = 2.0 * rng.random(N_vis) + 1.0     # (в оригинале эта строка была закомментирована)
ym_v  = 4.0 * rng.random(N_vis) - 1.0
M_v   = np.sum((ym_v >= A1(xm_v)) & (ym_v <= A2(xm_v)))

plt.plot(t, A1(t), 'b-')
plt.plot(t, A2(t), 'r-')
plt.plot(xm_v, ym_v, '.k', markersize=2)
plt.title(f'Пункт 5. Монте-Карло: N={N_vis}, M={M_v}, S={8.0*M_v/N_vis:.3f}')
plt.grid(True)

# Кривые сходимости Монте-Карло (2 подграфика)
fig5 = plt.figure(figsize=(10, 8))

ax1 = fig5.add_subplot(2, 1, 1)
ax1.plot(np.arange(1, len(S_means) + 1), S_means, 'bo-')
ax1.set_xlabel('Номер серии i')
ax1.set_ylabel('S_mean')
ax1.set_title('Среднее значение площади')
ax1.grid(True)

ax2 = fig5.add_subplot(2, 1, 2)
log_ranges = np.log10(np.maximum(S_ranges, 1e-300))   # защита от log10(0)
ax2.plot(np.arange(1, len(S_ranges) + 1), log_ranges, 'ro-')
ax2.set_xlabel('Номер серии i')
ax2.set_ylabel('log10(Delta)')
ax2.set_title('Разброс Delta (лог. масштаб)')
ax2.grid(True)

# ============================================================
#  ПУНКТ 6. Аналитический метод (формула Ньютона–Лейбница)
# ============================================================
print('\n=== ПУНКТ 6. Аналитический метод ===')

S_analytic = F_antider(3.0) - F_antider(1.0)

print('Первообразная: F(t) = sin(W2*t)/W2 + 2*t + cos(W1*t)/W1')
print(f'S = F(3) - F(1) = {S_analytic:.6f}')

# ============================================================
#  ТАБЛИЦА 3. Итоговые результаты
# ============================================================
print('\n=== Таблица 3. Итоговые результаты ===')
print(f'{"Метод":<28} {"Площадь S":<14}')
print(f'{"Аналитический":<28} {S_analytic:<14.5f}')
print(f'{"Средние прямоугольники":<28} {S_rect:<14.5f}')
print(f'{"Монте-Карло":<28} {S_MC:<14.5f}')

plt.show()