clear()
W1 = 1.2*2*%pi;
W2 = 0.6*2*%pi;

function y = A1(t)
    y = sin(t*W1)
endfunction

function y = A2(t)
    y = cos(t*W2) + 2
endfunction

// ПУНКТ 1. Графики функций A1(t) и A2(t)
t = 1.0:0.001:3.0;

scf(1); clf()
plot(t, A1(t), 'b-', t, A2(t), 'r-')
legend('A1(t) = sin(W1*t)', 'A2(t) = cos(W2*t) + 2')
xtitle('Пункт 1. Функции A1(t) и A2(t)')
xgrid()

// ПУНКТ 2. Сетка с шагом hp = 0.4, узлы — в серединах отрезков
hp = 0.4;
ap = 1.0 + hp/2; 
bp = 3.0 - hp/2; 
tp = ap:hp:bp;

Ap1 = A1(tp);
Ap2 = A2(tp);

printf('Пункт 2. Сетка (hp = %.1f):\n', hp)
printf('tp  = '); disp(tp)
printf('Ap1 = '); disp(Ap1)
printf('Ap2 = '); disp(Ap2)

// ПУНКТ 3. Гистограмма (Ap2-Ap1) + график A2(t)-A1(t)
scf(2); clf()
bar(tp, Ap2 - Ap1, 1.0)
plot2d(t, A2(t) - A1(t), style=[color("red")])
xgrid()
xtitle('Пункт 3. Гистограмма Ap2-Ap1 и график A2(t)-A1(t)')

// ПУНКТ 4. Метод средних прямоугольников
printf('\n=== ПУНКТ 4. Метод средних прямоугольников ===\n')
printf('%-4s %-10s %-14s %-14s\n', '№', 'hp', 'S', 'Delta')

hp_vals    = [];
S_vals     = [];
Delta_vals = [];
hp_cur     = 0.4;
S_prev     = %nan;
iter_num   = 1;
conv4      = %f;

while iter_num <= 20
    ap = 1.0 + hp_cur/2;
    bp = 3.0 - hp_cur/2;
    tp_cur = ap:hp_cur:bp;
    Ap1_cur = A1(tp_cur);
    Ap2_cur = A2(tp_cur);
    S_cur = hp_cur * sum(Ap2_cur - Ap1_cur);

    if iter_num == 1
        Delta_cur = %nan;
        printf('%-4d %-10.5f %-14.6f %-14s\n', iter_num, hp_cur, S_cur, '--')
    else
        Delta_cur = abs(S_cur - S_prev);
        printf('%-4d %-10.5f %-14.6f %-14.6f\n', iter_num, hp_cur, S_cur, Delta_cur)
    end

    hp_vals    = [hp_vals, hp_cur];
    S_vals     = [S_vals, S_cur];
    Delta_vals = [Delta_vals, Delta_cur];

    if iter_num > 1 && Delta_cur < 0.001
        conv4 = %t;
        printf('--> Сходимость достигнута: Delta = %.6f < 0.001\n', Delta_cur)
        break
    end

    S_prev  = S_cur;
    hp_cur  = hp_cur / 2;
    iter_num = iter_num + 1;
end

S_rect = S_vals($);

// Кривая сходимости в логарифмическом масштабе
scf(3); clf()
idx = 2:length(Delta_vals);
y_vals = log10(Delta_vals(2:$));
plot2d(idx, y_vals, style=1)  // соединяем точки линией
plot2d(idx, y_vals, style=-1) // отображаем сами точки
xlabel('Номер сетки i')
ylabel('log_{10}(\Delta)')
xtitle('Пункт 4. Кривая сходимости (лог. масштаб)')
xgrid()

// ПУНКТ 5. Метод Монте-Карло
printf('\n=== ПУНКТ 5. Метод Монте-Карло ===\n')

N0       = 100;
max_iters = 12;
N_vals   = [];
S_means  = [];
S_ranges = [];
N        = N0;
conv5    = %f;

for iter = 1:max_iters
    S_exps = zeros(1,5);
    for j = 1:5
        xm = 2*rand(1,N) + 1.0;   // [1, 3]
        ym = 4*rand(1,N) - 1.0;   // [-1, 3]
        M  = sum(ym >= A1(xm) & ym <= A2(xm));
        S_exps(j) = 8 * M / N;    // площадь прямоугольника = 2*4 = 8
    end

    S_mean  = mean(S_exps);
    S_range = max(S_exps) - min(S_exps);

    N_vals   = [N_vals, N];
    S_means  = [S_means, S_mean];
    S_ranges = [S_ranges, S_range];

    printf('N = %-8d  S_mean = %-12.6f  Delta = %-12.6f\n', N, S_mean, S_range)

    if S_range < 0.001 && iter > 1
        conv5 = %t;
        printf('--> Сходимость достигнута при N = %d\n', N)
        break
    end
    N = 2*N;
end

if ~conv5
    printf('--> Внимание: Delta < 0.001 не достигнуто за %d итераций.\n', max_iters)
    printf('    Для чистого Монте-Карло нужно N ~ 10^8 (нецелесообразно).\n')
end

S_MC = S_means($);

// Визуализация: 100 точек в прямоугольнике + линии A1, A2
scf(4); clf()
N_vis = 1000;
//xm_v = 2*rand(1,N_vis) + 1.0;
ym_v = 4*rand(1,N_vis) - 1.0;
M_v  = sum(ym_v >= A1(xm_v) & ym_v <= A2(xm_v));

plot(t, A1(t), 'b-', t, A2(t), 'r-')
plot(xm_v, ym_v, '.k')
xtitle(sprintf('Пункт 5. Монте-Карло: N=%d, M=%d, S=%.3f', N_vis, M_v, 8*M_v/N_vis))
xgrid()

// Кривые сходимости Монте-Карло
scf(5); clf()
subplot(2,1,1)
plot(1:length(S_means), S_means, 'bo-')
xlabel('Номер серии i')
ylabel('S_{mean}')
xtitle('Среднее значение площади')
xgrid()

subplot(2,1,2)
plot(1:length(S_ranges), log10(S_ranges), 'ro-')
xlabel('Номер серии i')
ylabel('log_{10}(\Delta)')
xtitle('Разброс \Delta (лог. масштаб)')
xgrid()

// ПУНКТ 6. Аналитический метод (формула Ньютона-Лейбница)
printf('\n=== ПУНКТ 6. Аналитический метод ===\n')

// S = ∫_1^3 (A2 - A1) dt = ∫_1^3 (cos(W2*t) + 2 - sin(W1*t)) dt
//   = [ sin(W2*t)/W2 + 2*t + cos(W1*t)/W1 ]_1^3
function y = F_antider(t)
    y = sin(t*W2)/W2 + 2*t + cos(t*W1)/W1
endfunction

S_analytic = F_antider(3) - F_antider(1);

printf('Первообразная: F(t) = sin(W2*t)/W2 + 2*t + cos(W1*t)/W1\n')
printf('S = F(3) - F(1) = %.6f\n', S_analytic)

// ТАБЛИЦА 3. Итоговые результаты
printf('\n=== Таблица 3. Итоговые результаты ===\n')
printf('%-28s %-14s\n', 'Метод', 'Площадь S')
printf('%-28s %-14.5f\n', 'Аналитический',            S_analytic)
printf('%-28s %-14.5f\n', 'Средние прямоугольники',   S_rect)
printf('%-28s %-14.5f\n', 'Монте-Карло',              S_MC)
