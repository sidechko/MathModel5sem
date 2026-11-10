clear()
x = [0.1, 0.2, 0.3, 0.5, 0.8, 1.1, 1.9, 2.2, 2.5];
y = [3.0, 2.0, 1.5, 2.0, 2.5, 1.9, 1.0, 1.2, 1.6];

// вычисление коэффициентов разных порядков
p3 = polyfit(x, y, 3);
p5 = polyfit(x, y, 5);
p7 = polyfit(x, y, 7);

xx = linspace(x(1), x($), 100); // равномерная сетка из 100 элементов

// расчёт ординат аппроксимирующих полиномов
Yap3 = horner(p3, xx);
Yap5 = horner(p5, xx);
Yap7 = horner(p7, xx);

// строим график
clf();
plot(x, y, 'o', xx, Yap3, xx, Yap5, xx, Yap7);
xgrid();
legend('Test', 'Yap3', 'Yap5', 'Yap7');
