x = [0.1, 0.2, 0.3, 0.5, 0.8, 1.1, 1.9, 2.2, 2.5];
y = [3.0, 2.0, 1.5, 2.0, 2.5, 1.9, 1.0, 1.2, 1.6];
xx = linspace(x(1), x($), 100); // равномерная сетка из 100 элементов

// линейная интерполяция
Ylin = interp1(x, y, xx, 'linear');

// сплайн-интерполяция
spline_coeff = splin(x, y);
Yspline = interp(xx, x, y, spline_coeff);

// ступенчатая интерполяция (ближайший сосед)
Ynearest = interp1(x, y, xx, 'nearest');

// строим график
clf();
plot(x, y, 'o', xx, Ylin, xx, Yspline, xx, Ynearest);
xgrid();
legend('Test', 'Linear', 'Spline', 'Nearest');

