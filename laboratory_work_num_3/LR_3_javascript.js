(function () {
  'use strict';

  // ==================== УТИЛИТЫ ====================
  const logEl = document.getElementById('log');
  function log(s) { logEl.textContent += s + '\n'; }

  function linspace(a, b, n) {
    const out = new Array(n);
    if (n === 1) { out[0] = a; return out; }
    const step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = a + step * i;
    out[n - 1] = b;
    return out;
  }

  // Детерминированный PRNG (mulberry32) — воспроизводимый, но не PCG64.
  function mulberry32(seed) {
    let s = seed >>> 0;
    return function () {
      s = (s + 0x6D2B79F5) >>> 0;
      let t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ==================== ФУНКЦИИ ====================
  const W1 = 1.2 * 2 * Math.PI;
  const W2 = 0.6 * 2 * Math.PI;

  const A1 = (t) => Math.sin(t * W1);
  const A2 = (t) => Math.cos(t * W2) + 2.0;
  const F_antider = (t) => Math.sin(t * W2) / W2 + 2.0 * t + Math.cos(t * W1) / W1;

  // ==================== ПУНКТ 1 ====================
  const t = linspace(1.0, 3.0, 2001);

  Plotly.newPlot('fig1', [
    { x: t, y: t.map(A1), name: 'A1(t) = sin(W1*t)',
      type: 'scatter', mode: 'lines', line: { color: 'blue' } },
    { x: t, y: t.map(A2), name: 'A2(t) = cos(W2*t) + 2',
      type: 'scatter', mode: 'lines', line: { color: 'red' } }
  ], {
    title: 'Пункт 1. Функции A1(t) и A2(t)',
    xaxis: { gridcolor: '#dcdcdc', zeroline: false },
    yaxis: { gridcolor: '#dcdcdc', zeroline: false },
    legend: { bgcolor: 'rgba(255,255,255,0.85)', bordercolor: '#bbb', borderwidth: 1 },
    margin: { l: 60, r: 30, t: 50, b: 45 }
  });

  // ==================== ПУНКТ 2 ====================
  const hp0 = 0.4;
  const ap0 = 1.0 + hp0 / 2.0;
  const bp0 = 3.0 - hp0 / 2.0;
  const nPts0 = Math.round((bp0 - ap0) / hp0) + 1;
  const tp0 = linspace(ap0, bp0, nPts0);
  const Ap1_0 = tp0.map(A1);
  const Ap2_0 = tp0.map(A2);

  log(`Пункт 2. Сетка (hp = ${hp0.toFixed(1)}):`);
  log(`tp  = [${tp0.map(v => v.toFixed(4)).join(', ')}]`);
  log(`Ap1 = [${Ap1_0.map(v => v.toFixed(4)).join(', ')}]`);
  log(`Ap2 = [${Ap2_0.map(v => v.toFixed(4)).join(', ')}]`);
  log('');

  // ==================== ПУНКТ 3 ====================
  const diffP = tp0.map((_, i) => Ap2_0[i] - Ap1_0[i]);
  const diffCurve = t.map(ti => A2(ti) - A1(ti));

  Plotly.newPlot('fig3', [
    { x: tp0, y: diffP, name: 'Ap2-Ap1', type: 'bar',
      width: hp0,
      marker: { color: '#1f77b4', line: { color: 'black', width: 1 } },
      opacity: 0.7 },
    { x: t, y: diffCurve, name: 'A2(t)-A1(t)',
      type: 'scatter', mode: 'lines', line: { color: 'red' } }
  ], {
    title: 'Пункт 3. Гистограмма Ap2-Ap1 и график A2(t)-A1(t)',
    xaxis: { gridcolor: '#dcdcdc', zeroline: false },
    yaxis: { gridcolor: '#dcdcdc', zeroline: false },
    legend: { bgcolor: 'rgba(255,255,255,0.85)', bordercolor: '#bbb', borderwidth: 1 },
    margin: { l: 60, r: 30, t: 50, b: 45 }
  });

  // ==================== ПУНКТ 4 ====================
  log('=== ПУНКТ 4. Метод средних прямоугольников ===');
  log('№    hp         S              Delta');

  const hp_vals = [], S_vals = [], Delta_vals = [];
  let hp_cur = 0.4;
  let S_prev = NaN;
  let iter_num = 1;
  let conv4 = false;

  while (iter_num <= 20) {
    const ap_cur = 1.0 + hp_cur / 2.0;
    const bp_cur = 3.0 - hp_cur / 2.0;
    const n = Math.round((bp_cur - ap_cur) / hp_cur) + 1;
    const tp_cur = linspace(ap_cur, bp_cur, n);

    let sum = 0;
    for (let k = 0; k < n; k++) {
      const tk = tp_cur[k];
      sum += A2(tk) - A1(tk);
    }
    const S_cur = hp_cur * sum;

    let Delta_cur;
    if (iter_num === 1) {
      Delta_cur = NaN;
      log(`${String(iter_num).padEnd(4)} ${hp_cur.toFixed(5).padEnd(10)} ${S_cur.toFixed(6).padEnd(14)} --`);
    } else {
      Delta_cur = Math.abs(S_cur - S_prev);
      log(`${String(iter_num).padEnd(4)} ${hp_cur.toFixed(5).padEnd(10)} ${S_cur.toFixed(6).padEnd(14)} ${Delta_cur.toFixed(6)}`);
    }

    hp_vals.push(hp_cur);
    S_vals.push(S_cur);
    Delta_vals.push(Delta_cur);

    if (iter_num > 1 && Delta_cur < 0.001) {
      conv4 = true;
      log(`--> Сходимость достигнута: Delta = ${Delta_cur.toFixed(6)} < 0.001`);
      break;
    }

    S_prev = S_cur;
    hp_cur = hp_cur / 2.0;
    iter_num++;
  }
  const S_rect = S_vals[S_vals.length - 1];
  log('');

  const idx4 = [], y4 = [];
  for (let i = 1; i < Delta_vals.length; i++) {
    idx4.push(i + 1);
    y4.push(Math.log10(Delta_vals[i]));
  }

  Plotly.newPlot('fig4', [
    { x: idx4, y: y4, name: 'log10(Delta)',
      type: 'scatter', mode: 'lines+markers',
      line: { color: 'blue' }, marker: { color: 'blue', size: 7 } }
  ], {
    title: 'Пункт 4. Кривая сходимости (лог. масштаб)',
    xaxis: { title: 'Номер сетки i', gridcolor: '#dcdcdc', zeroline: false },
    yaxis: { title: 'log10(Delta)', gridcolor: '#dcdcdc', zeroline: false },
    margin: { l: 60, r: 30, t: 50, b: 50 }
  });

  // ==================== ПУНКТ 5 ====================
  log('=== ПУНКТ 5. Метод Монте-Карло ===');

  const rand = mulberry32(42);   // детерминированный seed
  const N0 = 100;
  const max_iters = 12;
  const N_vals = [], S_means = [], S_ranges = [];

  let N = N0;
  let conv5 = false;

  for (let it = 1; it <= max_iters; it++) {
    const S_exps = new Array(5);
    for (let j = 0; j < 5; j++) {
      let M = 0;
      for (let k = 0; k < N; k++) {
        const xm = 2.0 * rand() + 1.0;   // [1, 3]
        const ym = 4.0 * rand() - 1.0;   // [-1, 3]
        if (ym >= A1(xm) && ym <= A2(xm)) M++;
      }
      S_exps[j] = 8.0 * M / N;           // площадь прямоугольника = 2·4 = 8
    }

    let sum = 0, mn = Infinity, mx = -Infinity;
    for (const v of S_exps) { sum += v; if (v < mn) mn = v; if (v > mx) mx = v; }
    const S_mean = sum / S_exps.length;
    const S_range = mx - mn;

    N_vals.push(N);
    S_means.push(S_mean);
    S_ranges.push(S_range);

    log(`N = ${String(N).padEnd(8)}  S_mean = ${S_mean.toFixed(6).padEnd(12)}  Delta = ${S_range.toFixed(6)}`);

    if (S_range < 0.001 && it > 1) {
      conv5 = true;
      log(`--> Сходимость достигнута при N = ${N}`);
      break;
    }
    N *= 2;
  }

  if (!conv5) {
    log(`--> Внимание: Delta < 0.001 не достигнуто за ${max_iters} итераций.`);
    log('    Для чистого Монте-Карло нужно N ~ 10^8 (нецелесообразно).');
  }
  const S_MC = S_means[S_means.length - 1];
  log('');

  // --- Визуализация Монте-Карло (1000 точек) ---
  const N_vis = 1000;
  const xm_v = new Array(N_vis);
  const ym_v = new Array(N_vis);
  let M_v = 0;
  for (let k = 0; k < N_vis; k++) {
    const x = 2.0 * rand() + 1.0;
    const y = 4.0 * rand() - 1.0;
    xm_v[k] = x; ym_v[k] = y;
    if (y >= A1(x) && y <= A2(x)) M_v++;
  }

  Plotly.newPlot('fig5', [
    { x: t, y: t.map(A1), name: 'A1(t)', type: 'scatter', mode: 'lines',
      line: { color: 'blue' } },
    { x: t, y: t.map(A2), name: 'A2(t)', type: 'scatter', mode: 'lines',
      line: { color: 'red' } },
    { x: xm_v, y: ym_v, name: 'точки', type: 'scattergl', mode: 'markers',
      marker: { color: 'black', size: 3 } }
  ], {
    title: `Пункт 5. Монте-Карло: N=${N_vis}, M=${M_v}, S=${(8.0 * M_v / N_vis).toFixed(3)}`,
    xaxis: { gridcolor: '#dcdcdc', zeroline: false },
    yaxis: { gridcolor: '#dcdcdc', zeroline: false },
    legend: { bgcolor: 'rgba(255,255,255,0.85)', bordercolor: '#bbb', borderwidth: 1 },
    margin: { l: 60, r: 30, t: 50, b: 45 }
  });

  // --- Кривые сходимости МК (2 подграфика) ---
  const itNums = N_vals.map((_, i) => i + 1);
  const logRanges = S_ranges.map(v => Math.log10(Math.max(v, 1e-300)));

  Plotly.newPlot('fig6', [
    { x: itNums, y: S_means, name: 'S_mean',
      type: 'scatter', mode: 'lines+markers',
      line: { color: 'blue' }, marker: { color: 'blue', size: 7 },
      xaxis: 'x', yaxis: 'y' },
    { x: itNums, y: logRanges, name: 'log10(Delta)',
      type: 'scatter', mode: 'lines+markers',
      line: { color: 'red' }, marker: { color: 'red', size: 7 },
      xaxis: 'x2', yaxis: 'y2' }
  ], {
    xaxis:  { anchor: 'y',  domain: [0, 1], title: 'Номер серии i',
              gridcolor: '#dcdcdc', zeroline: false },
    yaxis:  { anchor: 'x',  domain: [0.57, 1.0], title: 'S_mean',
              gridcolor: '#dcdcdc', zeroline: false },
    xaxis2: { anchor: 'y2', domain: [0, 1], title: 'Номер серии i',
              gridcolor: '#dcdcdc', zeroline: false },
    yaxis2: { anchor: 'x2', domain: [0.0, 0.43], title: 'log10(Delta)',
              gridcolor: '#dcdcdc', zeroline: false },
    annotations: [
      { text: 'Среднее значение площади',
        xref: 'paper', yref: 'paper', x: 0.5, y: 1.0,
        showarrow: false, font: { size: 13 } },
      { text: 'Разброс Delta (лог. масштаб)',
        xref: 'paper', yref: 'paper', x: 0.5, y: 0.43,
        showarrow: false, font: { size: 13 } }
    ],
    margin: { l: 65, r: 30, t: 60, b: 50 }
  });

  // ==================== ПУНКТ 6 ====================
  log('=== ПУНКТ 6. Аналитический метод ===');
  const S_analytic = F_antider(3.0) - F_antider(1.0);
  log('Первообразная: F(t) = sin(W2*t)/W2 + 2*t + cos(W1*t)/W1');
  log(`S = F(3) - F(1) = ${S_analytic.toFixed(6)}`);
  log('');

  // ==================== ТАБЛИЦА 3 ====================
  log('=== Таблица 3. Итоговые результаты ===');
  log('Метод                        Площадь S');
  log(`Аналитический                ${S_analytic.toFixed(5)}`);
  log(`Средние прямоугольники       ${S_rect.toFixed(5)}`);
  log(`Монте-Карло                  ${S_MC.toFixed(5)}`);

  // Адаптивный resize
  window.addEventListener('resize', () => {
    ['fig1', 'fig3', 'fig4', 'fig5', 'fig6'].forEach(id => Plotly.Plots.resize(id));
  });
})();