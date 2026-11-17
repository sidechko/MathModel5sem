(function () {
  'use strict';

  // ---------- 1. Аналог np.linspace(1.0, 3.0, 2001) ----------
  const N = 2001;
  const tMin = 1.0;
  const tMax = 3.0;

  const t = new Array(N);
  for (let i = 0; i < N; i++) {
    t[i] = tMin + (tMax - tMin) * i / (N - 1);
  }

  // ---------- 2. Параметры ----------
  const W1 = 1.2 * 2 * Math.PI;
  const W2 = 0.6 * 2 * Math.PI;
  const W0 = 20.0 * Math.PI;
  const K  = 10.0;

  // ---------- 3. Сигналы ----------
  const A    = new Array(N);
  const negA = new Array(N);
  const Famp = new Array(N);
  const Ffm  = new Array(N);

  for (let i = 0; i < N; i++) {
    const a = Math.sin(t[i] * W1) + Math.cos(t[i] * W2);
    A[i]    = a;
    negA[i] = -a;
    Famp[i] = a * Math.cos(t[i] * W0);
    Ffm[i]  = Math.cos((W0 + K * a) * t[i]);
  }

  // ---------- 4. Трассы ----------
  // legend: 'legend' | 'legend2' | 'legend3' — своя легенда на каждый сабплот
  const trace = (y, name, xaxis, yaxis, legend) => ({
    x: t,
    y: y,
    name: name,
    type: 'scatter',
    mode: 'lines',
    line: { width: 1 },
    xaxis: xaxis,
    yaxis: yaxis,
    legend: legend
  });

  const data = [
    // Подграфик 1: A, -A, Famp
    trace(A,    'A',     'x',  'y',  'legend'),
    trace(negA, '-A',    'x',  'y',  'legend'),
    trace(Famp, 'F AM',  'x',  'y',  'legend'),

    // Подграфик 2: A, Ffm
    trace(A,    'A',     'x2', 'y2', 'legend2'),
    trace(Ffm,  'F FM',  'x2', 'y2', 'legend2'),

    // Подграфик 3: A, -A, Famp, Ffm
    trace(A,    'A',     'x3', 'y3', 'legend3'),
    trace(negA, '-A',    'x3', 'y3', 'legend3'),
    trace(Famp, 'F AM',  'x3', 'y3', 'legend3'),
    trace(Ffm,  'F FM',  'x3', 'y3', 'legend3')
  ];

  // ---------- 5. Раскладка (3 сабплота друг под другом) ----------
  const legendStyle = (yPos) => ({
    x: 0.99, xanchor: 'right',
    y: yPos, yanchor: 'top',
    bgcolor: 'rgba(255,255,255,0.85)',
    bordercolor: '#bbbbbb', borderwidth: 1,
    font: { size: 11 }
  });

  const layout = {
    margin: { l: 60, r: 30, t: 55, b: 55 },
    showlegend: true,

    legend:  legendStyle(0.965),
    legend2: legendStyle(0.615),
    legend3: legendStyle(0.265),

    // по своей оси X на каждый сабплот (как в matplotlib)
    xaxis:  { anchor: 'y',  domain: [0, 1], range: [tMin, tMax] },
    xaxis2: { anchor: 'y2', domain: [0, 1], range: [tMin, tMax] },
    xaxis3: { anchor: 'y3', domain: [0, 1], range: [tMin, tMax] },

    yaxis:  { domain: [0.720, 0.965] },
    yaxis2: { domain: [0.370, 0.615] },
    yaxis3: { domain: [0.020, 0.265] },

    // заголовки сабплотов (аналог ax.set_title)
    annotations: [
      { text: 'F amp', xref: 'paper', yref: 'paper', x: 0.5, y: 1.00, showarrow: false, font: { size: 14 } },
      { text: 'F fm',  xref: 'paper', yref: 'paper', x: 0.5, y: 0.65, showarrow: false, font: { size: 14 } },
      { text: 'Mixed', xref: 'paper', yref: 'paper', x: 0.5, y: 0.30, showarrow: false, font: { size: 14 } }
    ]
  };

  // ---------- 6. Отрисовка ----------
  Plotly.newPlot('plot', data, layout, {
    responsive: true,
    displaylogo: false
  });

  // аналог plt.close('all') при повторном запуске — просто создаём заново
  window.addEventListener('resize', () => Plotly.Plots.resize('plot'));
})();