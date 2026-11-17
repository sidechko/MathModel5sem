(function () {
  'use strict';

  // ==================== ЧИСЛЕННЫЕ УТИЛИТЫ ====================

  function linspace(a, b, n) {
    const arr = new Array(n);
    if (n === 1) { arr[0] = a; return arr; }
    const step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) arr[i] = a + step * i;
    arr[n - 1] = b;
    return arr;
  }

  function transpose(M) {
    const m = M.length, n = M[0].length;
    const T = Array.from({ length: n }, () => new Array(m));
    for (let i = 0; i < m; i++)
      for (let j = 0; j < n; j++) T[j][i] = M[i][j];
    return T;
  }

  function matMul(A, B) {
    const m = A.length, n = B[0].length, p = B.length;
    const C = Array.from({ length: m }, () => new Array(n).fill(0));
    for (let i = 0; i < m; i++)
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k < p; k++) s += A[i][k] * B[k][j];
        C[i][j] = s;
      }
    return C;
  }

  function matVec(A, v) {
    const m = A.length, n = A[0].length;
    const res = new Array(m);
    for (let i = 0; i < m; i++) {
      let s = 0;
      for (let j = 0; j < n; j++) s += A[i][j] * v[j];
      res[i] = s;
    }
    return res;
  }

  // Метод Гаусса с частичным выбором ведущего элемента
  function solveLinearSystem(A, b) {
    const n = A.length;
    const M = A.map(row => row.slice());
    const v = b.slice();
    for (let col = 0; col < n; col++) {
      let maxRow = col, maxVal = Math.abs(M[col][col]);
      for (let r = col + 1; r < n; r++)
        if (Math.abs(M[r][col]) > maxVal) { maxVal = Math.abs(M[r][col]); maxRow = r; }
      if (maxRow !== col) {
        [M[col], M[maxRow]] = [M[maxRow], M[col]];
        [v[col], v[maxRow]] = [v[maxRow], v[col]];
      }
      for (let r = col + 1; r < n; r++) {
        const f = M[r][col] / M[col][col];
        if (f === 0) continue;
        for (let c = col; c < n; c++) M[r][c] -= f * M[col][c];
        v[r] -= f * v[col];
      }
    }
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let s = v[i];
      for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / M[i][i];
    }
    return x;
  }

  // Vandermonde (возрастающие степени): V[i][j] = t[i]^j
  function vander(t, degree) {
    const m = t.length, n = degree + 1;
    const V = Array.from({ length: m }, () => new Array(n));
    for (let i = 0; i < m; i++) {
      let p = 1;
      for (let j = 0; j < n; j++) { V[i][j] = p; p *= t[i]; }
    }
    return V;
  }

  // Возвращает объект: {lagrange: {x, y}} или {coef: [...]}
  function lsPolyFit(tNodes, yNodes, degree) {
    const m = tNodes.length, n = degree + 1;

    if (m === n) {
      // Квадратная система = интерполяция. Лагранж устойчивее.
      return { lagrange: { x: tNodes.slice(), y: yNodes.slice() } };
    }

    const V = vander(tNodes, degree);
    const Vt = transpose(V);

    if (m > n) {
      // Переопределённая: (V^T V) c = V^T y
      const VtV = matMul(Vt, V);
      const Vty = matVec(Vt, yNodes);
      return { coef: solveLinearSystem(VtV, Vty) };
    } else {
      // Недоопределённая: c = V^T (V V^T)^{-1} y  (minimum-norm, как np.lstsq)
      const VVt = matMul(V, Vt);
      const u = solveLinearSystem(VVt, yNodes);
      return { coef: matVec(Vt, u) };
    }
  }

  // Барицентрическая формула Лагранжа (устойчива для больших n)
  function makeLagrangeEval(xNodes, yNodes) {
    const n = xNodes.length;
    const w = new Array(n).fill(1);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        if (j !== i) w[i] /= (xNodes[i] - xNodes[j]);

    return function (x) {
      let num = 0, den = 0;
      for (let i = 0; i < n; i++) {
        const diff = x - xNodes[i];
        if (Math.abs(diff) < 1e-14) return yNodes[i];
        const t = w[i] / diff;
        num += t * yNodes[i];
        den += t;
      }
      return num / den;
    };
  }

  function polyEval(fitResult, tNew) {
    if (fitResult.lagrange) {
      const f = makeLagrangeEval(fitResult.lagrange.x, fitResult.lagrange.y);
      return tNew.map(f);
    }
    const c = fitResult.coef;
    return tNew.map(t => {
      let s = 0, p = 1;
      for (let k = 0; k < c.length; k++) { s += c[k] * p; p *= t; }
      return s;
    });
  }

  // Линейная интерполяция (аналог np.interp — с зажимом на краях)
  function interpLinear(tNodes, yNodes, tNew) {
    const n = tNodes.length;
    const t0 = tNodes[0], tN = tNodes[n - 1];
    return tNew.map(x => {
      if (x < t0) return yNodes[0];
      if (x > tN) return yNodes[n - 1];
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (tNodes[mid] <= x) lo = mid; else hi = mid;
      }
      if (x === tNodes[lo]) return yNodes[lo];
      const w = (x - tNodes[lo]) / (tNodes[hi] - tNodes[lo]);
      return yNodes[lo] * (1 - w) + yNodes[hi] * w;
    });
  }

  function interpNearest(tNodes, yNodes, tNew) {
    return tNew.map(x => {
      let best = 0, bestD = Math.abs(x - tNodes[0]);
      for (let i = 1; i < tNodes.length; i++) {
        const d = Math.abs(x - tNodes[i]);
        if (d < bestD) { bestD = d; best = i; }
      }
      return yNodes[best];
    });
  }

  // Кубический сплайн "not-a-knot"
  function cubicSplineNotAKnot(x, y, xNew) {
    const n = x.length;

    if (n === 2) {
      return xNew.map(xn => (xn < x[0] || xn > x[1]) ? NaN
        : y[0] + (y[1] - y[0]) * (xn - x[0]) / (x[1] - x[0]));
    }
    if (n === 3) {
      const fit = lsPolyFit(x, y, 2);
      return polyEval(fit, xNew).map((v, k) =>
        (xNew[k] < x[0] || xNew[k] > x[n - 1]) ? NaN : v);
    }

    const h = new Array(n - 1), d = new Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
      h[i] = x[i + 1] - x[i];
      d[i] = (y[i + 1] - y[i]) / h[i];
    }

    const A = Array.from({ length: n }, () => new Array(n).fill(0));
    const rhs = new Array(n).fill(0);

    for (let i = 1; i < n - 1; i++) {
      A[i][i - 1] = h[i - 1];
      A[i][i]     = 2 * (h[i - 1] + h[i]);
      A[i][i + 1] = h[i];
      rhs[i]      = 6 * (d[i] - d[i - 1]);
    }
    A[0][0] = h[1]; A[0][1] = -(h[0] + h[1]); A[0][2] = h[0];
    A[n - 1][n - 3] = -h[n - 2];
    A[n - 1][n - 2] =  h[n - 3] + h[n - 2];
    A[n - 1][n - 1] = -h[n - 3];

    const M = solveLinearSystem(A, rhs);

    return xNew.map(xn => {
      if (xn < x[0] || xn > x[n - 1]) return NaN;
      let i;
      if (xn === x[n - 1]) i = n - 2;
      else { i = 0; while (i < n - 2 && x[i + 1] <= xn) i++; }
      const hi = h[i], dx = xn - x[i];
      return M[i]     * Math.pow(x[i + 1] - xn, 3) / (6 * hi)
           + M[i + 1] * Math.pow(dx, 3)            / (6 * hi)
           + (y[i]     - M[i]     * hi * hi / 6) * (x[i + 1] - xn) / hi
           + (y[i + 1] - M[i + 1] * hi * hi / 6) * dx             / hi;
    });
  }

  // ==================== ДАННЫЕ ====================

  const W1 = 1.2 * 2 * Math.PI;
  const W2 = 0.6 * 2 * Math.PI;

  const t  = linspace(1.0, 3.0, 2001);
  const A  = t.map(ti => Math.sin(ti * W1) + Math.cos(ti * W2));
  const tt = linspace(t[0], t[t.length - 1], 100);

  function makeNodes(N) {
    const tn = linspace(1.0, 3.0, N);
    const An = tn.map(ti => Math.sin(ti * W1) + Math.cos(ti * W2));
    return { tn, An };
  }

  // ---- N = 6 ----
  const { tn: t6, An: A6 } = makeNodes(6);

  const ap3_6 = polyEval(lsPolyFit(t6, A6, 3), tt);
  const ap5_6 = polyEval(lsPolyFit(t6, A6, 5), tt);
  const ap7_6 = polyEval(lsPolyFit(t6, A6, 7), tt);
  const ap9_6 = polyEval(lsPolyFit(t6, A6, 9), tt);

  const in0_6 = interpLinear(t6, A6, tt);
  const in1_6 = cubicSplineNotAKnot(t6, A6, tt);
  const in2_6 = interpNearest(t6, A6, tt);
  const in3_6 = polyEval(lsPolyFit(t6, A6, t6.length - 1), tt);

  // ---- N = 12 ----
  const { tn: t12, An: A12 } = makeNodes(12);

  const ap3_12 = polyEval(lsPolyFit(t12, A12, 3), tt);
  const ap5_12 = polyEval(lsPolyFit(t12, A12, 5), tt);
  const ap7_12 = polyEval(lsPolyFit(t12, A12, 7), tt);
  const ap9_12 = polyEval(lsPolyFit(t12, A12, 9), tt);

  const in0_12 = interpLinear(t12, A12, tt);
  const in1_12 = cubicSplineNotAKnot(t12, A12, tt);
  const in2_12 = interpNearest(t12, A12, tt);

  // ---- N = 24 ----
  const { tn: t24, An: A24 } = makeNodes(24);

  const ap3_24 = polyEval(lsPolyFit(t24, A24, 3), tt);
  const ap5_24 = polyEval(lsPolyFit(t24, A24, 5), tt);
  const ap7_24 = polyEval(lsPolyFit(t24, A24, 7), tt);
  const ap9_24 = polyEval(lsPolyFit(t24, A24, 9), tt);

  const in0_24 = interpLinear(t24, A24, tt);
  const in1_24 = cubicSplineNotAKnot(t24, A24, tt);
  const in2_24 = interpNearest(t24, A24, tt);

  // ---- Экстраполяция ----
  const tExt = linspace(0.5, 3.5, 301);
  const AExt = tExt.map(ti => Math.sin(ti * W1) + Math.cos(ti * W2));

  const ap3_e = polyEval(lsPolyFit(t24, A24, 3), tExt);
  const ap5_e = polyEval(lsPolyFit(t24, A24, 5), tExt);
  const ap7_e = polyEval(lsPolyFit(t24, A24, 7), tExt);
  const ap9_e = polyEval(lsPolyFit(t24, A24, 9), tExt);

  // Полином степени 23 (интерполяционный) — как pL_24 в Python
  const in3_e = polyEval(lsPolyFit(t24, A24, t24.length - 1), tExt);

  // ==================== ГРАФИК ====================

  function lineTr(x, y, name, xaxis, yaxis, legend) {
    return { x, y, name, type: 'scatter', mode: 'lines',
             line: { width: 1.5 }, xaxis, yaxis, legend };
  }
  function markerTr(x, y, name, xaxis, yaxis, legend) {
    return { x, y, name, type: 'scatter', mode: 'markers',
             marker: { size: 7 }, xaxis, yaxis, legend };
  }

  const data = [];

  // (1) Пункт 4 — N=6, МНК
  data.push(lineTr(t,  A, 'A(t)',   'x',  'y',  'legend'));
  data.push(markerTr(t6, A6, 'A0(t0)','x', 'y',  'legend'));
  data.push(lineTr(tt, ap3_6, 'Yap3','x',  'y',  'legend'));
  data.push(lineTr(tt, ap5_6, 'Yap5','x',  'y',  'legend'));
  data.push(lineTr(tt, ap7_6, 'Yap7','x',  'y',  'legend'));
  data.push(lineTr(tt, ap9_6, 'Yap9','x',  'y',  'legend'));

  // (2) Пункт 5 — N=6, интерполяция
  data.push(lineTr(t,  A, 'A(t)',   'x2', 'y2', 'legend2'));
  data.push(markerTr(t6, A6, 'A0(t0)','x2', 'y2', 'legend2'));
  data.push(lineTr(tt, in0_6, 'Yint0 (лин.)',    'x2', 'y2', 'legend2'));
  data.push(lineTr(tt, in1_6, 'Yint1 (сплайн)',  'x2', 'y2', 'legend2'));
  data.push(lineTr(tt, in2_6, 'Yint2 (nearest)', 'x2', 'y2', 'legend2'));
  data.push(lineTr(tt, in3_6, 'Yint3 (полином)', 'x2', 'y2', 'legend2'));

  // (3) Пункт 6.1 — N=12
  data.push(lineTr(t,   A,   'A(t)',  'x3', 'y3', 'legend3'));
  data.push(markerTr(t12, A12, 'A0(t0)','x3', 'y3', 'legend3'));
  data.push(lineTr(tt, ap3_12, 'Yap3', 'x3', 'y3', 'legend3'));
  data.push(lineTr(tt, ap5_12, 'Yap5', 'x3', 'y3', 'legend3'));
  data.push(lineTr(tt, ap7_12, 'Yap7', 'x3', 'y3', 'legend3'));
  data.push(lineTr(tt, ap9_12, 'Yap9', 'x3', 'y3', 'legend3'));
  data.push(lineTr(tt, in0_12, 'Yint0','x3', 'y3', 'legend3'));
  data.push(lineTr(tt, in1_12, 'Yint1','x3', 'y3', 'legend3'));
  data.push(lineTr(tt, in2_12, 'Yint2','x3', 'y3', 'legend3'));

  // (4) Пункт 6.2 — N=24
  data.push(lineTr(t,   A,   'A(t)',  'x4', 'y4', 'legend4'));
  data.push(markerTr(t24, A24, 'A0(t0)','x4', 'y4', 'legend4'));
  data.push(lineTr(tt, ap3_24, 'Yap3', 'x4', 'y4', 'legend4'));
  data.push(lineTr(tt, ap5_24, 'Yap5', 'x4', 'y4', 'legend4'));
  data.push(lineTr(tt, ap7_24, 'Yap7', 'x4', 'y4', 'legend4'));
  data.push(lineTr(tt, ap9_24, 'Yap9', 'x4', 'y4', 'legend4'));
  data.push(lineTr(tt, in0_24, 'Yint0','x4', 'y4', 'legend4'));
  data.push(lineTr(tt, in1_24, 'Yint1','x4', 'y4', 'legend4'));
  data.push(lineTr(tt, in2_24, 'Yint2','x4', 'y4', 'legend4'));

  // (5) Экстраполяция
  data.push(lineTr(tExt, AExt, 'A(t)',  'x5', 'y5', 'legend5'));
  data.push(markerTr(t24, A24, 'A0(t0)','x5', 'y5', 'legend5'));
  data.push(lineTr(tExt, ap3_e, 'Yap3', 'x5', 'y5', 'legend5'));
  data.push(lineTr(tExt, ap5_e, 'Yap5', 'x5', 'y5', 'legend5'));
  data.push(lineTr(tExt, ap7_e, 'Yap7', 'x5', 'y5', 'legend5'));
  data.push(lineTr(tExt, ap9_e, 'Yap9', 'x5', 'y5', 'legend5'));
  data.push(lineTr(tExt, in3_e, 'Yint3 (полином)', 'x5', 'y5', 'legend5'));

  // ----- Layout -----
  const LEGEND_BASE = {
    bgcolor: 'rgba(255,255,255,0.88)',
    bordercolor: '#999',
    borderwidth: 1
  };
  const legendAt = (x, y, size) => ({
    ...LEGEND_BASE,
    x, xanchor: 'right', y, yanchor: 'top',
    font: { size: size || 9 }
  });

  const axX = (anchor, domain) => ({
    anchor, domain, showgrid: true, gridcolor: '#dcdcdc', zeroline: false
  });
  const axY = (anchor, domain) => ({
    anchor, domain, showgrid: true, gridcolor: '#dcdcdc', zeroline: false
  });

  const layout = {
    margin: { l: 60, r: 25, t: 55, b: 45 },
    showlegend: true,

    xaxis:  axX('y',  [0.00, 0.47]),   yaxis:  axY('x',  [0.71, 0.97]),
    xaxis2: axX('y2', [0.53, 1.00]),   yaxis2: axY('x2', [0.71, 0.97]),
    xaxis3: axX('y3', [0.00, 0.47]),   yaxis3: axY('x3', [0.37, 0.63]),
    xaxis4: axX('y4', [0.53, 1.00]),   yaxis4: axY('x4', [0.37, 0.63]),
    xaxis5: axX('y5', [0.00, 0.47]),   yaxis5: axY('x5', [0.03, 0.29]),

    legend:  legendAt(0.47, 0.965, 9),
    legend2: legendAt(1.00, 0.965, 9),
    legend3: legendAt(0.47, 0.625, 8),
    legend4: legendAt(1.00, 0.625, 8),
    legend5: legendAt(0.47, 0.285, 8),

    annotations: [
      { text: 'Пункт 4. МНК-аппроксимация (h=0.4)',
        xref: 'paper', yref: 'paper', x: 0.235, y: 0.975,
        showarrow: false, yanchor: 'bottom', font: { size: 13 } },
      { text: 'Пункт 5. Интерполяция (h=0.4)',
        xref: 'paper', yref: 'paper', x: 0.765, y: 0.975,
        showarrow: false, yanchor: 'bottom', font: { size: 13 } },
      { text: 'Пункт 6.1 МНК-аппроксимация (h=0.2)',
        xref: 'paper', yref: 'paper', x: 0.235, y: 0.635,
        showarrow: false, yanchor: 'bottom', font: { size: 13 } },
      { text: 'Пункт 6.2 МНК-аппроксимация (h=0.4)',
        xref: 'paper', yref: 'paper', x: 0.765, y: 0.635,
        showarrow: false, yanchor: 'bottom', font: { size: 13 } },
      { text: 'Экстраполяция на [0.5 3.5] (h=0.1)',
        xref: 'paper', yref: 'paper', x: 0.235, y: 0.295,
        showarrow: false, yanchor: 'bottom', font: { size: 13 } }
    ]
  };

  Plotly.newPlot('plot', data, layout, {
    responsive: true,
    displaylogo: false
  });
})();