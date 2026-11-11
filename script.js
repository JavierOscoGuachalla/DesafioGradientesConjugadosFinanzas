/*CORE-START*/
const TOL_SYM = 1e-9;
const clone = A => A.map(r => r.slice());
const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const norm = a => Math.sqrt(dot(a, a));
const matVec = (A, x) => A.map(r => dot(r, x));
const transpose = A => A[0].map((_, j) => A.map(r => r[j]));

function det(A) {
  const n = A.length, M = clone(A); let d = 1;
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    if (M[p][k] === 0) return 0;
    if (p !== k) { [M[p], M[k]] = [M[k], M[p]]; d = -d; }
    d *= M[k][k];
    for (let i = k + 1; i < n; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j < n; j++) M[i][j] -= f * M[k][j]; }
  }
  return d;
}
function chol(A) {
  const n = A.length, L = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
    if (i === j) { if (!(s > 0)) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
  }
  return L;
}
function solveLower(L, b) { const n = b.length, y = new Array(n); for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; } return y; }
function solveUpperT(L, y) { const n = y.length, x = new Array(n); for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; } return x; }
function cholSolve(L, b) { return solveUpperT(L, solveLower(L, b)); }
function gaussSolve(A, b) {
  const n = A.length, M = A.map((r, i) => [...r, b[i]]);
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    if (Math.abs(M[p][k]) < 1e-300) throw new Error('Matriz singular');
    [M[p], M[k]] = [M[k], M[p]];
    for (let i = k + 1; i < n; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j]; }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
  return x;
}
function luSolve(A, b) { // Doolittle con pivoteo parcial: PA = LU
  const n = A.length, U = clone(A), L = Array.from({ length: n }, () => new Array(n).fill(0)), perm = [...Array(n).keys()];
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(U[i][k]) > Math.abs(U[p][k])) p = i;
    if (Math.abs(U[p][k]) < 1e-300) throw new Error('Matriz singular');
    if (p !== k) { [U[p], U[k]] = [U[k], U[p]]; [L[p], L[k]] = [L[k], L[p]]; [perm[p], perm[k]] = [perm[k], perm[p]]; }
    for (let i = k + 1; i < n; i++) { const f = U[i][k] / U[k][k]; L[i][k] = f; for (let j = k; j < n; j++) U[i][j] -= f * U[k][j]; }
  }
  for (let i = 0; i < n; i++) L[i][i] = 1;
  const pb = perm.map(i => b[i]), y = new Array(n), x = new Array(n);
  for (let i = 0; i < n; i++) { let s = pb[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= U[i][k] * x[k]; x[i] = s / U[i][i]; }
  return x;
}
function directSolve(A, b, method) {
  if (method === 'lu') return luSolve(A, b);
  if (method === 'gauss') return gaussSolve(A, b);
  const L = chol(A); if (!L) throw new Error('Cholesky falló: A no es definida positiva');
  return cholSolve(L, b);
}
function eigSym(A0) { // Jacobi cíclico; devuelve autovalores ascendentes
  const n = A0.length, a = clone(A0);
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0, tot = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { tot += a[i][j] ** 2; if (i < j) off += a[i][j] ** 2; }
    if (off <= 1e-30 * tot) break;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      if (a[p][q] === 0) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
    }
  }
  return a.map((r, i) => r[i]).sort((x, y) => x - y);
}
const sym = A => A.map((r, i) => r.map((v, j) => (v + A[j][i]) / 2));
function cond2(A) { const e = eigSym(sym(A)); return e[0] > 0 ? e[e.length - 1] / e[0] : Infinity; }

function makePrecond(A, type, omega = 1) {
  const n = A.length, I = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => +(i === j)));
  if (type === 'none') return { type, name: 'Ninguno (M = I)', M: I, solve: r => r.slice() };
  if (type === 'jacobi') {
    const M = A.map((r, i) => r.map((_, j) => i === j ? A[i][i] : 0));
    return { type, name: 'Jacobi (M = diag A)', M, solve: r => r.map((v, i) => v / A[i][i]) };
  }
  if (type === 'ic0') {
    const L = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
      if (i !== j && Math.abs(A[i][j]) < 1e-14) continue; // patrón de A
      let s = A[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (!(s > 0)) return { error: 'IC(0): pivote no positivo (ruptura de la factorización)' }; L[i][i] = Math.sqrt(s); }
      else L[i][j] = s / L[j][j];
    }
    const M = L.map((r, i) => L.map((_, j) => dot(L[i], L[j])));
    return { type, name: 'Cholesky incompleto IC(0)', M, solve: r => cholSolve(L, r) };
  }
  if (type === 'ssor') {
    const w = omega, K = A.map((r, i) => r.map((v, j) => j < i ? w * v : (i === j ? v : 0))); // D + wL
    const c = 1 / (w * (2 - w));
    const KD = K.map((r, i) => r.map((v, j) => v / A[j][j])); // K D^-1
    const M = KD.map((r, i) => K.map((_, j) => c * dot(r, K[j])));
    const solve = r => {
      const rp = r.map(v => v * w * (2 - w));
      const v = solveLower(K, rp);            // K v = r'
      const u = v.map((x, i) => x * A[i][i]); // u = D v
      return solveUpperT(K, u);               // K^T z = u
    };
    return { type, name: `SSOR (ω = ${omega})`, M, solve };
  }
  throw new Error('Precondicionador desconocido');
}
function condPrecond(A, M) {
  const C = chol(M); if (!C) return NaN;
  const n = A.length, Y = Array.from({ length: n }, () => new Array(n));
  for (let j = 0; j < n; j++) { const col = solveLower(C, A.map(r => r[j])); for (let i = 0; i < n; i++) Y[i][j] = col[i]; }
  const Yt = transpose(Y), B = Array.from({ length: n }, () => new Array(n));
  for (let j = 0; j < n; j++) { const col = solveLower(C, Yt.map(r => r[j])); for (let i = 0; i < n; i++) B[i][j] = col[i]; }
  return cond2(B);
}
function pcg(A, b, x0, P, tol, maxIter) {
  let x = x0.slice(); let Ax = matVec(A, x); let r = b.map((v, i) => v - Ax[i]);
  const nb = norm(b), rel = () => nb === 0 ? norm(r) : norm(r) / nb;
  const hist = [rel()]; let converged = hist[0] < tol, breakdown = false, iters = 0;
  if (!converged) {
    let z = P.solve(r), p = z.slice(), rz = dot(r, z);
    for (let k = 1; k <= maxIter; k++) {
      const Ap = matVec(A, p), pAp = dot(p, Ap);
      if (!(pAp > 0)) { breakdown = true; break; }
      const alpha = rz / pAp;
      x = x.map((v, i) => v + alpha * p[i]); r = r.map((v, i) => v - alpha * Ap[i]);
      iters = k; hist.push(rel());
      if (hist[hist.length - 1] < tol) { converged = true; break; }
      z = P.solve(r); const rzN = dot(r, z), beta = rzN / rz;
      p = z.map((v, i) => v + beta * p[i]); rz = rzN;
    }
  }
  const Axf = matVec(A, x), rf = b.map((v, i) => v - Axf[i]);
  return { x, iters, hist, converged, breakdown, resRel: nb === 0 ? norm(rf) : norm(rf) / nb };
}
function validate(A, b) {
  const errors = [], n = A.length;
  if (n < 4) errors.push('n debe ser ≥ 4.');
  if (!A.every(r => r.length === n) || b.length !== n) { errors.push('Dimensiones inconsistentes (A debe ser n×n y b de longitud n).'); return { ok: false, errors }; }
  if (!A.every(r => r.every(Number.isFinite)) || !b.every(Number.isFinite)) { errors.push('Hay celdas vacías o no numéricas en A o b.'); return { ok: false, errors }; }
  let asym = 0, mx = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { asym = Math.max(asym, Math.abs(A[i][j] - A[j][i])); mx = Math.max(mx, Math.abs(A[i][j])); }
  const symOk = asym <= TOL_SYM * Math.max(1, mx);
  if (!symOk) errors.push(`A no es simétrica: max|aᵢⱼ − aⱼᵢ| = ${asym.toExponential(3)}.`);
  const As = sym(A), minors = []; for (let k = 1; k <= n; k++) minors.push(det(As.slice(0, k).map(r => r.slice(0, k))));
  const eig = eigSym(As), pd = minors.every(m => m > 0) && chol(As) !== null && eig[0] > 0;
  if (!pd) errors.push(`A no es definida positiva (λ_min = ${eig[0].toExponential(3)}; menores principales: ${minors.map(m => m.toExponential(2)).join(', ')}).`);
  if (b.every(v => v === 0)) errors.push('Aviso: b = 0 ⇒ solución trivial x = 0 (el residuo relativo usa ‖r‖₂).');
  const fatal = errors.filter(e => !e.startsWith('Aviso'));
  return { ok: fatal.length === 0, errors, asym, symOk, minors, eig, pd };
}
function validateConfig(r, n) { // valida y sanea los parámetros de ejecución
  const issues = [], c = { precond: r.precond, method: r.method, mirror: r.mirror };
  if (Number.isFinite(r.tol) && r.tol > 0 && r.tol < 1) c.tol = r.tol; else { issues.push({ id: 'cfgTol', msg: 'Tolerancia inválida: debe cumplirse 0 < tol < 1 (se usa 1e-6).' }); c.tol = 1e-6; }
  if (Number.isInteger(r.maxIter) && r.maxIter >= 1 && r.maxIter <= 100000) c.maxIter = r.maxIter; else { issues.push({ id: 'cfgMax', msg: 'N máx. inválido: entero entre 1 y 100000 (se usa 1000).' }); c.maxIter = 1000; }
  if (Number.isFinite(r.omega) && r.omega > 0 && r.omega < 2) c.omega = r.omega; else { issues.push({ id: 'cfgOmega', msg: 'ω inválido: debe cumplirse 0 < ω < 2 (se usa 1).' }); c.omega = 1; }
  c.x0 = Array.from({ length: n }, (_, i) => { if (Number.isFinite(r.x0[i])) return r.x0[i]; issues.push({ id: 'x0_' + i, msg: `x₀${i + 1} no es numérico (se usa 0).` }); return 0; });
  return { cfg: c, issues };
}
function solveScenario(A, b, cfg) {
  const val = validate(A, b);
  if (!val.ok) return { valid: false, val };
  try {
    const n = A.length, xref = directSolve(A, b, cfg.method);
    let P = makePrecond(A, cfg.precond, cfg.omega), precondNote = null;
    if (P.error) { precondNote = P.error + ' → se usa Jacobi.'; P = makePrecond(A, 'jacobi'); }
    const Id = makePrecond(A, 'none'), x0 = cfg.x0.slice(0, n);
    const gc = pcg(A, b, x0, Id, cfg.tol, cfg.maxIter), gcp = pcg(A, b, x0, P, cfg.tol, cfg.maxIter);
    const Ea = xref.map((v, i) => Math.abs(v - gcp.x[i])), Er = xref.map((v, i) => v === 0 ? null : Ea[i] / Math.abs(v));
    const Ep = Er.map(v => v === null ? null : v * 100), EaGC = xref.map((v, i) => Math.abs(v - gc.x[i]));
    return { valid: true, val, xref, gc, gcp, P, precondNote, Ea, Er, Ep, EaGC, kA: cond2(A), kMA: condPrecond(A, P.M) };
  } catch (e) { return { valid: false, val: { ok: false, errors: ['Error numérico: ' + e.message] } }; }
}

/* ---- Modelo de negocio: ajuste presupuestario ----
   J(x) = ½Σ wᵢ(xᵢ−dᵢ)² + ½Σ_{i<j} kᵢⱼ(xᵢ−xⱼ)² + ½ρ(Σxᵢ−B)²   ⇒   (W + L + ρ·11ᵀ)x = W·d + ρ·B·1 */
const NAMES = ['Desarrollo de producto', 'Infraestructura cloud', 'QA y pruebas', 'Seguridad informática', 'Capacitación', 'Soporte y mantenimiento'];
const MEANING = [
  'Planilla de ingenieros y líderes técnicos (mayor rubro).',
  'Servidores, servicios cloud y licencias; sensible a inflación y tipo de cambio.',
  'Automatización, entornos de prueba y testers; ligado a Desarrollo.',
  'Auditorías, herramientas y cumplimiento; ligado a Infraestructura y QA.',
  'Cursos y certificaciones del personal; rubro más postergable.',
  'Atención a clientes y mantenimiento correctivo; ligado a Desarrollo.'];
const KPAIRS = [[0, 1, .5], [0, 2, .8], [0, 5, .4], [1, 3, .5], [2, 3, .3], [4, 5, .2]];
function buildK(mult) {
  const n = NAMES.length, k = Array.from({ length: n }, () => new Array(n).fill(0.1 * mult));
  KPAIRS.forEach(([i, j, v]) => { k[i][j] = k[j][i] = v * mult; });
  for (let i = 0; i < n; i++) k[i][i] = 0; return k;
}
const DEFS = {
  ideal: { label: 'Ideal', w: [4, 3, 3, 3, 2, 3], d: [500, 200, 120, 90, 60, 130], rho: 2, B: 1000, kmult: 1,
    text: 'Año estable: las solicitudes suman 1 100 frente a un techo de 1 000 (déficit moderado de 9 %). Las prioridades wᵢ son comparables entre categorías (misma escala de medición), el techo se aplica con una penalización suave (ρ = 2) y los vínculos entre categorías son los de la operación normal. A queda bien escalada ⇒ κ(A) bajo y convergencia rápida.' },
  estres: { label: 'Bajo estrés', w: [4, 3, 3, 5, 1.5, 3], d: [560, 280, 130, 120, 70, 150], rho: 6, B: 900, kmult: 1.5,
    text: 'Shock financiero: la inflación y el tipo de cambio elevan las solicitudes (cloud +40 %, seguridad +33 %), los ingresos caen y el techo baja a 900. La dirección endurece el tope (ρ = 6), sube la prioridad de Seguridad (5) y baja la de Capacitación (1,5), y la crisis refuerza la interdependencia (kᵢⱼ × 1,5). Cambian A y b; se analiza cuánto se desplaza el equilibrio y si cambia el condicionamiento.' },
  mal: { label: 'Mal condicionado', w: [1, 10, 100, 1000, 10000, 3], d: [500, 200, 120, 90, 60, 130], rho: 5, B: 1000, kmult: 1,
    text: 'Mal escalamiento: las prioridades wᵢ abarcan cuatro órdenes de magnitud (1 a 10 000) porque cada rubro se pondera en unidades distintas (miles de Bs, Bs, USD). La matriz sigue siendo simétrica definida positiva pero con κ(A) ≫ 1. El mal escalamiento diagonal se corrige con Jacobi; el fuerte acoplamiento por el techo ρ entre categorías de baja prioridad se ataca con IC(0) o SSOR.' }
};
function buildAb(d) {
  const n = d.w.length, K = buildK(d.kmult), A = Array.from({ length: n }, () => new Array(n)), b = [];
  for (let i = 0; i < n; i++) {
    let s = 0; for (let j = 0; j < n; j++) if (j !== i) s += K[i][j];
    for (let j = 0; j < n; j++) A[i][j] = i === j ? d.w[i] + s + d.rho : d.rho - K[i][j];
    b.push(d.w[i] * d.d[i] + d.rho * d.B);
  }
  return { A, b };
}
const shares = x => { const s = x.reduce((a, v) => a + v, 0); return Math.abs(s) < 1e-300 ? x.map(() => NaN) : x.map(v => v / s); };
/*CORE-END*/

/* ================= CAPA DE ESTADO Y CONTROLADOR ================= */
const KEYS = ['ideal', 'estres', 'mal'], N = NAMES.length;
const COLORS = { gc: '#fbbf24', gcp: '#2dd4bf' };
const S = {}; let active = 'ideal', cfgIssues = [];
const cfg = { precond: 'jacobi', omega: 1, tol: 1e-6, maxIter: 1000, method: 'cholesky', x0: new Array(N).fill(0), mirror: true };
function initScenarios() { KEYS.forEach(k => { const ab = buildAb(DEFS[k]); S[k] = { A: ab.A, b: ab.b }; }); }
const $ = id => document.getElementById(id);
function f(v) {
  if (v === null || v === undefined) return '—';
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return '0';
  const a = Math.abs(v); return (a >= 1e5 || a < 1e-3) ? v.toExponential(4) : String(+v.toPrecision(7));
}
const fv = v => '[' + v.map(f).join(', ') + ']';
const pill = ok => `<span class="pill ${ok ? 'ok' : 'bad'}">${ok ? 'Sí' : 'No'}</span>`;
const esc = s => String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
const sub = i => '₁₂₃₄₅₆₇₈₉'[i];

function readCfg() {
  const raw = { precond: $('cfgPre').value, omega: parseFloat($('cfgOmega').value), tol: parseFloat($('cfgTol').value),
    maxIter: Number($('cfgMax').value), method: $('cfgDir').value, mirror: $('cfgMirror').checked,
    x0: [...document.querySelectorAll('.x0in')].map(e => e.value.trim() === '' ? NaN : parseFloat(e.value)) };
  if ($('cfgMax').value.trim() === '') raw.maxIter = NaN;
  const v = validateConfig(raw, N); Object.assign(cfg, v.cfg); cfgIssues = v.issues;
}
function buildInputs() {
  const A = S[active].A, b = S[active].b;
  $('tblA').innerHTML = '<tr><th></th>' + NAMES.map((_, j) => `<th>x${sub(j)}</th>`).join('') + '</tr>' +
    A.map((r, i) => `<tr><th class="l">x${sub(i)}</th>` + r.map((v, j) => `<td><input type="number" step="any" data-i="${i}" data-j="${j}" value="${v}"></td>`).join('') + '</tr>').join('');
  $('tblB').innerHTML = '<tr><th>i</th><th>bᵢ</th></tr>' + b.map((v, i) => `<tr><th>x${sub(i)}</th><td><input type="number" step="any" data-b="${i}" value="${v}"></td></tr>`).join('');
  if (!document.querySelector('.x0in')) $('x0box').innerHTML = NAMES.map((_, i) => `<div style="min-width:80px"><label>x₀${sub(i)}</label><input class="x0in" id="x0_${i}" type="number" step="any" value="0"></div>`).join('');
  $('entActivo').textContent = DEFS[active].label;
}
function onCellInput(e) {
  const t = e.target;
  if (t.dataset.i !== undefined) {
    const i = +t.dataset.i, j = +t.dataset.j, v = t.value.trim() === '' ? NaN : parseFloat(t.value);
    S[active].A[i][j] = v;
    if (cfg.mirror && i !== j) { S[active].A[j][i] = v; const m = document.querySelector(`input[data-i="${j}"][data-j="${i}"]`); if (m) m.value = t.value; }
  } else if (t.dataset.b !== undefined) { S[active].b[+t.dataset.b] = t.value.trim() === '' ? NaN : parseFloat(t.value); }
}
function markInputs() { // validación visual por celda
  document.querySelectorAll('#entrada input').forEach(e => e.classList.remove('bad', 'asym'));
  const A = S[active].A, b = S[active].b;
  for (let i = 0; i < N; i++) {
    const bi = document.querySelector(`input[data-b="${i}"]`); if (bi && !Number.isFinite(b[i])) bi.classList.add('bad');
    for (let j = 0; j < N; j++) {
      const c = document.querySelector(`input[data-i="${i}"][data-j="${j}"]`); if (!c) continue;
      if (!Number.isFinite(A[i][j]) || (i === j && A[i][j] <= 0)) c.classList.add('bad');
      else if (i !== j && Number.isFinite(A[j][i]) && Math.abs(A[i][j] - A[j][i]) > TOL_SYM * Math.max(1, Math.abs(A[i][j]))) c.classList.add('asym');
    }
  }
  cfgIssues.forEach(s => { const el = $(s.id); if (el) el.classList.add('bad'); });
}
function resultsAll() { const R = {}; KEYS.forEach(k => R[k] = solveScenario(S[k].A, S[k].b, cfg)); return R; }

/* ================= COMPONENTES (VISTA) ================= */
function drawConv(cv, series, tol) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight; cv.width = w * dpr; cv.height = h * dpr;
  const c = cv.getContext('2d'); c.scale(dpr, dpr); c.clearRect(0, 0, w, h);
  const L = 56, R = 14, T = 16, B = 36, pw = w - L - R, ph = h - T - B;
  const lg = v => Math.log10(Math.max(v, 1e-17));
  let maxLen = 2, lo = lg(tol) - 0.5, hi = -1;
  series.forEach(s => { maxLen = Math.max(maxLen, s.hist.length); s.hist.forEach(v => { lo = Math.min(lo, lg(v) - 0.3); hi = Math.max(hi, lg(v) + 0.3); }); });
  lo = Math.floor(lo); hi = Math.ceil(hi);
  const X = k => L + (k / (maxLen - 1)) * pw, Y = v => T + (hi - lg(v)) / (hi - lo) * ph;
  c.font = '11px system-ui'; c.fillStyle = '#93a4c3'; c.strokeStyle = '#243352'; c.lineWidth = 1;
  for (let e = lo; e <= hi; e++) { const y = T + (hi - e) / (hi - lo) * ph; c.beginPath(); c.moveTo(L, y); c.lineTo(w - R, y); c.stroke(); c.textAlign = 'right'; c.fillText('1e' + e, L - 6, y + 4); }
  const step = Math.max(1, Math.ceil((maxLen - 1) / 10)); c.textAlign = 'center';
  for (let k = 0; k < maxLen; k += step) c.fillText(k, X(k), h - B + 16);
  c.fillText('iteración k', L + pw / 2, h - 6);
  c.save(); c.translate(12, T + ph / 2); c.rotate(-Math.PI / 2); c.fillText('‖rₖ‖₂ / ‖b‖₂', 0, 0); c.restore();
  c.strokeStyle = '#f87171'; c.setLineDash([6, 4]); c.beginPath(); c.moveTo(L, Y(tol)); c.lineTo(w - R, Y(tol)); c.stroke(); c.setLineDash([]);
  c.fillStyle = '#f87171'; c.textAlign = 'left'; c.fillText('tolerancia', L + 6, Y(tol) - 5);
  series.forEach(s => {
    c.strokeStyle = s.color; c.lineWidth = 2.2; c.beginPath();
    s.hist.forEach((v, k) => { const x = X(k), y = Y(v); k ? c.lineTo(x, y) : c.moveTo(x, y); }); c.stroke();
    c.fillStyle = s.color; s.hist.forEach((v, k) => { c.beginPath(); c.arc(X(k), Y(v), 3, 0, 7); c.fill(); });
  });
}
function invalidBox(r) { return `<div class="banner err"><b>No calculable.</b> ${r.val.errors.map(e => '• ' + esc(e)).join('<br>')}</div>`; }
const stats = (x, k) => { const d = DEFS[k].d, tot = x.reduce((a, v) => a + v, 0); return { tot, B: DEFS[k].B, cuts: x.map((v, i) => (v - d[i]) / d[i]) }; };

function interpret(k, r) {
  if (!r.valid) return 'Sistema no válido: no hay interpretación.';
  const st = stats(r.xref, k), imin = st.cuts.reduce((m, v, i) => v < st.cuts[m] ? i : m, 0), neg = r.xref.map((v, i) => v < 0 ? NAMES[i] : null).filter(Boolean);
  let s = `Se asignan ${f(+st.tot.toFixed(1))} (techo ${st.B}); el mayor recorte relativo cae en ${NAMES[imin]} (${(st.cuts[imin] * 100).toFixed(1)} % frente a lo solicitado)` + (neg.length ? `, y ${neg.join(', ')} quedaría con asignación negativa (no factible)` : '') + '. ';
  s += `κ(A) = ${f(r.kA)}; con ${r.P.name}, κ(M⁻¹A) = ${f(r.kMA)}. GC: ${r.gc.iters} it., GCP: ${r.gcp.iters} it. `;
  if (k === 'ideal') s += 'El ajuste es suave y estable: todas las categorías reciben cerca de lo pedido.';
  else if (k === 'estres') s += 'El estrés desplaza el punto de equilibrio (recortes fuertes en rubros postergables) más de lo que altera la dificultad numérica.';
  else s += r.kA > 1e3 ? 'Sistema muy sensible a perturbaciones: pequeñas variaciones de b cambian la asignación; el precondicionador compensa el mal escalamiento entre rubros.' : 'El condicionamiento resultó moderado con los datos actuales.';
  return s;
}
function renderDatos() {
  let h = '<table><tr><th class="l">Escenario</th><th class="l">Prioridades w (1…6)</th><th class="l">Solicitudes d (miles de Bs)</th><th>Techo B</th><th>ρ</th><th>Factor de k</th></tr>';
  KEYS.forEach(k => { const d = DEFS[k]; h += `<tr><td class="l">${d.label}</td><td class="l m">${d.w.join(' / ')}</td><td class="l m">${d.d.join(' / ')} (Σ=${d.d.reduce((a, v) => a + v, 0)})</td><td class="m">${d.B}</td><td class="m">${d.rho}</td><td class="m">×${d.kmult}</td></tr>`; });
  h += '</table><p class="mut" style="margin-bottom:0">Vínculos kᵢⱼ base: Desarrollo–Cloud 0,5; Desarrollo–QA 0,8; Desarrollo–Soporte 0,4; Cloud–Seguridad 0,5; QA–Seguridad 0,3; Capacitación–Soporte 0,2; resto 0,1.</p>';
  $('datosEscenarios').innerHTML = h;
  $('tblVars').innerHTML = '<tr><th class="l">Variable</th><th class="l">Categoría</th><th class="l">Significado</th></tr>' + NAMES.map((n, i) => `<tr><td class="l">x${sub(i)}</td><td class="l"><b>${n}</b></td><td class="l">${MEANING[i]}</td></tr>`).join('');
}
function renderFormulacion(R) {
  const A = S[active].A, b = S[active].b; $('fmActivo').textContent = DEFS[active].label;
  let eq = ''; A.forEach((r, i) => { eq += r.map((v, j) => `${f(v).padStart(9)}·x${sub(j)}`).join(' + ') + `  =  ${f(b[i])}\n`; });
  $('fmSistema').textContent = eq;
  let h = '<table><tr><th class="l">Escenario</th><th>max|aᵢⱼ−aⱼᵢ|</th><th>Simétrica</th><th class="l">Autovalores λ (asc.)</th><th class="l">Menores principales Δ₁…Δₙ</th><th>Definida positiva</th></tr>';
  KEYS.forEach(k => { const v = validate(S[k].A, S[k].b);
    h += `<tr><td class="l">${DEFS[k].label}</td>` + (v.eig ? `<td class="m">${f(v.asym)}</td><td>${pill(v.symOk)}</td><td class="l m">${v.eig.map(f).join(', ')}</td><td class="l m">${v.minors.map(f).join(', ')}</td><td>${pill(v.pd)}</td>` : `<td colspan="5" class="l bad">${esc(v.errors.join(' '))}</td>`) + '</tr>'; });
  $('fmPruebas').innerHTML = h + '</table>';
}
function renderEntradaMsg(R) {
  markInputs();
  const r = R[active]; let m = '';
  if (cfgIssues.length) m += `<div class="banner wrn"><b>Parámetros de ejecución corregidos:</b><br>${cfgIssues.map(s => '• ' + esc(s.msg)).join('<br>')}</div>`;
  if (!r.valid) m += invalidBox(r);
  else {
    m += `<div class="banner good">✔ A es simétrica y definida positiva (λ_min = ${f(r.val.eig[0])}). Cálculo actualizado.</div>`;
    if (r.precondNote) m += `<div class="banner wrn">${esc(r.precondNote)}</div>`;
    if (r.val.errors.length) m += `<div class="banner wrn">${r.val.errors.map(esc).join('<br>')}</div>`;
  }
  $('entMsg').innerHTML = m;
}
function renderEscenarios(R) {
  KEYS.forEach(k => $('bt_' + k).classList.toggle('act', k === active));
  $('escJust').innerHTML = `<b>${DEFS[active].label}.</b> ${DEFS[active].text}`;
  const rows = [
    ['Solución numérica de referencia', r => fv(r.xref)], ['Solución GCP', r => fv(r.gcp.x)],
    ['Error absoluto (Eₐ) por variable', r => fv(r.Ea)], ['Error porcentual (Eₚ) por variable', r => fv(r.Ep.map(v => v === null ? NaN : v))],
    ['Iteraciones GC (sin precondicionar)', r => r.gc.iters + (r.gc.converged ? '' : ' (no conv.)')], ['Iteraciones GCP (precondicionado)', r => r.gcp.iters + (r.gcp.converged ? '' : ' (no conv.)')],
    ['Número de condición κ(A)', r => f(r.kA)], ['Número de condición κ(M⁻¹A)', r => f(r.kMA)],
    ['Residuo relativo final GCP', r => f(r.gcp.resRel)], ['¿Convergió el algoritmo?', r => pill(r.gcp.converged)],
    ['Interpretación financiera', null]
  ];
  let h = '<table><tr><th class="l">Indicador</th>' + KEYS.map(k => `<th>${DEFS[k].label}${k === active ? ' ◀' : ''}</th>`).join('') + '</tr>';
  rows.forEach(([lab, fn]) => {
    h += `<tr><td class="l">${lab}</td>` + KEYS.map(k => {
      const r = R[k]; if (!r.valid) return '<td class="bad l">Sistema no válido</td>';
      if (!fn) return `<td class="l" style="font-size:12.5px">${esc(interpret(k, r))}</td>`;
      return `<td class="m">${fn(r)}</td>`; }).join('') + '</tr>'; });
  $('tblComp').innerHTML = h + `</table><p class="mut">Precondicionador: ${esc(cfg.precond === 'ssor' ? 'SSOR (ω=' + cfg.omega + ')' : cfg.precond)} · tolerancia ${cfg.tol} · método directo: ${cfg.method}. Las columnas reflejan los datos actuales (incluidas ediciones). Unidades: miles de Bs.</p>`;
}
function renderRef(R) {
  const r = R[active]; if (!r.valid) { $('refBox').innerHTML = invalidBox(r); return; }
  const sh = shares(r.xref);
  $('refBox').innerHTML = `<div class="card"><div class="kpi"><div><b>${cfg.method === 'lu' ? 'LU (Doolittle)' : cfg.method === 'gauss' ? 'Gauss' : 'Cholesky'}</b><span>método directo</span></div>
    <div><b>${f(norm(matVec(S[active].A, r.xref).map((v, i) => v - S[active].b[i])) / (norm(S[active].b) || 1))}</b><span>residuo relativo de la referencia</span></div></div>
    <table><tr><th class="l">Variable</th><th class="l">Categoría</th><th>x_ref (miles de Bs)</th><th>Participación (%)</th></tr>` +
    r.xref.map((v, i) => `<tr><td class="l">x${sub(i)}</td><td class="l">${NAMES[i]}</td><td class="m">${f(v)}</td><td class="m">${f(sh[i] * 100)}</td></tr>`).join('') +
    `</table><p class="mut" style="margin-bottom:0">Vector de referencia: ${fv(r.xref)}</p></div>`;
}
function renderGCP(R) {
  const r = R[active]; if (!r.valid) { $('gcpBox').innerHTML = invalidBox(r); $('preCmp').innerHTML = ''; drawConv($('cvConv'), [], cfg.tol); return; }
  const g = (o, nm, col) => `<div class="card"><b><span class="sw" style="background:${col}"></span>${nm}</b>
    <div class="kpi"><div><b>${o.iters}</b><span>iteraciones</span></div><div><b>${f(o.resRel)}</b><span>residuo relativo final</span></div><div><b>${o.converged ? 'Sí' : 'No'}</b><span>convergió</span></div></div>
    <div class="vec">x = ${fv(o.x)}</div><p class="mut" style="margin:8px 0 0">${o.converged ? 'Alcanzó la tolerancia.' : o.breakdown ? 'Ruptura numérica (pᵀAp ≤ 0).' : 'Alcanzó N máx. sin cumplir la tolerancia.'}</p></div>`;
  const dI = r.gc.iters - r.gcp.iters;
  $('gcpBox').innerHTML = (r.precondNote ? `<div class="banner wrn">${esc(r.precondNote)}</div>` : '') +
    `<div class="grid2">${g(r.gc, 'GC simple (M = I)', COLORS.gc)}${g(r.gcp, 'GCP · ' + r.P.name, COLORS.gcp)}</div>
     <div class="kpi"><div><b>${f(r.kA)}</b><span>κ(A)</span></div><div><b>${f(r.kMA)}</b><span>κ(M⁻¹A)</span></div><div><b>${dI >= 0 ? '−' + dI : '+' + (-dI)}</b><span>variación de iteraciones con GCP</span></div></div>`;
  drawConv($('cvConv'), [{ hist: r.gc.hist, color: COLORS.gc }, { hist: r.gcp.hist, color: COLORS.gcp }], cfg.tol);
  let h = '<table><tr><th class="l">Precondicionador</th><th>Iteraciones</th><th>Residuo relativo</th><th>κ(M⁻¹A)</th><th>Convergió</th></tr>';
  [['none', 'Ninguno (GC)'], ['jacobi', 'Jacobi'], ['ic0', 'IC(0)'], ['ssor', `SSOR (ω=${cfg.omega})`]].forEach(([t, nm]) => {
    const P = makePrecond(S[active].A, t, cfg.omega);
    if (P.error) { h += `<tr><td class="l">${nm}</td><td colspan="4" class="l warn">${esc(P.error)}</td></tr>`; return; }
    const o = pcg(S[active].A, S[active].b, cfg.x0, P, cfg.tol, cfg.maxIter);
    h += `<tr><td class="l">${nm}</td><td class="m">${o.iters}</td><td class="m">${f(o.resRel)}</td><td class="m">${f(condPrecond(S[active].A, P.M))}</td><td>${pill(o.converged)}</td></tr>`; });
  const maxK = Math.min(12, Math.max(r.gc.hist.length, r.gcp.hist.length));
  let t2 = '<h3>Bitácora de las primeras iteraciones</h3><table><tr><th>k</th><th>‖r‖/‖b‖ GC</th><th>‖r‖/‖b‖ GCP</th></tr>';
  for (let k = 0; k < maxK; k++) t2 += `<tr><td>${k}</td><td class="m">${k < r.gc.hist.length ? f(r.gc.hist[k]) : '—'}</td><td class="m">${k < r.gcp.hist.length ? f(r.gcp.hist[k]) : '—'}</td></tr>`;
  $('preCmp').innerHTML = h + '</table><p class="mut">Con A densa, IC(0) coincide con la factorización de Cholesky completa (no hay relleno que descartar), por eso M⁻¹A ≈ I y converge en 1 iteración.</p>' + t2 + '</table>';
}
function renderErrores(R) {
  const r = R[active];
  if (!r.valid) { $('errBox').innerHTML = invalidBox(r); $('errExtra').innerHTML = ''; ['edaBox', 'portBox', 'sensBox'].forEach(i => $(i).innerHTML = '<span class="bad">No calculable</span>'); return; }
  $('errBox').innerHTML = '<table><tr><th class="l">Var.</th><th>x_ref</th><th>x_GCP</th><th>Eₐ</th><th>Eᵣ</th><th>Eₚ (%)</th><th>Eₐ del GC simple</th></tr>' +
    r.xref.map((v, i) => `<tr><td class="l">x${sub(i)}</td><td class="m">${f(v)}</td><td class="m">${f(r.gcp.x[i])}</td><td class="m">${f(r.Ea[i])}</td><td class="m">${r.Er[i] === null ? 'no definido (x_ref = 0)' : f(r.Er[i])}</td><td class="m">${r.Ep[i] === null ? '—' : f(r.Ep[i])}</td><td class="m">${f(r.EaGC[i])}</td></tr>`).join('') + '</table>';
  const nb = norm(S[active].b), erOk = r.Er.filter(v => v !== null);
  $('errExtra').innerHTML = `<div class="kpi"><div><b>${nb === 0 ? 'no def.' : f(r.gcp.resRel)}</b><span>residuo relativo global ‖b−Ax_GCP‖₂/‖b‖₂</span></div><div><b>${f(r.kA)}</b><span>κ(A)</span></div><div><b>${f(r.kMA)}</b><span>κ(M⁻¹A)</span></div>
    <div><b>${erOk.length ? f(Math.max(...erOk)) : '—'}</b><span>error relativo máximo</span></div></div>
    <p class="mut">Cota teórica: ‖x−x̂‖/‖x‖ ≤ κ(A)·‖r‖/‖b‖ = ${f(r.kA * r.gcp.resRel)}. Un residuo pequeño no implica error pequeño cuando κ(A) es grande.</p>`;
  const A = S[active].A, b = S[active].b, x = r.xref, d = DEFS[active].d, st = stats(x, active);
  $('edaBox').innerHTML = '<table><tr><th class="l">Categoría</th><th>Solicitado d</th><th>Asignado x</th><th>Δ vs solicitud</th></tr>' +
    x.map((v, i) => `<tr><td class="l">${NAMES[i]}</td><td class="m">${d[i]}</td><td class="m">${f(+v.toFixed(2))}</td><td class="m ${st.cuts[i] < 0 ? 'bad' : 'ok'}">${(st.cuts[i] * 100).toFixed(1)} %</td></tr>`).join('') +
    '</table><p class="mut" style="margin-bottom:0">«Solicitado» corresponde a los datos originales del escenario.</p>';
  const J = 0.5 * dot(x, matVec(A, x)) - dot(b, x);
  $('portBox').innerHTML = `<div class="kpi"><div><b>${f(+st.tot.toFixed(2))}</b><span>total asignado Σxᵢ</span></div><div><b>${st.B}</b><span>techo B del escenario</span></div><div><b>${f(+((st.tot - st.B) / st.B * 100).toFixed(2))} %</b><span>exceso sobre el techo (penalización suave)</span></div>
    <div><b>${f(J)}</b><span>J(x) = ½xᵀAx − bᵀx (mínimo, sin constante)</span></div><div><b>${x.filter(v => v < 0).length}</b><span>categorías con asignación negativa</span></div></div>`;
  let h = '<table><tr><th class="l">Perturbación</th>' + NAMES.map(n => `<th>Δx ${n.split(' ')[0]}</th>`).join('') + '</tr>';
  for (let i = 0; i < b.length; i++) { const b2 = b.slice(); b2[i] *= 1.1; let x2; try { x2 = directSolve(A, b2, 'cholesky'); } catch (e) { continue; }
    h += `<tr><td class="l">b${sub(i)} +10 %</td>` + x2.map((v, j) => `<td class="m">${f(+(v - x[j]).toFixed(3))}</td>`).join('') + '</tr>'; }
  $('sensBox').innerHTML = h + '</table><p class="mut">Subir la presión de una categoría le quita fondos a las demás por el techo común; el efecto se amplifica con κ(A) y con el acoplamiento.</p>';
}
function renderConc(R) {
  const ok = KEYS.every(k => R[k].valid);
  if (!ok) { $('concBox').innerHTML = '<span class="bad">Alguno de los sistemas no es válido; corrige los datos para ver las conclusiones.</span>'; return; }
  const li = KEYS.map(k => { const r = R[k], e = r.Ep.filter(v => v !== null); return `<li><b>${DEFS[k].label}:</b> κ(A) = ${f(r.kA)}, κ(M⁻¹A) = ${f(r.kMA)}; GC ${r.gc.iters} it. vs GCP ${r.gcp.iters} it.; Eₚ máx = ${f(Math.max(...e))} %; residuo ${f(r.gcp.resRel)}.</li>`; }).join('');
  $('concBox').innerHTML = `<p><b>Interpretación técnica.</b></p><ul>${li}</ul>
    <ul><li>En los tres escenarios el GCP alcanzó la tolerancia ${cfg.tol} y coincide con la solución directa. El precondicionador (${esc(R.mal.P.name)}) reduce el número de condición efectivo, sobre todo en el escenario mal condicionado, donde Jacobi corrige el escalamiento y IC(0)/SSOR reducen aún más las iteraciones (módulo 7).</li>
    <li>Un residuo pequeño no garantiza un error pequeño: la cota κ(A)·‖r‖/‖b‖ crece con el mal condicionamiento.</li>
    <li>Con A densa, IC(0) equivale a Cholesky completo; en sistemas grandes y dispersos IC(0) y SSOR dan un compromiso entre costo por iteración y mejora de convergencia.</li></ul>
    <p><b>Interpretación financiera.</b> En el escenario ideal el techo se cumple casi sin sacrificar ninguna categoría. El estrés desplaza el equilibrio: los recortes se concentran en los rubros de baja prioridad (p. ej. Capacitación) y casi no cambia la dificultad numérica. El escenario mal condicionado muestra el riesgo de modelo: si las prioridades se expresan en escalas muy distintas, pequeños errores en los datos producen grandes cambios en la asignación; conviene normalizar las unidades antes de resolver.</p>
    <p><b>Calidad del software.</b> Toda entrada se valida (dimensiones, valores numéricos, simetría, definición positiva y parámetros de ejecución) y un fallo en un componente no detiene a los demás; las pruebas del módulo 9 lo comprueban.</p>`;
}

/* ---- Registro de componentes con error boundary ---- */
const Components = [
  { name: 'Formulación', host: 'fmPruebas', render: renderFormulacion },
  { name: 'Entrada de datos', host: 'entMsg', render: renderEntradaMsg },
  { name: 'Tres escenarios', host: 'tblComp', render: renderEscenarios },
  { name: 'Solución de referencia', host: 'refBox', render: renderRef },
  { name: 'GCP vs GC', host: 'gcpBox', render: renderGCP },
  { name: 'Análisis de errores', host: 'errBox', render: renderErrores },
  { name: 'Conclusiones', host: 'concBox', render: renderConc }
];
function render() {
  readCfg(); let R;
  try { R = resultsAll(); } catch (e) { console.error(e); R = {}; KEYS.forEach(k => R[k] = { valid: false, val: { ok: false, errors: ['Error interno: ' + e.message] } }); }
  Components.forEach(c => { try { c.render(R); } catch (e) { console.error(c.name, e); const h = $(c.host); if (h) h.innerHTML = `<div class="banner err">El componente «${esc(c.name)}» falló: ${esc(e.message)}. Los demás componentes siguen funcionando.</div>`; } });
}
function load(k) { active = k; buildInputs(); render(); }

/* ================= MÓDULO QA ================= */
function runQA() {
  const tests = [], origActive = active;
  const add = (name, entrada, esperado, fn) => { let obt, ok; try { const r = fn(); ok = r.ok; obt = r.obt; } catch (e) { ok = false; obt = 'Excepción no controlada: ' + e.message; } tests.push({ name, entrada, esperado, obt, ok }); };
  const base = { precond: 'jacobi', omega: 1, tol: 1e-6, maxIter: 1000, method: 'cholesky', x0: new Array(N).fill(0) };
  const fresh = k => buildAb(DEFS[k]);
  add('1. Carga de escenarios', 'Ideal, Bajo estrés y Mal condicionado (datos originales)', 'Los 3 sistemas son válidos y GCP converge con residuo < 1e-6', () => {
    const rs = KEYS.map(k => { const a = fresh(k); return solveScenario(a.A, a.b, base); });
    return { ok: rs.every(r => r.valid && r.gcp.converged && r.gcp.resRel < 1e-6), obt: KEYS.map((k, i) => `${DEFS[k].label}: ${rs[i].valid ? (rs[i].gcp.converged ? 'OK, ' + rs[i].gcp.iters + ' it., res=' + f(rs[i].gcp.resRel) : 'no convergió') : 'inválido'}`).join(' | ') }; });
  add('2a. Matriz no simétrica', 'A ideal con a₁₂ = 9 y a₂₁ = 1', 'Error capturado: «A no es simétrica», sin colapsar', () => {
    const a = fresh('ideal'); a.A[0][1] = 9; a.A[1][0] = 1; const r = solveScenario(a.A, a.b, base);
    return { ok: !r.valid && r.val.errors.some(e => e.includes('simétrica')), obt: r.valid ? 'Aceptó la matriz (incorrecto)' : r.val.errors[0] }; });
  add('2b. Simétrica NO definida positiva', 'A ideal con a₁₂ = a₂₁ = 50', 'Error capturado: «no es definida positiva», sin colapsar', () => {
    const a = fresh('ideal'); a.A[0][1] = 50; a.A[1][0] = 50; const r = solveScenario(a.A, a.b, base);
    return { ok: !r.valid && r.val.errors.some(e => e.includes('definida positiva')), obt: r.valid ? 'Aceptó la matriz (incorrecto)' : r.val.errors[0].slice(0, 120) }; });
  add('2c. Celda vacía / no numérica', 'A ideal con a₃₃ = NaN', 'Notifica celdas no numéricas sin colapsar', () => {
    const a = fresh('ideal'); a.A[2][2] = NaN; const r = solveScenario(a.A, a.b, base);
    return { ok: !r.valid && r.val.errors[0].includes('no numéricas'), obt: r.valid ? 'Aceptó' : r.val.errors[0] }; });
  add('2d. Parámetros de ejecución inválidos', 'tol = −1, ω = 3, N máx. = 0, x₀₂ = NaN', 'validateConfig informa 4 problemas y devuelve valores por defecto seguros', () => {
    const v = validateConfig({ ...base, tol: -1, omega: 3, maxIter: 0, mirror: true, x0: [0, NaN, 0, 0, 0, 0] }, N);
    return { ok: v.issues.length === 4 && v.cfg.tol === 1e-6 && v.cfg.maxIter === 1000 && v.cfg.omega === 1 && v.cfg.x0[1] === 0, obt: `${v.issues.length} avisos; tol=${v.cfg.tol}, N máx.=${v.cfg.maxIter}, ω=${v.cfg.omega}` }; });
  add('3. Edición dinámica (núcleo)', 'Ideal: a₁₁ + 5; a₁₂ = a₂₁ = 2; b₂ + 100', 'La referencia cambia, A·x_ref = b y GCP coincide con la referencia (Eₚ < 0,01 %)', () => {
    const a0 = fresh('ideal'), r0 = solveScenario(a0.A, a0.b, base), a = fresh('ideal');
    a.A[0][0] += 5; a.A[0][1] = a.A[1][0] = 2; a.b[1] += 100; const r = solveScenario(a.A, a.b, base);
    const res = norm(matVec(a.A, r.xref).map((v, i) => v - a.b[i])) / norm(a.b), changed = norm(r.xref.map((v, i) => v - r0.xref[i])) > 1e-6, maxEp = Math.max(...r.Ep);
    return { ok: r.valid && changed && res < 1e-12 && maxEp < 0.01, obt: `x_ref cambió: ${changed ? 'sí' : 'no'}; residuo directo = ${f(res)}; Eₚ máx = ${f(maxEp)} %` }; });
  add('4. Comparativa de precondicionamiento', 'Escenario mal condicionado, tol 1e-6, Jacobi', 'Iteraciones GCP ≤ iteraciones GC y κ(M⁻¹A) < κ(A)', () => {
    const a = fresh('mal'), r = solveScenario(a.A, a.b, base), ic = pcg(a.A, a.b, base.x0, makePrecond(a.A, 'ic0'), base.tol, base.maxIter);
    return { ok: r.valid && r.gcp.converged && r.gcp.iters <= r.gc.iters && r.kMA < r.kA, obt: `GC = ${r.gc.iters} it. (κ(A) = ${f(r.kA)}); GCP-Jacobi = ${r.gcp.iters} it. (κ(M⁻¹A) = ${f(r.kMA)}); GCP-IC(0) = ${ic.iters} it.` }; });
  add('5. Consistencia de métodos directos', 'Cholesky vs LU vs Gauss en los 3 escenarios', 'Diferencia relativa entre soluciones < 1e-6', () => {
    let worst = 0; KEYS.forEach(k => { const a = fresh(k), c = directSolve(a.A, a.b, 'cholesky'); ['lu', 'gauss'].forEach(m => { const d = directSolve(a.A, a.b, m); worst = Math.max(worst, norm(c.map((v, i) => v - d[i])) / norm(c)); }); });
    return { ok: worst < 1e-6, obt: 'Máxima diferencia relativa = ' + f(worst) }; });
  add('6. Vector inicial x₀ ≠ 0', 'Ideal, x₀ = [100, −50, 30, 80, 10, −20]', 'Converge al mismo x que la referencia (dif. relativa < 1e-5)', () => {
    const a = fresh('ideal'), r2 = solveScenario(a.A, a.b, { ...base, x0: [100, -50, 30, 80, 10, -20] });
    const d = norm(r2.xref.map((v, i) => v - r2.gcp.x[i])) / norm(r2.xref);
    return { ok: r2.gcp.converged && d < 1e-5, obt: `converge en ${r2.gcp.iters} it., diferencia = ${f(d)}` }; });
  add('7. Coherencia financiera', 'Los 3 escenarios con datos originales', 'Σxᵢ dentro del 10 % del techo B y ninguna asignación negativa', () => {
    const o = KEYS.map(k => { const a = fresh(k), r = solveScenario(a.A, a.b, base), t = r.xref.reduce((s, v) => s + v, 0); return { k, dev: Math.abs(t - DEFS[k].B) / DEFS[k].B, neg: r.xref.some(v => v < 0) }; });
    return { ok: o.every(e => e.dev < 0.10 && !e.neg), obt: o.map(e => `${DEFS[e.k].label}: desv. ${(e.dev * 100).toFixed(2)} %${e.neg ? ' (hay negativos)' : ''}`).join(' | ') }; });
  add('8. Funcional UI: cambio de escenario', 'Clic en «Bajo estrés» (load) y regreso', 'El escenario activo y la tabla comparativa se actualizan', () => {
    load('estres'); const ok1 = $('entActivo').textContent === 'Bajo estrés' && $('tblComp').querySelectorAll('tr:first-child th').length === 4 && $('fmActivo').textContent === 'Bajo estrés';
    const obt = `entActivo = «${$('entActivo').textContent}», columnas de la tabla = ${$('tblComp').querySelectorAll('tr:first-child th').length}`; load(origActive); return { ok: ok1, obt }; });
  add('9. Funcional UI: edición de una celda', 'Escribir 2 en a₁₂ mediante un evento de entrada real', 'Se actualiza el estado, la celda espejo a₂₁ y se recalcula sin errores', () => {
    const inp = document.querySelector('input[data-i="0"][data-j="1"]'), mir = document.querySelector('input[data-i="1"][data-j="0"]');
    inp.value = '2'; inp.dispatchEvent(new Event('input', { bubbles: true }));
    const ok1 = S[active].A[0][1] === 2 && S[active].A[1][0] === 2 && mir.value === '2' && $('entMsg').textContent.includes('simétrica y definida positiva');
    const obt = `A₁₂ = ${S[active].A[0][1]}, A₂₁ = ${S[active].A[1][0]}, celda espejo = ${mir.value}`; $('btnReset').click(); return { ok: ok1, obt }; });
  add('10. Funcional UI: ruptura de simetría', 'Desactivar espejo y escribir 9 solo en a₁₂', 'El par se marca en amarillo y se muestra el aviso «no es simétrica» sin colapsar', () => {
    $('cfgMirror').checked = false; const inp = document.querySelector('input[data-i="0"][data-j="1"]');
    inp.value = '9'; inp.dispatchEvent(new Event('input', { bubbles: true }));
    const ok1 = $('entMsg').textContent.includes('no es simétrica') && inp.classList.contains('asym') && $('refBox').textContent.includes('No calculable');
    const obt = `aviso mostrado: ${$('entMsg').textContent.includes('no es simétrica') ? 'sí' : 'no'}; celda marcada: ${inp.classList.contains('asym') ? 'sí' : 'no'}; referencia: ${$('refBox').textContent.slice(0, 14)}…`;
    $('cfgMirror').checked = true; $('btnReset').click(); return { ok: ok1, obt }; });
  add('11. Aislamiento de fallos (error boundary)', 'Un componente lanza una excepción deliberada', 'Se muestra un aviso en ese componente y los demás siguen funcionando', () => {
    const c = { name: 'Prueba', host: 'refBox', render() { throw new Error('fallo simulado'); } }; Components.push(c); render(); Components.pop();
    const ok1 = $('refBox').textContent.includes('falló') && $('tblComp').querySelectorAll('tr').length > 3; const obt = `aviso: ${$('refBox').textContent.slice(0, 40)}…; tabla comparativa intacta: ${$('tblComp').querySelectorAll('tr').length > 3 ? 'sí' : 'no'}`;
    render(); return { ok: ok1, obt }; });
  load(origActive);
  const pass = tests.filter(t => t.ok).length;
  $('qaBox').innerHTML = '<table><tr><th class="l">Prueba</th><th class="l">Entrada</th><th class="l">Resultado esperado</th><th class="l">Resultado obtenido</th><th>Conclusión</th></tr>' +
    tests.map(t => `<tr><td class="l">${t.name}</td><td class="l">${esc(t.entrada)}</td><td class="l">${esc(t.esperado)}</td><td class="l m">${esc(t.obt)}</td><td>${t.ok ? '<span class="pill ok">APROBADA</span>' : '<span class="pill bad">FALLÓ</span>'}</td></tr>`).join('') + '</table>';
  $('qaSum').innerHTML = `<b class="${pass === tests.length ? 'ok' : 'bad'}">${pass}/${tests.length} pruebas aprobadas</b>`;
}

/* ================= ARRANQUE Y EVENTOS ================= */
document.addEventListener('DOMContentLoaded', () => {
  initScenarios(); renderDatos(); buildInputs(); render();
  document.body.addEventListener('input', e => {
    if (!e.target.matches('input,select')) return;
    if (e.target.dataset.i !== undefined || e.target.dataset.b !== undefined) { readCfg(); onCellInput(e); }
    render();
  });
  KEYS.forEach(k => $('bt_' + k).addEventListener('click', () => { load(k); document.getElementById('escenarios').scrollIntoView(); }));
  $('btAll').addEventListener('click', () => { render(); $('tblComp').scrollIntoView({ behavior: 'smooth', block: 'center' }); });
  $('btnRun').addEventListener('click', render);
  $('btnReset').addEventListener('click', () => { const ab = buildAb(DEFS[active]); S[active] = { A: ab.A, b: ab.b }; buildInputs(); render(); });
  $('btnQA').addEventListener('click', runQA);
  window.addEventListener('resize', () => { try { renderGCP(resultsAll()); } catch (e) { console.error(e); } });
  runQA(); // al final: las pruebas funcionales necesitan los eventos ya enlazados
});
