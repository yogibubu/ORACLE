"use strict";

const STORAGE_KEY = "wmsRotWavefunctionPayload";
const PLOTLY_DOWNLOAD_OPTIONS = { format: "png", scale: 5 };

let payload = null;
let displayLevels = [];
let selectedLevelIndex = 0;
let selectedM = 0;

const logFactorialCache = [0];
const levelsDownloadButton = document.getElementById("download-levels-png");
const waveDownloadButton = document.getElementById("download-wave-png");

function logFactorial(n) {
  if (n < 0) return NaN;
  for (let i = logFactorialCache.length; i <= n; i += 1) {
    logFactorialCache[i] = logFactorialCache[i - 1] + Math.log(i);
  }
  return logFactorialCache[n];
}

function withHighResDownloadConfig(config) {
  const baseConfig = (config && typeof config === "object") ? config : {};
  const existingOptions = (baseConfig.toImageButtonOptions && typeof baseConfig.toImageButtonOptions === "object")
    ? baseConfig.toImageButtonOptions
    : {};
  const existingScale = Number(existingOptions.scale);
  return {
    ...baseConfig,
    toImageButtonOptions: {
      format: "png",
      ...existingOptions,
      scale: Number.isFinite(existingScale) ? Math.max(existingScale, PLOTLY_DOWNLOAD_OPTIONS.scale) : PLOTLY_DOWNLOAD_OPTIONS.scale,
    },
  };
}

function sanitizePlotFilename(value, fallback = "wms_wavefunctions_plot") {
  const safe = String(value || "")
    .trim()
    .replace(/\s+/g, "_")
    .replace(/[^\w.-]+/g, "_")
    .replace(/^_+|_+$/g, "");
  return safe || fallback;
}

async function downloadPlotlyFigureAsPng(plotElement, filenameBase) {
  if (!plotElement) {
    throw new Error("Plot element not found.");
  }
  const width = Math.max(960, Math.round(plotElement.clientWidth || plotElement.offsetWidth || 960));
  const height = Math.max(540, Math.round(plotElement.clientHeight || plotElement.offsetHeight || 540));
  const exportOptions = {
    format: "png",
    filename: sanitizePlotFilename(filenameBase, "wms_wavefunctions_plot"),
    width,
    height,
    scale: PLOTLY_DOWNLOAD_OPTIONS.scale,
  };
  if (typeof Plotly.downloadImage === "function") {
    await Plotly.downloadImage(plotElement, exportOptions);
    return;
  }
  if (typeof Plotly.toImage === "function") {
    const dataUrl = await Plotly.toImage(plotElement, exportOptions);
    const link = document.createElement("a");
    link.href = dataUrl;
    link.download = `${exportOptions.filename}.png`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    return;
  }
  throw new Error("Plotly image export is not available.");
}

function setPlotDownloadButtonState(button, ready) {
  if (!button) return;
  button.disabled = !ready;
}

// Real-valued Wigner small-d using factorial-sum formula.
function wignerSmallD(J, M, K, theta) {
  const j = Math.trunc(Number(J));
  const mp = Math.trunc(Number(M));
  const m = Math.trunc(Number(K));
  if (!Number.isFinite(j) || !Number.isFinite(mp) || !Number.isFinite(m)) return 0.0;
  if (j < 0 || Math.abs(mp) > j || Math.abs(m) > j) return 0.0;

  const prefLog =
    0.5 * (
      logFactorial(j + m) + logFactorial(j - m) +
      logFactorial(j + mp) + logFactorial(j - mp)
    );

  const ch = Math.cos(theta / 2.0);
  const sh = Math.sin(theta / 2.0);

  const kMin = Math.max(0, m - mp);
  const kMax = Math.min(j + m, j - mp);
  let sum = 0.0;

  for (let k = kMin; k <= kMax; k += 1) {
    const a = j + m - k;
    const b = k;
    const c = mp - m + k;
    const d = j - mp - k;
    if (a < 0 || b < 0 || c < 0 || d < 0) continue;

    const denLog = logFactorial(a) + logFactorial(b) + logFactorial(c) + logFactorial(d);
    const coeff = Math.exp(prefLog - denLog);
    const pCos = 2 * j + m - mp - 2 * k;
    const pSin = mp - m + 2 * k;
    const phase = ((k + mp - m) % 2 === 0) ? 1.0 : -1.0;
    sum += phase * coeff * Math.pow(ch, pCos) * Math.pow(sh, pSin);
  }
  return sum;
}

function getPayloadFromStorage() {
  const raw = localStorage.getItem(STORAGE_KEY);
  if (!raw) return null;
  try {
    return JSON.parse(raw);
  } catch (err) {
    console.error("Invalid viewer payload:", err);
    return null;
  }
}

function levelEnergyCm(level) {
  return Number(level.energy_cm);
}

function formatQuantumNumber(value) {
  const num = Number(value);
  if (!Number.isFinite(num)) return "-";
  const rounded = Math.round(num);
  if (Math.abs(num - rounded) < 1e-9) return String(rounded);
  return num.toString();
}

function getQuadrupoleNucleusCount(meta) {
  const nuclei = Array.isArray(meta?.quadrupole?.nuclei) ? meta.quadrupole.nuclei : [];
  const explicitCount = Number(meta?.quadrupole?.count);
  if (Number.isFinite(explicitCount) && explicitCount >= 0) return explicitCount;
  return nuclei.length;
}

function getEffectiveNuclearSpin(level, meta) {
  const quadCount = getQuadrupoleNucleusCount(meta);
  const i12Dom = Number(level?.I12_dom);
  if (quadCount >= 2 && Number.isFinite(i12Dom) && i12Dom >= 0.0) {
    return {
      value: i12Dom,
      source: "level-I12_dom",
      label: "I12",
      html: "I<sub>12</sub>",
      description: "dominant coupled nuclear-spin label stored on this hyperfine level",
    };
  }

  const metaSpin = Number(meta?.quadrupole?.i_nuc);
  if (Number.isFinite(metaSpin) && metaSpin > 0.0) {
    return {
      value: metaSpin,
      source: "meta",
      label: "I",
      html: "I",
      description: "nuclear spin from quadrupole metadata",
    };
  }

  const inferred = inferINucFromLevels(displayLevels);
  if (Number.isFinite(inferred) && inferred > 0.0) {
    return {
      value: inferred,
      source: "inferred-from-F-manifold",
      label: quadCount >= 2 ? "I12" : "I",
      html: quadCount >= 2 ? "I<sub>12</sub>" : "I",
      description: quadCount >= 2
        ? "coupled nuclear spin inferred from the visible F/J manifold"
        : "nuclear spin inferred from the visible F/J manifold",
    };
  }

  return null;
}

function getSpinColumnLabel(meta) {
  const quadCount = getQuadrupoleNucleusCount(meta);
  if (quadCount >= 2) return "I12";
  if (quadCount === 1) return "I";
  const hasI12 = displayLevels.some((level) => Number.isFinite(Number(level?.I12_dom)));
  return hasI12 ? "I12" : "I / I12";
}

function getSpinColumnValue(level, meta) {
  const spinInfo = getEffectiveNuclearSpin(level, meta);
  if (!spinInfo || !Number.isFinite(Number(level?.F))) return "-";
  return formatQuantumNumber(spinInfo.value);
}

function levelLabel(level) {
  const bits = [`J=${level.J}`];
  if (Number.isFinite(Number(level.F))) bits.push(`F=${formatQuantumNumber(level.F)}`);
  if (Number.isFinite(Number(level.I12_dom))) bits.push(`I12=${formatQuantumNumber(level.I12_dom)}`);
  bits.push(`Ka=${level.Ka}`, `Kc=${level.Kc}`);
  if (Number.isFinite(Number(level.alpha))) bits.push(`alpha=${level.alpha}`);
  return bits.join(" ");
}

function sortLevels(levels) {
  return [...levels].sort((a, b) => {
    const e = levelEnergyCm(a) - levelEnergyCm(b);
    if (Math.abs(e) > 1e-14) return e;
    const j = Number(a.J) - Number(b.J);
    if (j !== 0) return j;
    const f = Number(a.F || 0) - Number(b.F || 0);
    if (f !== 0) return f;
    return Number(a.alpha || 0) - Number(b.alpha || 0);
  });
}

function setMetaLine() {
  const meta = payload?.meta || {};
  const mode = displayLevels.length && Number.isFinite(Number(displayLevels[0]?.F)) ? "hyperfine (F)" : "rotational";
  const generated = meta.generated_at_utc ? `, generated: ${meta.generated_at_utc}` : "";
  const trunc = meta.truncated_payload ? ", payload reduced (high J cut)" : "";
  const kInfo = getKAxisInfo(meta);
  const kTag = kInfo.abc ? `, K-axis: ${kInfo.abc}` : ", K-axis: none (spherical)";
  const hasHfs = Array.isArray(payload?.hfs_levels) && payload.hfs_levels.length > 0;
  const quadActive = Boolean(meta?.quadrupole?.active);
  const hfsWarn = quadActive && !hasHfs
    ? ", warning: quadrupole active but no F-level payload (recalculate in WMS-Rot and reload viewer)"
    : "";
  const metaWarn = (!quadActive && hasHfs)
    ? ", note: F-levels present but quadrupole metadata is incomplete (using F/J fallback for I or I12 when possible)"
    : "";
  const quadCount = getQuadrupoleNucleusCount(meta);
  const quadTag = quadCount > 0 ? `, quadrupole nuclei: ${quadCount}` : "";
  const spinTag = quadCount >= 2 ? ", cone model uses level-resolved I12 when available" : "";
  const line = `rotor: ${meta.rotor_type || "?"}, rep: ${meta.representation || "?"}, reduction: ${meta.reduction || "?"}, levels shown: ${mode}${kTag}${quadTag}${spinTag}${generated}${trunc}${hfsWarn}${metaWarn}`;
  document.getElementById("meta-line").textContent = line;
}

function setEmptyState(message) {
  document.getElementById("meta-line").textContent = message;
  document.getElementById("levels-table-body").innerHTML = "<tr><td colspan='4'>No levels available.</td></tr>";
  const spinHeader = document.getElementById("levels-spin-header");
  if (spinHeader) spinHeader.textContent = "I / I12";
  document.getElementById("selected-level").textContent = message;
  const mSelect = document.getElementById("m-select");
  if (mSelect) {
    mSelect.innerHTML = "<option value='0'>0</option>";
    mSelect.value = "0";
  }
  const mHelp = document.getElementById("m-help");
  if (mHelp) mHelp.textContent = "M selector enabled after loading levels.";
  selectedM = 0;
  Plotly.react("levels-plot", [], {
    xaxis: { visible: false },
    yaxis: { visible: false },
    annotations: [{ text: message, x: 0.5, y: 0.5, xref: "paper", yref: "paper", showarrow: false }],
    margin: { t: 20, r: 20, b: 40, l: 45 },
  }, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  Plotly.react("wave-3d", [], {
    scene: {},
    annotations: [{ text: message, x: 0.5, y: 0.5, xref: "paper", yref: "paper", showarrow: false }],
    margin: { t: 20, r: 20, b: 20, l: 20 },
  }, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadButtonState(levelsDownloadButton, false);
  setPlotDownloadButtonState(waveDownloadButton, false);
}

function syncLevelsTableHeader() {
  const spinHeader = document.getElementById("levels-spin-header");
  if (!spinHeader) return;
  spinHeader.textContent = getSpinColumnLabel(payload?.meta || {});
}

function renderLevelsTable() {
  syncLevelsTableHeader();
  const tbody = document.getElementById("levels-table-body");
  tbody.innerHTML = "";
  const meta = payload?.meta || {};
  displayLevels.forEach((level, idx) => {
    const tr = document.createElement("tr");
    tr.dataset.index = String(idx);
    if (idx === selectedLevelIndex) tr.classList.add("is-selected");
    const fVal = Number(level.F);
    const fTxt = Number.isFinite(fVal) ? formatQuantumNumber(fVal) : "-";
    const spinTxt = getSpinColumnValue(level, meta);
    tr.innerHTML = `<td>${levelLabel(level)}</td><td>${fTxt}</td><td>${spinTxt}</td><td>${levelEnergyCm(level).toFixed(7)}</td>`;
    tr.addEventListener("click", () => selectLevel(idx));
    tbody.appendChild(tr);
  });
}

function renderLevelPlot() {
  const grouped = new Map();
  displayLevels.forEach((level, idx) => {
    const J = Number(level.J);
    if (!grouped.has(J)) grouped.set(J, []);
    grouped.get(J).push({ idx, e: levelEnergyCm(level) });
  });

  const traces = [];
  const centersByIdx = new Map();
  const jVals = [...grouped.keys()].sort((a, b) => a - b);

  jVals.forEach((J) => {
    const rows = grouped.get(J).sort((a, b) => a.e - b.e);
    const n = rows.length;
    const x = [];
    const y = [];
    const custom = [];
    const text = [];
    rows.forEach((row, pos) => {
      const frac = n > 1 ? (pos / (n - 1)) - 0.5 : 0.0;
      const xc = J + 0.48 * frac;
      const half = 0.16;
      const label = levelLabel(displayLevels[row.idx]);
      x.push(xc - half, xc + half, null);
      y.push(row.e, row.e, null);
      custom.push(row.idx, row.idx, null);
      text.push(label, label, null);
      centersByIdx.set(row.idx, { x: xc, y: row.e });
    });
    const hue = (205 + (J * 17)) % 360;
    traces.push({
      type: "scatter",
      mode: "lines",
      x,
      y,
      customdata: custom,
      text,
      line: { width: 3.2, color: `hsl(${hue}, 58%, 38%)` },
      hovertemplate: "%{text}<br>E=%{y:.7f} cm^-1<extra></extra>",
      showlegend: false,
    });
  });

  const selectedCenter = centersByIdx.get(selectedLevelIndex);
  if (selectedCenter) {
    traces.push({
      type: "scatter",
      mode: "lines",
      x: [selectedCenter.x - 0.19, selectedCenter.x + 0.19],
      y: [selectedCenter.y, selectedCenter.y],
      line: { width: 6, color: "#f08c00" },
      hoverinfo: "skip",
      showlegend: false,
    });
  }

  const layout = {
    margin: { t: 20, r: 20, b: 50, l: 60 },
    xaxis: { title: "J (line-centered levels)", zeroline: false },
    yaxis: { title: "Rotational energy (cm^-1)" },
    paper_bgcolor: "rgba(255,255,255,0)",
    plot_bgcolor: "rgba(255,255,255,0.85)",
  };

  Plotly.react("levels-plot", traces, layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadButtonState(levelsDownloadButton, traces.length > 0);
  const plot = document.getElementById("levels-plot");
  if (plot.removeAllListeners) plot.removeAllListeners("plotly_click");
  plot.on("plotly_click", (ev) => {
    const idx = Number(ev.points?.[0]?.customdata);
    if (Number.isFinite(idx)) selectLevel(idx);
  });
}

function getKAxisInfo(meta) {
  const axisLabels = meta?.axis_labels || { x: "X", y: "Y", z: "Z" };
  const rotor = String(meta?.rotor_type || "").toLowerCase();
  const rep = String(meta?.representation || "Ir");
  const repMap = {
    ir: "A",
    il: "A",
    iir: "B",
    iil: "B",
    iiir: "C",
    iiil: "C",
  };
  let abc = meta?.k_axis_abc ? String(meta.k_axis_abc).toUpperCase() : null;
  if (!abc) {
    if (rotor === "linear" || rotor === "prolate") abc = "A";
    else if (rotor === "oblate") abc = "C";
    else if (rotor === "spherical") abc = null;
    else abc = repMap[rep.toLowerCase()] || String(axisLabels.z || "").toUpperCase();
  }
  let cart = null;
  if (abc) {
    for (const [cartAxis, abcAxis] of Object.entries(axisLabels)) {
      if (String(abcAxis).toUpperCase() === abc) {
        cart = cartAxis;
        break;
      }
    }
  }
  const text = abc
    ? `K is quantized on the ${abc} principal axis for this rotor/representation setup.`
    : "No unique K axis (spherical-top limit).";
  return { abc, cart, text };
}

function normalizeVec3(v) {
  const n = Math.hypot(v[0], v[1], v[2]);
  if (!Number.isFinite(n) || n <= 1e-15) return null;
  return [v[0] / n, v[1] / n, v[2] / n];
}

function matVec3(m, v) {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  ];
}

function matMul3(a, b) {
  const out = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    }
  }
  return out;
}

function dotVec3(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function crossVec3(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function quadForm3(m, v) {
  return dotVec3(v, matVec3(m, v));
}

function asAxisObj(v) {
  return { x: v[0], y: v[1], z: v[2] };
}

function axisObjToVec(a) {
  return [Number(a?.x || 0), Number(a?.y || 0), Number(a?.z || 0)];
}

// Dominant eigenvector of a real symmetric 3x3 matrix by power iteration.
function dominantAxisSymmetric3(matrix, useAbsSpectrum = false) {
  const mat = useAbsSpectrum ? matMul3(matrix, matrix) : matrix;
  let v = normalizeVec3([1.0, 0.7, 0.3]);
  if (!v) return null;
  for (let it = 0; it < 42; it += 1) {
    const w = matVec3(mat, v);
    const vn = normalizeVec3(w);
    if (!vn) return null;
    const d = Math.hypot(vn[0] - v[0], vn[1] - v[1], vn[2] - v[2]);
    const df = Math.hypot(vn[0] + v[0], vn[1] + v[1], vn[2] + v[2]);
    v = vn;
    if (Math.min(d, df) < 1e-11) break;
  }
  const av = matVec3(matrix, v);
  const lambda = v[0] * av[0] + v[1] * av[1] + v[2] * av[2];
  return { x: v[0], y: v[1], z: v[2], eigenvalue: lambda };
}

function labFrameFromProbability(prob, nTheta, nPhi) {
  let wSum = 0.0;
  let sxx = 0.0;
  let syy = 0.0;
  let szz = 0.0;
  let sxy = 0.0;
  let sxz = 0.0;
  let syz = 0.0;

  for (let it = 0; it <= nTheta; it += 1) {
    const theta = Math.PI * (it / nTheta);
    const st = Math.sin(theta);
    const ct = Math.cos(theta);
    if (Math.abs(st) < 1e-15) continue;
    for (let ip = 0; ip < nPhi; ip += 1) {
      const p = Number(prob[it][ip] || 0.0);
      if (!Number.isFinite(p) || p <= 0.0) continue;
      const phi = 2.0 * Math.PI * (ip / nPhi);
      const ux = st * Math.cos(phi);
      const uy = st * Math.sin(phi);
      const uz = ct;
      const w = p * st;
      wSum += w;
      sxx += w * ux * ux;
      syy += w * uy * uy;
      szz += w * uz * uz;
      sxy += w * ux * uy;
      sxz += w * ux * uz;
      syz += w * uy * uz;
    }
  }
  if (wSum <= 1e-15) return null;

  const moment = [
    [sxx / wSum, sxy / wSum, sxz / wSum],
    [sxy / wSum, syy / wSum, syz / wSum],
    [sxz / wSum, syz / wSum, szz / wSum],
  ];
  const tr = (moment[0][0] + moment[1][1] + moment[2][2]) / 3.0;
  const devNorm = Math.hypot(
    moment[0][0] - tr,
    moment[1][1] - tr,
    moment[2][2] - tr,
    moment[0][1],
    moment[0][2],
    moment[1][2]
  );
  if (devNorm < 1e-6) return null;
  const zPack = dominantAxisSymmetric3(moment, false);
  if (!zPack) return null;
  const zVec = normalizeVec3([zPack.x, zPack.y, zPack.z]);
  if (!zVec) return null;

  const helper = Math.abs(zVec[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
  const proj = dotVec3(helper, zVec);
  const e2 = normalizeVec3([
    helper[0] - proj * zVec[0],
    helper[1] - proj * zVec[1],
    helper[2] - proj * zVec[2],
  ]);
  if (!e2) return null;
  const e3 = normalizeVec3(crossVec3(zVec, e2));
  if (!e3) return null;

  const a = quadForm3(moment, e2);
  const b = dotVec3(e2, matVec3(moment, e3));
  const c = quadForm3(moment, e3);
  const disc = Math.sqrt(Math.max(0, (a - c) * (a - c) + 4.0 * b * b));
  const lambdaX = 0.5 * (a + c + disc);
  const lambdaY = 0.5 * (a + c - disc);

  let p = 1.0;
  let q = 0.0;
  if (Math.abs(b) > 1e-14 || Math.abs(lambdaX - a) > 1e-14) {
    p = b;
    q = lambdaX - a;
    const pqNorm = Math.hypot(p, q);
    if (pqNorm > 1e-15) {
      p /= pqNorm;
      q /= pqNorm;
    } else {
      p = 1.0;
      q = 0.0;
    }
  } else if (a < c) {
    p = 0.0;
    q = 1.0;
  }

  const xVec = normalizeVec3([
    p * e2[0] + q * e3[0],
    p * e2[1] + q * e3[1],
    p * e2[2] + q * e3[2],
  ]);
  if (!xVec) return null;
  let yVec = normalizeVec3(crossVec3(zVec, xVec));
  if (!yVec) return null;
  if (dotVec3(crossVec3(xVec, yVec), zVec) < 0.0) {
    yVec = [-yVec[0], -yVec[1], -yVec[2]];
  }

  return {
    xAxis: asAxisObj(xVec),
    yAxis: asAxisObj(yVec),
    zAxis: asAxisObj(zVec),
    eigenvalues: {
      x: lambdaX,
      y: lambdaY,
      z: zPack.eigenvalue,
    },
    moment,
  };
}

function hyperfineConeGeometry(level, state, meta) {
  const F = Number(level?.F);
  const J = Number(state?.J ?? level?.J);
  const spinInfo = getEffectiveNuclearSpin(level, meta);
  const iNuc = Number(spinInfo?.value);
  const iNucSource = spinInfo?.source || "missing";
  if (!Number.isFinite(F) || !Number.isFinite(J) || !Number.isFinite(iNuc)) return null;
  if (F <= 0.0 || iNuc <= 0.0) return null;

  const jAxis = rmsRotationAxis(state);
  if (!jAxis) return null;
  const j2 = J * (J + 1.0);
  const i2 = iNuc * (iNuc + 1.0);
  const ff = F * (F + 1.0);
  if (ff <= 1e-12) return null;

  if (F < Math.abs(J - iNuc) - 1e-9 || F > (J + iNuc) + 1e-9) return null;

  const cJ = (ff + j2 - i2) / (2.0 * ff);
  const cI = (ff + i2 - j2) / (2.0 * ff);
  if (!Number.isFinite(cJ) || !Number.isFinite(cI)) return null;

  const denomJF = 2.0 * Math.sqrt(Math.max(1e-16, ff * j2));
  const denomIF = 2.0 * Math.sqrt(Math.max(1e-16, ff * i2));
  let cosJF = (ff + j2 - i2) / denomJF;
  let cosIF = (ff + i2 - j2) / denomIF;
  cosJF = Math.max(-1.0, Math.min(1.0, cosJF));
  cosIF = Math.max(-1.0, Math.min(1.0, cosIF));
  const thetaJF = Math.acos(cosJF);
  const thetaIF = Math.acos(cosIF);

  const fSign = cJ >= 0.0 ? 1.0 : -1.0;
  const fAxis = {
    x: fSign * jAxis.x,
    y: fSign * jAxis.y,
    z: fSign * jAxis.z,
  };

  const fVec = axisObjToVec(fAxis);
  const helper = Math.abs(fVec[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
  const proj = dotVec3(helper, fVec);
  const e1 = normalizeVec3([
    helper[0] - proj * fVec[0],
    helper[1] - proj * fVec[1],
    helper[2] - proj * fVec[2],
  ]);
  if (!e1) return null;
  const e2 = normalizeVec3(crossVec3(fVec, e1));
  if (!e2) return null;

  const jRef = normalizeVec3([
    Math.cos(thetaJF) * fVec[0] + Math.sin(thetaJF) * e1[0],
    Math.cos(thetaJF) * fVec[1] + Math.sin(thetaJF) * e1[1],
    Math.cos(thetaJF) * fVec[2] + Math.sin(thetaJF) * e1[2],
  ]);
  const iRef = normalizeVec3([
    Math.cos(thetaIF) * fVec[0] - Math.sin(thetaIF) * e1[0],
    Math.cos(thetaIF) * fVec[1] - Math.sin(thetaIF) * e1[1],
    Math.cos(thetaIF) * fVec[2] - Math.sin(thetaIF) * e1[2],
  ]);
  if (!jRef || !iRef) return null;

  return {
    fAxis,
    jAxis,
    jRef: asAxisObj(jRef),
    iRef: asAxisObj(iRef),
    e1: asAxisObj(e1),
    e2: asAxisObj(e2),
    thetaJF,
    thetaIF,
    cosJF,
    cosIF,
    cJ,
    cI,
    F,
    J,
    iNuc,
    iNucSource,
    spinLabel: spinInfo?.label || "I",
    spinLabelHtml: spinInfo?.html || "I",
    spinDescription: spinInfo?.description || "nuclear spin",
  };
}

function inferINucFromLevels(levels) {
  if (!Array.isArray(levels) || !levels.length) return null;
  let maxDiff = 0.0;
  let found = false;
  levels.forEach((lvl) => {
    const f = Number(lvl?.F);
    const j = Number(lvl?.J);
    if (!Number.isFinite(f) || !Number.isFinite(j)) return;
    found = true;
    const d = Math.abs(f - j);
    if (d > maxDiff) maxDiff = d;
  });
  if (!found || maxDiff <= 1e-9) return null;
  return Math.round(maxDiff * 2.0) / 2.0;
}

function getSelectedFrameMode() {
  const labRadio = document.getElementById("frame-lab");
  if (labRadio && labRadio.checked) return "lab";
  return "body";
}

function bodyToLabVec(vBody, labFrame) {
  if (!labFrame) return [vBody[0], vBody[1], vBody[2]];
  const xLabBody = axisObjToVec(labFrame.xAxis);
  const yLabBody = axisObjToVec(labFrame.yAxis);
  const zLabBody = axisObjToVec(labFrame.zAxis);
  return [
    dotVec3(vBody, xLabBody),
    dotVec3(vBody, yLabBody),
    dotVec3(vBody, zLabBody),
  ];
}

function transformVecForFrame(vBody, frameMode, labFrame) {
  if (frameMode === "lab") return bodyToLabVec(vBody, labFrame);
  return [vBody[0], vBody[1], vBody[2]];
}

function transformAxisForFrame(axisObj, frameMode, labFrame) {
  const v = transformVecForFrame(axisObjToVec(axisObj), frameMode, labFrame);
  return asAxisObj(v);
}

function transformGridForFrame(xGrid, yGrid, zGrid, frameMode, labFrame) {
  if (frameMode !== "lab" || !labFrame) {
    return { x: xGrid, y: yGrid, z: zGrid };
  }
  const xOut = [];
  const yOut = [];
  const zOut = [];
  for (let i = 0; i < xGrid.length; i += 1) {
    const rowX = [];
    const rowY = [];
    const rowZ = [];
    for (let j = 0; j < xGrid[i].length; j += 1) {
      const v = bodyToLabVec([xGrid[i][j], yGrid[i][j], zGrid[i][j]], labFrame);
      rowX.push(v[0]);
      rowY.push(v[1]);
      rowZ.push(v[2]);
    }
    xOut.push(rowX);
    yOut.push(rowY);
    zOut.push(rowZ);
  }
  return { x: xOut, y: yOut, z: zOut };
}

function axisMagnitudeLabel(axis, axisLabels) {
  return `${axisLabels.x}:${Math.abs(axis.x).toFixed(3)}, ${axisLabels.y}:${Math.abs(axis.y).toFixed(3)}, ${axisLabels.z}:${Math.abs(axis.z).toFixed(3)}`;
}

function buildPrecessionConeTraces(axisObj, e1Obj, e2Obj, halfAngle, color, name, rMax = 1.3) {
  const n = axisObjToVec(axisObj);
  const e1 = axisObjToVec(e1Obj);
  const e2 = axisObjToVec(e2Obj);
  const cth = Math.cos(halfAngle);
  const sth = Math.sin(halfAngle);
  if (!Number.isFinite(cth) || !Number.isFinite(sth)) return [];
  if (Math.abs(sth) < 1e-4) return [];

  const nR = 14;
  const nPhi = 52;
  const traces = [];
  for (const sign of [1, -1]) {
    const x = [];
    const y = [];
    const z = [];
    const c = [];
    for (let ir = 0; ir <= nR; ir += 1) {
      const rho = (rMax * ir) / nR;
      const rowX = [];
      const rowY = [];
      const rowZ = [];
      const rowC = [];
      for (let ip = 0; ip <= nPhi; ip += 1) {
        const phi = (2.0 * Math.PI * ip) / nPhi;
        const circ = [
          Math.cos(phi) * e1[0] + Math.sin(phi) * e2[0],
          Math.cos(phi) * e1[1] + Math.sin(phi) * e2[1],
          Math.cos(phi) * e1[2] + Math.sin(phi) * e2[2],
        ];
        const dir = [
          sign * cth * n[0] + sth * circ[0],
          sign * cth * n[1] + sth * circ[1],
          sign * cth * n[2] + sth * circ[2],
        ];
        rowX.push(rho * dir[0]);
        rowY.push(rho * dir[1]);
        rowZ.push(rho * dir[2]);
        rowC.push(0.0);
      }
      x.push(rowX);
      y.push(rowY);
      z.push(rowZ);
      c.push(rowC);
    }
    traces.push({
      type: "surface",
      x,
      y,
      z,
      surfacecolor: c,
      cmin: 0,
      cmax: 1,
      colorscale: [[0, color], [1, color]],
      opacity: 0.22,
      showscale: false,
      hovertemplate: `${name} cone, half-angle=${(halfAngle * 180.0 / Math.PI).toFixed(2)} deg<extra></extra>`,
      name: `${name} cone`,
      legendgroup: `${name}-cone`,
      showlegend: sign > 0,
    });
  }
  return traces;
}

function syncMSelector(state) {
  const select = document.getElementById("m-select");
  const help = document.getElementById("m-help");
  const J = Math.max(0, Math.trunc(Number(state?.J || 0)));

  if (!Number.isInteger(selectedM) || Math.abs(selectedM) > J) selectedM = 0;

  select.innerHTML = "";
  for (let m = -J; m <= J; m += 1) {
    const opt = document.createElement("option");
    opt.value = String(m);
    opt.textContent = String(m);
    select.appendChild(opt);
  }

  select.value = String(selectedM);
  if (select.value !== String(selectedM)) {
    selectedM = 0;
    select.value = "0";
  }
  help.textContent = `M = projection of J on Z_lab (allowed: ${-J} to ${J})`;
}

function buildWaveSurface(state, gamma, mProj) {
  const J = Number(state.J);
  const kValues = state.k_values || [];
  const coeffRe = state.coeff_re || [];
  const coeffIm = state.coeff_im || [];
  const nTheta = 52;
  const nPhi = 104;

  const prob = Array.from({ length: nTheta + 1 }, () => Array(nPhi + 1).fill(0));
  let maxProb = 0;

  for (let it = 0; it <= nTheta; it += 1) {
    const theta = Math.PI * (it / nTheta);
    for (let ip = 0; ip <= nPhi; ip += 1) {
      const phi = 2.0 * Math.PI * (ip / nPhi);
      let re = 0.0;
      let im = 0.0;
      for (let i = 0; i < kValues.length; i += 1) {
        const K = Number(kValues[i]);
        const dMK = wignerSmallD(J, mProj, K, theta);
        const baseRe = dMK * Math.cos(K * phi);
        const baseIm = dMK * Math.sin(K * phi);
        const cr = Number(coeffRe[i] || 0.0);
        const ci = Number(coeffIm[i] || 0.0);
        re += cr * baseRe - ci * baseIm;
        im += cr * baseIm + ci * baseRe;
      }
      const p = re * re + im * im;
      prob[it][ip] = p;
      if (p > maxProb) maxProb = p;
    }
  }

  const x = [];
  const y = [];
  const z = [];
  const surfaceColor = [];
  if (maxProb <= 0) maxProb = 1.0;

  for (let it = 0; it <= nTheta; it += 1) {
    const theta = Math.PI * (it / nTheta);
    const rowX = [];
    const rowY = [];
    const rowZ = [];
    const rowC = [];
    for (let ip = 0; ip <= nPhi; ip += 1) {
      const phi = 2.0 * Math.PI * (ip / nPhi);
      const pNorm = prob[it][ip] / maxProb;
      const r = Math.pow(Math.max(0, pNorm), gamma);
      rowX.push(r * Math.sin(theta) * Math.cos(phi));
      rowY.push(r * Math.sin(theta) * Math.sin(phi));
      rowZ.push(r * Math.cos(theta));
      rowC.push(pNorm);
    }
    x.push(rowX);
    y.push(rowY);
    z.push(rowZ);
    surfaceColor.push(rowC);
  }
  const labFrame = labFrameFromProbability(prob, nTheta, nPhi);
  return {
    x,
    y,
    z,
    surfaceColor,
    maxProb,
    labFrame: labFrame || null,
  };
}

function rmsRotationAxis(state) {
  const J = Number(state.J);
  const denom = J * (J + 1);
  if (denom <= 1e-12) return null;
  const c = state.coord_expect || {};
  const vx = Math.sqrt(Math.max(0, Number(c.x || 0) / denom));
  const vy = Math.sqrt(Math.max(0, Number(c.y || 0) / denom));
  const vz = Math.sqrt(Math.max(0, Number(c.z || 0) / denom));
  const n = Math.hypot(vx, vy, vz);
  if (n <= 1e-12) return null;
  return { x: vx / n, y: vy / n, z: vz / n };
}

function renderSelectedInfo(level, state, mProj, surface, hyperfineGeom, frameMode) {
  const panel = document.getElementById("selected-level");
  const purity = Number(level.purity);
  const purityTxt = Number.isFinite(purity) ? `, purity=${purity.toFixed(4)}` : "";
  const axisLabels = payload?.meta?.axis_labels || { x: "X", y: "Y", z: "Z" };
  const plotAxisLabels = frameMode === "lab"
    ? { x: "X_lab", y: "Y_lab", z: "Z_lab" }
    : axisLabels;
  const frameTxt = frameMode === "lab"
    ? "Current 3D frame: laboratory (X_lab/Y_lab/Z_lab)."
    : "Current 3D frame: body-fixed molecular axes (A/B/C).";
  const kInfo = getKAxisInfo(payload?.meta || {});
  const coords = state?.coord_expect || {};
  const coordTxt = `<J^2> (${axisLabels.x}, ${axisLabels.y}, ${axisLabels.z}) = (${Number(coords.x || 0).toFixed(4)}, ${Number(coords.y || 0).toFixed(4)}, ${Number(coords.z || 0).toFixed(4)})`;
  const J = Number(state?.J || 0);
  const denom = J * (J + 1);
  let axisShare = "Second-moment fractions unavailable (J=0).";
  let interpTxt = `Displayed quantity: probability map for the direction of Z_lab in the molecular A/B/C frame, using M=${mProj}.`;
  let rmsDegTxt = "";
  let labAxisTxt = "Lab-frame principal axes unavailable (near-isotropic angular map).";
  let labAxisDegTxt = "";
  let quadAxisTxt = "Hyperfine cone model unavailable (requires level with F and usable I or I12 data).";
  let spinLevelTxt = "";
  let hfsWarnTxt = "";
  const hasHfs = Array.isArray(payload?.hfs_levels) && payload.hfs_levels.length > 0;
  const quadActive = Boolean(payload?.meta?.quadrupole?.active);
  if (quadActive && !hasHfs) {
    hfsWarnTxt = "Quadrupole constants are active in metadata, but this payload has no resolved F levels. Recompute spectrum in WMS-Rot and reopen/reload this viewer to refresh localStorage payload.";
  }
  if (denom > 1e-12) {
    const fx = Number(coords.x || 0) / denom;
    const fy = Number(coords.y || 0) / denom;
    const fz = Number(coords.z || 0) / denom;
    axisShare = `Second-moment fractions <J_${axisLabels.x}^2>, <J_${axisLabels.y}^2>, <J_${axisLabels.z}^2> over J(J+1): ${fx.toFixed(3)} / ${fy.toFixed(3)} / ${fz.toFixed(3)}.`;
    const compsDesc = [
      { axis: axisLabels.x, val: fx },
      { axis: axisLabels.y, val: fy },
      { axis: axisLabels.z, val: fz },
    ].sort((a, b) => b.val - a.val);
    interpTxt =
      `Displayed quantity: probability map for the direction of Z_lab in the molecular A/B/C frame, using M=${mProj}. ` +
      `A longer lobe along ${compsDesc[0].axis} means Z_lab is more often found near +/-${compsDesc[0].axis}. ` +
      `This is not atomic/electronic density and not a direct arrow of J.`;
    const nearTie01 = Math.abs(compsDesc[0].val - compsDesc[1].val) < 0.06;
    const nearTie12 = Math.abs(compsDesc[1].val - compsDesc[2].val) < 0.06;
    if (nearTie01 || nearTie12) {
      rmsDegTxt =
        "Note: two <J_i^2> components are nearly equal, so the orange RMS direction is not unique; the shown line is a plotting convention (square-root second moments).";
    }
  }
  interpTxt += ` ${frameTxt}`;
  if (surface?.labFrame) {
    if (frameMode === "body") {
      const lx = surface.labFrame.xAxis;
      const ly = surface.labFrame.yAxis;
      const lz = surface.labFrame.zAxis;
      labAxisTxt =
        `Black/gray lab triad from eigvecs(Q_lab): ` +
        `X_lab ${axisMagnitudeLabel(lx, axisLabels)}; ` +
        `Y_lab ${axisMagnitudeLabel(ly, axisLabels)}; ` +
        `Z_lab ${axisMagnitudeLabel(lz, axisLabels)}.`;
      const ex = Number(surface.labFrame.eigenvalues?.x);
      const ey = Number(surface.labFrame.eigenvalues?.y);
      if (Number.isFinite(ex) && Number.isFinite(ey) && Math.abs(ex - ey) < 0.02) {
        labAxisDegTxt = "Note: lambda_X and lambda_Y are close, so X_lab/Y_lab can rotate without changing probability significantly.";
      }
    } else {
      const bx = asAxisObj(bodyToLabVec([1, 0, 0], surface.labFrame));
      const by = asAxisObj(bodyToLabVec([0, 1, 0], surface.labFrame));
      const bz = asAxisObj(bodyToLabVec([0, 0, 1], surface.labFrame));
      labAxisTxt =
        `Black/gray body-axis triad in lab frame: ` +
        `${axisLabels.x} ${axisMagnitudeLabel(bx, plotAxisLabels)}; ` +
        `${axisLabels.y} ${axisMagnitudeLabel(by, plotAxisLabels)}; ` +
        `${axisLabels.z} ${axisMagnitudeLabel(bz, plotAxisLabels)}.`;
    }
  } else if (frameMode === "lab") {
    labAxisTxt = "Lab frame requested but unavailable: angular map is near-isotropic, so principal lab triad cannot be resolved numerically.";
  }
  if (Number.isFinite(Number(level?.F))) {
    const spinInfo = getEffectiveNuclearSpin(level, payload?.meta || {});
    if (spinInfo) {
      if (spinInfo.source === "level-I12_dom") {
        spinLevelTxt = `Dominant coupled nuclear-spin label for this selected F level: ${spinInfo.html}=${formatQuantumNumber(spinInfo.value)}.`;
      } else if (spinInfo.source === "meta") {
        spinLevelTxt = `Nuclear spin used in the cone model: ${spinInfo.html}=${formatQuantumNumber(spinInfo.value)}.`;
      } else {
        spinLevelTxt = `Fallback ${spinInfo.description}: ${spinInfo.html}=${formatQuantumNumber(spinInfo.value)}.`;
      }
    }
  }
  if (hyperfineGeom) {
    const tJ = (hyperfineGeom.thetaJF * 180.0 / Math.PI).toFixed(2);
    const tI = (hyperfineGeom.thetaIF * 180.0 / Math.PI).toFixed(2);
    const srcTxt = hyperfineGeom.iNucSource === "level-I12_dom"
      ? " [using dominant I12 stored on this hyperfine level]"
      : hyperfineGeom.iNucSource === "meta"
        ? ""
        : ` [${hyperfineGeom.spinLabel} inferred from the visible F/J manifold because explicit metadata is incomplete]`;
    const useLabProjection = frameMode === "lab" && Boolean(surface?.labFrame);
    const fAxisShown = useLabProjection
      ? transformAxisForFrame(hyperfineGeom.fAxis, frameMode, surface?.labFrame || null)
      : hyperfineGeom.fAxis;
    const fAxisLabels = useLabProjection ? plotAxisLabels : axisLabels;
    quadAxisTxt =
      `Hyperfine vector-coupling cones (semiclassical): ` +
      `F-axis ${axisMagnitudeLabel(fAxisShown, fAxisLabels)}, ` +
      `theta(J,F)=${tJ} deg, theta(${hyperfineGeom.spinLabel},F)=${tI} deg, ` +
      `c_J=${hyperfineGeom.cJ.toFixed(4)}, c_${hyperfineGeom.spinLabel}=${hyperfineGeom.cI.toFixed(4)} ` +
      `for (F=${formatQuantumNumber(hyperfineGeom.F)}, J=${formatQuantumNumber(hyperfineGeom.J)}, ${hyperfineGeom.spinLabel}=${formatQuantumNumber(hyperfineGeom.iNuc)}).${srcTxt}`;
  } else if (Number.isFinite(Number(level?.F)) && Number(state?.J) === 0) {
    quadAxisTxt =
      "F is defined for this level, but the cone construction is intentionally omitted at J=0 because the rotational direction from <J_i^2>/J(J+1) is undefined.";
  }
  panel.innerHTML = `
    <strong>${levelLabel(level)}</strong><br>
    E = ${levelEnergyCm(level).toFixed(7)} cm^-1 (${Number(level.energy_mhz || 0).toFixed(4)} MHz)${purityTxt}<br>
    M = ${mProj} (projection of J on Z_lab)<br>
    ${coordTxt}<br>
    ${axisShare}<br>
    ${interpTxt}<br>
    ${labAxisTxt}<br>
    ${labAxisDegTxt ? `${labAxisDegTxt}<br>` : ""}
    ${hfsWarnTxt ? `${hfsWarnTxt}<br>` : ""}
    ${spinLevelTxt ? `${spinLevelTxt}<br>` : ""}
    ${quadAxisTxt}<br>
    Note: J, I/I12, F graphics represent expectation-value/semi-classical geometry, not sharp simultaneous vector observables in a single quantum eigenstate.<br>
    ${rmsDegTxt ? `${rmsDegTxt}<br>` : ""}
    ${kInfo.text}
  `;
}

function renderWave3D() {
  const level = displayLevels[selectedLevelIndex];
  if (!level) return;
  const parentKey = level.parent_key || level.id;
  const state = payload?.states?.[parentKey];
  if (!state) {
    setEmptyState("Missing parent wavefunction coefficients.");
    return;
  }
  syncMSelector(state);
  selectedM = Number.parseInt(document.getElementById("m-select").value, 10);
  if (!Number.isFinite(selectedM)) selectedM = 0;

  const gamma = Number(document.getElementById("gamma-slider").value);
  document.getElementById("gamma-value").textContent = gamma.toFixed(2);
  const surface = buildWaveSurface(state, gamma, selectedM);
  const frameMode = getSelectedFrameMode();
  const axisLabels = payload?.meta?.axis_labels || { x: "X", y: "Y", z: "Z" };
  const plotAxisLabels = frameMode === "lab"
    ? { x: "X_lab", y: "Y_lab", z: "Z_lab" }
    : axisLabels;
  const hyperfineGeom = hyperfineConeGeometry(level, state, payload?.meta || {});
  const toggleFAxis = document.getElementById("toggle-f-axis");
  const toggleICone = document.getElementById("toggle-i-cone");
  const showFAxis = toggleFAxis ? Boolean(toggleFAxis.checked) : true;
  const showICone = toggleICone ? Boolean(toggleICone.checked) : true;
  renderSelectedInfo(level, state, selectedM, surface, hyperfineGeom, frameMode);
  const kInfo = getKAxisInfo(payload?.meta || {});
  const xLabel = frameMode === "body"
    ? `${axisLabels.x}${kInfo.abc === String(axisLabels.x).toUpperCase() ? " (K)" : ""}`
    : plotAxisLabels.x;
  const yLabel = frameMode === "body"
    ? `${axisLabels.y}${kInfo.abc === String(axisLabels.y).toUpperCase() ? " (K)" : ""}`
    : plotAxisLabels.y;
  const zLabel = frameMode === "body"
    ? `${axisLabels.z}${kInfo.abc === String(axisLabels.z).toUpperCase() ? " (K)" : ""}`
    : plotAxisLabels.z;
  const surfacePlot = transformGridForFrame(surface.x, surface.y, surface.z, frameMode, surface.labFrame || null);

  const traces = [
    {
      type: "surface",
      x: surfacePlot.x,
      y: surfacePlot.y,
      z: surfacePlot.z,
      surfacecolor: surface.surfaceColor,
      cmin: 0,
      cmax: 1,
      colorscale: "YlGnBu",
      opacity: 0.92,
      showscale: true,
      colorbar: { title: `|Psi|^2 (norm.), M=${selectedM}` },
      hovertemplate: "x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}<extra></extra>",
    },
    {
      type: "scatter3d",
      mode: "lines",
      x: [-1.35, 1.35],
      y: [0, 0],
      z: [0, 0],
      line: { color: "#d1495b", width: 6 },
      name: xLabel,
      hoverinfo: "skip",
    },
    {
      type: "scatter3d",
      mode: "lines",
      x: [0, 0],
      y: [-1.35, 1.35],
      z: [0, 0],
      line: { color: "#2f9e44", width: 6 },
      name: yLabel,
      hoverinfo: "skip",
    },
    {
      type: "scatter3d",
      mode: "lines",
      x: [0, 0],
      y: [0, 0],
      z: [-1.35, 1.35],
      line: { color: "#1c7ed6", width: 6 },
      name: zLabel,
      hoverinfo: "skip",
    },
    {
      type: "scatter3d",
      mode: "text",
      x: [1.44, 0, 0],
      y: [0, 1.44, 0],
      z: [0, 0, 1.44],
      text: [xLabel, yLabel, zLabel],
      textfont: { size: 16, color: "#102a43" },
      hoverinfo: "skip",
      showlegend: false,
    },
  ];

  if (kInfo.cart) {
    const kBodyBasis = {
      x: [1, 0, 0],
      y: [0, 1, 0],
      z: [0, 0, 1],
    };
    const kBody = kBodyBasis[kInfo.cart];
    const kPlot = transformVecForFrame(kBody, frameMode, surface.labFrame || null);
    traces.push({
      type: "scatter3d",
      mode: "lines",
      x: [-1.45 * kPlot[0], 1.45 * kPlot[0]],
      y: [-1.45 * kPlot[1], 1.45 * kPlot[1]],
      z: [-1.45 * kPlot[2], 1.45 * kPlot[2]],
      line: { color: "#0ea5a5", width: 10 },
      name: frameMode === "lab" ? `K-axis (${kInfo.abc}, in lab)` : `K-axis (${kInfo.abc})`,
      hovertemplate: `K-axis = ${kInfo.abc}<extra></extra>`,
    });
  }

  const rmsAxis = rmsRotationAxis(state);
  if (rmsAxis) {
    const rmsPlot = transformVecForFrame(axisObjToVec(rmsAxis), frameMode, surface.labFrame || null);
    traces.push({
      type: "scatter3d",
      mode: "lines",
      x: [-rmsPlot[0], rmsPlot[0]],
      y: [-rmsPlot[1], rmsPlot[1]],
      z: [-rmsPlot[2], rmsPlot[2]],
      line: { color: "#f08c00", width: 8 },
      name: "RMS second-moment axis",
      hovertemplate: "Convention from <J_i^2> second moments<extra></extra>",
    });
  }

  if (surface.labFrame) {
    if (frameMode === "body") {
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-surface.labFrame.xAxis.x, surface.labFrame.xAxis.x],
        y: [-surface.labFrame.xAxis.y, surface.labFrame.xAxis.y],
        z: [-surface.labFrame.xAxis.z, surface.labFrame.xAxis.z],
        line: { color: "#111111", width: 8 },
        name: "X_lab principal axis",
        hovertemplate: "From eigvecs(Q_lab): X_lab<extra></extra>",
      });
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-surface.labFrame.yAxis.x, surface.labFrame.yAxis.x],
        y: [-surface.labFrame.yAxis.y, surface.labFrame.yAxis.y],
        z: [-surface.labFrame.yAxis.z, surface.labFrame.yAxis.z],
        line: { color: "#4d4d4d", width: 8 },
        name: "Y_lab principal axis",
        hovertemplate: "From eigvecs(Q_lab): Y_lab<extra></extra>",
      });
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-surface.labFrame.zAxis.x, surface.labFrame.zAxis.x],
        y: [-surface.labFrame.zAxis.y, surface.labFrame.zAxis.y],
        z: [-surface.labFrame.zAxis.z, surface.labFrame.zAxis.z],
        line: { color: "#888888", width: 8 },
        name: "Z_lab principal axis",
        hovertemplate: "From eigvecs(Q_lab): Z_lab<extra></extra>",
      });
    } else {
      const bodyXLab = bodyToLabVec([1, 0, 0], surface.labFrame);
      const bodyYLab = bodyToLabVec([0, 1, 0], surface.labFrame);
      const bodyZLab = bodyToLabVec([0, 0, 1], surface.labFrame);
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-bodyXLab[0], bodyXLab[0]],
        y: [-bodyXLab[1], bodyXLab[1]],
        z: [-bodyXLab[2], bodyXLab[2]],
        line: { color: "#111111", width: 8 },
        name: `${axisLabels.x} body axis`,
        hovertemplate: `${axisLabels.x} body axis in lab frame<extra></extra>`,
      });
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-bodyYLab[0], bodyYLab[0]],
        y: [-bodyYLab[1], bodyYLab[1]],
        z: [-bodyYLab[2], bodyYLab[2]],
        line: { color: "#4d4d4d", width: 8 },
        name: `${axisLabels.y} body axis`,
        hovertemplate: `${axisLabels.y} body axis in lab frame<extra></extra>`,
      });
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-bodyZLab[0], bodyZLab[0]],
        y: [-bodyZLab[1], bodyZLab[1]],
        z: [-bodyZLab[2], bodyZLab[2]],
        line: { color: "#888888", width: 8 },
        name: `${axisLabels.z} body axis`,
        hovertemplate: `${axisLabels.z} body axis in lab frame<extra></extra>`,
      });
    }
  }

  if (hyperfineGeom) {
    const fAxisPlot = transformAxisForFrame(hyperfineGeom.fAxis, frameMode, surface.labFrame || null);
    const e1Plot = transformAxisForFrame(hyperfineGeom.e1, frameMode, surface.labFrame || null);
    const e2Plot = transformAxisForFrame(hyperfineGeom.e2, frameMode, surface.labFrame || null);
    const jRefPlot = transformAxisForFrame(hyperfineGeom.jRef, frameMode, surface.labFrame || null);
    const iRefPlot = transformAxisForFrame(hyperfineGeom.iRef, frameMode, surface.labFrame || null);
    const spinConeName = hyperfineGeom.spinLabel === "I12" ? "I12" : "I";
    if (showFAxis) {
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [-fAxisPlot.x, fAxisPlot.x],
        y: [-fAxisPlot.y, fAxisPlot.y],
        z: [-fAxisPlot.z, fAxisPlot.z],
        line: { color: "#0b3d91", width: 9 },
        name: "F-axis (effective)",
        hovertemplate: "Effective F-axis used for semiclassical cone model<extra></extra>",
      });
    }
    traces.push(...buildPrecessionConeTraces(fAxisPlot, e1Plot, e2Plot, hyperfineGeom.thetaJF, "#f08c00", "J", 1.3));
    if (showICone) {
      traces.push(...buildPrecessionConeTraces(fAxisPlot, e1Plot, e2Plot, hyperfineGeom.thetaIF, "#2a9d8f", spinConeName, 1.3));
    }

    traces.push({
      type: "scatter3d",
      mode: "lines",
      x: [0, 1.3 * jRefPlot.x],
      y: [0, 1.3 * jRefPlot.y],
      z: [0, 1.3 * jRefPlot.z],
      line: { color: "#f08c00", width: 6 },
      name: "J reference generator",
      hovertemplate: "One generator on J precession cone<extra></extra>",
    });
    if (showICone) {
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: [0, 1.3 * iRefPlot.x],
        y: [0, 1.3 * iRefPlot.y],
        z: [0, 1.3 * iRefPlot.z],
        line: { color: "#2a9d8f", width: 6 },
        name: `${spinConeName} reference generator`,
        hovertemplate: `One generator on ${spinConeName} precession cone<extra></extra>`,
      });
    }
  }

  const layout = {
    margin: { t: 20, r: 20, b: 20, l: 20 },
    scene: {
      aspectmode: "cube",
      xaxis: { title: plotAxisLabels.x, range: [-1.5, 1.5] },
      yaxis: { title: plotAxisLabels.y, range: [-1.5, 1.5] },
      zaxis: { title: plotAxisLabels.z, range: [-1.5, 1.5] },
      camera: { eye: { x: 1.45, y: 1.45, z: 1.2 } },
    },
    paper_bgcolor: "rgba(255,255,255,0)",
    legend: { orientation: "h", y: -0.05 },
  };
  Plotly.react("wave-3d", traces, layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadButtonState(waveDownloadButton, traces.length > 0);
}

function selectLevel(index) {
  if (!Number.isInteger(index) || index < 0 || index >= displayLevels.length) return;
  selectedLevelIndex = index;
  renderLevelsTable();
  renderWave3D();
}

function buildDisplayLevels() {
  const hfs = Array.isArray(payload?.hfs_levels) ? payload.hfs_levels : [];
  const rot = Array.isArray(payload?.rot_levels) ? payload.rot_levels : [];
  if (hfs.length > 0) return sortLevels(hfs);
  return sortLevels(rot);
}

function loadViewer() {
  payload = getPayloadFromStorage();
  if (!payload) {
    setEmptyState("No payload found. Run WMS-Rot and click Open Wavefunction Viewer.");
    return;
  }
  displayLevels = buildDisplayLevels();
  if (!displayLevels.length) {
    setEmptyState("Payload loaded but contains no levels.");
    return;
  }
  selectedLevelIndex = 0;
  syncLevelsTableHeader();
  selectedM = 0;
  setMetaLine();
  renderLevelPlot();
  renderLevelsTable();
  renderWave3D();
  if (window.MathJax && window.MathJax.typesetPromise) {
    window.MathJax.typesetPromise().catch(() => {});
  }
}

document.getElementById("reload-data").addEventListener("click", loadViewer);
document.getElementById("open-main").addEventListener("click", () => {
  window.open("../index.html", "_blank");
});
document.getElementById("gamma-slider").addEventListener("input", () => {
  if (!displayLevels.length) return;
  renderWave3D();
});
document.getElementById("m-select").addEventListener("change", (ev) => {
  if (!displayLevels.length) return;
  selectedM = Number.parseInt(ev.target.value, 10);
  if (!Number.isFinite(selectedM)) selectedM = 0;
  renderWave3D();
});
document.getElementById("toggle-f-axis").addEventListener("change", () => {
  if (!displayLevels.length) return;
  renderWave3D();
});
document.getElementById("toggle-i-cone").addEventListener("change", () => {
  if (!displayLevels.length) return;
  renderWave3D();
});
document.getElementById("frame-body").addEventListener("change", () => {
  if (!displayLevels.length) return;
  renderWave3D();
});
document.getElementById("frame-lab").addEventListener("change", () => {
  if (!displayLevels.length) return;
  renderWave3D();
});
if (levelsDownloadButton) {
  levelsDownloadButton.addEventListener("click", async () => {
    levelsDownloadButton.disabled = true;
    try {
      await downloadPlotlyFigureAsPng(
        document.getElementById("levels-plot"),
        `wms_wavefunctions_levels_${payload?.meta?.representation || "plot"}`
      );
    } catch (error) {
      console.error("Unable to download levels PNG:", error);
      alert(error.message || "Unable to download the high-resolution PNG.");
    } finally {
      setPlotDownloadButtonState(levelsDownloadButton, displayLevels.length > 0);
    }
  });
}
if (waveDownloadButton) {
  waveDownloadButton.addEventListener("click", async () => {
    waveDownloadButton.disabled = true;
    try {
      const level = displayLevels[selectedLevelIndex];
      const label = level ? levelLabel(level) : "map";
      await downloadPlotlyFigureAsPng(
        document.getElementById("wave-3d"),
        `wms_wavefunctions_${label}`
      );
    } catch (error) {
      console.error("Unable to download wavefunction PNG:", error);
      alert(error.message || "Unable to download the high-resolution PNG.");
    } finally {
      setPlotDownloadButtonState(waveDownloadButton, displayLevels.length > 0);
    }
  });
}

loadViewer();
