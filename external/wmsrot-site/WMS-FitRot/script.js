const uploadButton = document.getElementById('upload-button');
const uploadInput = document.getElementById('upload-input');
const uploadedTabs = document.getElementById('uploaded-tabs');
const uploadedPanels = document.getElementById('uploaded-panels');
const overlayList = document.getElementById('overlay-list');
const overlayPlot = document.getElementById('overlay-plot');
const overlayDownloadButton = document.getElementById('overlay-download-png');
const overlayToast = document.getElementById('overlay-toast');
const autofitRestoreInitialButton = document.getElementById('autofit-restore-initial');
const autofitRunButton = document.getElementById('autofit-run');
const autofitRunLocalButton = document.getElementById('autofit-run-local');
const autofitImportSpfitButton = document.getElementById('autofit-import-spfit');
const autofitImportInput = document.getElementById('autofit-import-input');
const autofitDownloadLogButton = document.getElementById('autofit-download-log');
const autofitUploadLogButton = document.getElementById('autofit-upload-log');
const autofitUploadLogInput = document.getElementById('autofit-upload-log-input');
const autofitCancelButton = document.getElementById('autofit-cancel');
const autofitTheorySelect = document.getElementById('autofit-theory-select');
const autofitExpSelect = document.getElementById('autofit-exp-select');
const autofitExpMode = document.getElementById('autofit-exp-mode');
const autofitUseInputAssignment = document.getElementById('autofit-use-input-assignment');
const autofitLabelList = document.getElementById('autofit-label-list');
const autofitControlList = document.getElementById('autofit-control-list');
const autofitManualList = document.getElementById('autofit-manual-list');
const autofitManualTitle = document.getElementById('autofit-manual-title');
const autofitManualHelp = document.getElementById('autofit-manual-help');
const autofitLabelInput = document.getElementById('autofit-label-input');
const autofitAddLabel = document.getElementById('autofit-add-label');
const autofitClearLabels = document.getElementById('autofit-clear-labels');
const autofitClickMode = document.getElementById('autofit-click-mode');
const autofitAssignmentMode = document.getElementById('autofit-assignment-mode');
const autofitClickModeField = document.getElementById('autofit-click-mode-field');
const autofitStatus = document.getElementById('autofit-status');
const autofitOutput = document.getElementById('autofit-output');
const autofitProgressBar = document.getElementById('autofit-progress-bar');
const autofitDownloads = document.getElementById('autofit-downloads');
const sensitivityPanel = document.getElementById('sensitivity-panel');
const sensitivityDashboard = document.getElementById('sensitivity-dashboard');
const sensitivityHideWeakToggle = document.getElementById('sensitivity-hide-weak');
const sensitivityParameterTable = document.getElementById('sensitivity-parameter-table');
const sensitivityCorrelationHeatmap = document.getElementById('sensitivity-correlation-heatmap');
const sensitivityRankingChart = document.getElementById('sensitivity-ranking-chart');
const sensitivityParameterDetail = document.getElementById('sensitivity-parameter-detail');
const sensitivityDiagnostics = document.getElementById('sensitivity-diagnostics');
const sensitivitySuggestions = document.getElementById('sensitivity-suggestions');
const autofitDiagnostics = document.getElementById('autofit-diagnostics');
const spfitImportReport = document.getElementById('spfit-import-report');
const autofitWarning = document.getElementById('autofit-warning');
const autofitControlBlock = document.getElementById('autofit-control-block');
const autofitManualBlock = document.getElementById('autofit-manual-block');
const autofitManualPending = document.getElementById('autofit-manual-pending');
const autofitParamBoundGrid = document.getElementById('autofit-param-bound-grid');
const autofitFitQuarticScale = document.getElementById('autofit-fit-s');
const autofitQuarticScaleSlider = document.getElementById('autofit-quartic-scale');
const autofitQuarticScaleValue = document.getElementById('autofit-quartic-scale-value');
const autofitQuarticScaleFactor = document.getElementById('autofit-quartic-scale-factor');
const overlayClickSummary = document.getElementById('overlay-click-summary');
const overlayAutoFigure = document.getElementById('overlay-figure-auto');
const overlayManualFigure = document.getElementById('overlay-figure-manual');
const overlayTargetFit = document.getElementById('overlay-target-fit');
const overlayTargetControl = document.getElementById('overlay-target-control');
const overlayManualStatus = document.getElementById('overlay-manual-status');
const overlayStepTheory = document.getElementById('overlay-step-theory');
const overlayStepExperimental = document.getElementById('overlay-step-experimental');
const overlayStepPair = document.getElementById('overlay-step-pair');
const autofitManualDisabledOptionGroups = Array.from(document.querySelectorAll('[data-manual-disabled-group]'));

const autofitOptions = {
  deltaAssign: document.getElementById('autofit-delta-assign'),
  deltaMatch: document.getElementById('autofit-delta-match'),
  kmax: document.getElementById('autofit-kmax'),
  maxCombos: document.getElementById('autofit-max-combos'),
  jmax: document.getElementById('autofit-jmax'),
  intensityCut: document.getElementById('autofit-intensity-cut'),
  span: document.getElementById('autofit-span'),
  finalSpan: document.getElementById('autofit-final-span'),
  mctrl: document.getElementById('autofit-mctrl'),
  lambdaRms: document.getElementById('autofit-lambda-rms'),
  lambdaBias: document.getElementById('autofit-lambda-bias'),
  kappaCol: document.getElementById('autofit-kappa-col'),
  usePowell: document.getElementById('autofit-use-powell'),
  sampleTemp: document.getElementById('autofit-sample-temp'),
  beamSize: document.getElementById('autofit-beam-size'),
  neighborK: document.getElementById('autofit-neighbor-k'),
  lsLoss: document.getElementById('autofit-ls-loss'),
  lsFScaleMode: document.getElementById('autofit-ls-f-scale-mode'),
  lsFScale: document.getElementById('autofit-ls-f-scale'),
  lsMaxNfev: document.getElementById('autofit-ls-max-nfev'),
};
const QUARTIC_PARAM_KEYS = ['DJ', 'DJK', 'DK', 'dJ', 'dK'];
const SEXTIC_PARAM_KEYS = ['HJ', 'HJK', 'HKJ', 'HK', 'h1', 'h2', 'h3'];
const DISTORTION_PARAM_KEYS = [...QUARTIC_PARAM_KEYS, ...SEXTIC_PARAM_KEYS];
const MODEL_PARAM_KEYS = ['A', 'B', 'C', ...DISTORTION_PARAM_KEYS];
const FIT_SCALE_PARAM_KEY = 's';
const BOUND_PARAM_KEYS = [...MODEL_PARAM_KEYS, FIT_SCALE_PARAM_KEY];
const SPFIT_DEFAULT_LINE_ERROR_MHZ = 0.05;
const SPFIT_DEFAULT_LINE_WEIGHT = 1;
const SPFIT_FIXED_PARAM_ERROR = 1e-37;
const SPFIT_PARAM_INCLUDE_EPS = 1e-12;
const SPFIT_HEADER_LINE_WINDOW = 50;
const SPFIT_DEFAULT_REPRESENTATION = 'Ir';
const SPFIT_FIT_ERROR_FLOORS_MHZ = Object.freeze({
  rotational: 1e-3,
  quartic: 1e-8,
  sextic: 1e-10,
});
const SPFIT_FIT_ERROR_RELATIVE = Object.freeze({
  rotational: 1e-3,
  quartic: 0.25,
  sextic: 0.5,
});
const FIT_STATUS_LOG_TYPE = 'wms-fitrot-status-log';
const FIT_STATUS_LOG_VERSION = 2;
const AUTOFIT_DIAGNOSTIC_HELP = {
  fittedLines: {
    title: 'Fitted Lines',
    bodyHtml: '<p>Number of assigned transitions used in the final fit residual vector. If <code>r_i = f_obs,i - f_calc,i</code>, this card reports the count <code>N = |{r_i}|</code>.</p>',
  },
  exportableLines: {
    title: 'Exportable Lines',
    bodyHtml: '<p>Number of fitted lines that pass the diagnostic gate and can be written to the SPFIT export. In automatic mode this is the count of lines classified as <code>solid</code>; in manual mode it is the count of non-rejected lines.</p>',
  },
  fitParameters: {
    title: 'Fit Parameters',
    bodyHtml: '<p>Dimension <code>p</code> of the fitted parameter vector <code>theta</code>. In this workflow <code>theta</code> can contain <code>A</code>, <code>B</code>, <code>C</code>, and any enabled distortion constants.</p>',
  },
  rmsLabels: {
    title: 'RMS(labels)',
    bodyHtml: '<p>Root-mean-square of the final fitted residuals, with <code>r_i = f_obs,i - f_calc,i</code> in MHz: <code>RMS = sqrt((1/N) sum_i r_i^2)</code>.</p>',
  },
  controlScore: {
    title: 'Control Score',
    bodyHtml: '<p>Ranking metric computed on control-line residuals after the local fit: <code>score = n - lambda_rms * MAD - lambda_bias * |bias| - kappa_col * collisions</code>.</p><p>Here <code>n</code> is the number of matched control lines, <code>MAD = median(|r_j - median(r)|)</code>, <code>bias = mean(r_j)</code>, and <code>collisions</code> counts reused experimental peaks.</p>',
  },
  peakQuality: {
    title: 'Peak Quality q',
    bodyHtml: '<p>Median peak-quality value across fitted assignments. If <code>q</code> is supplied in the peak list, that value is used directly.</p><p>Otherwise WMS-FitRot estimates <code>q = clip((I - median(I)) / (5 s) + 1, 0.1, 5)</code>, where <code>I</code> is peak intensity and <code>s</code> is a MAD-based intensity scale with standard-deviation fallback.</p>',
  },
  exportSigma: {
    title: 'Export Sigma',
    bodyHtml: `<p>Frequency uncertainty <code>sigma_export</code> assigned to every exported SPFIT observation in the generated <code>.lin</code> file. It is currently fixed at <code>${SPFIT_DEFAULT_LINE_ERROR_MHZ.toFixed(4)} MHz</code>.</p>`,
  },
  ambiguousReject: {
    title: 'Ambiguous / Reject',
    bodyHtml: '<p>Counts reported as <code>N_ambiguous / N_reject</code>. A line becomes ambiguous when a warning-level rule is triggered, such as elevated <code>|r|</code>, small alternative-gap margin, or moderate <code>q</code>.</p><p>A line is rejected when a hard gate is hit, such as missing residual, very large <code>|r|</code>, or <code>q &lt; 0.35</code>.</p>',
  },
};
const DISTORTION_BOUNDS_FLOORS_MHZ = {
  DJ: 0.05,
  DJK: 0.05,
  DK: 0.05,
  dJ: 0.01,
  dK: 0.01,
  HJ: 0.001,
  HJK: 0.001,
  HKJ: 0.001,
  HK: 0.001,
  h1: 0.001,
  h2: 0.001,
  h3: 0.001,
};
const PLOTLY_DOWNLOAD_OPTIONS = { format: 'png', scale: 5 };
const SENSITIVITY_TABLE_COLUMNS = [
  {
    key: 'name',
    label: 'Parameter',
    tooltip: 'Fitted parameter name. This is the stable ordering used to index the Jacobian, covariance, and correlation matrices.',
    type: 'string',
  },
  {
    key: 'value',
    label: 'Value',
    tooltip: 'Final fitted parameter value.',
    type: 'number',
  },
  {
    key: 'weighted_l2_sensitivity',
    label: 'Sensitivity',
    tooltip: 'Weighted local sensitivity ||J_col||. Larger values mean the residual vector changes more strongly when this parameter changes.',
    type: 'number',
  },
  {
    key: 'relative_weighted_l2_sensitivity',
    label: 'Relative Sensitivity',
    tooltip: 'Scale-aware local sensitivity |p| * ||J_col||. Useful when parameters have different natural magnitudes.',
    type: 'number',
  },
  {
    key: 'sigma',
    label: 'Sigma',
    tooltip: 'Formal 1σ uncertainty derived from the covariance matrix diagonal.',
    type: 'number',
  },
  {
    key: 'max_abs_correlation',
    label: 'Max Correlation',
    tooltip: 'Largest absolute off-diagonal correlation involving this parameter. Values near 1 indicate poor separability from another parameter.',
    type: 'number',
  },
  {
    key: 'status',
    label: 'Status',
    tooltip: 'Weak means low Jacobian sensitivity. Correlated means high local correlation with another fitted parameter. Good means neither warning triggered.',
    type: 'string',
  },
];
const SENSITIVITY_COND_THRESHOLDS = {
  good: 1e3,
  warn: 1e6,
};
const SENSITIVITY_DEFAULT_SORT = {
  key: 'weighted_l2_sensitivity',
  direction: 'desc',
};
const MODEL_PARAM_LABELS = {
  A: 'A',
  B: 'B',
  C: 'C',
  DJ: 'DJ',
  DJK: 'DJK',
  DK: 'DK',
  dJ: 'dJ / d1',
  dK: 'dK / d2',
  HJ: 'HJ',
  HJK: 'HJK',
  HKJ: 'HKJ',
  HK: 'HK',
  h1: 'h1',
  h2: 'h2',
  h3: 'h3',
  s: 's',
};
const autofitFitParamInputs = Object.fromEntries(
  DISTORTION_PARAM_KEYS.map((key) => [key, document.getElementById(`autofit-fit-${key}`)])
);
const autofitParamBoundInputs = new Map();
let cachedAutofitInputModelParamsText = null;
let cachedAutofitInputModelParams = null;
let quarticScaleRecalcTimer = null;

const spectra = new Map();
const overlaySettings = new Map();

const DEFAULT_COLORS = [
  '#1f77b4',
  '#ff7f0e',
  '#2ca02c',
  '#d62728',
  '#9467bd',
  '#8c564b',
  '#e377c2',
  '#7f7f7f',
  '#bcbd22',
  '#17becf'
];
let colorIndex = 0;
let spectrumCounter = 0;
let activeSpectrumId = null;
const queueKey = 'fitSpectraQueue';

let autofitInputText = null;
let autofitInputName = null;
const selectedFitLabels = [];
const selectedControlLabels = [];
let autofitWorker = null;
let autofitRunInProgress = false;
let overlayToastTimer = null;
const overlayHighlightsFit = new Map();
const overlayHighlightsControl = new Map();
const manualAssignments = [];
let pendingManualSelection = null;
let lastAutofitSummary = null;
let lastAutofitContext = null;
let lastAutofitDiagnostics = null;
let lastAutofitFitMetrics = null;
let pendingSpfitImportState = null;
let lastSpfitImportSummary = null;
const autofitSensitivityState = {
  hoveredParam: null,
  selectedParam: null,
  selectedParams: [],
  tableSortKey: SENSITIVITY_DEFAULT_SORT.key,
  tableSortDirection: SENSITIVITY_DEFAULT_SORT.direction,
  hideWeak: false,
};
const MANUAL_VISIBLE_OPTION_KEYS = new Set([
  'jmax',
  'intensityCut',
  'span',
  'finalSpan',
  'usePowell',
  'lsLoss',
  'lsFScaleMode',
  'lsFScale',
  'lsMaxNfev',
]);

function getNextColor() {
  const color = DEFAULT_COLORS[colorIndex % DEFAULT_COLORS.length];
  colorIndex += 1;
  return color;
}

function normalizeIntensityMode(value) {
  if (typeof value !== 'string') return 'absorbance';
  const normalized = value.trim().toLowerCase();
  return normalized === 'trasmittance' || normalized === 'transmittance'
    ? 'transmittance'
    : 'absorbance';
}

function applyIntensityMode(values, mode) {
  const sign = normalizeIntensityMode(mode) === 'transmittance' ? -1 : 1;
  return values.map((value) => {
    const num = Number(value);
    if (!Number.isFinite(num)) return Number.NaN;
    return sign * Math.abs(num);
  });
}

function computeMedian(values) {
  if (!Array.isArray(values) || !values.length) {
    return Number.NaN;
  }
  const sorted = values.slice().sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  if (sorted.length % 2 === 0) {
    return (sorted[mid - 1] + sorted[mid]) / 2;
  }
  return sorted[mid];
}

function computeMad(values, median) {
  if (!Array.isArray(values) || !values.length) {
    return Number.NaN;
  }
  const diffs = values.map((value) => Math.abs(value - median));
  return computeMedian(diffs);
}

function autoPeakThreshold(values) {
  const median = computeMedian(values);
  if (!Number.isFinite(median)) return Number.NaN;
  const mad = computeMad(values, median);
  if (!Number.isFinite(mad)) return Number.NaN;
  return median + 5 * mad;
}

function quadraticVertex(x0, y0, x1, y1, x2, y2) {
  // Fit in local coordinates around x1 to avoid catastrophic cancellation when
  // frequencies are large (GHz) but point spacing is tiny (kHz).
  const u0 = x0 - x1;
  const u2 = x2 - x1;
  if (![u0, u2, y0, y1, y2].every(Number.isFinite)) return null;
  if (u0 === 0 || u2 === 0 || u0 === u2) return null;

  const d0 = y0 - y1;
  const d2 = y2 - y1;
  const det = u0 * u2 * (u0 - u2);
  if (!Number.isFinite(det) || det === 0) return null;

  const a = ((d0 * u2) - (d2 * u0)) / det;
  const b = (((u0 * u0) * d2) - ((u2 * u2) * d0)) / det;
  if (!Number.isFinite(a) || !Number.isFinite(b) || a === 0) {
    return null;
  }

  const uv = -b / (2 * a);
  const xv = x1 + uv;
  const yv = y1 + (a * uv * uv) + (b * uv);
  if (!Number.isFinite(xv) || !Number.isFinite(yv)) return null;
  if (yv < 0) return null;
  const minX = Math.min(x0, x2);
  const maxX = Math.max(x0, x2);
  if (xv < minX || xv > maxX) return null;
  return { x: xv, y: yv };
}

function computeDatStickSpectrum(freqs, intensities, threshold = null, mode = 'absorbance') {
  const length = Math.min(Array.isArray(freqs) ? freqs.length : 0, Array.isArray(intensities) ? intensities.length : 0);
  if (length < 3) return { freqs: [], ints: [], baseInts: [] };

  const x = new Array(length);
  const y = new Array(length);
  const finiteValues = [];
  for (let i = 0; i < length; i++) {
    const xi = Number(freqs[i]);
    const yi = Number(intensities[i]);
    x[i] = xi;
    if (!Number.isFinite(yi)) {
      y[i] = Number.NaN;
      continue;
    }
    const absVal = Math.abs(yi);
    y[i] = absVal;
    if (Number.isFinite(absVal)) {
      finiteValues.push(absVal);
    }
  }

  if (finiteValues.length < 3) {
    return { freqs: [], ints: [], baseInts: [] };
  }

  let minHeight = null;
  if (threshold == null || !Number.isFinite(Number(threshold))) {
    const auto = autoPeakThreshold(finiteValues);
    if (Number.isFinite(auto)) minHeight = auto;
  } else {
    minHeight = Math.abs(Number(threshold));
  }

  const candidates = [];
  for (let i = 1; i < length - 1; i++) {
    const y0 = y[i - 1];
    const y1 = y[i];
    const y2 = y[i + 1];
    if (!Number.isFinite(y0) || !Number.isFinite(y1) || !Number.isFinite(y2)) continue;
    const dyLeft = y1 - y0;
    const dyRight = y2 - y1;
    if (!(dyLeft > 0 && dyRight <= 0)) continue;
    const d2 = y0 - 2 * y1 + y2;
    if (!(d2 < 0)) continue;
    if (minHeight != null && y1 < minHeight) continue;
    const xi = x[i];
    if (!Number.isFinite(xi)) continue;
    candidates.push({ index: i, x: xi, y: y1 });
  }

  if (!candidates.length) {
    return { freqs: [], ints: [], baseInts: [] };
  }

  const refined = [];
  for (const candidate of candidates) {
    const idx = candidate.index;
    if (idx <= 0 || idx >= length - 1) {
      refined.push({ x: candidate.x, y: candidate.y });
      continue;
    }
    const x0 = x[idx - 1];
    const x1 = x[idx];
    const x2 = x[idx + 1];
    const y0 = y[idx - 1];
    const y1 = y[idx];
    const y2 = y[idx + 1];
    if (![x0, x1, x2, y0, y1, y2].every(Number.isFinite)) {
      refined.push({ x: candidate.x, y: candidate.y });
      continue;
    }
    const vertex = quadraticVertex(x0, y0, x1, y1, x2, y2);
    refined.push(vertex || { x: candidate.x, y: candidate.y });
  }

  refined.sort((a, b) => a.x - b.x);
  const baseInts = refined.map((peak) => peak.y);
  const orientedPeaks = applyIntensityMode(baseInts, mode);
  return {
    freqs: refined.map((peak) => peak.x),
    ints: orientedPeaks,
    baseInts,
  };
}

function parseStickDatRows(text) {
  const rows = [];
  if (typeof text !== 'string' || !text.trim()) return rows;
  const numericTokenPattern = /^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?$/;
  const parseNumericToken = (token) => {
    if (!numericTokenPattern.test(String(token || '').trim())) return Number.NaN;
    const value = Number(token);
    return Number.isFinite(value) ? value : Number.NaN;
  };
  String(text).split(/\r?\n/).forEach((rawLine) => {
    const line = String(rawLine || '').trim();
    if (!line || line.startsWith('#') || line.startsWith('!') || line.startsWith(';')) return;
    const parts = line.split(/\s+/);
    if (parts.length < 2) return;
    const freq = parseNumericToken(parts[0]);
    const intensity = parseNumericToken(parts[1]);
    if (!Number.isFinite(freq) || !Number.isFinite(intensity)) return;
    let index = 2;
    let snr = Number.NaN;
    let q = Number.NaN;
    if (index < parts.length) {
      const parsed = parseNumericToken(parts[index]);
      if (Number.isFinite(parsed)) {
        snr = parsed;
        index += 1;
      }
    }
    if (index < parts.length) {
      const parsed = parseNumericToken(parts[index]);
      if (Number.isFinite(parsed)) {
        q = parsed;
        index += 1;
      }
    }
    let label = parts.slice(index).join(' ').trim();
    if ((label.startsWith('"') && label.endsWith('"')) || (label.startsWith("'") && label.endsWith("'"))) {
      label = label.slice(1, -1).trim();
    }
    rows.push({
      freq,
      intensity,
      snr,
      q,
      label: label || null,
    });
  });
  rows.sort((a, b) => a.freq - b.freq);
  return rows;
}

function parseDatSpectrum(text) {
  const pairs = parseStickDatRows(text);
  if (!pairs.length) return { freqs: [], intensities: [], labels: [] };
  pairs.sort((a, b) => a.freq - b.freq);
  return {
    freqs: pairs.map((pair) => pair.freq),
    intensities: pairs.map((pair) => pair.intensity),
    labels: pairs.map((pair) => pair.label),
  };
}

function normalizeValues(values) {
  const finite = values.filter((v) => Number.isFinite(v));
  if (!finite.length) return values.map(() => Number.NaN);
  const min = Math.min(...finite);
  const max = Math.max(...finite);
  const range = max - min;
  if (!Number.isFinite(range) || range === 0) {
    return values.map(() => 0);
  }
  return values.map((value) => {
    if (!Number.isFinite(value)) return Number.NaN;
    return (value - min) / range;
  });
}

function withHighResDownloadConfig(config) {
  const baseConfig = (config && typeof config === 'object') ? config : {};
  const existingOptions = (baseConfig.toImageButtonOptions && typeof baseConfig.toImageButtonOptions === 'object')
    ? baseConfig.toImageButtonOptions
    : {};
  const existingScale = Number(existingOptions.scale);
  return {
    ...baseConfig,
    toImageButtonOptions: {
      format: 'png',
      ...existingOptions,
      scale: Number.isFinite(existingScale) ? Math.max(existingScale, PLOTLY_DOWNLOAD_OPTIONS.scale) : PLOTLY_DOWNLOAD_OPTIONS.scale,
    },
  };
}

function sanitizePlotFilename(value, fallback = 'wms_fitrot_plot') {
  const safe = String(value || '')
    .trim()
    .replace(/\s+/g, '_')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return safe || fallback;
}

async function downloadPlotlyFigureAsPng(plotElement, filenameBase) {
  if (!plotElement) {
    throw new Error('Plot element not found.');
  }
  const width = Math.max(960, Math.round(plotElement.clientWidth || plotElement.offsetWidth || 960));
  const height = Math.max(540, Math.round(plotElement.clientHeight || plotElement.offsetHeight || 540));
  const exportOptions = {
    format: 'png',
    filename: sanitizePlotFilename(filenameBase, 'wms_fitrot_plot'),
    width,
    height,
    scale: PLOTLY_DOWNLOAD_OPTIONS.scale,
  };
  if (typeof Plotly.downloadImage === 'function') {
    await Plotly.downloadImage(plotElement, exportOptions);
    return;
  }
  if (typeof Plotly.toImage === 'function') {
    const dataUrl = await Plotly.toImage(plotElement, exportOptions);
    const link = document.createElement('a');
    link.href = dataUrl;
    link.download = `${exportOptions.filename}.png`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    return;
  }
  throw new Error('Plotly image export is not available.');
}

const plotDownloadControls = new WeakMap();

function registerPlotDownloadControl(plotElement, button, filenameResolver) {
  if (!plotElement || !button) return;
  plotDownloadControls.set(plotElement, { button, filenameResolver });
  plotElement.dataset.plotReady = 'false';
  button.disabled = true;
  button.addEventListener('click', async () => {
    if (button.disabled) return;
    button.disabled = true;
    try {
      const filename = typeof filenameResolver === 'function' ? filenameResolver() : filenameResolver;
      await downloadPlotlyFigureAsPng(plotElement, filename);
    } catch (error) {
      console.error('Unable to download Plotly PNG:', error);
      alert(error.message || 'Unable to download the high-resolution PNG.');
    } finally {
      button.disabled = plotElement.dataset.plotReady !== 'true';
    }
  });
}

function setPlotDownloadReady(plotElement, ready) {
  if (!plotElement) return;
  plotElement.dataset.plotReady = ready ? 'true' : 'false';
  const control = plotDownloadControls.get(plotElement);
  if (control?.button) {
    control.button.disabled = !ready;
  }
}

if (overlayPlot && overlayDownloadButton) {
  registerPlotDownloadControl(overlayPlot, overlayDownloadButton, () => 'wms_fitrot_overlay');
}

function buildStickTrace(freqs, intensities, color, name, labels = []) {
  const x = [];
  const y = [];
  const text = [];
  const customdata = [];
  const length = Math.min(freqs.length, intensities.length);
  for (let i = 0; i < length; i++) {
    const fx = Number(freqs[i]);
    const iy = Number(intensities[i]);
    if (!Number.isFinite(fx) || !Number.isFinite(iy)) continue;
    const label = labels[i] ? String(labels[i]) : '';
    x.push(fx, fx, null);
    y.push(0, iy, null);
    text.push(label, label, null);
    customdata.push(label, label, null);
  }
  return {
    x,
    y,
    type: 'scatter',
    mode: 'lines',
    name,
    text,
    customdata,
    hovertemplate: text.length ? '%{x}<br>%{text}<extra></extra>' : undefined,
    line: { color, width: 1.3 },
  };
}

function buildLineTrace(freqs, intensities, color, name, labels = []) {
  const hasLabels = Array.isArray(labels) && labels.length === freqs.length;
  return {
    x: freqs,
    y: intensities,
    type: 'scatter',
    mode: 'lines',
    name,
    text: hasLabels ? labels : undefined,
    hovertemplate: hasLabels ? '%{x}<br>%{text}<extra></extra>' : undefined,
    line: { color, width: 1.5 },
  };
}

function renderEmptyPlot(target, message) {
  const layout = {
    xaxis: { title: 'Frequency (MHz)' },
    yaxis: { title: 'Intensity' },
    uirevision: 'overlay-v1',
    annotations: [
      {
        text: message,
        xref: 'paper',
        yref: 'paper',
        x: 0.5,
        y: 0.5,
        showarrow: false,
        font: { size: 14, color: '#4a5568' },
      }
    ],
    margin: { l: 70, r: 30, t: 30, b: 60 },
  };
  Plotly.react(target, [], layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadReady(target, false);
}

function showOverlayToast(message) {
  if (!overlayToast) return;
  overlayToast.textContent = message;
  overlayToast.hidden = false;
  overlayToast.setAttribute('aria-hidden', 'false');
  overlayToast.classList.add('is-visible');
  clearTimeout(overlayToastTimer);
  overlayToastTimer = setTimeout(() => {
    overlayToast.classList.remove('is-visible');
    overlayToast.setAttribute('aria-hidden', 'true');
    setTimeout(() => {
      if (!overlayToast.classList.contains('is-visible')) {
        overlayToast.hidden = true;
      }
    }, 220);
  }, 2400);
}

function setAutofitWarningVisible(visible) {
  if (!autofitWarning) return;
  autofitWarning.hidden = !visible;
}


class DatModeDialog {
  constructor(root) {
    this.root = root;
    this.form = root?.querySelector('[data-role="dialog-form"]') || null;
    this.fileNameEl = root?.querySelector('[data-role="file-name"]') || null;
    this.confirmBtn = root?.querySelector('[data-role="dialog-confirm"]') || null;
    this.cancelBtn = root?.querySelector('[data-role="dialog-cancel"]') || null;
    this.backdrop = root?.querySelector('[data-role="dialog-backdrop"]') || null;
    this.activeResolver = null;
    this.boundKeydown = (event) => this.handleKeydown(event);
    this.registerEvents();
  }

  registerEvents() {
    if (this.confirmBtn) {
      this.confirmBtn.addEventListener('click', () => this.confirm());
    }
    if (this.cancelBtn) {
      this.cancelBtn.addEventListener('click', () => this.cancel());
    }
    if (this.backdrop) {
      this.backdrop.addEventListener('click', () => this.cancel());
    }
  }

  open({ fileName } = {}) {
    if (!this.root) return Promise.resolve('profile');
    if (this.fileNameEl) {
      this.fileNameEl.textContent = fileName || 'selected';
    }
    if (this.form) {
      const defaultInput = this.form.querySelector('input[value="profile"]');
      if (defaultInput) defaultInput.checked = true;
    }
    this.root.hidden = false;
    this.root.setAttribute('aria-hidden', 'false');
    document.addEventListener('keydown', this.boundKeydown, { capture: true });
    return new Promise((resolve) => {
      this.activeResolver = resolve;
    });
  }

  close() {
    if (!this.root) return;
    this.root.hidden = true;
    this.root.setAttribute('aria-hidden', 'true');
    document.removeEventListener('keydown', this.boundKeydown, { capture: true });
  }

  confirm() {
    this.resolve(this.getSelectedMode());
  }

  cancel() {
    this.resolve(null);
  }

  resolve(value) {
    const resolver = this.activeResolver;
    this.activeResolver = null;
    this.close();
    if (typeof resolver === 'function') resolver(value);
  }

  getSelectedMode() {
    if (!this.form) return 'profile';
    const checked = this.form.querySelector('input[name="dat-file-mode"]:checked');
    return checked?.value === 'stick' ? 'stick' : 'profile';
  }

  handleKeydown(event) {
    if (event.key === 'Escape') {
      event.preventDefault();
      this.cancel();
    }
  }
}

const datModeDialog = new DatModeDialog(document.getElementById('dat-mode-dialog'));

function getSpectrumDisplayData(spectrum, view) {
  if (!spectrum) return { freqs: [], intensities: [] };
  if (view === 'stick') {
    const stick = spectrum.stick || { freqs: [], baseInts: [], labels: [] };
    return {
      freqs: stick.freqs || [],
      intensities: stick.baseInts || [],
      labels: stick.labels || [],
    };
  }
  return {
    freqs: spectrum.freqs || [],
    intensities: spectrum.baseAbs || [],
    labels: spectrum.labels || [],
  };
}

function renderSpectrumPlot(spectrum, target) {
  const view = spectrum.previewView || (spectrum.mode === 'profile' ? 'profile' : 'stick');
  const dataSet = getSpectrumDisplayData(spectrum, view);
  if (!dataSet.freqs.length || !dataSet.intensities.length) {
    renderEmptyPlot(target, 'No data available for this spectrum.');
    return;
  }
  const trace = view === 'stick'
    ? buildStickTrace(dataSet.freqs, dataSet.intensities, '#2b6cb0', `${spectrum.name} (stick)`, dataSet.labels)
    : buildLineTrace(dataSet.freqs, dataSet.intensities, '#2b6cb0', `${spectrum.name} (profile)`, dataSet.labels);

  const layout = {
    xaxis: { title: 'Frequency (MHz)' },
    yaxis: { title: 'Intensity' },
    margin: { l: 70, r: 30, t: 30, b: 60 },
    showlegend: false,
  };

  Plotly.react(target, [trace], layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadReady(target, true);
}

function updateActiveTabDisplay() {
  const tabs = uploadedTabs?.querySelectorAll('[role=\"tab\"]') || [];
  const panels = uploadedPanels?.querySelectorAll('[role=\"tabpanel\"]') || [];
  tabs.forEach((tab) => {
    const isActive = tab.dataset.id === activeSpectrumId;
    tab.classList.toggle('is-active', isActive);
    tab.setAttribute('aria-selected', isActive ? 'true' : 'false');
    tab.tabIndex = isActive ? 0 : -1;
  });
  panels.forEach((panel) => {
    const isActive = panel.dataset.id === activeSpectrumId;
    panel.hidden = !isActive;
  });
}

function setActiveSpectrum(id) {
  activeSpectrumId = id;
  updateActiveTabDisplay();
}

function renderUploadedTabs() {
  if (!uploadedTabs || !uploadedPanels) return;
  uploadedTabs.innerHTML = '';
  uploadedPanels.innerHTML = '';

  if (!spectra.size) {
    const message = document.createElement('p');
    message.className = 'empty-message';
    message.textContent = 'No spectra loaded.';
    uploadedPanels.appendChild(message);
    activeSpectrumId = null;
    return;
  }

  const items = Array.from(spectra.values());
  if (!activeSpectrumId || !spectra.has(activeSpectrumId)) {
    activeSpectrumId = items[0]?.id || null;
  }

  items.forEach((spectrum) => {
    const tabButton = document.createElement('button');
    tabButton.type = 'button';
    tabButton.className = 'tab-button';
    tabButton.id = `tab-${spectrum.id}`;
    tabButton.dataset.id = spectrum.id;
    tabButton.setAttribute('role', 'tab');
    tabButton.setAttribute('aria-controls', `panel-${spectrum.id}`);
    tabButton.textContent = spectrum.name;
    tabButton.addEventListener('click', () => setActiveSpectrum(spectrum.id));

    const panel = document.createElement('div');
    panel.className = 'tab-panel';
    panel.id = `panel-${spectrum.id}`;
    panel.dataset.id = spectrum.id;
    panel.setAttribute('role', 'tabpanel');
    panel.setAttribute('aria-labelledby', tabButton.id);

    const card = document.createElement('div');
    card.className = 'spectrum-card';
    card.dataset.id = spectrum.id;

    const header = document.createElement('div');
    header.className = 'spectrum-card__header';

    const titleBlock = document.createElement('div');
    const title = document.createElement('h3');
    title.className = 'spectrum-card__title';
    title.textContent = spectrum.name;
    const meta = document.createElement('p');
    meta.className = 'spectrum-card__meta';
    const typeLabel = spectrum.mode === 'profile'
      ? 'Type: sampled continuous profile'
      : 'Type: stick spectrum';
    meta.textContent = spectrum.isTheoretical ? `${typeLabel} · Theoretical` : typeLabel;
    titleBlock.appendChild(title);
    titleBlock.appendChild(meta);

    const removeBtn = document.createElement('button');
    removeBtn.type = 'button';
    removeBtn.className = 'ghost';
    removeBtn.textContent = 'Remove';
    removeBtn.addEventListener('click', () => {
      if (spectrum.isTheoretical) return;
      const wasActive = spectrum.id === activeSpectrumId;
      spectra.delete(spectrum.id);
      overlaySettings.delete(spectrum.id);
      for (let i = manualAssignments.length - 1; i >= 0; i--) {
        const item = manualAssignments[i];
        if (item.expSpectrumId === spectrum.id || item.theorySpectrumId === spectrum.id) {
          manualAssignments.splice(i, 1);
        }
      }
      if (isManualAssignmentMode()) {
        syncFitLabelsFromManualAssignments();
      }
      if (pendingManualSelection && (pendingManualSelection.theorySpectrumId === spectrum.id)) {
        pendingManualSelection = null;
      }
      if (wasActive) {
        const remaining = Array.from(spectra.values());
        activeSpectrumId = remaining[0]?.id || null;
      }
      renderLabelLists();
      renderUploadedTabs();
      renderOverlayList();
      updateOverlayPlot();
    });

    header.appendChild(titleBlock);
    if (!spectrum.isTheoretical) {
      header.appendChild(removeBtn);
    } else {
      const badge = document.createElement('span');
      badge.className = 'theory-badge';
      badge.textContent = 'Theoretical';
      header.appendChild(badge);
    }

    const controls = document.createElement('div');
    controls.className = 'spectrum-card__controls';

    if (spectrum.mode === 'profile') {
      const profileLabel = document.createElement('label');
      const profileInput = document.createElement('input');
      profileInput.type = 'radio';
      profileInput.name = `preview-${spectrum.id}`;
      profileInput.value = 'profile';
      profileInput.checked = spectrum.previewView !== 'stick';
      profileLabel.appendChild(profileInput);
      profileLabel.append(' Profile');

      const stickLabel = document.createElement('label');
      const stickInput = document.createElement('input');
      stickInput.type = 'radio';
      stickInput.name = `preview-${spectrum.id}`;
      stickInput.value = 'stick';
      stickInput.checked = spectrum.previewView === 'stick';
      stickLabel.appendChild(stickInput);
      stickLabel.append(' Derived stick');

      profileInput.addEventListener('change', () => {
        spectrum.previewView = 'profile';
        renderSpectrumPlot(spectrum, plot);
      });
      stickInput.addEventListener('change', () => {
        spectrum.previewView = 'stick';
        renderSpectrumPlot(spectrum, plot);
      });

      controls.appendChild(profileLabel);
      controls.appendChild(stickLabel);
    } else {
      const label = document.createElement('span');
      label.textContent = 'View: stick';
      controls.appendChild(label);
    }

    const plot = document.createElement('div');
    plot.className = 'spectrum-card__plot';
    const plotActions = document.createElement('div');
    plotActions.className = 'plot-download-row';
    const downloadBtn = document.createElement('button');
    downloadBtn.type = 'button';
    downloadBtn.className = 'plot-download-button';
    downloadBtn.textContent = 'Download PNG (high-res)';
    downloadBtn.disabled = true;
    plotActions.appendChild(downloadBtn);
    registerPlotDownloadControl(
      plot,
      downloadBtn,
      () => `wms_fitrot_${spectrum.name}_${spectrum.previewView === 'stick' ? 'stick' : 'profile'}`
    );

    card.appendChild(header);
    card.appendChild(controls);
    card.appendChild(plot);
    card.appendChild(plotActions);
    panel.appendChild(card);

    uploadedTabs.appendChild(tabButton);
    uploadedPanels.appendChild(panel);
    renderSpectrumPlot(spectrum, plot);
  });

  updateActiveTabDisplay();
  updateAutofitSelectors();
}

function ensureOverlaySettings(spectrum) {
  if (overlaySettings.has(spectrum.id)) return overlaySettings.get(spectrum.id);
  const settings = {
    enabled: true,
    view: spectrum.mode === 'profile' ? 'profile' : 'stick',
    intensityMode: 'absorbance',
    normalize: false,
    scale: 1,
    color: spectrum.isTheoretical ? '#f59e0b' : getNextColor(),
  };
  overlaySettings.set(spectrum.id, settings);
  return settings;
}

const scheduleOverlayUpdate = (() => {
  let timer;
  return () => {
    clearTimeout(timer);
    timer = setTimeout(() => updateOverlayPlot(), 120);
  };
})();

function renderOverlayList() {
  overlayList.innerHTML = '';
  if (!spectra.size) {
    const message = document.createElement('p');
    message.className = 'empty-message';
    message.textContent = 'Upload at least one spectrum to enable the overlay.';
    overlayList.appendChild(message);
    return;
  }

  spectra.forEach((spectrum) => {
    const settings = ensureOverlaySettings(spectrum);
    const item = document.createElement('div');
    item.className = 'overlay-item';
    if (spectrum.isTheoretical) {
      item.classList.add('overlay-item--theory');
    }

    const nameBlock = document.createElement('div');
    nameBlock.className = 'overlay-item__name';
    const enabledInput = document.createElement('input');
    enabledInput.type = 'checkbox';
    enabledInput.checked = settings.enabled;
    enabledInput.addEventListener('change', () => {
      settings.enabled = enabledInput.checked;
      scheduleOverlayUpdate();
    });
    const nameSpan = document.createElement('span');
    nameSpan.textContent = spectrum.name;
    nameBlock.appendChild(enabledInput);
    nameBlock.appendChild(nameSpan);
    if (spectrum.isTheoretical) {
      const badge = document.createElement('span');
      badge.className = 'theory-badge';
      badge.textContent = 'Theoretical';
      nameBlock.appendChild(badge);
    }

    const viewLabel = document.createElement('label');
    viewLabel.textContent = 'View';
    const viewSelect = document.createElement('select');
    if (spectrum.mode === 'profile') {
      const optionProfile = document.createElement('option');
      optionProfile.value = 'profile';
      optionProfile.textContent = 'Profile';
      const optionStick = document.createElement('option');
      optionStick.value = 'stick';
      optionStick.textContent = 'Derived stick';
      viewSelect.appendChild(optionProfile);
      viewSelect.appendChild(optionStick);
      viewSelect.value = settings.view;
      viewSelect.addEventListener('change', () => {
        settings.view = viewSelect.value;
        scheduleOverlayUpdate();
      });
      viewLabel.appendChild(viewSelect);
    } else {
      viewSelect.disabled = true;
      const optionStick = document.createElement('option');
      optionStick.value = 'stick';
      optionStick.textContent = 'Stick';
      viewSelect.appendChild(optionStick);
      viewSelect.value = 'stick';
      viewLabel.appendChild(viewSelect);
      settings.view = 'stick';
    }

    const modeLabel = document.createElement('label');
    modeLabel.textContent = 'Mode';
    const modeSelect = document.createElement('select');
    const optionAbs = document.createElement('option');
    optionAbs.value = 'absorbance';
    optionAbs.textContent = 'Absorbance';
    const optionTrans = document.createElement('option');
    optionTrans.value = 'trasmittance';
    optionTrans.textContent = 'Transmittance';
    modeSelect.appendChild(optionAbs);
    modeSelect.appendChild(optionTrans);
    modeSelect.value = settings.intensityMode;
    modeSelect.addEventListener('change', () => {
      settings.intensityMode = modeSelect.value;
      scheduleOverlayUpdate();
    });
    modeLabel.appendChild(modeSelect);

    const normalizeLabel = document.createElement('label');
    const normalizeInput = document.createElement('input');
    normalizeInput.type = 'checkbox';
    normalizeInput.checked = settings.normalize;
    normalizeInput.addEventListener('change', () => {
      settings.normalize = normalizeInput.checked;
      scheduleOverlayUpdate();
    });
    normalizeLabel.appendChild(normalizeInput);
    normalizeLabel.append(' Normalize 0-1');

    const scaleLabel = document.createElement('label');
    scaleLabel.textContent = 'Scale';
    const scaleInput = document.createElement('input');
    scaleInput.type = 'number';
    scaleInput.step = '0.1';
    scaleInput.value = settings.scale;
    scaleInput.addEventListener('input', () => {
      const value = Number(scaleInput.value);
      settings.scale = Number.isFinite(value) ? value : 1;
      scheduleOverlayUpdate();
    });
    scaleLabel.appendChild(scaleInput);

    const colorInput = document.createElement('input');
    colorInput.type = 'color';
    colorInput.value = settings.color;
    if (spectrum.isTheoretical) {
      colorInput.disabled = true;
    }
    colorInput.addEventListener('input', () => {
      settings.color = colorInput.value;
      scheduleOverlayUpdate();
    });

    item.appendChild(nameBlock);
    item.appendChild(viewLabel);
    item.appendChild(modeLabel);
    item.appendChild(normalizeLabel);
    item.appendChild(scaleLabel);
    item.appendChild(colorInput);

    overlayList.appendChild(item);
  });
}

function updateOverlayPlot() {
  if (!spectra.size) {
    renderEmptyPlot(overlayPlot, 'Upload at least one spectrum to view the overlay.');
    return;
  }

  const traces = [];
  const plottedSpectrumMap = new Map();
  spectra.forEach((spectrum) => {
    const settings = overlaySettings.get(spectrum.id);
    if (!settings || !settings.enabled) return;

    const view = settings.view || (spectrum.mode === 'profile' ? 'profile' : 'stick');
    const dataSet = getSpectrumDisplayData(spectrum, view);
    if (!dataSet.freqs.length || !dataSet.intensities.length) return;

    let values = dataSet.intensities.map((value) => Math.abs(Number(value)));
    if (settings.normalize) {
      values = normalizeValues(values);
    }

    const oriented = applyIntensityMode(values, settings.intensityMode);
    const factor = Number.isFinite(Number(settings.scale)) ? Number(settings.scale) : 1;
    const scaled = oriented.map((value) => (Number.isFinite(value) ? value * factor : value));

    const labelSuffix = view === 'stick' ? 'stick' : 'profile';
    const traceName = `${spectrum.name} (${labelSuffix})`;

    const trace = view === 'stick'
      ? buildStickTrace(dataSet.freqs, scaled, settings.color, traceName, dataSet.labels)
      : buildLineTrace(dataSet.freqs, scaled, settings.color, traceName, dataSet.labels);

    plottedSpectrumMap.set(spectrum.id, {
      freqs: dataSet.freqs,
      labels: dataSet.labels || [],
      scaled,
    });
    trace.meta = { spectrumId: spectrum.id, view };
    traces.push(trace);

    const addHighlights = (highlightMap, color, lineColor) => {
      if (!highlightMap.size || !spectrum.isTheoretical) return;
      const highlightPoints = [];
      highlightMap.forEach((item) => {
        if (item.spectrumId !== spectrum.id) return;
        const targetX = Number(item.x);
        if (!Number.isFinite(targetX)) return;
        let bestIndex = -1;
        let bestDiff = Number.POSITIVE_INFINITY;
        const labels = dataSet.labels || [];
        for (let i = 0; i < dataSet.freqs.length; i++) {
          if (labels[i] !== item.label) continue;
          const fx = Number(dataSet.freqs[i]);
          if (!Number.isFinite(fx)) continue;
          const diff = Math.abs(fx - targetX);
          if (diff < bestDiff) {
            bestDiff = diff;
            bestIndex = i;
          }
        }
        if (bestIndex >= 0) {
          const yValue = scaled[bestIndex];
          if (Number.isFinite(yValue)) {
            highlightPoints.push({ x: dataSet.freqs[bestIndex], y: yValue });
          }
        }
      });
      if (highlightPoints.length) {
        traces.push({
          x: highlightPoints.map((p) => p.x),
          y: highlightPoints.map((p) => p.y),
          type: 'scatter',
          mode: 'markers',
          hoverinfo: 'skip',
          showlegend: false,
          marker: {
            size: 10,
            color,
            line: { color: lineColor, width: 2 },
          },
        });
      }
    };

    if (!isManualAssignmentMode()) {
      addHighlights(overlayHighlightsFit, '#f97316', '#b45309');
      addHighlights(overlayHighlightsControl, '#3b82f6', '#1d4ed8');
    }
  });

  if (isManualAssignmentMode() && (manualAssignments.length || pendingManualSelection?.label)) {
    const findPointByNearestFreq = (series, xTarget) => {
      if (!series || !Array.isArray(series.freqs) || !Array.isArray(series.scaled) || !series.freqs.length) {
        return null;
      }
      const target = Number(xTarget);
      if (!Number.isFinite(target)) return null;
      let bestIndex = -1;
      let bestDiff = Number.POSITIVE_INFINITY;
      for (let i = 0; i < series.freqs.length; i++) {
        const fx = Number(series.freqs[i]);
        const fy = Number(series.scaled[i]);
        if (!Number.isFinite(fx) || !Number.isFinite(fy)) continue;
        const diff = Math.abs(fx - target);
        if (diff < bestDiff) {
          bestDiff = diff;
          bestIndex = i;
        }
      }
      if (bestIndex < 0) return null;
      return {
        x: Number(series.freqs[bestIndex]),
        y: Number(series.scaled[bestIndex]),
      };
    };

    const findTheoryPoint = (series, item) => {
      if (!series || !Array.isArray(series.freqs) || !series.freqs.length) return null;
      let bestIndex = -1;
      let bestDiff = Number.POSITIVE_INFINITY;
      if (item?.label && Array.isArray(series.labels) && series.labels.length === series.freqs.length) {
        for (let i = 0; i < series.labels.length; i++) {
          if (series.labels[i] !== item.label) continue;
          const fx = Number(series.freqs[i]);
          const fy = Number(series.scaled[i]);
          if (!Number.isFinite(fx) || !Number.isFinite(fy)) continue;
          const diff = Math.abs(fx - Number(item.theoryX));
          if (diff < bestDiff) {
            bestDiff = diff;
            bestIndex = i;
          }
        }
      }
      if (bestIndex < 0) {
        return findPointByNearestFreq(series, item?.theoryX);
      }
      return {
        x: Number(series.freqs[bestIndex]),
        y: Number(series.scaled[bestIndex]),
      };
    };

    const lineX = [];
    const lineY = [];
    const lineText = [];
    const theoryMarkers = [];
    const expMarkers = [];

    manualAssignments.forEach((item) => {
      const theorySeries = plottedSpectrumMap.get(item.theorySpectrumId);
      const expSeries = plottedSpectrumMap.get(item.expSpectrumId);
      const theoryPoint = findTheoryPoint(theorySeries, item);
      const expPoint = findPointByNearestFreq(expSeries, item.expX);
      if (!theoryPoint || !expPoint) return;

      lineX.push(theoryPoint.x, expPoint.x, null);
      lineY.push(theoryPoint.y, expPoint.y, null);
      lineText.push(item.label, item.label, null);
      theoryMarkers.push({ x: theoryPoint.x, y: theoryPoint.y, label: item.label });
      expMarkers.push({ x: expPoint.x, y: expPoint.y, label: item.label });
    });

    if (lineX.length) {
      traces.push({
        x: lineX,
        y: lineY,
        text: lineText,
        type: 'scatter',
        mode: 'lines',
        showlegend: false,
        hovertemplate: 'Manual pair: %{text}<extra></extra>',
        line: { color: '#0f766e', width: 1.6, dash: 'dot' },
      });
    }
    if (theoryMarkers.length) {
      traces.push({
        x: theoryMarkers.map((p) => p.x),
        y: theoryMarkers.map((p) => p.y),
        text: theoryMarkers.map((p) => p.label),
        type: 'scatter',
        mode: 'markers',
        hovertemplate: 'Theory: %{text}<extra></extra>',
        showlegend: false,
        marker: {
          size: 11,
          color: '#f97316',
          line: { color: '#b45309', width: 2 },
        },
      });
    }
    if (expMarkers.length) {
      traces.push({
        x: expMarkers.map((p) => p.x),
        y: expMarkers.map((p) => p.y),
        text: expMarkers.map((p) => p.label),
        type: 'scatter',
        mode: 'markers',
        hovertemplate: 'Experimental: %{text}<extra></extra>',
        showlegend: false,
        marker: {
          size: 11,
          color: '#10b981',
          line: { color: '#047857', width: 2 },
        },
      });
    }

    if (pendingManualSelection?.label) {
      const pendingSeries = plottedSpectrumMap.get(pendingManualSelection.theorySpectrumId);
      const pendingPoint = findTheoryPoint(pendingSeries, pendingManualSelection);
      if (pendingPoint) {
        traces.push({
          x: [pendingPoint.x],
          y: [pendingPoint.y],
          text: [pendingManualSelection.label],
          type: 'scatter',
          mode: 'markers',
          hovertemplate: 'Pending theory: %{text}<extra></extra>',
          showlegend: false,
          marker: {
            size: 13,
            color: '#fb923c',
            line: { color: '#9a3412', width: 2 },
            symbol: 'diamond',
          },
        });
      }
    }
  }

  if (!traces.length) {
    renderEmptyPlot(overlayPlot, 'Select at least one spectrum for the overlay.');
    return;
  }

  const layout = {
    xaxis: { title: 'Frequency (MHz)' },
    yaxis: { title: 'Intensity' },
    legend: { orientation: 'h', y: -0.2 },
    margin: { l: 70, r: 30, t: 30, b: 80 },
    uirevision: 'overlay-v1',
  };

  Plotly.react(overlayPlot, traces, layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  setPlotDownloadReady(overlayPlot, true);
  bindOverlayClick();
}

async function handleFiles(files) {
  for (const file of files) {
    if (!file.name || !file.name.toLowerCase().endsWith('.dat')) {
      alert(`The file \"${file.name}\" does not have a .dat extension.`);
      continue;
    }
    const mode = await datModeDialog.open({ fileName: file.name });
    if (!mode) continue;
    try {
      const text = await file.text();
      addSpectrumFromText({ name: file.name, text, mode, setActive: true });
    } catch (error) {
      console.error('Unable to read file:', error);
      alert(`Unable to read the file \"${file.name}\".`);
    }
  }
}

uploadButton?.addEventListener('click', () => uploadInput?.click());
uploadInput?.addEventListener('change', (event) => {
  const files = Array.from(event.target.files || []);
  if (!files.length) return;
  handleFiles(files);
  uploadInput.value = '';
});

renderUploadedTabs();
renderOverlayList();
renderEmptyPlot(overlayPlot, 'Upload at least one spectrum to view the overlay.');

function updateAutofitSelectors() {
  if (!autofitTheorySelect || !autofitExpSelect) return;
  const currentTheory = autofitTheorySelect.value;
  const currentExp = autofitExpSelect.value;
  autofitTheorySelect.innerHTML = '';
  autofitExpSelect.innerHTML = '';
  const items = Array.from(spectra.values());
  if (!items.length) {
    const opt = document.createElement('option');
    opt.value = '';
    opt.textContent = 'No spectra loaded';
    autofitTheorySelect.appendChild(opt.cloneNode(true));
    autofitExpSelect.appendChild(opt);
    setQuarticScaleUiValue(1);
    updateRestoreInitialButtonState();
    return;
  }
  let theoryId = null;
  let firstExpId = null;
  items.forEach((spectrum) => {
    const opt1 = document.createElement('option');
    opt1.value = spectrum.id;
    opt1.textContent = spectrum.name;
    if (spectrum.isTheoretical) {
      opt1.textContent = `${spectrum.name} (theory)`;
      theoryId = spectrum.id;
    } else if (!firstExpId) {
      firstExpId = spectrum.id;
    }
    autofitTheorySelect.appendChild(opt1);
    const opt2 = document.createElement('option');
    opt2.value = spectrum.id;
    opt2.textContent = spectrum.name;
    autofitExpSelect.appendChild(opt2);
  });
  if (currentTheory && spectra.has(currentTheory)) {
    autofitTheorySelect.value = currentTheory;
  } else if (theoryId) {
    autofitTheorySelect.value = theoryId;
  }
  if (currentExp && spectra.has(currentExp)) {
    autofitExpSelect.value = currentExp;
  } else if (firstExpId) {
    autofitExpSelect.value = firstExpId;
  } else if (items[0]?.id) {
    autofitExpSelect.value = items[0].id;
  }
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();
  updateRestoreInitialButtonState();
}

function renderFitLabelList() {
  if (!autofitLabelList) return;
  autofitLabelList.innerHTML = '';
  if (!selectedFitLabels.length) {
    const empty = document.createElement('p');
    empty.className = 'empty-message';
    empty.textContent = 'No fit labels selected yet.';
    autofitLabelList.appendChild(empty);
    return;
  }
  selectedFitLabels.forEach((label) => {
    const chip = document.createElement('div');
    chip.className = 'label-chip';
    const span = document.createElement('span');
    span.textContent = label;
    const removeBtn = document.createElement('button');
    removeBtn.type = 'button';
    removeBtn.textContent = '×';
    removeBtn.addEventListener('click', () => {
      removeLabel(label, 'fit');
    });
    chip.appendChild(span);
    chip.appendChild(removeBtn);
    autofitLabelList.appendChild(chip);
  });
}

function renderControlLabelList() {
  if (!autofitControlList) return;
  autofitControlList.innerHTML = '';
  if (!selectedControlLabels.length) {
    const empty = document.createElement('p');
    empty.className = 'empty-message';
    empty.textContent = 'No control labels selected yet.';
    autofitControlList.appendChild(empty);
    return;
  }
  selectedControlLabels.forEach((label) => {
    const chip = document.createElement('div');
    chip.className = 'label-chip label-chip--control';
    const span = document.createElement('span');
    span.textContent = label;
    const removeBtn = document.createElement('button');
    removeBtn.type = 'button';
    removeBtn.textContent = '×';
    removeBtn.addEventListener('click', () => {
      removeLabel(label, 'control');
    });
    chip.appendChild(span);
    chip.appendChild(removeBtn);
    autofitControlList.appendChild(chip);
  });
}

function isManualAssignmentMode() {
  return (autofitAssignmentMode?.value || 'auto') === 'manual';
}

function isInputAssignmentModeEnabled() {
  return Boolean(autofitUseInputAssignment?.checked);
}

function isDirectAssignmentMode() {
  return isManualAssignmentMode() || isInputAssignmentModeEnabled();
}

function normalizeLabelValue(label) {
  let trimmed = String(label || '').trim();
  if ((trimmed.startsWith('"') && trimmed.endsWith('"')) || (trimmed.startsWith("'") && trimmed.endsWith("'"))) {
    trimmed = trimmed.slice(1, -1);
  }
  return trimmed;
}

function collectInputAssignmentPreview(expSpectrum = getSpectrumById(autofitExpSelect?.value)) {
  const preview = {
    entries: [],
    totalLines: 0,
    labeledCount: 0,
    unlabeledCount: 0,
    duplicateLabels: [],
  };
  if (!expSpectrum || expSpectrum.mode !== 'stick') return preview;
  const seen = new Set();
  const duplicateLabels = new Set();
  const freqs = Array.isArray(expSpectrum.freqs) ? expSpectrum.freqs : [];
  const labels = Array.isArray(expSpectrum.labels) ? expSpectrum.labels : [];
  const length = Math.min(freqs.length, labels.length);
  preview.totalLines = length;
  for (let i = 0; i < length; i += 1) {
    const expX = Number(freqs[i]);
    const label = normalizeLabelValue(labels[i]);
    if (!Number.isFinite(expX)) continue;
    if (!label) {
      preview.unlabeledCount += 1;
      continue;
    }
    if (seen.has(label)) {
      duplicateLabels.add(label);
      continue;
    }
    seen.add(label);
    preview.entries.push({ label, expX });
  }
  preview.labeledCount = preview.entries.length;
  preview.duplicateLabels = Array.from(duplicateLabels);
  return preview;
}

function buildInputAssignmentsFromSpectrum(theorySpectrum, expSpectrum) {
  const expMode = autofitExpMode?.value || 'profile';
  if (expMode !== 'stick') {
    throw new Error('Use input file assignation requires Experimental mode = Stick list.');
  }
  if (!theorySpectrum) {
    throw new Error('Select a theoretical spectrum before enabling Use input file assignation.');
  }
  if (!expSpectrum || expSpectrum.mode !== 'stick') {
    throw new Error('Use input file assignation requires an experimental stick spectrum.');
  }
  const preview = collectInputAssignmentPreview(expSpectrum);
  if (preview.duplicateLabels.length) {
    throw new Error(`Experimental stick file contains duplicate labels: ${preview.duplicateLabels.slice(0, 5).join(', ')}${preview.duplicateLabels.length > 5 ? ', …' : ''}.`);
  }
  if (preview.entries.length < 3) {
    throw new Error('Use input file assignation requires at least 3 labeled stick lines in the experimental file.');
  }
  const labelToFreq = buildTheoryLabelFrequencyMap(theorySpectrum);
  const missing = preview.entries
    .map((entry) => entry.label)
    .filter((label) => !labelToFreq.has(label));
  if (missing.length) {
    throw new Error(`These experimental labels were not found in the selected theoretical spectrum: ${missing.slice(0, 5).join(', ')}${missing.length > 5 ? ', …' : ''}.`);
  }
  return preview.entries.map((entry) => ({
    label: entry.label,
    theoryX: labelToFreq.get(entry.label),
    theorySpectrumId: theorySpectrum?.id || null,
    expX: entry.expX,
    expSpectrumId: expSpectrum?.id || null,
  }));
}

function syncFitLabelsFromManualAssignments() {
  selectedFitLabels.length = 0;
  manualAssignments.forEach((item) => {
    selectedFitLabels.push(item.label);
  });
}

function findManualAssignmentByExperimentalPoint(spectrumId, expX) {
  const target = Number(expX);
  if (!Number.isFinite(target)) return null;
  let best = null;
  let bestDiff = Number.POSITIVE_INFINITY;
  manualAssignments.forEach((item) => {
    if (item.expSpectrumId !== spectrumId) return;
    const current = Number(item.expX);
    if (!Number.isFinite(current)) return;
    const diff = Math.abs(current - target);
    if (diff < bestDiff) {
      bestDiff = diff;
      best = item;
    }
  });
  if (!best) return null;
  const toleranceMHz = Math.max(1e-6, Math.abs(target) * 1e-8);
  return bestDiff <= toleranceMHz ? best : null;
}

function removeManualAssignment(label, { skipRender = false, clearPending = true } = {}) {
  const idx = manualAssignments.findIndex((item) => item.label === label);
  if (idx >= 0) {
    manualAssignments.splice(idx, 1);
  }
  if (clearPending && pendingManualSelection?.label === label) {
    pendingManualSelection = null;
  }
  syncFitLabelsFromManualAssignments();
  overlayHighlightsFit.delete(label);
  if (!skipRender) {
    renderLabelLists();
    updateOverlayPlot();
  }
}

function clearManualAssignments({ skipRender = false, clearFitLabels = true } = {}) {
  manualAssignments.length = 0;
  pendingManualSelection = null;
  if (clearFitLabels) {
    selectedFitLabels.length = 0;
  }
  overlayHighlightsFit.clear();
  if (!skipRender) {
    renderLabelLists();
    updateOverlayPlot();
  }
}

function renderManualAssignmentList() {
  if (!autofitManualList) return;
  autofitManualList.innerHTML = '';
  const inputAssignmentMode = isInputAssignmentModeEnabled();
  const preview = inputAssignmentMode ? collectInputAssignmentPreview() : null;
  const displayedAssignments = inputAssignmentMode
    ? (preview?.entries || [])
    : manualAssignments;

  if (autofitManualTitle) {
    autofitManualTitle.textContent = inputAssignmentMode
      ? 'Input-file assignments'
      : 'Manual assignments (orange -> green)';
  }
  if (autofitManualHelp) {
    autofitManualHelp.innerHTML = inputAssignmentMode
      ? '<em>Represents direct theory-to-experiment pairs read from the selected experimental stick file. The third column must contain WMS-Rot labels; assignment search is skipped.</em>'
      : '<em>Represents explicit theory-to-experiment line pairs. Use this when you already trust assignments and want deterministic fitting.</em>';
  }

  if (!displayedAssignments.length) {
    const empty = document.createElement('p');
    empty.className = 'empty-message';
    empty.textContent = inputAssignmentMode
      ? 'No labeled lines were found in the selected experimental stick file.'
      : 'No manual pairs selected yet.';
    autofitManualList.appendChild(empty);
  } else {
    displayedAssignments.forEach((item) => {
      const chip = document.createElement('div');
      chip.className = 'label-chip label-chip--manual';
      const span = document.createElement('span');
      span.textContent = `${item.label} -> ${item.expX.toFixed(6)} MHz`;
      chip.appendChild(span);
      if (!inputAssignmentMode) {
        const removeBtn = document.createElement('button');
        removeBtn.type = 'button';
        removeBtn.textContent = '×';
        removeBtn.addEventListener('click', () => removeManualAssignment(item.label));
        chip.appendChild(removeBtn);
      }
      autofitManualList.appendChild(chip);
    });
  }
  if (autofitManualPending) {
    if (inputAssignmentMode) {
      if (preview?.duplicateLabels?.length) {
        autofitManualPending.textContent = `Duplicate labels detected in the experimental .dat: ${preview.duplicateLabels.slice(0, 5).join(', ')}${preview.duplicateLabels.length > 5 ? ', …' : ''}.`;
      } else if ((autofitExpMode?.value || 'profile') !== 'stick') {
        autofitManualPending.textContent = 'Set Experimental mode to Stick list to use input file assignation.';
      } else if (preview?.entries?.length) {
        const unlabeledSuffix = preview.unlabeledCount
          ? ` ${preview.unlabeledCount} unlabeled line${preview.unlabeledCount === 1 ? '' : 's'} ignored.`
          : '';
        autofitManualPending.textContent = `${preview.entries.length} labeled line${preview.entries.length === 1 ? '' : 's'} will be used directly from the experimental .dat.${unlabeledSuffix} Overlay clicks are disabled while this option is active.`;
      } else {
        autofitManualPending.textContent = 'Load a labeled stick `.dat` file to use input file assignation.';
      }
    } else if (pendingManualSelection?.label) {
      autofitManualPending.textContent = `Pending: ${pendingManualSelection.label}. Now click the experimental line (green).`;
    } else {
      autofitManualPending.textContent = 'Click one theoretical line (orange), then one experimental line (green) to create a 1:1 pair.';
    }
  }
}

function renderLabelLists() {
  renderFitLabelList();
  renderControlLabelList();
  renderManualAssignmentList();
  updateControlUiState();
}

function updateControlUiState() {
  const manualMode = isManualAssignmentMode();
  const inputAssignmentMode = isInputAssignmentModeEnabled();
  const directAssignmentMode = isDirectAssignmentMode();
  const hasControl = selectedControlLabels.length > 0;

  const setOptionVisibility = (input, visible) => {
    if (!input) return;
    const label = input.closest('label');
    if (!label) return;
    label.hidden = !visible;
    const hint = label.previousElementSibling;
    if (hint && hint.tagName === 'SMALL') {
      hint.hidden = !visible;
    }
  };

  Object.entries(autofitOptions).forEach(([key, input]) => {
    const visible = !directAssignmentMode || MANUAL_VISIBLE_OPTION_KEYS.has(key);
    setOptionVisibility(input, visible);
    if (input) {
      input.disabled = directAssignmentMode ? !MANUAL_VISIBLE_OPTION_KEYS.has(key) : false;
    }
  });

  if (autofitAssignmentMode) {
    autofitAssignmentMode.disabled = inputAssignmentMode;
  }
  if (autofitClickMode) {
    if (directAssignmentMode) {
      autofitClickMode.value = 'fit';
    }
    autofitClickMode.disabled = directAssignmentMode;
  }
  if (autofitLabelInput) {
    autofitLabelInput.disabled = directAssignmentMode;
  }
  if (autofitAddLabel) {
    autofitAddLabel.disabled = directAssignmentMode;
  }
  if (autofitControlBlock) {
    autofitControlBlock.classList.toggle('is-disabled', directAssignmentMode);
    autofitControlBlock.hidden = directAssignmentMode;
  }
  if (autofitManualBlock) {
    autofitManualBlock.hidden = !directAssignmentMode;
  }
  if (autofitOptions.mctrl) {
    autofitOptions.mctrl.disabled = directAssignmentMode ? true : hasControl;
  }
  autofitManualDisabledOptionGroups.forEach((group) => {
    group.classList.toggle('is-disabled-in-mode', directAssignmentMode);
    group.querySelectorAll('.autofit-mode-badge, .autofit-mode-note-block').forEach((element) => {
      element.hidden = !directAssignmentMode;
    });
  });
  syncLeastSquaresUiState();
  updateOverlayInteractionUi();
}

function syncLeastSquaresUiState() {
  const mode = String(autofitOptions.lsFScaleMode?.value || 'auto');
  if (!autofitOptions.lsFScale) return;
  const allowEdit = mode === 'custom' && !autofitOptions.lsFScaleMode?.disabled;
  autofitOptions.lsFScale.disabled = !allowEdit;
}

function getClickModeValue() {
  const value = autofitClickMode?.value || 'fit';
  return value === 'control' ? 'control' : 'fit';
}

function setClickModeValue(mode) {
  if (!autofitClickMode) return;
  autofitClickMode.value = mode === 'control' ? 'control' : 'fit';
}

function updateOverlayInteractionUi() {
  const manualMode = isManualAssignmentMode();
  const inputAssignmentMode = isInputAssignmentModeEnabled();
  const clickMode = getClickModeValue();
  const pendingLabel = normalizeLabelValue(pendingManualSelection?.label);

  if (overlayAutoFigure) {
    overlayAutoFigure.hidden = manualMode || inputAssignmentMode;
  }
  if (overlayManualFigure) {
    overlayManualFigure.hidden = !manualMode || inputAssignmentMode;
  }

  const syncTargetButton = (button, selected, disabled = false) => {
    if (!button) return;
    button.classList.toggle('is-selected', selected);
    button.setAttribute('aria-pressed', selected ? 'true' : 'false');
    button.disabled = disabled;
  };
  syncTargetButton(overlayTargetFit, clickMode === 'fit', manualMode || inputAssignmentMode);
  syncTargetButton(overlayTargetControl, clickMode === 'control', manualMode || inputAssignmentMode);

  [overlayStepTheory, overlayStepExperimental, overlayStepPair].forEach((element) => {
    element?.classList.remove('is-active', 'is-complete');
  });

  if (inputAssignmentMode) {
    const preview = collectInputAssignmentPreview();
    if (overlayClickSummary) {
      overlayClickSummary.innerHTML = preview.entries.length
        ? `<strong>Input-file assignment mode.</strong> ${preview.entries.length} labeled stick line${preview.entries.length === 1 ? '' : 's'} will be read directly from the experimental <code>.dat</code> file and used as explicit assignments. Overlay clicks are ignored until this option is disabled.`
        : '<strong>Input-file assignment mode.</strong> Load a labeled experimental stick <code>.dat</code> and keep Experimental mode set to Stick list. WMS-FitRot will then bypass assignment search and fit those labeled lines directly.';
    }
    if (overlayManualStatus) {
      overlayManualStatus.textContent = preview.entries.length
        ? `${preview.entries.length} labeled line${preview.entries.length === 1 ? '' : 's'} detected in the experimental stick file.`
        : 'No usable labeled stick lines detected yet.';
    }
    return;
  }

  if (manualMode) {
    if (overlayClickSummary) {
      overlayClickSummary.innerHTML = pendingLabel
        ? `<strong>Manual mode.</strong> The theoretical line <code>${pendingLabel}</code> is pending. Complete the pair by clicking the experimental counterpart on another trace in the overlay plot.`
        : '<strong>Manual mode.</strong> First click an orange theoretical line in the overlay plot, then click the experimental counterpart on another trace. Each accepted pair is stored below as an explicit theory-to-experiment assignment.';
    }
    if (pendingLabel) {
      overlayStepTheory?.classList.add('is-complete');
      overlayStepExperimental?.classList.add('is-active');
    } else {
      overlayStepTheory?.classList.add('is-active');
    }
    if (manualAssignments.length) {
      overlayStepPair?.classList.add('is-complete');
    }
    if (overlayManualStatus) {
      overlayManualStatus.textContent = pendingLabel
        ? `Pending theoretical line: ${pendingLabel}. The next accepted click must be an experimental line on a non-theoretical trace.`
        : manualAssignments.length
          ? `${manualAssignments.length} manual pair${manualAssignments.length === 1 ? '' : 's'} stored. Click a new orange theoretical line to start another pair, or click an existing paired line to remove it.`
          : 'Manual mode expects one theory click first, followed by one experimental click.';
    }
    return;
  }

  if (overlayClickSummary) {
    overlayClickSummary.innerHTML = clickMode === 'control'
      ? '<strong>Automatic mode, control selection.</strong> Click an orange theoretical peak in the overlay plot to add or remove it from the control list. Control labels are optional; if you do not provide them explicitly, WMS-FitRot can assign them automatically and use them for ranking and validation. Clicks on non-theoretical traces are ignored.'
      : '<strong>Automatic mode, fit selection.</strong> Click an orange theoretical peak in the overlay plot to add or remove it from the fit list. Fit labels define the residual vector used in optimization; clicks on non-theoretical traces are ignored.';
  }
}

function addLabel(label, mode = 'fit') {
  if (isDirectAssignmentMode()) return false;
  const trimmed = normalizeLabelValue(label);
  if (!trimmed) return false;
  const list = mode === 'control' ? selectedControlLabels : selectedFitLabels;
  const other = mode === 'control' ? selectedFitLabels : selectedControlLabels;
  if (other.includes(trimmed)) return false;
  if (list.includes(trimmed)) return false;
  list.push(trimmed);
  renderLabelLists();
  return true;
}

function removeLabel(label, mode = 'fit') {
  if (mode === 'fit' && isManualAssignmentMode()) {
    removeManualAssignment(label);
    return;
  }
  const list = mode === 'control' ? selectedControlLabels : selectedFitLabels;
  const idx = list.indexOf(label);
  if (idx >= 0) list.splice(idx, 1);
  if (mode === 'control') {
    overlayHighlightsControl.delete(label);
  } else {
    overlayHighlightsFit.delete(label);
  }
  renderLabelLists();
  updateOverlayPlot();
}

let overlayClickBound = false;

function handleManualOverlaySelection(point, spectrumId) {
  const theoryId = autofitTheorySelect?.value;
  if (!theoryId) {
    showOverlayToast('Select a theoretical spectrum first.');
    return;
  }

  if (spectrumId === theoryId) {
    const label = normalizeLabelValue(point.customdata || point.text);
    if (!label) {
      showOverlayToast('Clicked point has no theoretical label');
      return;
    }
    const existingPair = manualAssignments.find((item) => (
      item.label === label && item.theorySpectrumId === spectrumId
    ));
    if (existingPair) {
      removeManualAssignment(label);
      showOverlayToast('Manual pair removed');
      return;
    }
    if (pendingManualSelection?.label === label && pendingManualSelection?.theorySpectrumId === spectrumId) {
      pendingManualSelection = null;
      renderLabelLists();
      updateOverlayPlot();
      showOverlayToast('Pending manual selection cleared');
      return;
    }
    pendingManualSelection = {
      label,
      theoryX: Number(point.x),
      theorySpectrumId: spectrumId,
    };
    showOverlayToast('Theoretical line selected. Now click an experimental line.');
    renderLabelLists();
    updateOverlayPlot();
    return;
  }

  if (spectrumId && spectrumId !== theoryId) {
    const existingPair = findManualAssignmentByExperimentalPoint(spectrumId, point.x);
    if (existingPair) {
      removeManualAssignment(existingPair.label);
      showOverlayToast('Manual pair removed');
      return;
    }
    if (!pendingManualSelection?.label) {
      if (autofitExpSelect && autofitExpSelect.value !== spectrumId && spectra.has(spectrumId)) {
        autofitExpSelect.value = spectrumId;
      }
      showOverlayToast('Select a theoretical line first.');
      return;
    }
    if (autofitExpSelect && autofitExpSelect.value !== spectrumId && spectra.has(spectrumId)) {
      autofitExpSelect.value = spectrumId;
    }
    const expX = Number(point.x);
    if (!Number.isFinite(expX)) return;
    const pending = pendingManualSelection;
    if (!pending?.label) return;

    removeManualAssignment(pending.label, { skipRender: true, clearPending: false });
    manualAssignments.push({
      label: pending.label,
      theoryX: pending.theoryX,
      theorySpectrumId: pending.theorySpectrumId,
      expX,
      expSpectrumId: spectrumId,
    });
    const controlIdx = selectedControlLabels.indexOf(pending.label);
    if (controlIdx >= 0) {
      selectedControlLabels.splice(controlIdx, 1);
      overlayHighlightsControl.delete(pending.label);
    }
    syncFitLabelsFromManualAssignments();
    overlayHighlightsFit.set(pending.label, {
      label: pending.label,
      x: pending.theoryX,
      spectrumId: pending.theorySpectrumId,
    });
    pendingManualSelection = null;
    renderLabelLists();
    updateOverlayPlot();
    showOverlayToast('Manual pair added');
    return;
  }

  if (pendingManualSelection?.label) {
    const pending = pendingManualSelection;
    const fallbackExpId = autofitExpSelect?.value;
    const expX = Number(point?.x);
    if (pending?.label && fallbackExpId && fallbackExpId !== theoryId && spectra.has(fallbackExpId) && Number.isFinite(expX)) {
      removeManualAssignment(pending.label, { skipRender: true, clearPending: false });
      manualAssignments.push({
        label: pending.label,
        theoryX: pending.theoryX,
        theorySpectrumId: pending.theorySpectrumId,
        expX,
        expSpectrumId: fallbackExpId,
      });
      syncFitLabelsFromManualAssignments();
      pendingManualSelection = null;
      renderLabelLists();
      updateOverlayPlot();
      showOverlayToast('Manual pair added');
      return;
    }
  }
  showOverlayToast('Click an experimental trace (not the theoretical one).');
}

function bindOverlayClick() {
  if (overlayClickBound || !overlayPlot || typeof overlayPlot.on !== 'function') return;
  overlayPlot.on('plotly_click', (event) => {
    const point = event?.points?.[0];
    if (!point || !point.data) return;
    const traceFromCurve = Number.isInteger(point.curveNumber) ? overlayPlot?.data?.[point.curveNumber] : null;
    const spectrumId = point.data?.meta?.spectrumId || traceFromCurve?.meta?.spectrumId || null;
    if (!autofitTheorySelect) return;
    if (isInputAssignmentModeEnabled()) {
      showOverlayToast('Disable Use input file assignation to edit assignments manually.');
      return;
    }
    if (isManualAssignmentMode()) {
      handleManualOverlaySelection(point, spectrumId);
      return;
    }
    if (!spectrumId) return;
    if (autofitTheorySelect.value !== spectrumId) return;
    const label = normalizeLabelValue(point.customdata || point.text);
    if (!label) return;
    if (selectedFitLabels.includes(label)) {
      removeLabel(label, 'fit');
      showOverlayToast('Removed from fit list');
      return;
    }
    if (selectedControlLabels.includes(label)) {
      removeLabel(label, 'control');
      showOverlayToast('Removed from control list');
      return;
    }
    const mode = getClickModeValue();
    const added = addLabel(label, mode);
    if (added) {
      const target = mode === 'control' ? overlayHighlightsControl : overlayHighlightsFit;
      target.set(label, { label, x: point.x, spectrumId });
      updateOverlayPlot();
      showOverlayToast(mode === 'control' ? 'Added to control list' : 'Added to fit list');
    }
  });
  overlayClickBound = true;
}

function initWorker() {
  if (autofitWorker) return autofitWorker;
  autofitWorker = new Worker('autofit-worker.js');
  autofitWorker.addEventListener('message', (event) => {
    const { type, payload } = event.data || {};
    if (type === 'progress') {
      if (payload?.message) setAutofitStatus(payload.message);
      if (autofitProgressBar && Number.isFinite(payload?.progress)) {
        setAutofitProgress(payload.progress);
      }
      return;
    }
    if (type === 'result') {
      autofitRunInProgress = false;
      const outputText = payload?.output || '';
      const fitMetrics = payload?.fitMetrics && typeof payload.fitMetrics === 'object'
        ? payload.fitMetrics
        : null;
      if (autofitOutput) autofitOutput.textContent = outputText;
      lastAutofitSummary = parseAutofitSummary(outputText);
      applyAutomaticControlLines(lastAutofitSummary);
      const fittedSuffix = Number.isFinite(lastAutofitSummary?.fittedLines)
        ? ` Fitted lines: ${lastAutofitSummary.fittedLines}.`
        : '';
      setAutofitStatus(`Completed.${fittedSuffix}`);
      setAutofitWarningVisible(false);
      void window.WMSPwaEnhancements?.notifyBackground?.({
        title: 'WMS-FitRot',
        body: Number.isFinite(lastAutofitSummary?.fittedLines)
          ? `Autofit completed (${lastAutofitSummary.fittedLines} fitted lines).`
          : 'Autofit completed.',
        tag: 'wms-fitrot-autofit'
      });
      if (autofitRunButton) autofitRunButton.disabled = false;
      if (autofitCancelButton) autofitCancelButton.disabled = true;
      updateRestoreInitialButtonState();
      const finalModelParams = parseFinalModelParams(outputText);
      const finalQuarticScaleState = parseFinalQuarticScaleState(outputText);
      const finalConstants = finalModelParams || parseFinalConstants(outputText);
      const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
      const expSpectrum = getSpectrumById(autofitExpSelect?.value);
      const exportModelParams = normalizeModelParams({
        ...(getSpectrumModelParams(theorySpectrum) || {}),
        ...(finalModelParams || finalConstants || {}),
      });
      const diagnostics = buildAutofitDiagnostics({
        outputText,
        expSpectrum,
      });
      renderSensitivityDashboard(fitMetrics);
      renderAutofitDiagnostics(diagnostics);
      void populateSpfitSeedDownload({
        outputText,
        modelParams: exportModelParams,
        quarticScaleState: finalQuarticScaleState,
        theorySpectrum,
        expSpectrum,
        diagnostics,
      });
      if (finalConstants) {
        if (theorySpectrum) {
          if (finalQuarticScaleState) {
            theorySpectrum.quarticScale = finalQuarticScaleState.scale;
            theorySpectrum.quarticReferenceParams = finalQuarticScaleState.referenceParams;
          }
          theorySpectrum.fittedConstants = {
            A: finalConstants.A,
            B: finalConstants.B,
            C: finalConstants.C,
          };
          const mergedModelParams = normalizeModelParams({
            ...(getSpectrumModelParams(theorySpectrum) || {}),
            ...(finalModelParams || finalConstants),
          });
          if (mergedModelParams) {
            theorySpectrum.modelParams = mergedModelParams;
          }
          syncQuarticScaleUiFromSpectrum();
          updateAutofitParamBoundPlaceholders();
        }
        if (!recalcTheoreticalSpectrum(finalModelParams || finalConstants)) {
          setAutofitProgress(1);
        }
      } else {
        setAutofitProgress(1);
        setAutofitStatus(`Completed.${fittedSuffix} Unable to parse fitted constants for theory refresh.`);
      }
      return;
    }
    if (type === 'recalc') {
      const fittedSuffix = Number.isFinite(lastAutofitSummary?.fittedLines)
        ? ` Fitted lines: ${lastAutofitSummary.fittedLines}.`
        : '';
      try {
        const data = JSON.parse(payload?.output || '{}');
        applyRecalculatedSpectrum(data);
        setAutofitProgress(1);
        const importReport = finalizePendingSpfitImport();
        if (importReport) {
          const rmsText = importReport.residualComparison && Number.isFinite(importReport.residualComparison.afterRmsMHz)
            ? ` Residual RMS ${formatDiagnosticNumber(importReport.residualComparison.beforeRmsMHz, 4)} -> ${formatDiagnosticNumber(importReport.residualComparison.afterRmsMHz, 4)} MHz.`
            : '';
          setAutofitStatus(`SPFIT refinement imported.${rmsText}`);
        } else {
          setAutofitStatus(`Completed. Theory spectrum refreshed.${fittedSuffix}`);
        }
      } catch (error) {
        console.error('Unable to apply recalculated spectrum:', error);
        pendingSpfitImportState = null;
        setAutofitProgress(1);
        setAutofitStatus(`Completed. Theory refresh failed.${fittedSuffix}`);
      }
      return;
    }
    if (type === 'error') {
      autofitRunInProgress = false;
      lastAutofitSummary = null;
      pendingSpfitImportState = null;
      renderAutofitDiagnostics(null);
      setAutofitStatus('Autofit failed. See output.');
      setAutofitWarningVisible(false);
      void window.WMSPwaEnhancements?.notifyBackground?.({
        title: 'WMS-FitRot',
        body: 'Autofit failed.',
        tag: 'wms-fitrot-autofit'
      });
      if (autofitOutput) autofitOutput.textContent = payload?.error || 'Unknown error';
      if (autofitRunButton) autofitRunButton.disabled = false;
      if (autofitCancelButton) autofitCancelButton.disabled = true;
      updateRestoreInitialButtonState();
    }
  });
  return autofitWorker;
}

function parseFinalConstants(output) {
  const finalModel = parseFinalModelParams(output);
  if (finalModel) {
    return { A: finalModel.A, B: finalModel.B, C: finalModel.C };
  }
  if (!output) return null;
  const lines = String(output).split(/\r?\n/);
  const idx = lines.findIndex((line) => line.trim() === 'Final constants:');
  if (idx < 0 || idx + 1 >= lines.length) return null;
  const line = lines[idx + 1];
  const match = line.match(/A=([+\-]?\d+(?:\.\d+)?(?:[Ee][+\-]?\d+)?)(?:\(\d+\)|\s+\[[^\]]+\])?\s+B=([+\-]?\d+(?:\.\d+)?(?:[Ee][+\-]?\d+)?)(?:\(\d+\)|\s+\[[^\]]+\])?\s+C=([+\-]?\d+(?:\.\d+)?(?:[Ee][+\-]?\d+)?)(?:\(\d+\)|\s+\[[^\]]+\])?/);
  if (!match) return null;
  const A = Number(match[1]);
  const B = Number(match[2]);
  const C = Number(match[3]);
  if (![A, B, C].every(Number.isFinite)) return null;
  return { A, B, C };
}

function parseModelParamsFromInputText(text) {
  const values = parseAutofitRotationalSection(text);
  if (!values) return null;

  const read = (...aliases) => {
    for (const alias of aliases) {
      if (!values.has(alias)) continue;
      const value = Number(values.get(alias));
      if (Number.isFinite(value)) return value;
    }
    return undefined;
  };

  return normalizeModelParams({
    A: read('A_MHz', 'A'),
    B: read('B_MHz', 'B'),
    C: read('C_MHz', 'C'),
    DJ: read('DELTA J_MHz', 'DJ_MHz', 'DJ'),
    DJK: read('DELTA JK_MHz', 'DJK_MHz', 'DJK'),
    DK: read('DELTA K_MHz', 'DK_MHz', 'DK'),
    dJ: read('delta J_MHz', 'd1_MHz', 'dJ_MHz', 'd1', 'dJ'),
    dK: read('delta K_MHz', 'd2_MHz', 'dK_MHz', 'd2', 'dK'),
    HJ: read('PHI N_MHz', 'H N_MHz', 'HJ_MHz', 'HJ'),
    HJK: read('PHI NK_MHz', 'H NK_MHz', 'HJK_MHz', 'HJK'),
    HKJ: read('PHI KN_MHz', 'H KN_MHz', 'HKJ_MHz', 'HKJ'),
    HK: read('PHI K_MHz', 'H K_MHz', 'HK_MHz', 'HK'),
    h1: read('phi N_MHz', 'h1_MHz', 'h1'),
    h2: read('phi NK_MHz', 'h2_MHz', 'h2'),
    h3: read('phi K_MHz', 'h3_MHz', 'h3'),
  });
}

function getAutofitInputModelParams() {
  const source = typeof autofitInputText === 'string' ? autofitInputText : '';
  if (source === cachedAutofitInputModelParamsText) {
    return cachedAutofitInputModelParams;
  }
  cachedAutofitInputModelParamsText = source;
  cachedAutofitInputModelParams = parseModelParamsFromInputText(source);
  return cachedAutofitInputModelParams;
}

function normalizeQuarticReferenceParams(params) {
  if (!params || typeof params !== 'object') return null;
  const normalized = {};
  QUARTIC_PARAM_KEYS.forEach((key) => {
    const value = Number(params[key]);
    normalized[key] = Number.isFinite(value) ? value : 0;
  });
  return normalized;
}

function normalizeQuarticScale(value, fallback = 1) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    return Number.isFinite(Number(fallback)) && Number(fallback) > 0 ? Number(fallback) : 1;
  }
  return parsed;
}

function computeQuarticScaleFactor(scale) {
  const safeScale = normalizeQuarticScale(scale, 1);
  return 1 / (safeScale * safeScale);
}

function computeScaledQuarticParams(referenceParams, scale) {
  const reference = normalizeQuarticReferenceParams(referenceParams);
  if (!reference) return null;
  const factor = computeQuarticScaleFactor(scale);
  const scaled = {};
  QUARTIC_PARAM_KEYS.forEach((key) => {
    scaled[key] = reference[key] * factor;
  });
  return scaled;
}

function normalizeModelParams(params) {
  if (!params || typeof params !== 'object') return null;
  const normalized = {};
  MODEL_PARAM_KEYS.forEach((key) => {
    const value = Number(params[key]);
    if (!Number.isFinite(value)) return;
    normalized[key] = value;
  });
  if (!['A', 'B', 'C'].every((key) => Number.isFinite(normalized[key]))) {
    return null;
  }
  DISTORTION_PARAM_KEYS.forEach((key) => {
    if (!Number.isFinite(normalized[key])) normalized[key] = 0;
  });
  return normalized;
}

function parseFinalModelParams(output) {
  if (!output) return null;
  const line = String(output)
    .split(/\r?\n/)
    .find((entry) => String(entry || '').startsWith('FINAL_MODEL|'));
  if (!line) return null;
  const params = {};
  line.split('|').slice(1).forEach((segment) => {
    const idx = segment.indexOf('=');
    if (idx <= 0) return;
    const key = segment.slice(0, idx).trim();
    if (!MODEL_PARAM_KEYS.includes(key)) return;
    const value = Number(segment.slice(idx + 1).trim());
    if (!Number.isFinite(value)) return;
    params[key] = value;
  });
  return normalizeModelParams(params);
}

function parseFinalQuarticScaleState(output) {
  if (!output) return null;
  const line = String(output)
    .split(/\r?\n/)
    .find((entry) => String(entry || '').startsWith('FINAL_QUARTIC_SCALE|'));
  if (!line) return null;
  const result = { scale: 1, referenceParams: {} };
  line.split('|').slice(1).forEach((segment) => {
    const idx = segment.indexOf('=');
    if (idx <= 0) return;
    const key = segment.slice(0, idx).trim();
    const value = Number(segment.slice(idx + 1).trim());
    if (!Number.isFinite(value)) return;
    if (key === FIT_SCALE_PARAM_KEY) {
      result.scale = normalizeQuarticScale(value, 1);
      return;
    }
    if (!key.startsWith('ref_')) return;
    const param = key.slice(4);
    if (!QUARTIC_PARAM_KEYS.includes(param)) return;
    result.referenceParams[param] = value;
  });
  return {
    scale: normalizeQuarticScale(result.scale, 1),
    referenceParams: normalizeQuarticReferenceParams(result.referenceParams),
  };
}

function getSpectrumQuarticScale(spectrum) {
  return normalizeQuarticScale(spectrum?.quarticScale, 1);
}

function getSpectrumQuarticReferenceParams(spectrum, effectiveModelParams = null) {
  if (!spectrum) return null;
  const explicit = normalizeQuarticReferenceParams(spectrum.quarticReferenceParams);
  if (explicit) return explicit;
  if (spectrum.isTheoretical) {
    const fromInput = normalizeQuarticReferenceParams(getAutofitInputModelParams());
    if (fromInput) return fromInput;
  }
  const effective = normalizeModelParams(effectiveModelParams || spectrum.modelParams);
  if (!effective) return null;
  const scaleSquared = getSpectrumQuarticScale(spectrum) ** 2;
  const derived = {};
  QUARTIC_PARAM_KEYS.forEach((key) => {
    derived[key] = Number(effective[key]) * scaleSquared;
  });
  return normalizeQuarticReferenceParams(derived);
}

function getSpectrumModelParams(spectrum) {
  if (!spectrum) return null;
  const baseModelState = normalizeModelParams(spectrum.modelParams)
    || (spectrum.isTheoretical ? getAutofitInputModelParams() : null)
    || normalizeModelParams({
    A: spectrum?.fittedConstants?.A,
    B: spectrum?.fittedConstants?.B,
    C: spectrum?.fittedConstants?.C,
  });
  if (!baseModelState) return null;
  const quarticReference = getSpectrumQuarticReferenceParams(spectrum, baseModelState);
  if (!quarticReference) return baseModelState;
  return normalizeModelParams({
    ...baseModelState,
    ...computeScaledQuarticParams(quarticReference, getSpectrumQuarticScale(spectrum)),
  });
}

function getSelectedOptionalFitParams() {
  const selected = DISTORTION_PARAM_KEYS.filter((key) => autofitFitParamInputs[key]?.checked);
  if (autofitFitQuarticScale?.checked) {
    selected.unshift(FIT_SCALE_PARAM_KEY);
  }
  return selected;
}

function getModelParamDisplayLabel(key) {
  return MODEL_PARAM_LABELS[key] || key;
}

function getQuarticScaleState(spectrum) {
  return {
    scale: getSpectrumQuarticScale(spectrum),
    referenceParams: getSpectrumQuarticReferenceParams(spectrum),
  };
}

function syncQuarticScaleFactorLabel(scale) {
  if (!autofitQuarticScaleFactor) return;
  autofitQuarticScaleFactor.textContent = formatAutoBoundValue(computeQuarticScaleFactor(scale));
}

function setQuarticScaleUiValue(value) {
  const normalized = normalizeQuarticScale(value, 1);
  if (autofitQuarticScaleSlider) {
    const min = Number(autofitQuarticScaleSlider.min || 0.2);
    const max = Number(autofitQuarticScaleSlider.max || 5);
    autofitQuarticScaleSlider.value = String(Math.min(max, Math.max(min, normalized)));
  }
  if (autofitQuarticScaleValue) {
    autofitQuarticScaleValue.value = String(normalized);
  }
  syncQuarticScaleFactorLabel(normalized);
  return normalized;
}

function syncQuarticScaleUiFromSpectrum() {
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  return setQuarticScaleUiValue(getSpectrumQuarticScale(theorySpectrum));
}

function syncQuarticFitSelectionUi(changedKey = null) {
  if (autofitFitQuarticScale?.checked) {
    QUARTIC_PARAM_KEYS.forEach((key) => {
      if (autofitFitParamInputs[key]) {
        autofitFitParamInputs[key].checked = false;
      }
    });
    return;
  }
  if (changedKey && QUARTIC_PARAM_KEYS.includes(changedKey) && autofitFitQuarticScale) {
    autofitFitQuarticScale.checked = false;
  }
}

function scheduleQuarticScaleRecalc() {
  clearTimeout(quarticScaleRecalcTimer);
  quarticScaleRecalcTimer = setTimeout(() => {
    if (autofitRunInProgress) return;
    const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
    if (!theorySpectrum) return;
    const effectiveModelParams = getSpectrumModelParams(theorySpectrum);
    if (!effectiveModelParams) return;
    theorySpectrum.modelParams = effectiveModelParams;
    recalcTheoreticalSpectrum(effectiveModelParams);
  }, 180);
}

function commitQuarticScaleFromUi({ scheduleRecalc = true } = {}) {
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const uiValue = normalizeQuarticScale(
    autofitQuarticScaleValue?.value ?? autofitQuarticScaleSlider?.value ?? 1,
    getSpectrumQuarticScale(theorySpectrum)
  );
  setQuarticScaleUiValue(uiValue);
  if (!theorySpectrum) {
    updateAutofitParamBoundPlaceholders();
    return;
  }
  theorySpectrum.quarticReferenceParams = getSpectrumQuarticReferenceParams(theorySpectrum);
  theorySpectrum.quarticScale = uiValue;
  theorySpectrum.modelParams = getSpectrumModelParams(theorySpectrum);
  updateAutofitParamBoundPlaceholders();
  if (scheduleRecalc) {
    scheduleQuarticScaleRecalc();
  }
}

function renderAutofitParamBounds() {
  if (!autofitParamBoundGrid) return;
  autofitParamBoundGrid.innerHTML = '';
  autofitParamBoundInputs.clear();

  const header = document.createElement('div');
  header.className = 'autofit-bound-grid__header';
  ['Parameter', 'Min', 'Max'].forEach((label) => {
    const cell = document.createElement('span');
    cell.textContent = label;
    header.appendChild(cell);
  });
  autofitParamBoundGrid.appendChild(header);

  BOUND_PARAM_KEYS.forEach((key) => {
    const row = document.createElement('div');
    row.className = 'autofit-bound-grid__row';

    const param = document.createElement('span');
    param.className = 'autofit-bound-grid__param';
    param.textContent = getModelParamDisplayLabel(key);
    row.appendChild(param);

    const minInput = document.createElement('input');
    minInput.type = 'text';
    minInput.inputMode = 'decimal';
    minInput.placeholder = 'auto';
    minInput.setAttribute('aria-label', `${getModelParamDisplayLabel(key)} minimum bound`);
    row.appendChild(minInput);

    const maxInput = document.createElement('input');
    maxInput.type = 'text';
    maxInput.inputMode = 'decimal';
    maxInput.placeholder = 'auto';
    maxInput.setAttribute('aria-label', `${getModelParamDisplayLabel(key)} maximum bound`);
    row.appendChild(maxInput);

    autofitParamBoundInputs.set(key, { minInput, maxInput });
    autofitParamBoundGrid.appendChild(row);
  });

  updateAutofitParamBoundPlaceholders();
}

function formatAutoBoundValue(value) {
  const num = Number(value);
  if (!Number.isFinite(num)) return 'n/a';
  const abs = Math.abs(num);
  if ((abs > 0 && abs < 1e-3) || abs >= 1e6) {
    return num.toExponential(4);
  }
  return Number.parseFloat(num.toPrecision(9)).toString();
}

function computeAutomaticBoundsForParam(name, value, span) {
  const val = Number(value);
  const spanAbs = Math.abs(Number(span));
  if (!Number.isFinite(val) || !Number.isFinite(spanAbs)) return null;

  let lower;
  let upper;
  if (name === 'A' || name === 'B' || name === 'C' || name === FIT_SCALE_PARAM_KEY) {
    lower = Math.max(1e-6, val * (1 - spanAbs));
    upper = Math.max(lower * 1.001, val * (1 + spanAbs));
  } else {
    const baseFloor = DISTORTION_BOUNDS_FLOORS_MHZ[name] || 0.001;
    const spanScale = spanAbs > 0 ? Math.max(0.25, spanAbs / 0.3) : 0.25;
    let window = Math.max(Math.abs(val) * spanAbs, baseFloor * spanScale);
    if (!Number.isFinite(window) || window <= 0) {
      window = baseFloor * spanScale;
    }
    lower = val - window;
    upper = val + window;
    if (lower === upper) {
      upper = lower + Math.max(1e-12, baseFloor * 1e-6);
    }
  }

  if (!Number.isFinite(lower) || !Number.isFinite(upper)) return null;
  return { lower, upper };
}

function buildAutoBoundPlaceholder(primary, secondary) {
  const hasPrimary = Number.isFinite(Number(primary));
  const hasSecondary = Number.isFinite(Number(secondary));
  if (!hasPrimary && !hasSecondary) return 'auto';
  if (!hasSecondary) return `auto (${formatAutoBoundValue(primary)})`;
  if (!hasPrimary) return `auto (${formatAutoBoundValue(secondary)})`;
  const first = formatAutoBoundValue(primary);
  const second = formatAutoBoundValue(secondary);
  if (first === second) return `auto (${first})`;
  return `auto (${first} / ${second})`;
}

function buildAutoBoundTitle(kind, placeholder) {
  if (!placeholder || placeholder === 'auto') {
    return `Automatic ${kind} bound derived from the current theory model and Span / Final span.`;
  }
  return `Automatic ${kind} bound from Span / Final span: ${placeholder}`;
}

function updateAutofitParamBoundPlaceholders() {
  if (!autofitParamBoundInputs.size) return;
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const modelParams = getSpectrumModelParams(theorySpectrum);
  const quarticScale = getSpectrumQuarticScale(theorySpectrum);
  const localSpan = Number(autofitOptions.span?.value);
  const finalSpan = Number(autofitOptions.finalSpan?.value);

  BOUND_PARAM_KEYS.forEach((key) => {
    const inputs = autofitParamBoundInputs.get(key);
    if (!inputs) return;
    const value = key === FIT_SCALE_PARAM_KEY ? quarticScale : modelParams?.[key];
    const localBounds = computeAutomaticBoundsForParam(key, value, localSpan);
    const finalBounds = computeAutomaticBoundsForParam(key, value, finalSpan);
    const minPlaceholder = buildAutoBoundPlaceholder(localBounds?.lower, finalBounds?.lower);
    const maxPlaceholder = buildAutoBoundPlaceholder(localBounds?.upper, finalBounds?.upper);
    inputs.minInput.placeholder = minPlaceholder;
    inputs.maxInput.placeholder = maxPlaceholder;
    inputs.minInput.title = buildAutoBoundTitle('min', minPlaceholder);
    inputs.maxInput.title = buildAutoBoundTitle('max', maxPlaceholder);
  });
}

function parseOptionalNumericInput(input) {
  if (!input) return null;
  const raw = String(input.value || '').trim();
  if (!raw) return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
}

function getCustomParamBounds() {
  const bounds = {};
  BOUND_PARAM_KEYS.forEach((key) => {
    const inputs = autofitParamBoundInputs.get(key);
    if (!inputs) return;
    const min = parseOptionalNumericInput(inputs.minInput);
    const max = parseOptionalNumericInput(inputs.maxInput);
    if (min == null && max == null) return;
    if (min != null && max != null && max <= min) {
      throw new Error(`Invalid custom bounds for ${getModelParamDisplayLabel(key)}: max must be greater than min.`);
    }
    bounds[key] = { min, max };
  });
  return Object.keys(bounds).length ? bounds : null;
}

function buildModelParamCliArgs(modelParams) {
  const normalized = normalizeModelParams(modelParams);
  if (!normalized) return [];
  return MODEL_PARAM_KEYS.flatMap((key) => [`--${key} ${normalized[key]}`]);
}

function buildQuarticScaleCliArgs(quarticScaleState) {
  if (!quarticScaleState || typeof quarticScaleState !== 'object') return [];
  const args = [];
  const scale = normalizeQuarticScale(quarticScaleState.scale, 1);
  args.push(`--quartic-scale ${scale}`);
  const reference = normalizeQuarticReferenceParams(quarticScaleState.referenceParams);
  if (!reference) return args;
  QUARTIC_PARAM_KEYS.forEach((key) => {
    args.push(`--quartic-ref-${key} ${reference[key]}`);
  });
  return args;
}

function buildParamBoundCliArgs(customParamBounds) {
  if (!customParamBounds || typeof customParamBounds !== 'object') return [];
  const args = [];
  BOUND_PARAM_KEYS.forEach((key) => {
    const entry = customParamBounds[key];
    if (!entry || typeof entry !== 'object') return;
    const minValue = Number(entry.min);
    const maxValue = Number(entry.max);
    const hasMin = entry.min != null && Number.isFinite(minValue);
    const hasMax = entry.max != null && Number.isFinite(maxValue);
    if (!hasMin && !hasMax) return;
    args.push(`--param-bound ${key} ${hasMin ? minValue : 'nan'} ${hasMax ? maxValue : 'nan'}`);
  });
  return args;
}

function parseAutofitSummary(output) {
  const summary = {
    controlMode: null,
    controlLines: [],
    fittedLines: null,
  };
  if (!output) return summary;
  const seenControl = new Set();
  String(output).split(/\r?\n/).forEach((line) => {
    const text = String(line || '').trim();
    if (!text) return;
    if (text.startsWith('CTRL_MODE|')) {
      const parts = text.split('|');
      if (parts[1]) summary.controlMode = String(parts[1]).trim().toUpperCase();
      return;
    }
    if (text.startsWith('CTRL_LINE|')) {
      const parts = text.split('|');
      if (parts.length < 5) return;
      const mode = String(parts[1] || '').trim().toUpperCase();
      const label = normalizeLabelValue(parts[2]);
      if (!label || seenControl.has(label)) return;
      seenControl.add(label);
      const fObs = Number(parts[3]);
      const residual = Number(parts[4]);
      summary.controlLines.push({
        mode,
        label,
        fObs: Number.isFinite(fObs) ? fObs : null,
        residual: Number.isFinite(residual) ? residual : null,
      });
      if (!summary.controlMode && mode) summary.controlMode = mode;
      return;
    }
    if (text.startsWith('FITTED_LINES|')) {
      const match = text.match(/n\s*=\s*(\d+)/i);
      if (match) {
        const value = Number.parseInt(match[1], 10);
        summary.fittedLines = Number.isFinite(value) ? value : null;
      }
    }
  });
  return summary;
}

function parseAutofitFitParamNames(output) {
  if (!output) return [];
  const line = String(output)
    .split(/\r?\n/)
    .find((entry) => String(entry || '').startsWith('FIT_PARAMS|'));
  if (!line) return [];
  return line
    .slice('FIT_PARAMS|'.length)
    .split(',')
    .map((entry) => entry.trim())
    .filter((entry) => MODEL_PARAM_KEYS.includes(entry) || entry === FIT_SCALE_PARAM_KEY);
}

function parseTaggedMetricsLine(output, prefix) {
  if (!output) return null;
  const line = String(output)
    .split(/\r?\n/)
    .find((entry) => String(entry || '').startsWith(prefix));
  if (!line) return null;
  const result = {};
  line.slice(prefix.length).split('|').forEach((segment, index) => {
    const trimmed = segment.trim();
    if (!trimmed) return;
    const eqIdx = trimmed.indexOf('=');
    if (eqIdx < 0) {
      if (index === 0) result.mode = trimmed;
      return;
    }
    const key = trimmed.slice(0, eqIdx).trim();
    const rawValue = trimmed.slice(eqIdx + 1).trim();
    const numericValue = Number(rawValue);
    result[key] = Number.isFinite(numericValue) ? numericValue : rawValue;
  });
  return result;
}

function parseAutofitRmsMetrics(output) {
  const result = {
    initialRms: null,
    finalRms: null,
    finalScore: null,
    finalScoreText: null,
  };
  if (!output) return result;
  const matches = [];
  String(output).split(/\r?\n/).forEach((line) => {
    const match = String(line || '').match(/^RMS\(labels\)\s*=\s*([+\-]?\d+(?:\.\d+)?(?:[Ee][+\-]?\d+)?)\s*MHz(?:\s+Score\s*=\s*(.+))?$/);
    if (!match) return;
    const rms = Number(match[1]);
    const scoreText = match[2] ? match[2].trim() : null;
    let score = null;
    if (scoreText && !/^n\/a/i.test(scoreText)) {
      const parsed = Number(scoreText);
      if (Number.isFinite(parsed)) score = parsed;
    }
    matches.push({
      rms: Number.isFinite(rms) ? rms : null,
      score,
      scoreText,
    });
  });
  if (matches[0]) {
    result.initialRms = matches[0].rms;
  }
  if (matches.length) {
    const last = matches[matches.length - 1];
    result.finalRms = last.rms;
    result.finalScore = last.score;
    result.finalScoreText = last.scoreText;
  }
  return result;
}

function parseAutofitOptimizerDiagnostics(output) {
  return {
    powell: parseTaggedMetricsLine(output, 'OPT_POWELL|'),
    leastSquares: parseTaggedMetricsLine(output, 'OPT_LS|'),
  };
}

function parseAutofitControlMetrics(output) {
  return parseTaggedMetricsLine(output, 'CTRL_METRICS|');
}

function parseAutofitFinalAssignments(output) {
  if (!output) return [];
  const assignments = [];
  const seen = new Set();
  let inBlock = false;
  const stopPattern = /^(Fitted lines:|FITTED_LINES\||Control scoring source:|Control metrics:|CTRL_MODE\||CTRL_LINE\||CTRL_METRICS\||Final constants:|FINAL_MODEL\||FINAL_SIGMA\||FINAL_FORMATTED\||FINAL_UNCERTAINTY_STATUS\||FINAL_QUARTIC_SCALE\||FIT_PARAMS\||OPT_POWELL\||OPT_LS\|)/;
  String(output).split(/\r?\n/).forEach((rawLine) => {
    const line = String(rawLine || '');
    const trimmed = line.trim();
    if (!inBlock) {
      if (trimmed === 'Final assignments (label -> f_obs, residual):') {
        inBlock = true;
      }
      return;
    }
    if (!trimmed || stopPattern.test(trimmed)) return;
    const parts = line.split('|');
    if (parts.length < 3) return;
    const label = normalizeLabelValue(parts[0]);
    const fObs = Number(parts[1]);
    const residual = Number(parts[2]);
    if (!label || seen.has(label) || !Number.isFinite(fObs)) return;
    seen.add(label);
    assignments.push({
      label,
      fObs,
      residual: Number.isFinite(residual) ? residual : null,
    });
  });
  return assignments;
}

function relativeToleranceMHz(freqMHz, deltaPercent) {
  const freq = Math.abs(Number(freqMHz));
  const pct = Math.abs(Number(deltaPercent));
  if (!Number.isFinite(freq) || !Number.isFinite(pct)) return 1e-9;
  return Math.max(1e-9, freq * (pct * 0.01));
}

function parsePeakListText(text) {
  const peaks = parseStickDatRows(text).map((row) => ({
    freq: Number(row.freq),
    intensity: Number(row.intensity),
    snr: Number(row.snr),
    q: Number(row.q),
  }));
  if (!peaks.length) return peaks;
  peaks.sort((a, b) => a.freq - b.freq);
  const allQMissing = peaks.every((peak) => !Number.isFinite(peak.q));
  if (allQMissing) {
    const intensities = peaks.map((peak) => Number(peak.intensity)).filter(Number.isFinite);
    const median = computeMedian(intensities);
    const mad = computeMad(intensities, median);
    const std = intensities.length ? Math.sqrt(intensities.reduce((sum, value) => sum + ((value - median) ** 2), 0) / intensities.length) : 1;
    const scale = Number.isFinite(mad) && mad > 0 ? mad : (Number.isFinite(std) && std > 0 ? std : 1);
    peaks.forEach((peak) => {
      peak.q = Math.min(5, Math.max(0.1, ((peak.intensity - median) / (5 * scale)) + 1));
    });
  }
  const allSnrMissing = peaks.every((peak) => !Number.isFinite(peak.snr));
  if (allSnrMissing) {
    peaks.forEach((peak) => {
      peak.snr = Number.isFinite(peak.q) ? peak.q : 1;
    });
  }
  return peaks.map((peak) => ({
    freq: Number(peak.freq),
    intensity: Number(peak.intensity),
    snr: Number.isFinite(peak.snr) ? Number(peak.snr) : 1,
    q: Number.isFinite(peak.q) ? Number(peak.q) : 1,
  }));
}

function getAutofitExperimentalPeaks(expSpectrum, expMode) {
  if (!expSpectrum) return [];
  const mode = expMode || autofitExpMode?.value || 'profile';
  if (mode === 'stick') {
    const parsed = parsePeakListText(expSpectrum.rawText || '');
    if (parsed.length) return parsed;
    const freqs = expSpectrum.stick?.freqs || expSpectrum.freqs || [];
    const intensities = expSpectrum.stick?.baseInts || expSpectrum.baseAbs || [];
    return freqs.map((freq, index) => ({
      freq: Number(freq),
      intensity: Number(intensities[index]),
      snr: 1,
      q: 1,
    })).filter((peak) => Number.isFinite(peak.freq) && Number.isFinite(peak.intensity));
  }
  if (expSpectrum.mode === 'profile') {
    const stick = expSpectrum.stick || computeDatStickSpectrum(expSpectrum.freqs || [], expSpectrum.baseAbs || [], null, 'absorbance');
    const freqs = stick.freqs || [];
    const intensities = stick.baseInts || stick.ints || [];
    return freqs.map((freq, index) => ({
      freq: Number(freq),
      intensity: Number(intensities[index]),
      snr: Number(intensities[index]),
      q: 1,
    })).filter((peak) => Number.isFinite(peak.freq) && Number.isFinite(peak.intensity));
  }
  return parsePeakListText(expSpectrum.rawText || '');
}

function findNearestPeak(peaks, targetFreq) {
  const target = Number(targetFreq);
  if (!Number.isFinite(target) || !Array.isArray(peaks) || !peaks.length) return null;
  let best = null;
  let bestDiff = Number.POSITIVE_INFINITY;
  peaks.forEach((peak) => {
    const freq = Number(peak?.freq);
    if (!Number.isFinite(freq)) return;
    const diff = Math.abs(freq - target);
    if (diff < bestDiff) {
      bestDiff = diff;
      best = peak;
    }
  });
  return best ? { peak: best, diff: bestDiff } : null;
}

function formatDiagnosticNumber(value, digits = 4) {
  const num = Number(value);
  if (!Number.isFinite(num)) return 'n/a';
  const abs = Math.abs(num);
  if (abs === 0) return '0';
  if (abs >= 1000 || abs < 1e-3) {
    return num.toExponential(Math.max(1, digits - 1)).replace('e', 'E');
  }
  return Number.parseFloat(num.toPrecision(Math.max(2, digits))).toString();
}

function getSensitivityStatus(metrics) {
  if (metrics?.weakly_sensitive) {
    return { key: 'weak', label: 'weak', className: 'sensitivity-status sensitivity-status--weak' };
  }
  if (metrics?.poorly_identifiable) {
    return { key: 'correlated', label: 'correlated', className: 'sensitivity-status sensitivity-status--correlated' };
  }
  return { key: 'good', label: 'good', className: 'sensitivity-status sensitivity-status--good' };
}

function getFitMetricsParameterNames(fitMetrics) {
  if (!fitMetrics || typeof fitMetrics !== 'object') return [];
  const perParameter = fitMetrics.per_parameter;
  if (!perParameter || typeof perParameter !== 'object') return [];
  return Object.keys(perParameter);
}

function getFitMetricForParameter(fitMetrics, name) {
  if (!fitMetrics || typeof fitMetrics !== 'object') return null;
  if (!name) return null;
  return fitMetrics?.per_parameter?.[name] || null;
}

function getSensitivityVisibleParameterNames(fitMetrics) {
  const names = getFitMetricsParameterNames(fitMetrics);
  if (!autofitSensitivityState.hideWeak) return names;
  return names.filter((name) => !getFitMetricForParameter(fitMetrics, name)?.weakly_sensitive);
}

function getSensitivityTableRows(fitMetrics) {
  return getFitMetricsParameterNames(fitMetrics).map((name) => {
    const metrics = getFitMetricForParameter(fitMetrics, name) || {};
    const status = getSensitivityStatus(metrics);
    return {
      name,
      value: Number(metrics.value),
      sigma: Number(metrics.sigma),
      weighted_l2_sensitivity: Number(metrics.weighted_l2_sensitivity),
      relative_weighted_l2_sensitivity: Number(metrics.relative_weighted_l2_sensitivity),
      max_abs_correlation: Number(metrics.max_abs_correlation),
      most_correlated_with: metrics.most_correlated_with || null,
      weakly_sensitive: Boolean(metrics.weakly_sensitive),
      poorly_identifiable: Boolean(metrics.poorly_identifiable),
      status: status.label,
      statusKey: status.key,
      statusClassName: status.className,
    };
  });
}

function getSensitivitySortValue(row, key) {
  if (key === 'status') {
    return { good: 0, correlated: 1, weak: 2 }[row.statusKey] ?? 99;
  }
  return row[key];
}

function compareSensitivityTableRows(left, right, key, direction) {
  const leftValue = getSensitivitySortValue(left, key);
  const rightValue = getSensitivitySortValue(right, key);
  const dir = direction === 'asc' ? 1 : -1;
  if (typeof leftValue === 'string' || typeof rightValue === 'string') {
    return dir * String(leftValue || '').localeCompare(String(rightValue || ''), undefined, { numeric: true });
  }
  const leftNum = Number(leftValue);
  const rightNum = Number(rightValue);
  if (!Number.isFinite(leftNum) && !Number.isFinite(rightNum)) return left.name.localeCompare(right.name);
  if (!Number.isFinite(leftNum)) return 1;
  if (!Number.isFinite(rightNum)) return -1;
  if (leftNum === rightNum) return left.name.localeCompare(right.name);
  return dir * (leftNum - rightNum);
}

function toggleSensitivityTableSort(key) {
  if (autofitSensitivityState.tableSortKey === key) {
    autofitSensitivityState.tableSortDirection = autofitSensitivityState.tableSortDirection === 'asc' ? 'desc' : 'asc';
  } else {
    autofitSensitivityState.tableSortKey = key;
    autofitSensitivityState.tableSortDirection = key === 'name' || key === 'status' ? 'asc' : 'desc';
  }
  renderParameterTable();
  syncSensitivityHighlights();
}

function ensureSensitivitySelection(fitMetrics) {
  const allNames = getFitMetricsParameterNames(fitMetrics);
  if (!allNames.length) {
    autofitSensitivityState.selectedParam = null;
    autofitSensitivityState.selectedParams = [];
    autofitSensitivityState.hoveredParam = null;
    return;
  }
  if (!allNames.includes(autofitSensitivityState.selectedParam)) {
    const visibleNames = getSensitivityVisibleParameterNames(fitMetrics);
    autofitSensitivityState.selectedParam = visibleNames[0] || allNames[0];
  }
  autofitSensitivityState.selectedParams = Array.isArray(autofitSensitivityState.selectedParams)
    ? autofitSensitivityState.selectedParams.filter((name, index, values) => allNames.includes(name) && values.indexOf(name) === index)
    : [];
  if (!autofitSensitivityState.selectedParams.length && autofitSensitivityState.selectedParam) {
    autofitSensitivityState.selectedParams = [autofitSensitivityState.selectedParam];
  }
  if (autofitSensitivityState.hoveredParam && !allNames.includes(autofitSensitivityState.hoveredParam)) {
    autofitSensitivityState.hoveredParam = null;
  }
}

function setSensitivityHoveredParam(name) {
  const next = name || null;
  if (autofitSensitivityState.hoveredParam === next) return;
  autofitSensitivityState.hoveredParam = next;
  syncSensitivityHighlights();
}

function selectSensitivityParameter(name, selectedParams = null) {
  const uniqueParams = Array.isArray(selectedParams)
    ? selectedParams.filter((value, index, values) => value && values.indexOf(value) === index)
    : (name ? [name] : []);
  autofitSensitivityState.selectedParam = name || uniqueParams[0] || null;
  autofitSensitivityState.selectedParams = uniqueParams.length
    ? uniqueParams
    : (autofitSensitivityState.selectedParam ? [autofitSensitivityState.selectedParam] : []);
  renderParameterDetail();
  syncSensitivityHighlights();
}

function renderSensitivityEmptyState(container, message) {
  if (!container) return;
  container.innerHTML = '';
  const empty = document.createElement('p');
  empty.className = 'sensitivity-empty';
  empty.textContent = message;
  container.appendChild(empty);
}

function getSensitivityHighlightedParams() {
  const highlighted = new Set();
  if (autofitSensitivityState.hoveredParam) {
    highlighted.add(autofitSensitivityState.hoveredParam);
  }
  (autofitSensitivityState.selectedParams || []).forEach((name) => {
    if (name) highlighted.add(name);
  });
  return highlighted;
}

function getSensitivitySelectedParams() {
  return new Set((autofitSensitivityState.selectedParams || []).filter(Boolean));
}

function getSensitivityBarColors(names) {
  const highlighted = getSensitivityHighlightedParams();
  const selected = getSensitivitySelectedParams();
  return names.map((name) => {
    if (autofitSensitivityState.hoveredParam === name) return '#f59e0b';
    if (selected.has(name)) return '#0f766e';
    if (highlighted.has(name)) return '#38bdf8';
    return '#2563eb';
  });
}

function getSensitivityHeatmapTickText(names) {
  const highlighted = getSensitivityHighlightedParams();
  const selected = getSensitivitySelectedParams();
  return names.map((name) => {
    if (autofitSensitivityState.hoveredParam === name) return `<b>${name}</b>`;
    if (selected.has(name)) return `<span style="color:#0f766e"><b>${name}</b></span>`;
    if (highlighted.has(name)) return `<span style="color:#f59e0b">${name}</span>`;
    return name;
  });
}

function buildSensitivityHeatmapShapes(names) {
  const shapes = [];
  if (!Array.isArray(names) || !names.length) return shapes;
  const highlighted = getSensitivityHighlightedParams();
  const selected = getSensitivitySelectedParams();
  const count = names.length;
  names.forEach((name, index) => {
    if (!highlighted.has(name)) return;
    const isHovered = autofitSensitivityState.hoveredParam === name;
    const isSelected = selected.has(name);
    const fill = isHovered ? 'rgba(245, 158, 11, 0.10)' : 'rgba(15, 118, 110, 0.08)';
    const line = isHovered ? 'rgba(245, 158, 11, 0.55)' : 'rgba(15, 118, 110, 0.55)';
    shapes.push({
      type: 'rect',
      xref: 'x',
      yref: 'y',
      x0: -0.5,
      x1: count - 0.5,
      y0: index - 0.5,
      y1: index + 0.5,
      fillcolor: fill,
      line: { width: isSelected ? 1.5 : 1, color: line },
      layer: 'above',
    });
    shapes.push({
      type: 'rect',
      xref: 'x',
      yref: 'y',
      x0: index - 0.5,
      x1: index + 0.5,
      y0: -0.5,
      y1: count - 0.5,
      fillcolor: fill,
      line: { width: isSelected ? 1.5 : 1, color: line },
      layer: 'above',
    });
  });
  const locked = autofitSensitivityState.selectedParams || [];
  if (locked.length === 2) {
    const firstIndex = names.indexOf(locked[0]);
    const secondIndex = names.indexOf(locked[1]);
    if (firstIndex >= 0 && secondIndex >= 0) {
      [[firstIndex, secondIndex], [secondIndex, firstIndex]].forEach(([xIndex, yIndex]) => {
        shapes.push({
          type: 'rect',
          xref: 'x',
          yref: 'y',
          x0: xIndex - 0.5,
          x1: xIndex + 0.5,
          y0: yIndex - 0.5,
          y1: yIndex + 0.5,
          fillcolor: 'rgba(0, 0, 0, 0)',
          line: { width: 2.2, color: 'rgba(15, 23, 42, 0.55)' },
          layer: 'above',
        });
      });
    }
  }
  return shapes;
}

function syncSensitivityHighlights() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!fitMetrics) return;

  if (sensitivityParameterTable) {
    const highlighted = getSensitivityHighlightedParams();
    const selected = getSensitivitySelectedParams();
    sensitivityParameterTable.querySelectorAll('[data-param-name]').forEach((row) => {
      const name = row.getAttribute('data-param-name');
      row.classList.toggle('is-highlighted', highlighted.has(name));
      row.classList.toggle('is-selected', selected.has(name));
    });
  }

  const visibleNames = getSensitivityVisibleParameterNames(fitMetrics);
  if (sensitivityRankingChart && sensitivityRankingChart.data && sensitivityRankingChart.data.length && visibleNames.length) {
    const rankingNames = (fitMetrics.ranking_by_weighted_sensitivity || [])
      .map((entry) => entry?.name)
      .filter((name) => visibleNames.includes(name));
    Plotly.restyle(sensitivityRankingChart, {
      'marker.color': [getSensitivityBarColors(rankingNames)],
    });
  }

  if (sensitivityCorrelationHeatmap && sensitivityCorrelationHeatmap.data && sensitivityCorrelationHeatmap.data.length && visibleNames.length) {
    Plotly.relayout(sensitivityCorrelationHeatmap, {
      shapes: buildSensitivityHeatmapShapes(visibleNames),
      'xaxis.ticktext': getSensitivityHeatmapTickText(visibleNames),
      'yaxis.ticktext': getSensitivityHeatmapTickText(visibleNames),
    });
  }
}

function resizeSensitivityPlots() {
  if (typeof Plotly === 'undefined' || !Plotly?.Plots?.resize) return;
  [sensitivityRankingChart, sensitivityCorrelationHeatmap].forEach((plotEl) => {
    if (!plotEl) return;
    try {
      Plotly.Plots.resize(plotEl);
    } catch {
      // Ignore resize failures while the panel is transitioning.
    }
  });
}

function scheduleSensitivityPlotResize() {
  requestAnimationFrame(() => {
    requestAnimationFrame(() => {
      resizeSensitivityPlots();
    });
  });
}

function renderParameterTable() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivityParameterTable) return;
  if (!fitMetrics) {
    renderSensitivityEmptyState(sensitivityParameterTable, 'Run a fit to inspect parameter sensitivities.');
    return;
  }

  const visibleNames = new Set(getSensitivityVisibleParameterNames(fitMetrics));
  const rows = getSensitivityTableRows(fitMetrics)
    .filter((row) => visibleNames.has(row.name))
    .sort((left, right) => compareSensitivityTableRows(
      left,
      right,
      autofitSensitivityState.tableSortKey,
      autofitSensitivityState.tableSortDirection,
    ));

  if (!rows.length) {
    renderSensitivityEmptyState(sensitivityParameterTable, 'No parameters remain after hiding weak entries.');
    return;
  }

  sensitivityParameterTable.innerHTML = '';
  const wrap = document.createElement('div');
  wrap.className = 'sensitivity-table-wrap';

  const table = document.createElement('table');
  table.className = 'sensitivity-table';
  const thead = document.createElement('thead');
  const headRow = document.createElement('tr');
  SENSITIVITY_TABLE_COLUMNS.forEach((column) => {
    const th = document.createElement('th');
    const isActive = autofitSensitivityState.tableSortKey === column.key;
    th.setAttribute('aria-sort', isActive
      ? (autofitSensitivityState.tableSortDirection === 'asc' ? 'ascending' : 'descending')
      : 'none');
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'sensitivity-table__sort';
    button.title = column.tooltip;
    button.textContent = column.label;
    if (isActive) {
      button.dataset.direction = autofitSensitivityState.tableSortDirection;
    }
    button.addEventListener('click', () => toggleSensitivityTableSort(column.key));
    th.appendChild(button);
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);
  table.appendChild(thead);

  const tbody = document.createElement('tbody');
  rows.forEach((row) => {
    const tr = document.createElement('tr');
    tr.setAttribute('data-param-name', row.name);
    tr.addEventListener('mouseenter', () => setSensitivityHoveredParam(row.name));
    tr.addEventListener('mouseleave', () => setSensitivityHoveredParam(null));
    tr.addEventListener('click', () => selectSensitivityParameter(row.name));

    const nameCell = document.createElement('td');
    nameCell.textContent = row.name;
    nameCell.title = `Parameter ${row.name}`;
    tr.appendChild(nameCell);

    const valueCell = document.createElement('td');
    valueCell.textContent = formatDiagnosticNumber(row.value, 6);
    valueCell.title = `Final fitted value for ${row.name}`;
    tr.appendChild(valueCell);

    const sensitivityCell = document.createElement('td');
    sensitivityCell.textContent = formatDiagnosticNumber(row.weighted_l2_sensitivity, 4);
    sensitivityCell.title = `Weighted local sensitivity ||J_col|| for ${row.name}`;
    tr.appendChild(sensitivityCell);

    const relativeCell = document.createElement('td');
    relativeCell.textContent = formatDiagnosticNumber(row.relative_weighted_l2_sensitivity, 4);
    relativeCell.title = `Relative sensitivity |p| * ||J_col|| for ${row.name}`;
    tr.appendChild(relativeCell);

    const sigmaCell = document.createElement('td');
    sigmaCell.textContent = formatDiagnosticNumber(row.sigma, 4);
    sigmaCell.title = `Formal 1σ uncertainty for ${row.name}`;
    tr.appendChild(sigmaCell);

    const corrCell = document.createElement('td');
    corrCell.textContent = formatDiagnosticNumber(row.max_abs_correlation, 4);
    corrCell.title = row.most_correlated_with
      ? `${row.name} is most correlated with ${row.most_correlated_with}`
      : `No dominant correlation reported for ${row.name}`;
    tr.appendChild(corrCell);

    const statusCell = document.createElement('td');
    const badge = document.createElement('span');
    badge.className = row.statusClassName;
    badge.textContent = row.status;
    badge.title = row.status === 'weak'
      ? `Low Jacobian sensitivity for ${row.name}`
      : row.status === 'correlated'
        ? `${row.name} is strongly correlated with another fitted parameter`
        : `${row.name} looks locally well constrained`;
    statusCell.appendChild(badge);
    tr.appendChild(statusCell);

    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
  wrap.appendChild(table);
  sensitivityParameterTable.appendChild(wrap);
}

function renderSensitivityChart() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivityRankingChart) return;
  if (!fitMetrics || typeof Plotly === 'undefined') {
    renderSensitivityEmptyState(sensitivityRankingChart, 'Run a fit to populate the sensitivity ranking.');
    return;
  }

  const visibleNames = new Set(getSensitivityVisibleParameterNames(fitMetrics));
  const ranking = (fitMetrics.ranking_by_weighted_sensitivity || [])
    .filter((entry) => entry && visibleNames.has(entry.name));

  if (!ranking.length) {
    renderSensitivityEmptyState(sensitivityRankingChart, 'No ranked parameters remain after hiding weak entries.');
    return;
  }

  const names = ranking.map((entry) => entry.name);
  const sensitivityValues = ranking.map((entry) => Number(entry.weighted_l2_sensitivity));
  const customdata = ranking.map((entry) => {
    const parameterMetrics = getFitMetricForParameter(fitMetrics, entry.name) || {};
    return [
      entry.name,
      Number(entry.weighted_l2_sensitivity),
      Number(entry.relative_weighted_l2_sensitivity),
      Number(entry.sigma),
      parameterMetrics?.most_correlated_with || 'n/a',
    ];
  });

  const trace = {
    type: 'bar',
    orientation: 'h',
    x: sensitivityValues,
    y: names,
    marker: {
      color: getSensitivityBarColors(names),
      line: {
        color: 'rgba(15, 23, 42, 0.18)',
        width: 1,
      },
    },
    customdata,
    hovertemplate:
      '<b>%{customdata[0]}</b><br>' +
      'Sensitivity: %{customdata[1]:.6g}<br>' +
      'Relative sensitivity: %{customdata[2]:.6g}<br>' +
      'Sigma: %{customdata[3]:.6g}<extra></extra>',
  };
  const layout = {
    margin: { l: 96, r: 24, t: 12, b: 42 },
    height: Math.max(280, names.length * 28 + 90),
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#ffffff',
    showlegend: false,
    xaxis: {
      title: { text: 'Weighted L2 sensitivity' },
      gridcolor: '#e2e8f0',
      zerolinecolor: '#cbd5e1',
    },
    yaxis: {
      automargin: true,
      autorange: 'reversed',
    },
    hoverlabel: {
      bgcolor: '#0f172a',
      bordercolor: '#0f172a',
      font: { color: '#f8fafc' },
    },
  };

  if (typeof Plotly.purge === 'function') {
    Plotly.purge(sensitivityRankingChart);
  }
  sensitivityRankingChart.innerHTML = '';
  Plotly.react(sensitivityRankingChart, [trace], layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  if (typeof sensitivityRankingChart.removeAllListeners === 'function') {
    sensitivityRankingChart.removeAllListeners('plotly_hover');
    sensitivityRankingChart.removeAllListeners('plotly_unhover');
    sensitivityRankingChart.removeAllListeners('plotly_click');
  }
  sensitivityRankingChart.on('plotly_hover', (event) => {
    const name = event?.points?.[0]?.y;
    setSensitivityHoveredParam(name);
  });
  sensitivityRankingChart.on('plotly_unhover', () => setSensitivityHoveredParam(null));
  sensitivityRankingChart.on('plotly_click', (event) => {
    const name = event?.points?.[0]?.y;
    if (name) selectSensitivityParameter(name);
  });
}

function renderCorrelationHeatmap() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivityCorrelationHeatmap) return;
  if (!fitMetrics || typeof Plotly === 'undefined') {
    renderSensitivityEmptyState(sensitivityCorrelationHeatmap, 'Run a fit to inspect parameter correlations.');
    return;
  }

  const visibleNames = getSensitivityVisibleParameterNames(fitMetrics);
  if (!visibleNames.length) {
    renderSensitivityEmptyState(sensitivityCorrelationHeatmap, 'No parameters remain after hiding weak entries.');
    return;
  }

  const allNames = getFitMetricsParameterNames(fitMetrics);
  const rawMatrix = Array.isArray(fitMetrics.correlation) ? fitMetrics.correlation : [];
  const indexMap = new Map(allNames.map((name, index) => [name, index]));
  const z = visibleNames.map((rowName) => {
    const rowIndex = indexMap.get(rowName);
    return visibleNames.map((colName) => {
      const colIndex = indexMap.get(colName);
      const value = rawMatrix?.[rowIndex]?.[colIndex];
      return Number.isFinite(Number(value)) ? Number(value) : null;
    });
  });
  const customdata = visibleNames.map((rowName) => (
    visibleNames.map((colName) => [rowName, colName])
  ));
  const indices = visibleNames.map((_, index) => index);
  const trace = {
    type: 'heatmap',
    x: indices,
    y: indices,
    z,
    zmin: -1,
    zmax: 1,
    colorscale: [
      [0, '#2563eb'],
      [0.5, '#ffffff'],
      [1, '#dc2626'],
    ],
    customdata,
    xgap: 1,
    ygap: 1,
    hovertemplate:
      '<b>%{customdata[0]}</b> vs <b>%{customdata[1]}</b><br>' +
      'Correlation: %{z:.4f}<extra></extra>',
    colorbar: {
      title: 'corr',
      tickvals: [-1, -0.5, 0, 0.5, 1],
    },
  };
  const tickText = getSensitivityHeatmapTickText(visibleNames);
  const layout = {
    margin: { l: 90, r: 30, t: 12, b: 84 },
    height: Math.max(320, visibleNames.length * 30 + 120),
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#ffffff',
    xaxis: {
      tickmode: 'array',
      tickvals: indices,
      ticktext: tickText,
      side: 'bottom',
      tickangle: -35,
      automargin: true,
      constrain: 'domain',
      scaleanchor: 'y',
    },
    yaxis: {
      tickmode: 'array',
      tickvals: indices,
      ticktext: tickText,
      autorange: 'reversed',
      automargin: true,
      constrain: 'domain',
    },
    shapes: buildSensitivityHeatmapShapes(visibleNames),
    hoverlabel: {
      bgcolor: '#0f172a',
      bordercolor: '#0f172a',
      font: { color: '#f8fafc' },
    },
  };

  if (typeof Plotly.purge === 'function') {
    Plotly.purge(sensitivityCorrelationHeatmap);
  }
  sensitivityCorrelationHeatmap.innerHTML = '';
  Plotly.react(sensitivityCorrelationHeatmap, [trace], layout, withHighResDownloadConfig({ responsive: true, displaylogo: false }));
  if (typeof sensitivityCorrelationHeatmap.removeAllListeners === 'function') {
    sensitivityCorrelationHeatmap.removeAllListeners('plotly_click');
    sensitivityCorrelationHeatmap.removeAllListeners('plotly_hover');
    sensitivityCorrelationHeatmap.removeAllListeners('plotly_unhover');
  }
  sensitivityCorrelationHeatmap.on('plotly_hover', (event) => {
    const point = event?.points?.[0];
    const rowName = point?.customdata?.[0] || null;
    setSensitivityHoveredParam(rowName);
  });
  sensitivityCorrelationHeatmap.on('plotly_unhover', () => setSensitivityHoveredParam(null));
  sensitivityCorrelationHeatmap.on('plotly_click', (event) => {
    const point = event?.points?.[0];
    const rowName = point?.customdata?.[0] || null;
    const colName = point?.customdata?.[1] || null;
    if (!rowName || !colName) return;
    selectSensitivityParameter(rowName, rowName === colName ? [rowName] : [rowName, colName]);
  });
}

function renderParameterDetail() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivityParameterDetail) return;
  if (!fitMetrics) {
    renderSensitivityEmptyState(sensitivityParameterDetail, 'Select a parameter after running a fit.');
    return;
  }
  ensureSensitivitySelection(fitMetrics);
  const name = autofitSensitivityState.selectedParam;
  const metrics = getFitMetricForParameter(fitMetrics, name);
  if (!name || !metrics) {
    renderSensitivityEmptyState(sensitivityParameterDetail, 'Select a parameter to inspect its local fit diagnostics.');
    return;
  }

  sensitivityParameterDetail.innerHTML = '';
  const list = document.createElement('dl');
  list.className = 'sensitivity-detail-list';
  [
    ['Parameter', name],
    ['Value', formatDiagnosticNumber(metrics.value, 6)],
    ['Sigma', formatDiagnosticNumber(metrics.sigma, 4)],
    ['Sensitivity', formatDiagnosticNumber(metrics.weighted_l2_sensitivity, 4)],
    ['Relative sensitivity', formatDiagnosticNumber(metrics.relative_weighted_l2_sensitivity, 4)],
    ['Most correlated parameter', metrics.most_correlated_with || 'n/a'],
    ['Correlation value', formatDiagnosticNumber(metrics.max_abs_correlation, 4)],
  ].forEach(([label, value]) => {
    const dt = document.createElement('dt');
    dt.textContent = label;
    const dd = document.createElement('dd');
    dd.textContent = value;
    list.appendChild(dt);
    list.appendChild(dd);
  });
  sensitivityParameterDetail.appendChild(list);
}

function getConditionBadgeConfig(conditionNumber) {
  const value = Number(conditionNumber);
  if (Number.isFinite(value) && value < SENSITIVITY_COND_THRESHOLDS.good) {
    return { label: 'stable', className: 'sensitivity-diagnostic-badge sensitivity-diagnostic-badge--good' };
  }
  if (Number.isFinite(value) && value < SENSITIVITY_COND_THRESHOLDS.warn) {
    return { label: 'moderate', className: 'sensitivity-diagnostic-badge sensitivity-diagnostic-badge--warn' };
  }
  return { label: 'ill-conditioned', className: 'sensitivity-diagnostic-badge sensitivity-diagnostic-badge--bad' };
}

function renderDiagnostics() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivityDiagnostics) return;
  if (!fitMetrics) {
    renderSensitivityEmptyState(sensitivityDiagnostics, 'Fit stability metrics will appear after the next run.');
    return;
  }

  const conditionNumber = Number(fitMetrics.condition_number);
  const sigma2 = Number(fitMetrics.sigma2);
  const dof = Number(fitMetrics.dof);
  const badge = getConditionBadgeConfig(conditionNumber);

  sensitivityDiagnostics.innerHTML = '';
  const summary = document.createElement('div');
  summary.className = 'sensitivity-diagnostics-summary';
  const badgeEl = document.createElement('span');
  badgeEl.className = badge.className;
  badgeEl.textContent = badge.label;
  summary.appendChild(badgeEl);
  sensitivityDiagnostics.appendChild(summary);

  const stats = document.createElement('div');
  stats.className = 'sensitivity-diagnostics-grid';
  [
    ['Condition number', formatDiagnosticNumber(conditionNumber, 4)],
    ['Sigma²', formatDiagnosticNumber(sigma2, 4)],
    ['DOF', Number.isFinite(dof) ? String(dof) : 'n/a'],
  ].forEach(([label, value]) => {
    const card = document.createElement('div');
    card.className = 'sensitivity-diagnostic-stat';
    const key = document.createElement('span');
    key.className = 'sensitivity-diagnostic-stat__label';
    key.textContent = label;
    const val = document.createElement('strong');
    val.className = 'sensitivity-diagnostic-stat__value';
    val.textContent = value;
    card.appendChild(key);
    card.appendChild(val);
    stats.appendChild(card);
  });
  sensitivityDiagnostics.appendChild(stats);
}

function buildSensitivitySuggestions(fitMetrics) {
  const suggestions = [];
  getFitMetricsParameterNames(fitMetrics).forEach((name) => {
    const metrics = getFitMetricForParameter(fitMetrics, name);
    if (metrics?.weakly_sensitive) {
      suggestions.push(`Consider fixing ${name}`);
    }
  });

  const highCorrelationPairs = Array.isArray(fitMetrics?.high_correlation_pairs)
    ? fitMetrics.high_correlation_pairs
    : [];
  highCorrelationPairs.forEach((pair) => {
    const correlation = Math.abs(Number(pair?.correlation));
    if (correlation > 0.95 && pair?.param_i && pair?.param_j) {
      suggestions.push(`Strong correlation: ${pair.param_i} ↔ ${pair.param_j}`);
    }
  });

  const conditionNumber = Number(fitMetrics?.condition_number);
  if (Number.isFinite(conditionNumber) && conditionNumber >= SENSITIVITY_COND_THRESHOLDS.warn) {
    suggestions.push('Fit is ill-conditioned');
  }

  return suggestions.filter((item, index, values) => item && values.indexOf(item) === index);
}

function renderSuggestions() {
  const fitMetrics = lastAutofitFitMetrics;
  if (!sensitivitySuggestions) return;
  if (!fitMetrics) {
    renderSensitivityEmptyState(sensitivitySuggestions, 'Automatic suggestions will appear after the next fit.');
    return;
  }

  const items = buildSensitivitySuggestions(fitMetrics);
  if (!items.length) {
    renderSensitivityEmptyState(sensitivitySuggestions, 'No immediate sensitivity warnings were detected for this fit.');
    return;
  }

  sensitivitySuggestions.innerHTML = '';
  const list = document.createElement('ul');
  list.className = 'sensitivity-suggestions-list';
  items.forEach((item) => {
    const li = document.createElement('li');
    li.textContent = item;
    list.appendChild(li);
  });
  sensitivitySuggestions.appendChild(list);
}

function renderSensitivityDashboard(fitMetrics) {
  lastAutofitFitMetrics = fitMetrics && typeof fitMetrics === 'object'
    ? cloneJsonValue(fitMetrics)
    : null;

  if (sensitivityHideWeakToggle) {
    sensitivityHideWeakToggle.checked = autofitSensitivityState.hideWeak;
  }
  if (!sensitivityDashboard || !sensitivityPanel) return;
  if (!lastAutofitFitMetrics || !getFitMetricsParameterNames(lastAutofitFitMetrics).length) {
    sensitivityPanel.hidden = true;
    renderParameterTable();
    renderSensitivityChart();
    renderCorrelationHeatmap();
    renderParameterDetail();
    renderDiagnostics();
    renderSuggestions();
    return;
  }

  sensitivityPanel.hidden = false;
  if (!sensitivityPanel.open) {
    return;
  }
  ensureSensitivitySelection(lastAutofitFitMetrics);
  renderParameterTable();
  renderSensitivityChart();
  renderCorrelationHeatmap();
  renderParameterDetail();
  renderDiagnostics();
  renderSuggestions();
  syncSensitivityHighlights();
  scheduleSensitivityPlotResize();
}

function getReadyStateLabel(state) {
  if (state === 'ready') return 'Ready For SPFIT';
  if (state === 'review') return 'Review Before SPFIT';
  return 'Not Ready For SPFIT';
}

function getReadyStateClassName(state) {
  if (state === 'ready') return 'autofit-ready-badge autofit-ready-badge--ready';
  if (state === 'review') return 'autofit-ready-badge autofit-ready-badge--review';
  return 'autofit-ready-badge autofit-ready-badge--not-ready';
}

function buildAssignmentDiagnostics(assignments, expPeaks, options = {}) {
  const deltaAssignPct = Number(options.deltaAssignPct);
  const manualMode = Boolean(options.manualMode);
  const exportSigma = Number.isFinite(Number(options.exportSigmaMHz))
    ? Number(options.exportSigmaMHz)
    : SPFIT_DEFAULT_LINE_ERROR_MHZ;
  const residualAbsValues = assignments
    .map((item) => Math.abs(Number(item?.residual)))
    .filter(Number.isFinite);
  const medianResidual = Number.isFinite(computeMedian(residualAbsValues)) ? computeMedian(residualAbsValues) : exportSigma;
  const madResidual = Number.isFinite(computeMad(residualAbsValues, medianResidual)) ? computeMad(residualAbsValues, medianResidual) : exportSigma;
  const residualWarnThreshold = Math.max(exportSigma * 2.5, medianResidual + 3 * Math.max(madResidual, exportSigma));
  const residualRejectThreshold = Math.max(exportSigma * 5, medianResidual + 6 * Math.max(madResidual, exportSigma));
  const alternativeGapThreshold = Math.max(exportSigma * 2, medianResidual * 2, 0.05);

  return assignments.map((assignment) => {
    const observed = Number(assignment?.fObs);
    const residual = Number(assignment?.residual);
    const predicted = Number.isFinite(observed) && Number.isFinite(residual) ? observed - residual : Number.NaN;
    const candidateWindow = Number.isFinite(predicted)
      ? Math.min(relativeToleranceMHz(predicted, deltaAssignPct), Math.max(0.35, exportSigma * 8))
      : Math.max(0.35, exportSigma * 8);
    const candidates = Array.isArray(expPeaks)
      ? expPeaks
        .filter((peak) => Number.isFinite(Number(peak?.freq)) && Math.abs(Number(peak.freq) - predicted) <= candidateWindow)
        .map((peak) => ({
          ...peak,
          diffToPredicted: Math.abs(Number(peak.freq) - predicted),
          diffToObserved: Math.abs(Number(peak.freq) - observed),
        }))
        .sort((a, b) => a.diffToPredicted - b.diffToPredicted)
      : [];
    const nearestObserved = findNearestPeak(candidates.length ? candidates : expPeaks, observed);
    const assignedPeak = nearestObserved?.peak || null;
    const peakQuality = assignedPeak && Number.isFinite(Number(assignedPeak.q)) ? Number(assignedPeak.q) : 1;
    const secondCandidate = candidates.find((candidate) => {
      if (!assignedPeak) return true;
      return Math.abs(Number(candidate.freq) - Number(assignedPeak.freq)) > 1e-6;
    }) || null;
    const alternativeGap = secondCandidate && Number.isFinite(residual)
      ? Math.abs(Number(secondCandidate.freq) - predicted) - Math.abs(residual)
      : null;
    const reasons = [];
    let status = 'solid';
    if (!Number.isFinite(residual)) {
      status = 'reject';
      reasons.push('missing final residual');
    } else if (Math.abs(residual) > residualRejectThreshold) {
      status = 'reject';
      reasons.push(`residual ${formatDiagnosticNumber(Math.abs(residual), 4)} MHz exceeds reject threshold`);
    } else if (peakQuality < 0.35) {
      status = 'reject';
      reasons.push(`peak quality q=${formatDiagnosticNumber(peakQuality, 3)} is too low`);
    } else {
      if (alternativeGap != null && alternativeGap <= alternativeGapThreshold) {
        status = 'ambiguous';
        reasons.push(`small nearest-peak margin (${formatDiagnosticNumber(alternativeGap, 3)} MHz)`);
      }
      if (Math.abs(residual) > residualWarnThreshold) {
        status = 'ambiguous';
        reasons.push(`residual ${formatDiagnosticNumber(Math.abs(residual), 4)} MHz is elevated`);
      }
      if (peakQuality < 0.85) {
        status = 'ambiguous';
        reasons.push(`peak quality q=${formatDiagnosticNumber(peakQuality, 3)} is moderate`);
      }
      if (!reasons.length) {
        reasons.push('clean residual and no close alternative peak');
      }
    }
    const exportable = manualMode ? status !== 'reject' : status === 'solid';
    return {
      label: assignment.label,
      fObs: observed,
      residual,
      predictedFreq: predicted,
      candidateWindow,
      candidateCount: candidates.length,
      peakQuality,
      alternativeGap,
      status,
      exportable,
      reason: reasons.join('; '),
    };
  });
}

function buildAutofitDiagnostics({ outputText, expSpectrum } = {}) {
  const assignments = parseAutofitFinalAssignments(outputText);
  if (!assignments.length) return null;
  const context = lastAutofitContext || {};
  const expMode = context.expMode || autofitExpMode?.value || 'profile';
  const expPeaks = getAutofitExperimentalPeaks(expSpectrum, expMode);
  const rmsMetrics = parseAutofitRmsMetrics(outputText);
  const controlMetrics = parseAutofitControlMetrics(outputText);
  const optimizer = parseAutofitOptimizerDiagnostics(outputText);
  const fitParamNames = parseAutofitFitParamNames(outputText);
  const assignmentDiagnostics = buildAssignmentDiagnostics(assignments, expPeaks, {
    deltaAssignPct: context?.optionValues?.deltaAssign ?? autofitOptions.deltaAssign?.value ?? 0.1,
    manualMode: Boolean(context.manualMode),
    exportSigmaMHz: SPFIT_DEFAULT_LINE_ERROR_MHZ,
  });
  const exportableAssignments = assignmentDiagnostics
    .filter((item) => item.exportable)
    .map((item) => ({
      label: item.label,
      fObs: item.fObs,
      residual: item.residual,
    }));
  const severe = [];
  const warnings = [];
  const fitParamCount = fitParamNames.length || 3;
  const ambiguousCount = assignmentDiagnostics.filter((item) => item.status === 'ambiguous').length;
  const rejectCount = assignmentDiagnostics.filter((item) => item.status === 'reject').length;
  const exportableCount = exportableAssignments.length;
  const medianPeakQuality = computeMedian(
    assignmentDiagnostics
      .map((item) => Number(item.peakQuality))
      .filter(Number.isFinite)
  );

  if (!optimizer.leastSquares || Number(optimizer.leastSquares.success) !== 1) {
    severe.push('Final least-squares did not report success.');
  }
  if (Number.isFinite(Number(rmsMetrics.finalRms))) {
    const finalRms = Number(rmsMetrics.finalRms);
    if (finalRms > Math.max(0.25, SPFIT_DEFAULT_LINE_ERROR_MHZ * 6)) {
      severe.push(`Final RMS(labels) is high: ${formatDiagnosticNumber(finalRms, 4)} MHz.`);
    } else if (finalRms > Math.max(0.12, SPFIT_DEFAULT_LINE_ERROR_MHZ * 3)) {
      warnings.push(`Final RMS(labels) is not yet tight: ${formatDiagnosticNumber(finalRms, 4)} MHz.`);
    }
  } else {
    severe.push('Final RMS(labels) could not be parsed.');
  }
  if (exportableCount < Math.max(3, fitParamCount + 1)) {
    severe.push(`Too few exportable lines (${exportableCount}) for ${fitParamCount} fitted parameters.`);
  } else if (exportableCount < fitParamCount + 3) {
    warnings.push(`Line leverage is still limited: ${exportableCount} exportable lines for ${fitParamCount} fitted parameters.`);
  }
  if (ambiguousCount > 0) {
    warnings.push(`${ambiguousCount} fitted line${ambiguousCount === 1 ? '' : 's'} remain ambiguous.`);
  }
  if (rejectCount > 0) {
    warnings.push(`${rejectCount} fitted line${rejectCount === 1 ? '' : 's'} were rejected from automatic SPFIT export.`);
  }
  const lsOptimality = Number(optimizer.leastSquares?.optimality);
  if (Number.isFinite(lsOptimality)) {
    if (lsOptimality > 0.05) {
      severe.push(`Least-squares optimality is still high (${formatDiagnosticNumber(lsOptimality, 3)}).`);
    } else if (lsOptimality > 0.005) {
      warnings.push(`Least-squares optimality is moderate (${formatDiagnosticNumber(lsOptimality, 3)}).`);
    }
  }
  if (!context.manualMode) {
    const controlCount = Number(controlMetrics?.n);
    if (!controlMetrics || !Number.isFinite(controlCount) || controlCount <= 0) {
      warnings.push('No usable control-line metrics were produced.');
    } else {
      const controlMad = Number(controlMetrics.mad);
      const controlBias = Math.abs(Number(controlMetrics.bias));
      if (Number.isFinite(controlCount) && controlCount < 2) {
        warnings.push(`Only ${controlCount} control line contributed to the ranking.`);
      }
      if (Number.isFinite(controlMad) && controlMad > Math.max(0.1, SPFIT_DEFAULT_LINE_ERROR_MHZ * 3)) {
        warnings.push(`Control MAD is broad: ${formatDiagnosticNumber(controlMad, 4)} MHz.`);
      }
      if (Number.isFinite(controlBias) && controlBias > 0.1) {
        warnings.push(`Control bias is non-negligible: ${formatDiagnosticNumber(controlBias, 4)} MHz.`);
      }
    }
  }
  if (fitParamCount > 3 && exportableCount < fitParamCount * 2) {
    warnings.push('Optional fitted parameters may be weakly constrained by the current accepted line set.');
  }
  if (context?.optionValues?.usePowell && optimizer.powell && Number(optimizer.powell.success) !== 1) {
    warnings.push('Powell pre-refinement did not converge cleanly before least-squares.');
  }

  const readyState = severe.length ? 'not-ready' : (warnings.length ? 'review' : 'ready');
  return {
    readyState,
    readyLabel: getReadyStateLabel(readyState),
    rmsMetrics,
    controlMetrics,
    optimizer,
    fitParamNames,
    fitParamCount,
    assignmentDiagnostics,
    exportableAssignments,
    assignmentCount: assignments.length,
    exportableCount,
    ambiguousCount,
    rejectCount,
    medianPeakQuality: Number.isFinite(medianPeakQuality) ? medianPeakQuality : null,
    exportSigmaMHz: SPFIT_DEFAULT_LINE_ERROR_MHZ,
    severe,
    warnings,
    manualMode: Boolean(context.manualMode),
  };
}

function createHelpDisclosure({ title, bodyHtml, ariaLabel } = {}) {
  const details = document.createElement('details');
  details.className = 'help-disclosure';

  const summary = document.createElement('summary');
  summary.className = 'help-button';
  summary.textContent = '?';
  summary.setAttribute('aria-label', ariaLabel || (title ? `Show help for ${title}` : 'Show help'));
  details.appendChild(summary);

  const popover = document.createElement('div');
  popover.className = 'help-popover';
  if (title) {
    const heading = document.createElement('h4');
    heading.textContent = title;
    popover.appendChild(heading);
  }
  if (bodyHtml) {
    const body = document.createElement('div');
    body.className = 'help-popover__body';
    body.innerHTML = bodyHtml;
    popover.appendChild(body);
  }
  details.appendChild(popover);
  return details;
}

function appendDiagnosticStat(container, label, value, help = null) {
  const card = document.createElement('div');
  card.className = 'autofit-diagnostic-stat';
  const top = document.createElement('div');
  top.className = 'autofit-diagnostic-stat__top';
  const labelEl = document.createElement('span');
  labelEl.className = 'autofit-diagnostic-stat__label';
  labelEl.textContent = label;
  top.appendChild(labelEl);
  if (help) {
    top.appendChild(createHelpDisclosure({
      ...help,
      ariaLabel: `Show help for ${label}`,
    }));
  }
  const valueEl = document.createElement('span');
  valueEl.className = 'autofit-diagnostic-stat__value';
  valueEl.textContent = value;
  card.appendChild(top);
  card.appendChild(valueEl);
  container.appendChild(card);
}

function appendDiagnosticList(container, title, items, variant = 'warning') {
  if (!Array.isArray(items) || !items.length) return;
  const group = document.createElement('div');
  group.className = 'autofit-diagnostic-group';
  const heading = document.createElement('h4');
  heading.textContent = title;
  const list = document.createElement('ul');
  list.className = `autofit-diagnostic-list autofit-diagnostic-list--${variant}`;
  items.forEach((item) => {
    const li = document.createElement('li');
    li.textContent = item;
    list.appendChild(li);
  });
  group.appendChild(heading);
  group.appendChild(list);
  container.appendChild(group);
}

function renderAutofitDiagnostics(diagnostics) {
  lastAutofitDiagnostics = diagnostics || null;
  if (!autofitDiagnostics) return;
  autofitDiagnostics.innerHTML = '';
  if (!diagnostics) {
    autofitDiagnostics.hidden = true;
    return;
  }
  autofitDiagnostics.hidden = false;

  const header = document.createElement('div');
  header.className = 'autofit-diagnostics__header';
  const titleBlock = document.createElement('div');
  const title = document.createElement('h3');
  title.className = 'autofit-diagnostics__title';
  title.textContent = 'Fit Diagnostics';
  const subtitle = document.createElement('p');
  subtitle.className = 'autofit-diagnostics__subtitle';
  subtitle.textContent = diagnostics.manualMode
    ? 'Manual mode exports non-rejected lines; automatic mode exports only solid lines.'
    : 'Only solid lines are promoted automatically to the SPFIT seed package.';
  titleBlock.appendChild(title);
  titleBlock.appendChild(subtitle);
  const badge = document.createElement('span');
  badge.className = getReadyStateClassName(diagnostics.readyState);
  badge.textContent = diagnostics.readyLabel;
  header.appendChild(titleBlock);
  header.appendChild(badge);
  autofitDiagnostics.appendChild(header);

  const stats = document.createElement('div');
  stats.className = 'autofit-diagnostics__stats';
  appendDiagnosticStat(stats, 'Fitted Lines', String(diagnostics.assignmentCount), AUTOFIT_DIAGNOSTIC_HELP.fittedLines);
  appendDiagnosticStat(stats, 'Exportable Lines', String(diagnostics.exportableCount), AUTOFIT_DIAGNOSTIC_HELP.exportableLines);
  appendDiagnosticStat(stats, 'Fit Parameters', String(diagnostics.fitParamCount), AUTOFIT_DIAGNOSTIC_HELP.fitParameters);
  appendDiagnosticStat(
    stats,
    'RMS(labels)',
    `${formatDiagnosticNumber(diagnostics.rmsMetrics.finalRms, 4)} MHz`,
    AUTOFIT_DIAGNOSTIC_HELP.rmsLabels
  );
  appendDiagnosticStat(
    stats,
    'Control Score',
    diagnostics.controlMetrics && Number.isFinite(Number(diagnostics.controlMetrics.score))
      ? formatDiagnosticNumber(diagnostics.controlMetrics.score, 4)
      : 'n/a',
    AUTOFIT_DIAGNOSTIC_HELP.controlScore
  );
  appendDiagnosticStat(
    stats,
    'Peak Quality q',
    Number.isFinite(Number(diagnostics.medianPeakQuality))
      ? formatDiagnosticNumber(diagnostics.medianPeakQuality, 3)
      : 'n/a',
    AUTOFIT_DIAGNOSTIC_HELP.peakQuality
  );
  appendDiagnosticStat(
    stats,
    'Export Sigma',
    `${formatDiagnosticNumber(diagnostics.exportSigmaMHz, 4)} MHz`,
    AUTOFIT_DIAGNOSTIC_HELP.exportSigma
  );
  appendDiagnosticStat(
    stats,
    'Ambiguous / Reject',
    `${diagnostics.ambiguousCount} / ${diagnostics.rejectCount}`,
    AUTOFIT_DIAGNOSTIC_HELP.ambiguousReject
  );
  autofitDiagnostics.appendChild(stats);

  appendDiagnosticList(autofitDiagnostics, 'Critical Issues', diagnostics.severe, 'critical');
  appendDiagnosticList(autofitDiagnostics, 'Warnings', diagnostics.warnings, 'warning');

  const tableGroup = document.createElement('div');
  tableGroup.className = 'autofit-diagnostic-group';
  const tableHeading = document.createElement('h4');
  tableHeading.textContent = 'Assignment Confidence';
  tableGroup.appendChild(tableHeading);
  const tableWrap = document.createElement('div');
  tableWrap.className = 'autofit-assignments-wrap';
  const table = document.createElement('table');
  table.className = 'autofit-assignments-table';
  const thead = document.createElement('thead');
  const headRow = document.createElement('tr');
  ['Status', 'Export', 'Label', 'Residual', 'q', 'Alt gap', 'Candidates', 'Reason'].forEach((label) => {
    const th = document.createElement('th');
    th.textContent = label;
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);
  table.appendChild(thead);
  const tbody = document.createElement('tbody');
  diagnostics.assignmentDiagnostics.forEach((item) => {
    const row = document.createElement('tr');

    const statusCell = document.createElement('td');
    const status = document.createElement('span');
    status.className = `autofit-assignment-status autofit-assignment-status--${item.status}`;
    status.textContent = item.status;
    statusCell.appendChild(status);
    row.appendChild(statusCell);

    const exportCell = document.createElement('td');
    exportCell.className = `autofit-export-flag ${item.exportable ? 'autofit-export-flag--yes' : 'autofit-export-flag--no'}`;
    exportCell.textContent = item.exportable ? 'Yes' : 'No';
    row.appendChild(exportCell);

    const labelCell = document.createElement('td');
    labelCell.textContent = item.label;
    row.appendChild(labelCell);

    const residualCell = document.createElement('td');
    residualCell.textContent = `${formatDiagnosticNumber(Math.abs(item.residual), 4)} MHz`;
    row.appendChild(residualCell);

    const qCell = document.createElement('td');
    qCell.textContent = formatDiagnosticNumber(item.peakQuality, 3);
    row.appendChild(qCell);

    const gapCell = document.createElement('td');
    gapCell.textContent = item.alternativeGap == null ? 'n/a' : `${formatDiagnosticNumber(item.alternativeGap, 3)} MHz`;
    row.appendChild(gapCell);

    const candidateCell = document.createElement('td');
    candidateCell.textContent = String(item.candidateCount);
    row.appendChild(candidateCell);

    const reasonCell = document.createElement('td');
    reasonCell.textContent = item.reason;
    row.appendChild(reasonCell);

    tbody.appendChild(row);
  });
  table.appendChild(tbody);
  tableWrap.appendChild(table);
  tableGroup.appendChild(tableWrap);
  autofitDiagnostics.appendChild(tableGroup);
}

function parseAutofitInputMetadata(text) {
  const emptyMeta = {
    raw: null,
    reduction: null,
    reductionSource: 'missing',
    rotorType: null,
    representation: null,
    jMax: null,
    sD2Id: null,
    hfs: {
      nuclei: [],
    },
  };
  if (typeof text !== 'string' || !text.trim()) return emptyMeta;

  const raw = parseAutofitRotationalSection(text);
  if (!raw) return emptyMeta;

  const meta = buildAutofitInputMetadataFromRaw(raw);
  if (meta.reduction) return meta;

  const inferredText = String(text);
  if (/\b(delta J_MHz|delta K_MHz|DELTA J_MHz|DELTA K_MHz|PHI N_MHz|phi N_MHz)\b/.test(inferredText)) {
    meta.reduction = 'A';
    meta.reductionSource = 'inferred';
    return meta;
  }
  if (/\b(DJ_MHz|DJK_MHz|DK_MHz|d1_MHz|d2_MHz|H N_MHz|h1_MHz)\b/.test(inferredText)) {
    meta.reduction = 'S';
    meta.reductionSource = 'inferred';
  }
  return meta;
}

function parseAutofitRotationalSection(text) {
  if (typeof text !== 'string' || !text.trim()) return null;
  const raw = new Map();
  let inRotational = false;
  String(text).split(/\r?\n/).forEach((rawLine) => {
    const line = String(rawLine || '').trim();
    if (!line) return;
    if (line.startsWith('#')) {
      inRotational = line.toUpperCase() === '#ROTATIONAL';
      return;
    }
    if (!inRotational) return;
    const idx = line.indexOf('=');
    if (idx < 0) return;
    const key = line.slice(0, idx).trim().replace(/\s+/g, ' ');
    const value = line.slice(idx + 1).trim();
    if (!key) return;
    raw.set(key, value);
  });
  return raw.size ? raw : null;
}

function readAutofitInputNumberBySuffix(readNumber, baseKey, suffixes) {
  for (let i = 0; i < suffixes.length; i += 1) {
    const suffix = suffixes[i];
    const aliases = suffix
      ? [`${baseKey}${suffix}_MHz`]
      : [`${baseKey}_MHz`];
    const value = readNumber(`${baseKey}${suffix}`, aliases);
    if (Number.isFinite(value)) return value;
  }
  return 0;
}

function parseNonNegativeInteger(value) {
  const parsed = Number.parseInt(value, 10);
  return Number.isInteger(parsed) && parsed >= 0 ? parsed : null;
}

function normalizeSpfitRepresentationToken(value) {
  const raw = String(value || '').trim().replace(/\s+/g, '');
  if (!raw) return null;
  const lower = raw.toLowerCase();
  if (lower === 'ir') return 'Ir';
  if (lower === 'iir') return 'IIr';
  if (lower === 'iiir') return 'IIIr';
  if (lower === 'il') return 'Il';
  if (lower === 'iil') return 'IIl';
  if (lower === 'iiil') return 'IIIl';
  return raw;
}

function buildAutofitInputMetadataFromRaw(raw) {
  const readRaw = (key, alternates = []) => {
    if (!raw || !(raw instanceof Map)) return null;
    if (raw.has(key)) return raw.get(key);
    for (let i = 0; i < alternates.length; i += 1) {
      if (raw.has(alternates[i])) return raw.get(alternates[i]);
    }
    return null;
  };
  const readNumber = (key, alternates = []) => {
    const value = readRaw(key, alternates);
    if (value == null || value === '') return null;
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : null;
  };

  const rawRotorType = readRaw('rotor_type', ['rotor type']);
  let rotorType = null;
  if (rawRotorType) {
    const rotorLower = String(rawRotorType).trim().toLowerCase();
    if (rotorLower.includes('asymmetric')) rotorType = 'asymmetric';
    else if (rotorLower.includes('linear')) rotorType = 'linear';
    else if (rotorLower.includes('spherical')) rotorType = 'spherical';
    else if (rotorLower.includes('oblate')) rotorType = 'oblate';
    else if (rotorLower.includes('prolate')) rotorType = 'prolate';
  }

  const reductionRaw = readRaw('Watson Reduction');
  let reduction = null;
  if (reductionRaw) {
    const candidate = String(reductionRaw).trim().toUpperCase();
    if (candidate.startsWith('A')) reduction = 'A';
    else if (candidate.startsWith('S')) reduction = 'S';
  }

  return {
    raw: Object.fromEntries(raw.entries()),
    reduction,
    reductionSource: reduction ? 'explicit' : 'missing',
    rotorType,
    representation: normalizeSpfitRepresentationToken(readRaw('representation')),
    jMax: parseNonNegativeInteger(readRaw('J_MAX', ['J_max'])),
    sD2Id: parseNonNegativeInteger(readRaw('S_D2_ID', ['s_d2_id', 'S_d2_ID'])),
    hfs: {
      nuclei: [1, 2].map((index) => {
        const spin = index === 1
          ? readNumber('I_NUC_1', ['I_nuc_1', 'I_NUC', 'I_nuc'])
          : readNumber('I_NUC_2', ['I_nuc_2']);
        const suffixes = index === 1 ? ['_1', ''] : ['_2'];
        const chi = {
          aa: readAutofitInputNumberBySuffix(readNumber, 'chi_aa', suffixes),
          bb: readAutofitInputNumberBySuffix(readNumber, 'chi_bb', suffixes),
          cc: readAutofitInputNumberBySuffix(readNumber, 'chi_cc', suffixes),
          ab: readAutofitInputNumberBySuffix(readNumber, 'chi_ab', suffixes),
          ac: readAutofitInputNumberBySuffix(readNumber, 'chi_ac', suffixes),
          bc: readAutofitInputNumberBySuffix(readNumber, 'chi_bc', suffixes),
        };
        const spinValue = Number.isFinite(spin) ? Number(spin) : 0;
        const active = spinValue > 0.5 && Object.values(chi).some((value) => Math.abs(Number(value) || 0) > SPFIT_PARAM_INCLUDE_EPS);
        return {
          index,
          spin: spinValue,
          chi,
          active,
        };
      }),
    },
  };
}

function getSpfitRepresentationHeaderSign(representation) {
  const token = String(representation || '').trim().toLowerCase();
  return token.endsWith('l') ? -1 : 1;
}

function getSpfitRepresentationAxisCode(representation) {
  const token = String(representation || '').trim().toLowerCase();
  if (token.startsWith('iii')) return 3;
  if (token.startsWith('ii')) return 2;
  return 1;
}

function getSpfitParamGroup(key) {
  if (key === 'A' || key === 'B' || key === 'C') return 'rotational';
  if (QUARTIC_PARAM_KEYS.includes(key)) return 'quartic';
  if (SEXTIC_PARAM_KEYS.includes(key)) return 'sextic';
  return 'quartic';
}

function computeSpfitFloatingParamError(key, value, initialValue = null) {
  const numericValue = Number(value);
  const numericInitial = Number(initialValue);
  const group = getSpfitParamGroup(key);
  const floor = SPFIT_FIT_ERROR_FLOORS_MHZ[group] ?? 1e-6;
  const relative = SPFIT_FIT_ERROR_RELATIVE[group] ?? 0.1;
  const delta = Number.isFinite(numericInitial)
    ? Math.abs(numericValue - numericInitial)
    : 0;
  return Math.max(floor, delta, Math.abs(numericValue) * relative);
}

function syncAutofitJMaxFromInputText() {
  if (!autofitOptions.jmax) return;
  const inputMeta = parseAutofitInputMetadata(autofitInputText);
  if (!Number.isInteger(inputMeta.jMax) || inputMeta.jMax < 0) return;
  autofitOptions.jmax.value = String(inputMeta.jMax);
}

function sanitizeSpfitStem(name) {
  const raw = String(name || 'wms_fitrot_seed')
    .replace(/\.[^.]+$/, '')
    .trim()
    .toLowerCase();
  const cleaned = raw
    .replace(/[^a-z0-9]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || 'wms_fitrot_seed';
}

function createAutofitDownloadMessage(message) {
  const element = document.createElement('div');
  element.className = 'autofit-download-note';
  element.textContent = message;
  return element;
}

function formatSpfitScientific(value, decimals = 12) {
  const safeValue = Number.isFinite(Number(value)) ? Number(value) : 0;
  const text = safeValue.toExponential(decimals).replace('e', 'E');
  return text.replace(/E([+-])(\d+)$/, (_, sign, exp) => `E${sign}${exp.padStart(3, '0')}`);
}

function getSpfitParamSpec(key, reduction) {
  const upperReduction = reduction === 'A' ? 'A' : 'S';
  const specs = {
    A: { id: 10000, label: 'A' },
    B: { id: 20000, label: 'B' },
    C: { id: 30000, label: 'C' },
    DJ: { id: 200, label: upperReduction === 'A' ? 'DELTA J' : 'DJ' },
    DJK: { id: 1100, label: upperReduction === 'A' ? 'DELTA JK' : 'DJK' },
    DK: { id: 2000, label: upperReduction === 'A' ? 'DELTA K' : 'DK' },
    dJ: { id: 40100, label: upperReduction === 'A' ? 'delta J' : 'd1' },
    dK: { id: upperReduction === 'A' ? 41000 : 50000, label: upperReduction === 'A' ? 'delta K' : 'd2' },
    HJ: { id: 300, label: upperReduction === 'A' ? 'PHI N' : 'H N' },
    HJK: { id: 1200, label: upperReduction === 'A' ? 'PHI NK' : 'H NK' },
    HKJ: { id: 2100, label: upperReduction === 'A' ? 'PHI KN' : 'H KN' },
    HK: { id: 3000, label: upperReduction === 'A' ? 'PHI K' : 'H K' },
    h1: { id: upperReduction === 'A' ? 40200 : 50100, label: upperReduction === 'A' ? 'phi N' : 'h 1' },
    h2: { id: upperReduction === 'A' ? 41100 : 51000, label: upperReduction === 'A' ? 'phi NK' : 'h 2' },
    h3: { id: upperReduction === 'A' ? 42000 : 60000, label: upperReduction === 'A' ? 'phi K' : 'h 3' },
  };
  return specs[key] || null;
}

function buildSpfitQuadrupoleParameterEntries(inputMeta) {
  const nuclei = Array.isArray(inputMeta?.hfs?.nuclei) ? inputMeta.hfs.nuclei : [];
  const entries = [];
  nuclei.forEach((nucleus) => {
    if (!nucleus?.active) return;
    const chi = nucleus.chi || {};
    const prefix = 110000000 * Number(nucleus.index || 1);
    const chiAa = Number(chi.aa) || 0;
    const chiBb = Number(chi.bb) || 0;
    const chiCc = Number(chi.cc) || 0;
    const chiAb = Number(chi.ab) || 0;
    const chiAc = Number(chi.ac) || 0;
    const chiBc = Number(chi.bc) || 0;
    const etaTerm = 0.25 * (chiBb - chiCc);
    entries.push({
      key: `chi_aa_${nucleus.index}`,
      id: prefix + 10000,
      label: `N${nucleus.index} 3/2*Chi_aa`,
      value: 1.5 * chiAa,
      error: SPFIT_FIXED_PARAM_ERROR,
      floating: false,
    });
    if (Math.abs(etaTerm) > SPFIT_PARAM_INCLUDE_EPS) {
      entries.push({
        key: `chi_eta_${nucleus.index}`,
        id: prefix + 40000,
        label: `N${nucleus.index} 1/4*(Chi_bb-Chi_cc)`,
        value: etaTerm,
        error: SPFIT_FIXED_PARAM_ERROR,
        floating: false,
      });
    }
    if (Math.abs(chiAb) > SPFIT_PARAM_INCLUDE_EPS) {
      entries.push({
        key: `chi_ab_${nucleus.index}`,
        id: prefix + 610000,
        label: `N${nucleus.index} Chi_ab`,
        value: chiAb,
        error: SPFIT_FIXED_PARAM_ERROR,
        floating: false,
      });
    }
    if (Math.abs(chiAc) > SPFIT_PARAM_INCLUDE_EPS) {
      entries.push({
        key: `chi_ac_${nucleus.index}`,
        id: prefix + 410000,
        label: `N${nucleus.index} Chi_ac`,
        value: chiAc,
        error: SPFIT_FIXED_PARAM_ERROR,
        floating: false,
      });
    }
    if (Math.abs(chiBc) > SPFIT_PARAM_INCLUDE_EPS) {
      entries.push({
        key: `chi_bc_${nucleus.index}`,
        id: prefix + 210000,
        label: `N${nucleus.index} Chi_bc`,
        value: chiBc,
        error: SPFIT_FIXED_PARAM_ERROR,
        floating: false,
      });
    }
  });
  return entries;
}

function buildSpfitParameterEntries(modelParams, {
  reduction,
  fitParamNames = [],
  initialModelParams = null,
  extraEntries = [],
} = {}) {
  const normalized = normalizeModelParams(modelParams);
  if (!normalized) return [];
  const fitParamSet = new Set(Array.isArray(fitParamNames) ? fitParamNames : []);
  const paramEntries = MODEL_PARAM_KEYS
    .filter((key) => {
      const value = Number(normalized[key]);
      if (!Number.isFinite(value)) return false;
      if (key === 'A' || key === 'B' || key === 'C') return true;
      if (fitParamSet.has(key)) return true;
      return Math.abs(value) > SPFIT_PARAM_INCLUDE_EPS;
    })
    .map((key) => {
      const spec = getSpfitParamSpec(key, reduction);
      return spec ? {
        key,
        id: spec.id,
        label: spec.label,
        value: Number(normalized[key]),
        error: fitParamSet.has(key) || key === 'A' || key === 'B' || key === 'C'
          ? computeSpfitFloatingParamError(key, normalized[key], initialModelParams ? initialModelParams[key] : null)
          : SPFIT_FIXED_PARAM_ERROR,
        floating: fitParamSet.has(key) || key === 'A' || key === 'B' || key === 'C',
      } : null;
    })
    .filter(Boolean);
  return [...paramEntries, ...(Array.isArray(extraEntries) ? extraEntries : [])];
}

function parseSpfitAssignmentLabel(label) {
  const normalized = normalizeLabelValue(label);
  if (!normalized || /\balpha\s*=/i.test(normalized)) return null;
  const match = normalized.match(
    /^J\s*=\s*(-?\d+)(?:\s+F\s*=\s*([+\-]?\d+(?:\.\d+)?))?(?:\s+I12\s*=\s*([+\-]?\d+(?:\.\d+)?))?\s+Ka\s*=\s*(-?\d+)\s+Kc\s*=\s*(-?\d+)\s*->\s*J\s*=\s*(-?\d+)(?:\s+F\s*=\s*([+\-]?\d+(?:\.\d+)?))?(?:\s+I12\s*=\s*([+\-]?\d+(?:\.\d+)?))?\s+Ka\s*=\s*(-?\d+)\s+Kc\s*=\s*(-?\d+)(?:\s*\(([^)]+)\))?$/i
  );
  if (!match) return null;

  const lowerJ = Number.parseInt(match[1], 10);
  const lowerF = match[2] == null ? null : Number(match[2]);
  const lowerI12 = match[3] == null ? null : Number(match[3]);
  const lowerKa = Number.parseInt(match[4], 10);
  const lowerKc = Number.parseInt(match[5], 10);
  const upperJ = Number.parseInt(match[6], 10);
  const upperF = match[7] == null ? null : Number(match[7]);
  const upperI12 = match[8] == null ? null : Number(match[8]);
  const upperKa = Number.parseInt(match[9], 10);
  const upperKc = Number.parseInt(match[10], 10);
  const rotationalNumbers = [lowerJ, lowerKa, lowerKc, upperJ, upperKa, upperKc];
  if (!rotationalNumbers.every(Number.isInteger)) return null;

  const hasHyperfine = Number.isFinite(lowerF) || Number.isFinite(upperF);
  if (hasHyperfine && !(Number.isFinite(lowerF) && Number.isFinite(upperF))) return null;
  const hasI12 = Number.isFinite(lowerI12) || Number.isFinite(upperI12);
  if (hasI12 && !(Number.isFinite(lowerI12) && Number.isFinite(upperI12))) return null;

  return {
    rawLabel: normalized,
    lower: { J: lowerJ, F: lowerF, I12: lowerI12, Ka: lowerKa, Kc: lowerKc },
    upper: { J: upperJ, F: upperF, I12: upperI12, Ka: upperKa, Kc: upperKc },
    hasHyperfine,
    hasI12,
    branch: match[11] ? match[11].trim() : null,
  };
}

function buildSpfitLineRecords(assignments) {
  const records = [];
  let maxQuantum = 0;
  let lineMode = null;
  let ignoredI12 = false;
  const qnCollisions = new Map();
  assignments.forEach((assignment) => {
    const parsed = parseSpfitAssignmentLabel(assignment?.label);
    if (!parsed) {
      throw new Error('SPFIT export supports only asymmetric-top labels written as J[/F][/I12], Ka, Kc -> J[/F][/I12], Ka, Kc.');
    }
    const currentMode = parsed.hasHyperfine ? 'hyperfine' : 'rotational';
    if (!lineMode) {
      lineMode = currentMode;
    } else if (lineMode !== currentMode) {
      throw new Error('SPFIT export does not support mixing pure rotational and hyperfine assignments in the same .lin file.');
    }
    if (parsed.hasHyperfine && (!Number.isInteger(parsed.upper.F) || !Number.isInteger(parsed.lower.F))) {
      throw new Error('SPFIT export currently supports only integer F quantum numbers. Half-integer F labels need explicit Pickett quantum-number mapping before export.');
    }

    const quanta = parsed.hasHyperfine
      ? [
          parsed.upper.J, parsed.upper.Ka, parsed.upper.Kc, parsed.upper.F,
          parsed.lower.J, parsed.lower.Ka, parsed.lower.Kc, parsed.lower.F,
        ]
      : [
          parsed.upper.J, parsed.upper.Ka, parsed.upper.Kc,
          parsed.lower.J, parsed.lower.Ka, parsed.lower.Kc,
        ];
    if (parsed.hasI12) {
      ignoredI12 = true;
      const qnKey = quanta.join(',');
      const seenLabel = qnCollisions.get(qnKey);
      if (seenLabel && seenLabel !== parsed.rawLabel) {
        throw new Error('SPFIT export found multiple assignments that differ only by I12. This exporter would collapse them onto the same SPFIT quantum numbers.');
      }
      qnCollisions.set(qnKey, parsed.rawLabel);
    }
    maxQuantum = Math.max(maxQuantum, ...quanta.map((value) => Math.abs(Number(value) || 0)));
    records.push({
      ...assignment,
      ...parsed,
      quanta,
    });
  });
  return {
    records,
    maxQuantum,
    lineMode: lineMode || 'rotational',
    ignoredI12,
  };
}

function buildSpfitLinText(records, { lineErrorMHz = SPFIT_DEFAULT_LINE_ERROR_MHZ, weight = SPFIT_DEFAULT_LINE_WEIGHT } = {}) {
  return records.map((record) => {
    const qnFields = Array.from({ length: 12 }, (_, index) => {
      const value = index < record.quanta.length ? record.quanta[index] : null;
      return value == null ? '   ' : String(value).padStart(3, ' ');
    }).join('');
    const freqText = record.fObs.toFixed(6).padStart(15, ' ');
    const errText = Number(lineErrorMHz).toFixed(4).padStart(10, ' ');
    const weightText = formatSpfitScientific(weight, 3).padStart(10, ' ');
    const comment = record.hasHyperfine
      ? [
          `Ju=${record.upper.J}`,
          `Ka_u=${record.upper.Ka}`,
          `Kc_u=${record.upper.Kc}`,
          `F_u=${record.upper.F}`,
          Number.isFinite(record.upper.I12) ? `I12_u=${record.upper.I12}` : null,
          '<-',
          `Jl=${record.lower.J}`,
          `Ka_l=${record.lower.Ka}`,
          `Kc_l=${record.lower.Kc}`,
          `F_l=${record.lower.F}`,
          Number.isFinite(record.lower.I12) ? `I12_l=${record.lower.I12}` : null,
        ].filter(Boolean).join(' ')
      : `Ju=${record.upper.J} Ka_u=${record.upper.Ka} Kc_u=${record.upper.Kc} <- Jl=${record.lower.J} Ka_l=${record.lower.Ka} Kc_l=${record.lower.Kc}`;
    return `${qnFields}${freqText} ${errText} ${weightText} / ${comment} /`;
  }).join('\n') + '\n';
}

function formatSpfitTitleTimestamp(date = new Date()) {
  const weekdays = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'];
  const months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const hh = String(date.getHours()).padStart(2, '0');
  const mm = String(date.getMinutes()).padStart(2, '0');
  const ss = String(date.getSeconds()).padStart(2, '0');
  return `${weekdays[date.getDay()]} ${months[date.getMonth()]} ${String(date.getDate()).padStart(2, ' ')} ${hh}:${mm}:${ss} ${date.getFullYear()}`;
}

function buildSpfitOptionLine({ reduction, representation, nuclei, jMax }) {
  const safeReduction = reduction === 'A' ? 'a' : 's';
  const spinTag = (Array.isArray(nuclei) && nuclei.length)
    ? nuclei.map((nucleus) => String(Math.round((2 * Number(nucleus.spin)) + 1))).join('')
    : '1';
  const sign = getSpfitRepresentationHeaderSign(representation);
  const axisCode = getSpfitRepresentationAxisCode(representation);
  const safeJMax = Math.max(1, Number(jMax) || 1);
  return `${safeReduction}   ${spinTag}  ${sign}  0  ${String(safeJMax).padStart(2, ' ')}  0  ${axisCode}  1  1  0  0  0`;
}

function buildSpfitParText({ title, reduction, representation, paramEntries, lineCount, jMax, nuclei }) {
  const requestedLines = Math.max(1, Number(lineCount) || 1);
  const header = `${String(paramEntries.length).padStart(4, ' ')} ${String(requestedLines).padStart(4, ' ')} ${String(SPFIT_HEADER_LINE_WINDOW).padStart(4, ' ')}  0  0.0000E+000  1.0000E+006  -1.0000E+000  1.0000000000`;
  const optionLine = buildSpfitOptionLine({
    reduction,
    representation,
    nuclei,
    jMax: Math.max(1, Number(jMax) || 1),
  });
  const paramLines = paramEntries.map((entry) => {
    const idText = String(entry.id).padStart(13, ' ');
    const valueText = formatSpfitScientific(entry.value, 12).padStart(20, ' ');
    const errorText = formatSpfitScientific(entry.error, 8).padStart(15, ' ');
    return `${idText} ${valueText} ${errorText} / ${entry.label} /`;
  });
  const titleLine = `${String(title || 'wms_fitrot_seed').padEnd(55, ' ')}${formatSpfitTitleTimestamp()}`;
  return [titleLine, header, optionLine, ...paramLines, ''].join('\n');
}

function buildSpfitReadmeText({
  stem,
  reduction,
  representation,
  assignments,
  paramEntries,
  lineErrorMHz,
  weight,
  diagnostics = null,
  inputMeta = null,
  lineMode = 'rotational',
  warnings = [],
}) {
  const activeNuclei = Array.isArray(inputMeta?.hfs?.nuclei)
    ? inputMeta.hfs.nuclei.filter((nucleus) => nucleus.active)
    : [];
  const lines = [
    'WMS-FitRot -> SPFIT seed export',
    '',
    `Base name: ${stem}`,
    `Reduction: ${reduction}`,
    `Representation: ${representation || SPFIT_DEFAULT_REPRESENTATION}`,
    `Export mode: ${lineMode === 'hyperfine' ? 'hyperfine (J, Ka, Kc, F)' : 'pure rotational (J, Ka, Kc)'}`,
    `Lines exported: ${assignments.length}`,
    `Parameters exported: ${paramEntries.length}`,
    activeNuclei.length ? `Quadrupole nuclei exported: ${activeNuclei.map((nucleus) => `N${nucleus.index} (I=${String(nucleus.spin)})`).join(', ')}` : null,
    diagnostics?.readyLabel ? `Diagnostic state: ${diagnostics.readyLabel}` : null,
    '',
    'Files:',
    `- ${stem}.par : initial SPFIT parameter seed`,
    `- ${stem}.lin : assigned observed transitions`,
    '- session.json : provenance and export metadata',
    '- README.txt : export notes',
    '',
    'Current export rules:',
    '- The .par option card is derived from WMS reduction, representation, active quadrupole spins, and J_MAX.',
    `- All .lin uncertainties are set to ${Number(lineErrorMHz).toFixed(4)} MHz.`,
    `- All .lin weights are set to ${Number(weight).toFixed(3)}.`,
    '- Quantum numbers in .lin are written as upper state first, then lower state.',
    '- Parameters actively fitted by WMS-FitRot are exported floating; context-only parameters are exported fixed.',
    '- Hyperfine quadrupole constants are written as Pickett parameters only when F-resolved assignments are exported.',
    '- Review the .par file before publication-grade refinement.',
    diagnostics ? '- Only lines approved by the M2 diagnostic gate are exported automatically.' : null,
    '',
    'Suggested next step:',
    `- Run: spfit ${stem}`,
  ].filter(Boolean);
  if (warnings.length) {
    lines.push('', 'Warnings:');
    warnings.forEach((warning) => {
      lines.push(`- ${warning}`);
    });
  }
  return `${lines.join('\n')}\n`;
}

function buildSpfitSessionPayload({
  stem,
  reduction,
  reductionSource,
  representation,
  inputMeta,
  theorySpectrum,
  expSpectrum,
  assignments,
  records,
  lineMode,
  paramEntries,
  modelParams,
  quarticScaleState,
  fitParamNames,
  diagnostics = null,
  warnings = [],
}) {
  return {
    exportType: 'spfit-seed',
    exportVersion: 2,
    exportedAt: new Date().toISOString(),
    stem,
    reduction,
    reductionSource,
    representation,
    inputMeta: cloneJsonValue(inputMeta),
    lineMode,
    theorySpectrum: theorySpectrum ? {
      id: theorySpectrum.id,
      name: theorySpectrum.name,
    } : null,
    experimentalSpectrum: expSpectrum ? {
      id: expSpectrum.id,
      name: expSpectrum.name,
      mode: expSpectrum.mode,
    } : null,
    finalModelParams: normalizeModelParams(modelParams),
    finalQuarticScaleState: quarticScaleState || null,
    fitParamNames,
    exportedParameters: paramEntries.map((entry) => ({
      key: entry.key,
      id: entry.id,
      label: entry.label,
      value: entry.value,
      error: entry.error,
      floating: Boolean(entry.floating),
    })),
    assignments,
    spfitAssignments: records.map((record) => ({
      label: record.label,
      fObsMHz: record.fObs,
      residualMHz: record.residual,
      lower: record.lower,
      upper: record.upper,
      quanta: Array.isArray(record.quanta) ? record.quanta.slice() : [],
      hasHyperfine: Boolean(record.hasHyperfine),
      hasI12: Boolean(record.hasI12),
      branch: record.branch,
    })),
    diagnostics: cloneJsonValue(diagnostics),
    context: cloneJsonValue(lastAutofitContext),
    summary: cloneJsonValue(lastAutofitSummary),
    fitMetrics: cloneJsonValue(lastAutofitFitMetrics),
    warnings,
  };
}

async function createSpfitSeedExport({
  outputText,
  modelParams,
  quarticScaleState,
  theorySpectrum,
  expSpectrum,
  diagnostics = null,
}) {
  const assignments = diagnostics
    ? diagnostics.exportableAssignments
    : parseAutofitFinalAssignments(outputText);
  if (!assignments.length) {
    throw new Error('No approved assignments are currently available for SPFIT export.');
  }
  const inputMeta = parseAutofitInputMetadata(autofitInputText);
  const reduction = inputMeta.reduction || 'S';
  const reductionSource = inputMeta.reduction ? inputMeta.reductionSource : 'fallback';
  const representation = inputMeta.representation || SPFIT_DEFAULT_REPRESENTATION;
  const fitParamNames = parseAutofitFitParamNames(outputText);
  const initialModelParams = getAutofitInputModelParams();
  const { records, maxQuantum, lineMode, ignoredI12 } = buildSpfitLineRecords(assignments);
  const activeQuadrupoleNuclei = lineMode === 'hyperfine' && Array.isArray(inputMeta?.hfs?.nuclei)
    ? inputMeta.hfs.nuclei.filter((nucleus) => nucleus.active)
    : [];
  const quadrupoleEntries = lineMode === 'hyperfine'
    ? buildSpfitQuadrupoleParameterEntries(inputMeta)
    : [];
  const paramEntries = buildSpfitParameterEntries(modelParams, {
    reduction,
    fitParamNames,
    initialModelParams,
    extraEntries: quadrupoleEntries,
  });
  if (!paramEntries.length) {
    throw new Error('Unable to derive SPFIT parameters from the final model.');
  }
  if (ignoredI12 && activeQuadrupoleNuclei.length > 1) {
    throw new Error('SPFIT export with I12-resolved labels and more than one active quadrupole nucleus is not yet implemented reliably. The current exporter stops instead of writing ambiguous quantum numbers.');
  }
  const stem = sanitizeSpfitStem(theorySpectrum?.name || expSpectrum?.name || 'wms_fitrot_seed');
  const warnings = [];
  if (reductionSource !== 'explicit') {
    warnings.push(`Watson reduction was ${reductionSource === 'fallback' ? 'not found' : 'inferred'} from the WMS input; verify the ${stem}.par option card.`);
  }
  if (!inputMeta.representation) {
    warnings.push(`Representation was not found in the WMS input; the ${stem}.par option card defaults to ${SPFIT_DEFAULT_REPRESENTATION}.`);
  }
  if (lineMode === 'hyperfine' && !activeQuadrupoleNuclei.length) {
    warnings.push(`F-resolved assignments were exported, but no active quadrupole tensor was found in input.txt; verify the hyperfine block in ${stem}.par.`);
  }
  if (ignoredI12) {
    warnings.push('I12 labels were present in the WMS assignments and were kept only in comments; SPFIT quantum numbers were exported as J, Ka, Kc, F.');
  }
  if (diagnostics?.readyState && diagnostics.readyState !== 'ready') {
    warnings.push(`Diagnostic state is ${diagnostics.readyLabel}; inspect the session before using this seed for final refinement.`);
  }
  const zip = new JSZip();
  zip.file(`${stem}.par`, buildSpfitParText({
    title: stem,
    reduction,
    representation,
    paramEntries,
    lineCount: records.length,
    jMax: Math.max(
      maxQuantum,
      Number.isInteger(inputMeta.jMax) ? inputMeta.jMax : 0
    ),
    nuclei: activeQuadrupoleNuclei,
  }));
  zip.file(`${stem}.lin`, buildSpfitLinText(records, {
    lineErrorMHz: SPFIT_DEFAULT_LINE_ERROR_MHZ,
    weight: SPFIT_DEFAULT_LINE_WEIGHT,
  }));
  zip.file('README.txt', buildSpfitReadmeText({
    stem,
    reduction,
    representation,
    assignments,
    paramEntries,
    lineErrorMHz: SPFIT_DEFAULT_LINE_ERROR_MHZ,
    weight: SPFIT_DEFAULT_LINE_WEIGHT,
    diagnostics,
    inputMeta,
    lineMode,
    warnings,
  }));
  zip.file('session.json', `${JSON.stringify(buildSpfitSessionPayload({
    stem,
    reduction,
    reductionSource,
    representation,
    inputMeta,
    theorySpectrum,
    expSpectrum,
    assignments,
    records,
    lineMode,
    paramEntries,
    modelParams,
    quarticScaleState,
    fitParamNames,
    diagnostics,
    warnings,
  }), null, 2)}\n`);
  const blob = await zip.generateAsync({ type: 'blob' });
  return {
    stem,
    blob,
    warnings,
    lineCount: records.length,
    paramCount: paramEntries.length,
  };
}

function cloneJsonValue(value) {
  if (value == null) return value;
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (error) {
    console.warn('Unable to clone value for session export:', error);
    return null;
  }
}

function collectAutofitOptionInputState() {
  const state = {};
  Object.entries(autofitOptions).forEach(([key, input]) => {
    if (!input) return;
    state[key] = input.type === 'checkbox'
      ? Boolean(input.checked)
      : String(input.value ?? '');
  });
  return state;
}

function applyAutofitOptionInputState(state) {
  if (!state || typeof state !== 'object') return;
  Object.entries(autofitOptions).forEach(([key, input]) => {
    if (!input || !(key in state)) return;
    if (input.type === 'checkbox') {
      input.checked = Boolean(state[key]);
      return;
    }
    input.value = String(state[key] ?? '');
  });
}

function collectAutofitParamBoundInputState() {
  const state = {};
  BOUND_PARAM_KEYS.forEach((key) => {
    const inputs = autofitParamBoundInputs.get(key);
    if (!inputs) return;
    state[key] = {
      min: String(inputs.minInput?.value ?? ''),
      max: String(inputs.maxInput?.value ?? ''),
    };
  });
  return state;
}

function applyAutofitParamBoundInputState(state) {
  BOUND_PARAM_KEYS.forEach((key) => {
    const inputs = autofitParamBoundInputs.get(key);
    if (!inputs) return;
    const current = state && typeof state === 'object' ? state[key] : null;
    inputs.minInput.value = current?.min != null ? String(current.min) : '';
    inputs.maxInput.value = current?.max != null ? String(current.max) : '';
  });
}

function collectOptionalFitParamState() {
  const state = {
    [FIT_SCALE_PARAM_KEY]: Boolean(autofitFitQuarticScale?.checked),
  };
  DISTORTION_PARAM_KEYS.forEach((key) => {
    state[key] = Boolean(autofitFitParamInputs[key]?.checked);
  });
  return state;
}

function applyOptionalFitParamState(state) {
  if (autofitFitQuarticScale) {
    autofitFitQuarticScale.checked = Boolean(state?.[FIT_SCALE_PARAM_KEY]);
  }
  DISTORTION_PARAM_KEYS.forEach((key) => {
    if (autofitFitParamInputs[key]) {
      autofitFitParamInputs[key].checked = Boolean(state?.[key]);
    }
  });
  syncQuarticFitSelectionUi();
}

function updateRestoreInitialButtonState() {
  if (!autofitRestoreInitialButton) return;
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const initialState = normalizeSerializedSpectrumState(theorySpectrum?.initialState);
  autofitRestoreInitialButton.disabled = autofitRunInProgress || !theorySpectrum?.isTheoretical || !initialState?.text;
}

function syncInitialTheoryStateFromInputText() {
  const inputModelParams = getAutofitInputModelParams();
  spectra.forEach((spectrum) => {
    if (!spectrum?.isTheoretical) return;
    const initialState = normalizeSerializedSpectrumState(spectrum.initialState)
      || buildSpectrumInitialStateSnapshot(spectrum);
    if (!initialState) return;
    if (inputModelParams && !initialState.modelParams) {
      initialState.modelParams = cloneJsonValue(inputModelParams);
    }
    if (inputModelParams && !initialState.quarticReferenceParams) {
      initialState.quarticReferenceParams = normalizeQuarticReferenceParams(inputModelParams);
    }
    spectrum.initialState = initialState;
  });
  updateRestoreInitialButtonState();
}

function applySpectrumSerializedState(spectrum, serializedState) {
  if (!spectrum) return false;
  const state = normalizeSerializedSpectrumState(serializedState);
  if (!state?.text) return false;
  const parsed = parseDatSpectrum(state.text);
  if (!parsed.freqs.length) return false;

  const baseAbs = parsed.intensities.map((value) => Math.abs(Number(value)));
  spectrum.mode = state.mode;
  spectrum.freqs = cloneSpectrumFrequencies(parsed.freqs);
  spectrum.baseAbs = cloneSpectrumFrequencies(baseAbs);
  spectrum.labels = cloneSpectrumLabels(parsed.labels);
  spectrum.rawText = state.text;
  spectrum.previewView = state.mode === 'profile'
    ? (state.previewView === 'stick' ? 'stick' : 'profile')
    : 'stick';
  spectrum.modelParams = state.modelParams ? cloneJsonValue(state.modelParams) : null;
  spectrum.fittedConstants = state.fittedConstants ? cloneJsonValue(state.fittedConstants) : null;
  spectrum.quarticScale = normalizeQuarticScale(state.quarticScale, 1);
  spectrum.quarticReferenceParams = normalizeQuarticReferenceParams(state.quarticReferenceParams);

  if (spectrum.mode === 'profile') {
    const stickData = computeDatStickSpectrum(parsed.freqs, parsed.intensities, null, 'absorbance');
    spectrum.stick = {
      freqs: stickData.freqs,
      baseInts: stickData.baseInts,
      labels: [],
    };
  } else {
    spectrum.stick = {
      freqs: cloneSpectrumFrequencies(parsed.freqs),
      baseInts: cloneSpectrumFrequencies(baseAbs),
      labels: cloneSpectrumLabels(parsed.labels),
    };
  }
  updateHighlightsForSpectrum(spectrum);
  return true;
}

function buildFitStatusSpectrumEntry(spectrum) {
  if (!spectrum) return null;
  const settings = ensureOverlaySettings(spectrum);
  const initialState = normalizeSerializedSpectrumState(spectrum.initialState)
    || (spectrum.isTheoretical ? buildSpectrumInitialStateSnapshot(spectrum, {
      modelParamsOverride: getAutofitInputModelParams() || spectrum.modelParams,
      quarticReferenceParamsOverride: normalizeQuarticReferenceParams(getAutofitInputModelParams() || spectrum.quarticReferenceParams),
    }) : null);
  return {
    id: spectrum.id,
    name: spectrum.name,
    role: spectrum.isTheoretical ? 'theory' : 'experimental',
    mode: spectrum.mode === 'profile' ? 'profile' : 'stick',
    previewView: spectrum.previewView === 'profile' ? 'profile' : 'stick',
    text: serializeSpectrumText(spectrum),
    modelParams: normalizeModelParams(spectrum.modelParams),
    fittedConstants: spectrum.fittedConstants ? cloneJsonValue(spectrum.fittedConstants) : null,
    quarticScale: getSpectrumQuarticScale(spectrum),
    quarticReferenceParams: normalizeQuarticReferenceParams(spectrum.quarticReferenceParams),
    overlay: cloneJsonValue({
      enabled: settings.enabled,
      view: settings.view,
      intensityMode: settings.intensityMode,
      normalize: settings.normalize,
      scale: settings.scale,
      color: settings.color,
    }),
    initialState: initialState ? cloneJsonValue(initialState) : null,
  };
}

function buildFitStatusLogPayload() {
  return {
    type: FIT_STATUS_LOG_TYPE,
    version: FIT_STATUS_LOG_VERSION,
    savedAt: new Date().toISOString(),
    input: {
      name: autofitInputName || null,
      text: typeof autofitInputText === 'string' ? autofitInputText : '',
    },
    spectra: Array.from(spectra.values())
      .map((spectrum) => buildFitStatusSpectrumEntry(spectrum))
      .filter(Boolean),
    ui: {
      activeSpectrumId,
      theorySpectrumId: autofitTheorySelect?.value || null,
      experimentalSpectrumId: autofitExpSelect?.value || null,
      expMode: autofitExpMode?.value || 'profile',
      useInputAssignment: Boolean(autofitUseInputAssignment?.checked),
      assignmentMode: autofitAssignmentMode?.value || 'auto',
      clickMode: autofitClickMode?.value || 'fit',
      optionInputs: collectAutofitOptionInputState(),
      fitParamState: collectOptionalFitParamState(),
      customParamBounds: collectAutofitParamBoundInputState(),
      selectedFitLabels: selectedFitLabels.slice(),
      selectedControlLabels: selectedControlLabels.slice(),
      manualAssignments: cloneJsonValue(manualAssignments),
      pendingManualSelection: cloneJsonValue(pendingManualSelection),
      statusText: autofitStatus?.textContent || '',
      outputText: autofitOutput?.textContent || '',
      warningVisible: autofitWarning ? !autofitWarning.hidden : false,
    },
    lastAutofitSummary: cloneJsonValue(lastAutofitSummary),
    lastAutofitContext: cloneJsonValue(lastAutofitContext),
    fitMetrics: cloneJsonValue(lastAutofitFitMetrics),
    diagnostics: cloneJsonValue(lastAutofitDiagnostics),
    spfitImportSummary: cloneJsonValue(lastSpfitImportSummary),
  };
}

function clearCurrentFitStatusState() {
  if (autofitWorker) {
    autofitWorker.terminate();
    autofitWorker = null;
  }
  autofitRunInProgress = false;
  spectra.clear();
  overlaySettings.clear();
  overlayHighlightsFit.clear();
  overlayHighlightsControl.clear();
  selectedFitLabels.length = 0;
  selectedControlLabels.length = 0;
  manualAssignments.length = 0;
  pendingManualSelection = null;
  lastAutofitSummary = null;
  lastAutofitContext = null;
  lastAutofitFitMetrics = null;
  activeSpectrumId = null;
  colorIndex = 0;
  spectrumCounter = 0;
  autofitInputText = null;
  autofitInputName = null;
  if (autofitUseInputAssignment) {
    autofitUseInputAssignment.checked = false;
  }
  if (autofitOutput) {
    autofitOutput.textContent = '';
  }
  if (autofitDownloads) {
    autofitDownloads.innerHTML = '';
  }
  if (autofitRunButton) autofitRunButton.disabled = false;
  if (autofitCancelButton) autofitCancelButton.disabled = true;
  updateRestoreInitialButtonState();
  setAutofitWarningVisible(false);
  resetAutofitProgress();
  renderLabelLists();
  renderUploadedTabs();
  renderOverlayList();
  updateOverlayPlot();
  updateRestoreInitialButtonState();
}

function restoreFitStatusSelections(uiState = {}) {
  selectedFitLabels.length = 0;
  selectedControlLabels.length = 0;
  manualAssignments.length = 0;
  pendingManualSelection = null;
  overlayHighlightsFit.clear();
  overlayHighlightsControl.clear();

  if (Array.isArray(uiState.selectedControlLabels)) {
    uiState.selectedControlLabels.forEach((label) => {
      const normalized = normalizeLabelValue(label);
      if (!normalized || selectedControlLabels.includes(normalized)) return;
      selectedControlLabels.push(normalized);
    });
  }

  if (Array.isArray(uiState.manualAssignments)) {
    uiState.manualAssignments.forEach((item) => {
      const label = normalizeLabelValue(item?.label);
      const theorySpectrumId = spectra.has(item?.theorySpectrumId) ? item.theorySpectrumId : autofitTheorySelect?.value;
      const expSpectrumId = spectra.has(item?.expSpectrumId) ? item.expSpectrumId : autofitExpSelect?.value;
      const expX = Number(item?.expX);
      const theoryX = Number(item?.theoryX);
      if (!label || !theorySpectrumId || !expSpectrumId || !Number.isFinite(expX)) return;
      manualAssignments.push({
        label,
        theorySpectrumId,
        expSpectrumId,
        theoryX: Number.isFinite(theoryX) ? theoryX : null,
        expX,
      });
    });
  }

  if (!isManualAssignmentMode() && Array.isArray(uiState.selectedFitLabels)) {
    uiState.selectedFitLabels.forEach((label) => {
      const normalized = normalizeLabelValue(label);
      if (!normalized || selectedFitLabels.includes(normalized)) return;
      selectedFitLabels.push(normalized);
    });
  } else {
    syncFitLabelsFromManualAssignments();
  }

  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const labelToFreq = buildTheoryLabelFrequencyMap(theorySpectrum);
  selectedFitLabels.forEach((label) => {
    const freq = labelToFreq.get(label);
    if (Number.isFinite(freq) && theorySpectrum) {
      overlayHighlightsFit.set(label, { label, x: freq, spectrumId: theorySpectrum.id });
    }
  });
  selectedControlLabels.forEach((label) => {
    const freq = labelToFreq.get(label);
    if (Number.isFinite(freq) && theorySpectrum) {
      overlayHighlightsControl.set(label, { label, x: freq, spectrumId: theorySpectrum.id });
    }
  });
  manualAssignments.forEach((item) => {
    const freq = labelToFreq.get(item.label);
    if (Number.isFinite(freq)) {
      item.theoryX = freq;
    }
  });

  const pending = uiState.pendingManualSelection;
  if (pending?.label && theorySpectrum) {
    const normalizedLabel = normalizeLabelValue(pending.label);
    const theoryX = labelToFreq.get(normalizedLabel);
    if (normalizedLabel) {
      pendingManualSelection = {
        label: normalizedLabel,
        theorySpectrumId: theorySpectrum.id,
        theoryX: Number.isFinite(theoryX) ? theoryX : Number(pending?.theoryX),
      };
    }
  }
  renderLabelLists();
  updateOverlayPlot();
}

function restoreFitStatusLogPayload(payload) {
  if (!payload || typeof payload !== 'object') {
    throw new Error('Fit status log is empty or invalid.');
  }
  if (payload.type !== FIT_STATUS_LOG_TYPE) {
    throw new Error('Selected file is not a WMS-FitRot Fit Status LOG.');
  }
  if (autofitRunInProgress) {
    throw new Error('Cancel the current autofit before loading a Fit Status LOG.');
  }

  clearCurrentFitStatusState();

  autofitInputText = typeof payload?.input?.text === 'string' ? payload.input.text : null;
  autofitInputName = typeof payload?.input?.name === 'string' ? payload.input.name : null;

  const spectrumEntries = Array.isArray(payload.spectra) ? payload.spectra : [];
  spectrumEntries.forEach((entry, index) => {
    addSpectrumFromText({
      id: entry?.id || null,
      name: entry?.name || `Spectrum ${index + 1}`,
      text: typeof entry?.text === 'string' ? entry.text : '',
      mode: entry?.mode === 'profile' ? 'profile' : 'stick',
      role: entry?.role === 'theory' ? 'theory' : null,
      previewView: entry?.previewView,
      overlayState: entry?.overlay,
      modelParams: entry?.modelParams,
      fittedConstants: entry?.fittedConstants,
      quarticScale: entry?.quarticScale,
      quarticReferenceParams: entry?.quarticReferenceParams,
      initialState: entry?.initialState,
      setActive: index === 0,
      skipRender: true,
    });
  });

  renderUploadedTabs();
  renderOverlayList();
  updateOverlayPlot();
  syncAutofitJMaxFromInputText();
  syncInitialTheoryStateFromInputText();

  const uiState = payload.ui && typeof payload.ui === 'object' ? payload.ui : {};
  if (autofitTheorySelect && spectra.has(uiState.theorySpectrumId)) {
    autofitTheorySelect.value = uiState.theorySpectrumId;
  }
  if (autofitExpSelect && spectra.has(uiState.experimentalSpectrumId)) {
    autofitExpSelect.value = uiState.experimentalSpectrumId;
  }
  if (autofitExpMode) {
    autofitExpMode.value = uiState.expMode === 'stick' ? 'stick' : 'profile';
  }
  if (autofitUseInputAssignment) {
    autofitUseInputAssignment.checked = Boolean(uiState.useInputAssignment);
  }
  if (autofitAssignmentMode) {
    autofitAssignmentMode.value = uiState.assignmentMode === 'manual' ? 'manual' : 'auto';
  }
  if (autofitClickMode) {
    autofitClickMode.value = uiState.clickMode === 'control' ? 'control' : 'fit';
  }

  applyAutofitOptionInputState(uiState.optionInputs);
  applyOptionalFitParamState(uiState.fitParamState);
  applyAutofitParamBoundInputState(uiState.customParamBounds);
  syncLeastSquaresUiState();
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();
  restoreFitStatusSelections(uiState);

  activeSpectrumId = spectra.has(uiState.activeSpectrumId) ? uiState.activeSpectrumId : activeSpectrumId;
  updateActiveTabDisplay();
  lastAutofitSummary = cloneJsonValue(payload.lastAutofitSummary);
  lastAutofitContext = cloneJsonValue(payload.lastAutofitContext);
  lastAutofitFitMetrics = cloneJsonValue(payload.fitMetrics);
  renderSensitivityDashboard(lastAutofitFitMetrics);
  renderAutofitDiagnostics(payload.diagnostics || null);
  renderSpfitImportReport(payload.spfitImportSummary || null);
  setAutofitWarningVisible(Boolean(uiState.warningVisible));
  setAutofitStatus(uiState.statusText || 'Fit status log restored.');
  if (autofitOutput) {
    autofitOutput.textContent = uiState.outputText || '';
  }
  updateRestoreInitialButtonState();
  showOverlayToast('Fit status log restored');
}

function downloadFitStatusLog() {
  if (autofitRunInProgress) {
    alert('Wait for the current autofit to finish, or cancel it, before downloading the Fit Status LOG.');
    return;
  }
  const payload = buildFitStatusLogPayload();
  const theorySpectrum = getSpectrumById(payload?.ui?.theorySpectrumId);
  const baseName = sanitizePlotFilename(
    theorySpectrum?.name || autofitInputName || 'wms_fitrot_fit_status',
    'wms_fitrot_fit_status'
  );
  const blob = new Blob([`${JSON.stringify(payload, null, 2)}\n`], { type: 'application/json' });
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = `${baseName}_fit_status_log.json`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(link.href), 1000);
}

async function handleFitStatusLogUpload(file) {
  if (!file) return;
  const text = await file.text();
  let payload;
  try {
    payload = JSON.parse(text);
  } catch (error) {
    throw new Error('Fit Status LOG is not valid JSON.');
  }
  restoreFitStatusLogPayload(payload);
}

function restoreInitialTheoryParameters() {
  if (autofitRunInProgress) {
    alert('Cancel the current autofit before restoring the initial parameters.');
    return;
  }
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const initialState = normalizeSerializedSpectrumState(theorySpectrum?.initialState);
  if (!theorySpectrum?.isTheoretical || !initialState?.text) {
    alert('No restorable initial theoretical state is available for the selected spectrum.');
    return;
  }
  const confirmed = window.confirm(
    'This will restore the original theoretical parameters and spectrum as initially received from WMS-Rot, replacing the current fitted state. Continue?'
  );
  if (!confirmed) return;

  if (!applySpectrumSerializedState(theorySpectrum, initialState)) {
    alert('Unable to restore the initial theoretical spectrum.');
    return;
  }
  theorySpectrum.initialState = cloneJsonValue(initialState);
  resetAutofitProgress();
  lastAutofitSummary = null;
  lastAutofitContext = null;
  if (autofitOutput) {
    autofitOutput.textContent = '';
  }
  setAutofitWarningVisible(false);
  setAutofitStatus(`Initial theoretical parameters restored for ${theorySpectrum.name}.`);
  renderUploadedTabs();
  renderOverlayList();
  updateOverlayPlot();
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();
  updateRestoreInitialButtonState();
  showOverlayToast('Initial theory parameters restored');
}

function parseSpfitNumber(value) {
  if (value == null) return Number.NaN;
  const normalized = String(value).trim().replace(/D/gi, 'E');
  const parsed = Number(normalized);
  return Number.isFinite(parsed) ? parsed : Number.NaN;
}

function inferSpfitReductionFromText(text, fallback = null) {
  if (typeof text !== 'string' || !text.trim()) return fallback;
  const lines = String(text).split(/\r?\n/);
  if (lines.length >= 3) {
    const optionLine = String(lines[2] || '').trim();
    const reductionChar = optionLine.charAt(0).toUpperCase();
    if (reductionChar === 'A' || reductionChar === 'S') {
      return reductionChar;
    }
  }
  if (/\b(DELTA J|DELTA N|DELTA JK|DELTA NK|PHI N|PHI NK|PHI KN|PHI K|delta J|delta N|delta K|phi N|phi NK|phi K)\b/.test(text)) {
    return 'A';
  }
  if (/\b(DJ|DJK|DK|d1|d2|H N|H NK|H KN|H K|h 1|h 2|h 3)\b/.test(text)) {
    return 'S';
  }
  return fallback;
}

function canonicalizeSpfitParamLabel(label) {
  if (label == null) return '';
  return String(label)
    .replace(/^\s*[-+]\s*/, '')
    .replace(/\(.*?\)/g, ' ')
    .replace(/\[.*?\]/g, ' ')
    .replace(/[.,]/g, ' ')
    .replace(/\s+/g, ' ')
    .trim();
}

function getSpfitParamKeyById(id) {
  const numericId = Math.abs(Number.parseInt(id, 10));
  if (!Number.isFinite(numericId)) return null;
  const idMap = {
    10000: 'A',
    20000: 'B',
    30000: 'C',
    200: 'DJ',
    1100: 'DJK',
    2000: 'DK',
    40100: 'dJ',
    41000: 'dK',
    50000: 'dK',
    300: 'HJ',
    1200: 'HJK',
    2100: 'HKJ',
    3000: 'HK',
    40200: 'h1',
    41100: 'h2',
    42000: 'h3',
    50100: 'h1',
    51000: 'h2',
    60000: 'h3',
  };
  return idMap[numericId] || null;
}

function getSpfitParamKeyByLabel(label) {
  const normalized = canonicalizeSpfitParamLabel(label);
  if (!normalized) return null;
  const upper = normalized.toUpperCase();
  if (upper === 'A') return 'A';
  if (upper === 'B') return 'B';
  if (upper === 'C') return 'C';
  if (upper === 'A-B') return 'A_MINUS_B';
  if (upper === 'DELTA J' || upper === 'DELTA N' || upper === 'DJ') return 'DJ';
  if (upper === 'DELTA JK' || upper === 'DELTA NK' || upper === 'DJK') return 'DJK';
  if (upper === 'DELTA K' || upper === 'DK') return 'DK';
  if (normalized === 'delta J' || normalized === 'delta N' || upper === 'D1' || upper === 'DJ / D1') return 'dJ';
  if (normalized === 'delta K' || upper === 'D2' || upper === 'DK / D2') return 'dK';
  if (upper === 'PHI N' || upper === 'H N' || upper === 'HJ') return 'HJ';
  if (upper === 'PHI NK' || upper === 'H NK' || upper === 'HJK') return 'HJK';
  if (upper === 'PHI KN' || upper === 'H KN' || upper === 'HKJ') return 'HKJ';
  if (upper === 'PHI K' || upper === 'H K' || upper === 'HK') return 'HK';
  if (normalized === 'phi N' || upper === 'H1' || upper === 'H 1') return 'h1';
  if (normalized === 'phi NK' || upper === 'H2' || upper === 'H 2') return 'h2';
  if (normalized === 'phi K' || upper === 'H3' || upper === 'H 3') return 'h3';
  return null;
}

function parseSpfitParameterLine(rawLine) {
  const line = String(rawLine || '').trim();
  if (!line || /^[-=*]{3,}$/.test(line)) return null;

  const slashMatch = line.match(/^\s*(?:(\d+)\s+(\d+)\s+)?([+-]?\d+)\s+([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)\s+([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)\s*\/\s*([^/]+?)\s*\/?/i);
  if (slashMatch) {
    return {
      index: slashMatch[1] ? Number.parseInt(slashMatch[1], 10) : null,
      group: slashMatch[2] ? Number.parseInt(slashMatch[2], 10) : null,
      id: Number.parseInt(slashMatch[3], 10),
      value: parseSpfitNumber(slashMatch[4]),
      error: parseSpfitNumber(slashMatch[5]),
      label: canonicalizeSpfitParamLabel(slashMatch[6]),
      rawLine: String(rawLine || ''),
    };
  }

  const plainMatch = line.match(/^\s*(\d+)\s+(\d+)\s+([+-]?\d+)\s+([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)\s+([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)\s+(.+)$/i);
  if (plainMatch) {
    return {
      index: Number.parseInt(plainMatch[1], 10),
      group: Number.parseInt(plainMatch[2], 10),
      id: Number.parseInt(plainMatch[3], 10),
      value: parseSpfitNumber(plainMatch[4]),
      error: parseSpfitNumber(plainMatch[5]),
      label: canonicalizeSpfitParamLabel(plainMatch[6]),
      rawLine: String(rawLine || ''),
    };
  }

  return null;
}

function parseSpfitParameterText(text, { sourceName = 'SPFIT file', fallbackReduction = null } = {}) {
  const reduction = inferSpfitReductionFromText(text, fallbackReduction);
  const entries = [];
  const rawModelParams = {};
  const parameterErrors = {};
  const unsupportedLabels = [];
  let aMinusB = null;

  String(text || '').split(/\r?\n/).forEach((line) => {
    const entry = parseSpfitParameterLine(line);
    if (!entry) return;
    const key = getSpfitParamKeyByLabel(entry.label) || getSpfitParamKeyById(entry.id);
    const enriched = {
      ...entry,
      key,
      sourceName,
    };
    entries.push(enriched);
    if (!Number.isFinite(entry.value)) return;
    if (key === 'A_MINUS_B') {
      aMinusB = entry.value;
      return;
    }
    if (!MODEL_PARAM_KEYS.includes(key)) {
      if (entry.label) unsupportedLabels.push(entry.label);
      return;
    }
    rawModelParams[key] = entry.value;
    if (Number.isFinite(entry.error)) {
      parameterErrors[key] = Math.abs(entry.error);
    }
  });

  if (!Number.isFinite(rawModelParams.A) && Number.isFinite(aMinusB) && Number.isFinite(rawModelParams.B)) {
    rawModelParams.A = aMinusB + rawModelParams.B;
  }

  return {
    sourceName,
    reduction,
    entries,
    rawModelParams,
    parameterErrors,
    unsupportedLabels: Array.from(new Set(unsupportedLabels)),
  };
}

function parseSpfitNormalizedDiagonal(text, paramOrder = []) {
  if (typeof text !== 'string' || !text.trim()) return [];
  const lines = String(text).split(/\r?\n/);
  const start = lines.findIndex((line) => String(line || '').trim() === 'NORMALIZED DIAGONAL:');
  if (start < 0) return [];
  const result = [];
  for (let i = start + 1; i < lines.length; i++) {
    const line = String(lines[i] || '').trim();
    if (!line) {
      if (result.length) break;
      continue;
    }
    if (/^(MARQUARDT PARAMETER|WARNING:|MICROWAVE AVG|MICROWAVE RMS|END OF ITERATION|FIT COMPLETE)/i.test(line)) {
      break;
    }
    const matches = Array.from(line.matchAll(/(\d+)\s+([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)/gi));
    if (!matches.length) continue;
    matches.forEach((match) => {
      const index = Number.parseInt(match[1], 10);
      const value = parseSpfitNumber(match[2]);
      if (!Number.isFinite(index) || !Number.isFinite(value)) return;
      result.push({
        index,
        value,
        key: paramOrder[index - 1] || null,
      });
    });
  }
  return result;
}

function parseSpfitCorrelations(text, paramOrder = []) {
  if (typeof text !== 'string' || !text.trim()) return [];
  const lines = String(text).split(/\r?\n/);
  const start = lines.findIndex((line) => /FIT COMPLETE/i.test(String(line || '')));
  if (start < 0) return [];
  const correlations = [];
  for (let i = start + 1; i < lines.length; i++) {
    const line = String(lines[i] || '');
    const matches = Array.from(line.matchAll(/(\d{1,3})\s*(\d{1,3})\s*([+\-]?\d+\.\d+(?:[DE][+\-]?\d+)?)/gi));
    if (!matches.length) continue;
    matches.forEach((match) => {
      const leftIndex = Number.parseInt(match[1], 10);
      const rightIndex = Number.parseInt(match[2], 10);
      const value = parseSpfitNumber(match[3]);
      if (!Number.isFinite(leftIndex) || !Number.isFinite(rightIndex) || !Number.isFinite(value)) return;
      if (leftIndex === rightIndex) return;
      correlations.push({
        leftIndex,
        rightIndex,
        value,
        leftKey: paramOrder[leftIndex - 1] || null,
        rightKey: paramOrder[rightIndex - 1] || null,
      });
    });
  }
  return correlations;
}

function parseSpfitFitMetrics(text, paramOrder = []) {
  if (typeof text !== 'string' || !text.trim()) {
    return {
      microwaveRmsMHz: null,
      normalizedDiagonal: [],
      correlations: [],
    };
  }
  const rmsMatch = String(text).match(/MICROWAVE RMS\s*=\s*([+\-]?\d+(?:\.\d+)?(?:[DE][+\-]?\d+)?)\s*MHz/i);
  return {
    microwaveRmsMHz: rmsMatch ? parseSpfitNumber(rmsMatch[1]) : null,
    normalizedDiagonal: parseSpfitNormalizedDiagonal(text, paramOrder),
    correlations: parseSpfitCorrelations(text, paramOrder),
  };
}

function getSpfitImportSourceRank(name) {
  const lower = String(name || '').toLowerCase();
  if (lower.endsWith('.var')) return 400;
  if (lower.endsWith('.fit')) return 300;
  if (lower.endsWith('.out')) return 250;
  if (lower.endsWith('.par')) return 200;
  if (lower.endsWith('.txt')) return 100;
  return 0;
}

function chooseBestSpfitParameterSource(parsedFiles) {
  const candidates = (Array.isArray(parsedFiles) ? parsedFiles : [])
    .filter((item) => item && item.entries && item.entries.length);
  if (!candidates.length) return null;
  return candidates
    .slice()
    .sort((left, right) => {
      const rankDelta = getSpfitImportSourceRank(right.sourceName) - getSpfitImportSourceRank(left.sourceName);
      if (rankDelta !== 0) return rankDelta;
      const valueCountDelta = Object.keys(right.rawModelParams || {}).length - Object.keys(left.rawModelParams || {}).length;
      if (valueCountDelta !== 0) return valueCountDelta;
      return (right.entries?.length || 0) - (left.entries?.length || 0);
    })[0];
}

function chooseBestSpfitMetricsSource(fileMetrics) {
  const candidates = (Array.isArray(fileMetrics) ? fileMetrics : [])
    .filter((item) => item && (
      Number.isFinite(item.metrics?.microwaveRmsMHz)
      || (item.metrics?.normalizedDiagonal?.length || 0) > 0
      || (item.metrics?.correlations?.length || 0) > 0
    ));
  if (!candidates.length) return null;
  return candidates
    .slice()
    .sort((left, right) => getSpfitImportSourceRank(right.sourceName) - getSpfitImportSourceRank(left.sourceName))[0];
}

function normalizeSpfitBaselineAssignments(source) {
  if (!Array.isArray(source)) return [];
  return source.map((entry) => {
    const label = normalizeLabelValue(entry?.label);
    const fObs = Number(entry?.fObsMHz ?? entry?.fObs);
    const residual = Number(entry?.residualMHz ?? entry?.residual);
    if (!label || !Number.isFinite(fObs)) return null;
    return {
      label,
      fObs,
      residual: Number.isFinite(residual) ? residual : null,
    };
  }).filter(Boolean);
}

function getSpfitImportBaselineAssignments(session) {
  const sessionAssignments = normalizeSpfitBaselineAssignments(session?.spfitAssignments || session?.assignments);
  if (sessionAssignments.length) return sessionAssignments;
  const diagnosticsAssignments = normalizeSpfitBaselineAssignments(lastAutofitDiagnostics?.exportableAssignments);
  if (diagnosticsAssignments.length) return diagnosticsAssignments;
  return [];
}

function resolveSpectrumByReference(reference, predicate = null) {
  if (!reference || typeof reference !== 'object') return null;
  const items = Array.from(spectra.values()).filter((item) => (typeof predicate === 'function' ? predicate(item) : true));
  if (!items.length) return null;
  if (reference.id) {
    const byId = items.find((item) => item.id === reference.id);
    if (byId) return byId;
  }
  if (reference.name) {
    const target = String(reference.name).trim().toLowerCase();
    const byName = items.find((item) => String(item?.name || '').trim().toLowerCase() === target);
    if (byName) return byName;
  }
  return null;
}

function restoreImportedSessionSelections(session, theorySpectrum, expSpectrum) {
  const context = session?.context;
  if (!context || typeof context !== 'object') return;
  if (autofitUseInputAssignment) {
    autofitUseInputAssignment.checked = Boolean(context.useInputAssignments);
  }
  if (autofitAssignmentMode) {
    autofitAssignmentMode.value = context.manualMode && !context.useInputAssignments ? 'manual' : 'auto';
  }
  if (Array.isArray(context.selectedFitLabels)) {
    selectedFitLabels.length = 0;
    context.selectedFitLabels.forEach((label) => {
      const normalized = normalizeLabelValue(label);
      if (!normalized || selectedFitLabels.includes(normalized)) return;
      selectedFitLabels.push(normalized);
    });
  }
  if (Array.isArray(context.selectedControlLabels)) {
    selectedControlLabels.length = 0;
    context.selectedControlLabels.forEach((label) => {
      const normalized = normalizeLabelValue(label);
      if (!normalized || selectedControlLabels.includes(normalized)) return;
      selectedControlLabels.push(normalized);
    });
  }
  pendingManualSelection = null;
  manualAssignments.length = 0;
  if (Array.isArray(context.manualAssignments) && theorySpectrum && expSpectrum) {
    context.manualAssignments.forEach((item) => {
      const label = normalizeLabelValue(item?.label);
      const theoryX = Number(item?.theoryX);
      const expX = Number(item?.expX);
      if (!label || !Number.isFinite(theoryX) || !Number.isFinite(expX)) return;
      manualAssignments.push({
        label,
        theoryX,
        theorySpectrumId: theorySpectrum.id,
        expX,
        expSpectrumId: expSpectrum.id,
      });
    });
  }
  overlayHighlightsFit.clear();
  overlayHighlightsControl.clear();
  const labels = theorySpectrum?.labels || [];
  const freqs = theorySpectrum?.freqs || [];
  const labelToFreq = new Map();
  for (let i = 0; i < labels.length; i++) {
    const label = normalizeLabelValue(labels[i]);
    const freq = Number(freqs[i]);
    if (!label || !Number.isFinite(freq) || labelToFreq.has(label)) continue;
    labelToFreq.set(label, freq);
  }
  selectedFitLabels.forEach((label) => {
    const freq = labelToFreq.get(label);
    if (Number.isFinite(freq) && theorySpectrum) {
      overlayHighlightsFit.set(label, { label, x: freq, spectrumId: theorySpectrum.id });
    }
  });
  selectedControlLabels.forEach((label) => {
    const freq = labelToFreq.get(label);
    if (Number.isFinite(freq) && theorySpectrum) {
      overlayHighlightsControl.set(label, { label, x: freq, spectrumId: theorySpectrum.id });
    }
  });
  renderLabelLists();
  updateOverlayPlot();
}

function buildSpfitResidualComparison(assignments, theorySpectrum) {
  const baseline = normalizeSpfitBaselineAssignments(assignments);
  if (!baseline.length || !theorySpectrum) return null;
  const labels = theorySpectrum.labels || [];
  const freqs = theorySpectrum.freqs || [];
  const labelToFreq = new Map();
  for (let i = 0; i < labels.length; i++) {
    const label = normalizeLabelValue(labels[i]);
    const freq = Number(freqs[i]);
    if (!label || !Number.isFinite(freq) || labelToFreq.has(label)) continue;
    labelToFreq.set(label, freq);
  }

  const rows = [];
  let beforeSumSq = 0;
  let afterSumSq = 0;
  let beforeCount = 0;
  let afterCount = 0;
  let improvedCount = 0;
  let worsenedCount = 0;
  let unchangedCount = 0;
  let missingCount = 0;

  baseline.forEach((item) => {
    const predictedFreq = labelToFreq.get(item.label);
    const beforeResidual = Number(item.residual);
    const afterResidual = Number.isFinite(predictedFreq) ? item.fObs - predictedFreq : Number.NaN;
    const beforeAbs = Math.abs(beforeResidual);
    const afterAbs = Math.abs(afterResidual);
    let trend = 'same';
    if (!Number.isFinite(predictedFreq) || !Number.isFinite(afterResidual)) {
      missingCount += 1;
    } else if (Number.isFinite(beforeResidual)) {
      const delta = afterAbs - beforeAbs;
      if (delta < -1e-9) {
        trend = 'better';
        improvedCount += 1;
      } else if (delta > 1e-9) {
        trend = 'worse';
        worsenedCount += 1;
      } else {
        unchangedCount += 1;
      }
    }
    if (Number.isFinite(beforeResidual)) {
      beforeSumSq += beforeResidual ** 2;
      beforeCount += 1;
    }
    if (Number.isFinite(afterResidual)) {
      afterSumSq += afterResidual ** 2;
      afterCount += 1;
    }
    rows.push({
      label: item.label,
      fObs: item.fObs,
      predictedFreq,
      beforeResidual,
      afterResidual,
      trend,
    });
  });

  return {
    rows,
    assignmentCount: baseline.length,
    matchedCount: afterCount,
    missingCount,
    improvedCount,
    worsenedCount,
    unchangedCount,
    beforeRmsMHz: beforeCount ? Math.sqrt(beforeSumSq / beforeCount) : null,
    afterRmsMHz: afterCount ? Math.sqrt(afterSumSq / afterCount) : null,
  };
}

function buildSpfitParameterComparisonRows(beforeModelParams, afterModelParams, parameterErrors = {}) {
  const beforeNormalized = normalizeModelParams(beforeModelParams) || normalizeModelParams({
    ...(getAutofitInputModelParams() || {}),
    ...(beforeModelParams || {}),
  }) || {};
  const afterNormalized = normalizeModelParams(afterModelParams) || {};
  return MODEL_PARAM_KEYS
    .filter((key) => {
      const beforeValue = Number(beforeNormalized[key]);
      const afterValue = Number(afterNormalized[key]);
      const delta = Number.isFinite(beforeValue) && Number.isFinite(afterValue) ? afterValue - beforeValue : Number.NaN;
      if (['A', 'B', 'C'].includes(key)) return true;
      if (Number.isFinite(Number(parameterErrors[key]))) return true;
      return Number.isFinite(delta) && Math.abs(delta) > 1e-15;
    })
    .map((key) => {
      const beforeValue = Number(beforeNormalized[key]);
      const afterValue = Number(afterNormalized[key]);
      const delta = Number.isFinite(beforeValue) && Number.isFinite(afterValue) ? afterValue - beforeValue : Number.NaN;
      return {
        key,
        label: getModelParamDisplayLabel(key),
        beforeValue,
        afterValue,
        delta,
        error: Number(parameterErrors[key]),
      };
    });
}

function summarizeTopCorrelations(correlations, limit = 5) {
  if (!Array.isArray(correlations) || !correlations.length) return [];
  const unique = new Map();
  correlations.forEach((entry) => {
    const left = Number(entry?.leftIndex);
    const right = Number(entry?.rightIndex);
    const value = Number(entry?.value);
    if (!Number.isFinite(left) || !Number.isFinite(right) || !Number.isFinite(value)) return;
    const first = Math.min(left, right);
    const second = Math.max(left, right);
    const key = `${first}:${second}`;
    if (!unique.has(key) || Math.abs(value) > Math.abs(unique.get(key).value)) {
      unique.set(key, {
        leftIndex: first,
        rightIndex: second,
        leftKey: entry.leftKey || null,
        rightKey: entry.rightKey || null,
        value,
      });
    }
  });
  return Array.from(unique.values())
    .sort((left, right) => Math.abs(right.value) - Math.abs(left.value))
    .slice(0, Math.max(1, limit));
}

function getSpfitImportBadgeClassName(appliedState) {
  return appliedState === 'partial'
    ? 'spfit-import-badge spfit-import-badge--partial'
    : 'spfit-import-badge spfit-import-badge--applied';
}

function appendSpfitImportList(container, title, items, variant = 'warning') {
  if (!Array.isArray(items) || !items.length) return;
  const group = document.createElement('div');
  group.className = 'spfit-import-group';
  const heading = document.createElement('h4');
  heading.textContent = title;
  const list = document.createElement('ul');
  list.className = `spfit-import-list spfit-import-list--${variant}`;
  items.forEach((item) => {
    const li = document.createElement('li');
    li.textContent = item;
    list.appendChild(li);
  });
  group.appendChild(heading);
  group.appendChild(list);
  container.appendChild(group);
}

function renderSpfitImportReport(report) {
  lastSpfitImportSummary = report || null;
  if (!spfitImportReport) return;
  spfitImportReport.innerHTML = '';
  if (!report) {
    spfitImportReport.hidden = true;
    return;
  }
  spfitImportReport.hidden = false;

  const header = document.createElement('div');
  header.className = 'spfit-import-report__header';
  const titleBlock = document.createElement('div');
  const title = document.createElement('h3');
  title.className = 'spfit-import-report__title';
  title.textContent = 'SPFIT Import Summary';
  const subtitle = document.createElement('p');
  subtitle.className = 'spfit-import-report__subtitle';
  subtitle.textContent = `Imported from ${report.sourceLabel}. ${report.fileCount} file${report.fileCount === 1 ? '' : 's'} processed.`;
  titleBlock.appendChild(title);
  titleBlock.appendChild(subtitle);
  const badge = document.createElement('span');
  badge.className = getSpfitImportBadgeClassName(report.appliedState);
  badge.textContent = report.appliedLabel;
  header.appendChild(titleBlock);
  header.appendChild(badge);
  spfitImportReport.appendChild(header);

  const stats = document.createElement('div');
  stats.className = 'autofit-diagnostics__stats';
  appendDiagnosticStat(stats, 'Theory Spectrum', report.theorySpectrumName || 'n/a');
  appendDiagnosticStat(stats, 'Parameter Source', report.parameterSource || 'n/a');
  appendDiagnosticStat(stats, 'Reduction', report.reduction || 'n/a');
  appendDiagnosticStat(stats, 'Imported Params', String(report.importedParamCount));
  appendDiagnosticStat(
    stats,
    'SPFIT RMS',
    Number.isFinite(report.microwaveRmsMHz) ? `${formatDiagnosticNumber(report.microwaveRmsMHz, 4)} MHz` : 'n/a'
  );
  appendDiagnosticStat(
    stats,
    'Residual RMS',
    report.residualComparison && Number.isFinite(report.residualComparison.afterRmsMHz)
      ? `${formatDiagnosticNumber(report.residualComparison.beforeRmsMHz, 4)} -> ${formatDiagnosticNumber(report.residualComparison.afterRmsMHz, 4)} MHz`
      : 'n/a'
  );
  appendDiagnosticStat(
    stats,
    'Improved / Worse',
    report.residualComparison
      ? `${report.residualComparison.improvedCount} / ${report.residualComparison.worsenedCount}`
      : 'n/a'
  );
  appendDiagnosticStat(
    stats,
    'Correlations',
    report.topCorrelations.length ? String(report.topCorrelations.length) : 'n/a'
  );
  spfitImportReport.appendChild(stats);

  const contextGroup = document.createElement('div');
  contextGroup.className = 'spfit-import-group';
  const contextTitle = document.createElement('h4');
  contextTitle.textContent = 'Applied Context';
  const contextText = document.createElement('p');
  contextText.className = 'spfit-import-text';
  contextText.textContent = report.sessionLoaded
    ? 'Session metadata restored fit/control selections and provided the baseline assignment set used for residual comparison.'
    : 'No session.json was available; the comparison uses only the spectra currently loaded in WMS-FitRot.';
  contextGroup.appendChild(contextTitle);
  contextGroup.appendChild(contextText);
  spfitImportReport.appendChild(contextGroup);

  appendSpfitImportList(spfitImportReport, 'Warnings', report.warnings, 'warning');

  if (report.topCorrelations.length) {
    const corrGroup = document.createElement('div');
    corrGroup.className = 'spfit-import-group';
    const corrTitle = document.createElement('h4');
    corrTitle.textContent = 'Strongest Parsed Correlations';
    corrGroup.appendChild(corrTitle);
    const wrap = document.createElement('div');
    wrap.className = 'spfit-import-table-wrap';
    const table = document.createElement('table');
    table.className = 'spfit-import-table';
    const head = document.createElement('thead');
    const headRow = document.createElement('tr');
    ['Left', 'Right', 'Correlation'].forEach((label) => {
      const th = document.createElement('th');
      th.textContent = label;
      headRow.appendChild(th);
    });
    head.appendChild(headRow);
    table.appendChild(head);
    const body = document.createElement('tbody');
    report.topCorrelations.forEach((entry) => {
      const row = document.createElement('tr');
      const left = document.createElement('td');
      left.textContent = entry.leftKey || `#${entry.leftIndex}`;
      const right = document.createElement('td');
      right.textContent = entry.rightKey || `#${entry.rightIndex}`;
      const value = document.createElement('td');
      value.textContent = formatDiagnosticNumber(entry.value, 4);
      row.appendChild(left);
      row.appendChild(right);
      row.appendChild(value);
      body.appendChild(row);
    });
    table.appendChild(body);
    wrap.appendChild(table);
    corrGroup.appendChild(wrap);
    spfitImportReport.appendChild(corrGroup);
  }

  if (report.parameterRows.length) {
    const paramGroup = document.createElement('div');
    paramGroup.className = 'spfit-import-group';
    const paramTitle = document.createElement('h4');
    paramTitle.textContent = 'Parameter Changes';
    paramGroup.appendChild(paramTitle);
    const wrap = document.createElement('div');
    wrap.className = 'spfit-import-table-wrap';
    const table = document.createElement('table');
    table.className = 'spfit-import-table';
    const head = document.createElement('thead');
    const headRow = document.createElement('tr');
    ['Parameter', 'Before', 'After', 'Delta', 'Error'].forEach((label) => {
      const th = document.createElement('th');
      th.textContent = label;
      headRow.appendChild(th);
    });
    head.appendChild(headRow);
    table.appendChild(head);
    const body = document.createElement('tbody');
    report.parameterRows.forEach((entry) => {
      const row = document.createElement('tr');
      const label = document.createElement('td');
      label.textContent = entry.label;
      const before = document.createElement('td');
      before.textContent = formatDiagnosticNumber(entry.beforeValue, 6);
      const after = document.createElement('td');
      after.textContent = formatDiagnosticNumber(entry.afterValue, 6);
      const delta = document.createElement('td');
      delta.textContent = formatDiagnosticNumber(entry.delta, 6);
      const deltaAbs = Math.abs(Number(entry.delta));
      delta.className = !Number.isFinite(deltaAbs) || deltaAbs <= 1e-15
        ? 'spfit-import-delta--same'
        : (entry.delta > 0 ? 'spfit-import-delta--positive' : 'spfit-import-delta--negative');
      const error = document.createElement('td');
      error.textContent = Number.isFinite(entry.error) ? formatDiagnosticNumber(entry.error, 4) : 'n/a';
      row.appendChild(label);
      row.appendChild(before);
      row.appendChild(after);
      row.appendChild(delta);
      row.appendChild(error);
      body.appendChild(row);
    });
    table.appendChild(body);
    wrap.appendChild(table);
    paramGroup.appendChild(wrap);
    spfitImportReport.appendChild(paramGroup);
  }

  if (report.residualComparison?.rows?.length) {
    const residualGroup = document.createElement('div');
    residualGroup.className = 'spfit-import-group';
    const residualTitle = document.createElement('h4');
    residualTitle.textContent = 'Residual Comparison On Exported Lines';
    residualGroup.appendChild(residualTitle);
    const wrap = document.createElement('div');
    wrap.className = 'spfit-import-table-wrap';
    const table = document.createElement('table');
    table.className = 'spfit-import-table';
    const head = document.createElement('thead');
    const headRow = document.createElement('tr');
    ['Label', 'Before', 'After', 'Trend'].forEach((label) => {
      const th = document.createElement('th');
      th.textContent = label;
      headRow.appendChild(th);
    });
    head.appendChild(headRow);
    table.appendChild(head);
    const body = document.createElement('tbody');
    report.residualComparison.rows.forEach((entry) => {
      const row = document.createElement('tr');
      const label = document.createElement('td');
      label.textContent = entry.label;
      const before = document.createElement('td');
      before.textContent = Number.isFinite(entry.beforeResidual)
        ? `${formatDiagnosticNumber(entry.beforeResidual, 4)} MHz`
        : 'n/a';
      const after = document.createElement('td');
      after.textContent = Number.isFinite(entry.afterResidual)
        ? `${formatDiagnosticNumber(entry.afterResidual, 4)} MHz`
        : 'missing';
      const trend = document.createElement('td');
      trend.textContent = entry.trend;
      trend.className = `spfit-import-delta--${entry.trend === 'better' ? 'better' : entry.trend === 'worse' ? 'worse' : 'same'}`;
      row.appendChild(label);
      row.appendChild(before);
      row.appendChild(after);
      row.appendChild(trend);
      body.appendChild(row);
    });
    table.appendChild(body);
    wrap.appendChild(table);
    residualGroup.appendChild(wrap);
    spfitImportReport.appendChild(residualGroup);
  }
}

async function parseSpfitImportFiles(files) {
  const fileEntries = await Promise.all((Array.isArray(files) ? files : []).map(async (file) => ({
    name: file.name,
    text: await file.text(),
  })));

  let session = null;
  const warnings = [];
  fileEntries.forEach((entry) => {
    const lower = String(entry.name || '').toLowerCase();
    if (!(lower === 'session.json' || lower.endsWith('/session.json') || lower.endsWith('.json'))) return;
    try {
      const parsed = JSON.parse(entry.text);
      if (parsed && typeof parsed === 'object' && (parsed.exportType === 'spfit-seed' || parsed.spfitAssignments || parsed.finalModelParams)) {
        session = parsed;
      }
    } catch (error) {
      warnings.push(`Unable to parse ${entry.name} as JSON.`);
    }
  });

  const parsedFiles = [];
  const fileMetrics = [];
  fileEntries.forEach((entry) => {
    const lower = String(entry.name || '').toLowerCase();
    if (lower.endsWith('.json')) return;
    const parsed = parseSpfitParameterText(entry.text, {
      sourceName: entry.name,
      fallbackReduction: session?.reduction || null,
    });
    if (parsed.entries.length) {
      parsedFiles.push(parsed);
    }
    fileMetrics.push({
      sourceName: entry.name,
      metrics: parseSpfitFitMetrics(entry.text, parsed.entries.map((item) => item.key).filter(Boolean)),
    });
  });

  const bestParamSource = chooseBestSpfitParameterSource(parsedFiles);
  const bestMetricsSource = chooseBestSpfitMetricsSource(fileMetrics);
  if (!bestParamSource && !session?.finalModelParams) {
    throw new Error('No usable SPFIT parameters were found. Import at least a .var or .par file.');
  }

  const reduction = bestParamSource?.reduction || session?.reduction || null;
  const theorySpectrum = resolveSpectrumByReference(session?.theorySpectrum, (item) => item.isTheoretical)
    || getSpectrumById(autofitTheorySelect?.value)
    || Array.from(spectra.values()).find((item) => item.isTheoretical)
    || null;
  const expSpectrum = resolveSpectrumByReference(session?.experimentalSpectrum, (item) => !item.isTheoretical)
    || getSpectrumById(autofitExpSelect?.value)
    || Array.from(spectra.values()).find((item) => !item.isTheoretical)
    || null;

  const beforeModelParams = normalizeModelParams(
    (theorySpectrum && getSpectrumModelParams(theorySpectrum))
    || session?.finalModelParams
    || getAutofitInputModelParams()
  );
  const importedRawParams = {
    ...(session?.finalModelParams || {}),
    ...(bestParamSource?.rawModelParams || {}),
  };
  const importedModelParams = normalizeModelParams({
    ...(beforeModelParams || {}),
    ...importedRawParams,
  });
  if (!importedModelParams) {
    throw new Error('Imported SPFIT parameters are incomplete for A, B, C.');
  }

  const parameterErrors = {
    ...(bestParamSource?.parameterErrors || {}),
  };
  const metrics = bestMetricsSource?.metrics || {
    microwaveRmsMHz: null,
    normalizedDiagonal: [],
    correlations: [],
  };

  if (bestParamSource?.unsupportedLabels?.length) {
    warnings.push(`Ignored unsupported SPFIT labels: ${bestParamSource.unsupportedLabels.slice(0, 6).join(', ')}${bestParamSource.unsupportedLabels.length > 6 ? ', …' : ''}.`);
  }
  if (!session) {
    warnings.push('session.json not provided; only the currently loaded WMS-FitRot context can be restored.');
  }
  if (theorySpectrum == null) {
    warnings.push('No theoretical spectrum matched the imported session; the refined constants cannot be applied until one is loaded.');
  }

  return {
    session,
    theorySpectrum,
    expSpectrum,
    reduction,
    beforeModelParams,
    importedModelParams,
    parameterErrors,
    parameterSource: bestParamSource?.sourceName || (session ? 'session.json' : 'n/a'),
    importedParamCount: Object.keys(bestParamSource?.rawModelParams || session?.finalModelParams || {}).length,
    microwaveRmsMHz: metrics.microwaveRmsMHz,
    normalizedDiagonal: metrics.normalizedDiagonal,
    correlations: metrics.correlations,
    baselineAssignments: getSpfitImportBaselineAssignments(session),
    warnings,
    fileNames: fileEntries.map((item) => item.name),
  };
}

function finalizePendingSpfitImport() {
  if (!pendingSpfitImportState) return null;
  const state = pendingSpfitImportState;
  pendingSpfitImportState = null;
  const theorySpectrum = getSpectrumById(state.theorySpectrumId)
    || getSpectrumById(autofitTheorySelect?.value)
    || null;
  const residualComparison = buildSpfitResidualComparison(state.baselineAssignments, theorySpectrum);
  const topCorrelations = summarizeTopCorrelations(state.correlations, 5);
  const report = {
    appliedState: residualComparison && residualComparison.missingCount > 0 ? 'partial' : 'applied',
    appliedLabel: residualComparison && residualComparison.missingCount > 0 ? 'Imported With Gaps' : 'SPFIT Imported',
    sourceLabel: state.fileNames.join(', '),
    fileCount: state.fileNames.length,
    sessionLoaded: Boolean(state.session),
    theorySpectrumName: theorySpectrum?.name || state.theorySpectrumName || null,
    reduction: state.reduction || null,
    parameterSource: state.parameterSource || null,
    importedParamCount: state.importedParamCount || 0,
    microwaveRmsMHz: state.microwaveRmsMHz,
    topCorrelations,
    warnings: state.warnings || [],
    parameterRows: buildSpfitParameterComparisonRows(state.beforeModelParams, state.importedModelParams, state.parameterErrors),
    residualComparison,
  };
  renderSpfitImportReport(report);
  return report;
}

async function handleSpfitImportFiles(files) {
  const selectedFiles = Array.isArray(files) ? files : [];
  if (!selectedFiles.length) return;
  pendingSpfitImportState = null;
  renderSpfitImportReport(null);
  setAutofitWarningVisible(false);
  setAutofitStatus('Parsing SPFIT import package…');
  setAutofitProgress(0.08);

  const parsed = await parseSpfitImportFiles(selectedFiles);
  if (!parsed.theorySpectrum) {
    renderSpfitImportReport({
      appliedState: 'partial',
      appliedLabel: 'Import Parsed',
      sourceLabel: parsed.fileNames.join(', '),
      fileCount: parsed.fileNames.length,
      sessionLoaded: Boolean(parsed.session),
      theorySpectrumName: null,
      reduction: parsed.reduction,
      parameterSource: parsed.parameterSource,
      importedParamCount: parsed.importedParamCount,
      microwaveRmsMHz: parsed.microwaveRmsMHz,
      topCorrelations: summarizeTopCorrelations(parsed.correlations, 5),
      warnings: parsed.warnings,
      parameterRows: buildSpfitParameterComparisonRows(parsed.beforeModelParams, parsed.importedModelParams, parsed.parameterErrors),
      residualComparison: null,
    });
    setAutofitProgress(1);
    setAutofitStatus('SPFIT files parsed, but no theoretical spectrum is loaded to apply the refinement.');
    return;
  }

  if (autofitTheorySelect && spectra.has(parsed.theorySpectrum.id)) {
    autofitTheorySelect.value = parsed.theorySpectrum.id;
  }
  if (parsed.expSpectrum && autofitExpSelect && spectra.has(parsed.expSpectrum.id)) {
    autofitExpSelect.value = parsed.expSpectrum.id;
  }
  if (parsed.session) {
    lastAutofitContext = cloneJsonValue(parsed.session.context) || lastAutofitContext;
    lastAutofitSummary = cloneJsonValue(parsed.session.summary) || lastAutofitSummary;
    restoreImportedSessionSelections(parsed.session, parsed.theorySpectrum, parsed.expSpectrum);
  }

  parsed.theorySpectrum.fittedConstants = {
    A: parsed.importedModelParams.A,
    B: parsed.importedModelParams.B,
    C: parsed.importedModelParams.C,
  };
  parsed.theorySpectrum.modelParams = parsed.importedModelParams;
  parsed.theorySpectrum.quarticScale = 1;
  parsed.theorySpectrum.quarticReferenceParams = normalizeQuarticReferenceParams(parsed.importedModelParams);
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();

  pendingSpfitImportState = {
    ...parsed,
    theorySpectrumId: parsed.theorySpectrum.id,
    theorySpectrumName: parsed.theorySpectrum.name,
  };

  setAutofitProgress(0.35);
  if (!recalcTheoreticalSpectrum(parsed.importedModelParams)) {
    const report = finalizePendingSpfitImport();
    setAutofitProgress(1);
    setAutofitStatus(report
      ? 'SPFIT parameters applied, but the theory spectrum could not be recalculated automatically.'
      : 'SPFIT parameters parsed, but the theory spectrum could not be recalculated automatically.');
    return;
  }
  setAutofitStatus(`Applying refined constants from ${parsed.parameterSource}…`);
}

async function populateSpfitSeedDownload({
  outputText,
  modelParams,
  quarticScaleState,
  theorySpectrum,
  expSpectrum,
  diagnostics = null,
}) {
  if (!autofitDownloads) return;
  try {
    const result = await createSpfitSeedExport({
      outputText,
      modelParams,
      quarticScaleState,
      theorySpectrum,
      expSpectrum,
      diagnostics,
    });
    autofitDownloads.appendChild(
      createDownloadLinkFromBlob(
        `${result.stem}_spfit_seed.zip`,
        result.blob,
        `Download SPFIT seed ZIP (${result.lineCount} lines)`
      )
    );
    autofitDownloads.appendChild(
      createAutofitDownloadMessage(
        diagnostics?.readyLabel
          ? `SPFIT seed ready: ${result.paramCount} parameters, ${result.lineCount} approved lines. ${diagnostics.readyLabel}.`
          : `SPFIT seed ready: ${result.paramCount} parameters, ${result.lineCount} lines.`
      )
    );
    result.warnings.forEach((warning) => {
      autofitDownloads.appendChild(createAutofitDownloadMessage(`Warning: ${warning}`));
    });
  } catch (error) {
    console.warn('Unable to prepare SPFIT export:', error);
    autofitDownloads.appendChild(
      createAutofitDownloadMessage(`SPFIT export unavailable: ${String(error?.message || error)}`)
    );
  }
}

function applyAutomaticControlLines(summary) {
  if (!summary) return false;
  if (isDirectAssignmentMode()) return false;
  if (selectedControlLabels.length > 0) return false;
  if (String(summary.controlMode || '').toUpperCase() !== 'AUTO') return false;
  if (!Array.isArray(summary.controlLines) || !summary.controlLines.length) return false;

  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  if (!theorySpectrum) return false;
  const labelToFreq = new Map();
  const labels = theorySpectrum.labels || [];
  const freqs = theorySpectrum.freqs || [];
  for (let i = 0; i < labels.length; i++) {
    const label = normalizeLabelValue(labels[i]);
    const fx = Number(freqs[i]);
    if (!label || !Number.isFinite(fx) || labelToFreq.has(label)) continue;
    labelToFreq.set(label, fx);
  }

  selectedControlLabels.length = 0;
  overlayHighlightsControl.clear();
  summary.controlLines.forEach((item) => {
    const label = normalizeLabelValue(item?.label);
    if (!label || selectedControlLabels.includes(label)) return;
    selectedControlLabels.push(label);
    const theoryX = labelToFreq.get(label);
    const x = Number.isFinite(theoryX) ? theoryX : item?.fObs;
    if (Number.isFinite(x)) {
      overlayHighlightsControl.set(label, { label, x, spectrumId: theorySpectrum.id });
    }
  });
  if (!selectedControlLabels.length) return false;
  renderLabelLists();
  updateOverlayPlot();
  showOverlayToast(`Automatic control lines selected: ${selectedControlLabels.length}`);
  return true;
}

function recalcTheoreticalSpectrum(constants) {
  if (!autofitInputText) return false;
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  if (!theorySpectrum) {
    setAutofitStatus('Completed. Theory spectrum not found for refresh.');
    return false;
  }
  const modelParams = normalizeModelParams({
    ...(getSpectrumModelParams(theorySpectrum) || {}),
    ...(constants && typeof constants === 'object' ? constants : {}),
  });
  const worker = initWorker();
  worker.postMessage({
    type: 'recalc',
    payload: {
      inputText: autofitInputText,
      A: constants?.A,
      B: constants?.B,
      C: constants?.C,
      currentModelParams: modelParams,
      options: {
        jmax: autofitOptions.jmax?.value || '10',
        intensityCut: autofitOptions.intensityCut?.value || '1e-9',
      },
    },
  });
  return true;
}

function applyRecalculatedSpectrum(data) {
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  if (!theorySpectrum) return;
  const freqs = Array.isArray(data?.freqs) ? data.freqs : [];
  const ints = Array.isArray(data?.ints) ? data.ints : [];
  const labels = Array.isArray(data?.labels) ? data.labels : [];
  if (!freqs.length || !ints.length) return;

  const baseAbs = ints.map((value) => Math.abs(Number(value)));
  theorySpectrum.freqs = freqs;
  theorySpectrum.baseAbs = baseAbs;
  theorySpectrum.labels = labels;
  theorySpectrum.mode = 'stick';
  theorySpectrum.stick = {
    freqs,
    baseInts: baseAbs,
    labels,
  };
  theorySpectrum.previewView = 'stick';
  theorySpectrum.rawText = formatStickDat(freqs, baseAbs, labels);

  const settings = ensureOverlaySettings(theorySpectrum);
  settings.normalize = true;
  settings.scale = 1;

  updateHighlightsForSpectrum(theorySpectrum);
  renderUploadedTabs();
  renderOverlayList();
  updateOverlayPlot();
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();
  updateRestoreInitialButtonState();
  showOverlayToast('Theory spectrum updated with fitted constants');
}

function updateHighlightsForSpectrum(spectrum) {
  if (!spectrum || !spectrum.isTheoretical) return;
  const labelToFreq = new Map();
  const labels = spectrum.labels || [];
  for (let i = 0; i < labels.length; i++) {
    const lbl = labels[i];
    if (lbl == null) continue;
    labelToFreq.set(lbl, spectrum.freqs[i]);
  }
  const updateMap = (map) => {
    map.forEach((item, key) => {
      if (item.spectrumId !== spectrum.id) return;
      const fx = labelToFreq.get(item.label);
      if (Number.isFinite(fx)) {
        item.x = fx;
      }
    });
  };
  updateMap(overlayHighlightsFit);
  updateMap(overlayHighlightsControl);
  manualAssignments.forEach((item) => {
    if (item.theorySpectrumId !== spectrum.id) return;
    const fx = labelToFreq.get(item.label);
    if (Number.isFinite(fx)) {
      item.theoryX = fx;
    }
  });
  if (pendingManualSelection?.theorySpectrumId === spectrum.id) {
    const fx = labelToFreq.get(pendingManualSelection.label);
    if (Number.isFinite(fx)) {
      pendingManualSelection.theoryX = fx;
    }
  }
}

function getSpectrumById(id) {
  if (!id || !spectra.has(id)) return null;
  return spectra.get(id);
}

function cloneSpectrumFrequencies(values) {
  return Array.isArray(values)
    ? values
      .map((value) => Number(value))
      .filter(Number.isFinite)
    : [];
}

function cloneSpectrumLabels(values) {
  return Array.isArray(values)
    ? values.map((value) => (value == null ? null : String(value)))
    : [];
}

function normalizeSerializedSpectrumState(state) {
  if (!state || typeof state !== 'object') return null;
  const mode = state.mode === 'profile' ? 'profile' : 'stick';
  const previewView = state.previewView === 'profile' ? 'profile' : 'stick';
  const text = typeof state.text === 'string' ? state.text : '';
  const fittedConstants = {};
  ['A', 'B', 'C'].forEach((key) => {
    const value = Number(state?.fittedConstants?.[key]);
    if (Number.isFinite(value)) {
      fittedConstants[key] = value;
    }
  });
  return {
    mode,
    previewView,
    text,
    modelParams: normalizeModelParams(state.modelParams),
    fittedConstants: Object.keys(fittedConstants).length ? fittedConstants : null,
    quarticScale: normalizeQuarticScale(state.quarticScale, 1),
    quarticReferenceParams: normalizeQuarticReferenceParams(state.quarticReferenceParams),
  };
}

function formatStickDat(freqs, intensities, labels = []) {
  const lines = [];
  for (let i = 0; i < freqs.length; i++) {
    const f = Number(freqs[i]);
    const v = Number(intensities[i]);
    if (!Number.isFinite(f) || !Number.isFinite(v)) continue;
    const label = labels[i] == null ? '' : String(labels[i]).trim();
    lines.push(label ? `${f} ${v} ${label}` : `${f} ${v}`);
  }
  return lines.join('\n');
}

function serializeSpectrumText(spectrum) {
  if (!spectrum) return '';
  if (spectrum.mode === 'stick') {
    return formatStickDat(spectrum.freqs || [], spectrum.baseAbs || [], spectrum.labels || []);
  }
  return typeof spectrum.rawText === 'string' ? spectrum.rawText : '';
}

function buildSpectrumInitialStateSnapshot(spectrum, {
  textOverride = null,
  modelParamsOverride = null,
  quarticReferenceParamsOverride = null,
  quarticScaleOverride = null,
} = {}) {
  if (!spectrum) return null;
  const text = typeof textOverride === 'string' ? textOverride : serializeSpectrumText(spectrum);
  return {
    mode: spectrum.mode === 'profile' ? 'profile' : 'stick',
    previewView: spectrum.previewView === 'profile' ? 'profile' : 'stick',
    text,
    modelParams: normalizeModelParams(modelParamsOverride ?? spectrum.modelParams),
    fittedConstants: spectrum?.fittedConstants ? cloneJsonValue(spectrum.fittedConstants) : null,
    quarticScale: normalizeQuarticScale(
      quarticScaleOverride ?? spectrum.quarticScale,
      1
    ),
    quarticReferenceParams: normalizeQuarticReferenceParams(
      quarticReferenceParamsOverride ?? spectrum.quarticReferenceParams
    ),
  };
}

function syncSpectrumCounterWithId(id) {
  const match = /^spectrum-(\d+)$/.exec(String(id || ''));
  if (!match) return;
  const value = Number.parseInt(match[1], 10);
  if (Number.isFinite(value)) {
    spectrumCounter = Math.max(spectrumCounter, value);
  }
}

function applyStoredOverlaySettings(spectrum, state) {
  if (!spectrum || !state || typeof state !== 'object') return;
  const current = ensureOverlaySettings(spectrum);
  current.enabled = state.enabled !== false;
  current.view = spectrum.mode === 'profile' && state.view === 'profile' ? 'profile' : 'stick';
  current.intensityMode = normalizeIntensityMode(state.intensityMode);
  current.normalize = Boolean(state.normalize);
  const scale = Number(state.scale);
  current.scale = Number.isFinite(scale) ? scale : 1;
  current.color = typeof state.color === 'string' && state.color.trim()
    ? state.color
    : current.color;
}

function buildTheoryLabelFrequencyMap(theorySpectrum) {
  const labelToFreq = new Map();
  const labels = theorySpectrum?.labels || [];
  const freqs = theorySpectrum?.freqs || [];
  for (let i = 0; i < labels.length; i++) {
    const label = normalizeLabelValue(labels[i]);
    const freq = Number(freqs[i]);
    if (!label || !Number.isFinite(freq) || labelToFreq.has(label)) continue;
    labelToFreq.set(label, freq);
  }
  return labelToFreq;
}

function formatAutofitLabelText(labels, theorySpectrum) {
  const labelToFreq = buildTheoryLabelFrequencyMap(theorySpectrum);
  return (Array.isArray(labels) ? labels : [])
    .map((value) => {
      const label = normalizeLabelValue(value);
      if (!label) return null;
      const freq = labelToFreq.get(label);
      return Number.isFinite(freq) ? `${label} | ${freq.toFixed(9)}` : label;
    })
    .filter(Boolean)
    .join('\n');
}

function formatDirectAssignmentsText(assignments = manualAssignments) {
  return (Array.isArray(assignments) ? assignments : [])
    .map((item) => {
      const label = normalizeLabelValue(item?.label);
      const theoryX = Number(item?.theoryX);
      const expX = Number(item?.expX);
      if (!label || !Number.isFinite(expX)) return null;
      return Number.isFinite(theoryX)
        ? `${label} | ${theoryX.toFixed(9)} | ${expX.toFixed(9)}`
        : `${label} | ${expX.toFixed(9)}`;
    })
    .filter(Boolean)
    .join('\n');
}

function inferManualRequiredJMax(assignments = manualAssignments) {
  let maxJ = Number.NEGATIVE_INFINITY;
  const pattern = /J\s*=\s*(-?\d+)/g;
  (Array.isArray(assignments) ? assignments : []).forEach((item) => {
    const label = String(item?.label || '');
    pattern.lastIndex = 0;
    let match = pattern.exec(label);
    while (match) {
      const value = Number.parseInt(match[1], 10);
      if (Number.isFinite(value)) {
        maxJ = Math.max(maxJ, value);
      }
      match = pattern.exec(label);
    }
  });
  return Number.isFinite(maxJ) ? Math.max(0, maxJ) : null;
}

function getAutofitOptionValues(manualMode = false, assignments = manualAssignments) {
  const values = {
    deltaAssign: String(autofitOptions.deltaAssign?.value || '0.1'),
    deltaMatch: String(autofitOptions.deltaMatch?.value || '0.1'),
    kmax: String(autofitOptions.kmax?.value || '6'),
    maxCombos: String(autofitOptions.maxCombos?.value || '5000'),
    jmax: String(autofitOptions.jmax?.value || '30'),
    intensityCut: String(autofitOptions.intensityCut?.value || '1e-9'),
    span: String(autofitOptions.span?.value || '0.3'),
    finalSpan: String(autofitOptions.finalSpan?.value || '0.1'),
    mctrl: String(autofitOptions.mctrl?.value || '50'),
    lambdaRms: String(autofitOptions.lambdaRms?.value || '1.0'),
    lambdaBias: String(autofitOptions.lambdaBias?.value || '0.3'),
    kappaCol: String(autofitOptions.kappaCol?.value || '0.5'),
    usePowell: autofitOptions.usePowell ? Boolean(autofitOptions.usePowell.checked) : true,
    sampleTemp: String(autofitOptions.sampleTemp?.value || '1.0'),
    beamSize: String(autofitOptions.beamSize?.value || '200'),
    neighborK: String(autofitOptions.neighborK?.value || '3'),
    lsLoss: String(autofitOptions.lsLoss?.value || 'soft_l1'),
    lsFScaleMode: String(autofitOptions.lsFScaleMode?.value || 'auto'),
    lsFScale: String(autofitOptions.lsFScale?.value || '1.0'),
    lsMaxNfev: String(autofitOptions.lsMaxNfev?.value || '200'),
  };
  const notes = [];

  if (manualMode) {
    const requiredJ = inferManualRequiredJMax(assignments);
    const currentJ = Number(values.jmax);
    if (Number.isFinite(requiredJ) && (!Number.isFinite(currentJ) || currentJ < requiredJ)) {
      values.jmax = String(requiredJ);
      notes.push(`J max auto-set to ${requiredJ}`);
    }
    const cut = Number(values.intensityCut);
    if (!Number.isFinite(cut) || cut > 1e-9) {
      values.intensityCut = '1e-9';
      notes.push('Intensity cut auto-set to 1e-9');
    }
  }

  return { values, notes };
}

function setAutofitStatus(message) {
  if (!autofitStatus) return;
  autofitStatus.textContent = message;
}

function setAutofitProgress(progress) {
  if (!autofitProgressBar) return;
  const clamped = Math.min(100, Math.max(0, Number(progress) * 100));
  autofitProgressBar.style.width = `${clamped}%`;
}

function resetAutofitProgress() {
  setAutofitProgress(0);
  if (autofitDownloads) autofitDownloads.innerHTML = '';
  renderSensitivityDashboard(null);
  pendingSpfitImportState = null;
  renderAutofitDiagnostics(null);
  renderSpfitImportReport(null);
}

function buildCliArgs(expMode, {
  hasControlLabels = false,
  manualMode = false,
  useInputAssignments = false,
  optionValues = null,
  currentModelParams = null,
  fitParams = null,
  customParamBounds = null,
  quarticScaleState = null,
} = {}) {
  const opts = optionValues || getAutofitOptionValues(manualMode).values;
  const normalizedModelParams = normalizeModelParams(currentModelParams);
  const selectedFitParams = Array.isArray(fitParams) ? fitParams : [];
  const args = [
    'python3 fitting_script.py',
    '--input-file input.txt',
    '--labels-file labels.txt',
    ...(manualMode ? ['--manual-assignments-file manual_assignments.txt'] : []),
    ...(useInputAssignments ? ['--use-input-assignments'] : []),
    ...(!manualMode && hasControlLabels ? ['--control-labels-file control_labels.txt'] : []),
    expMode === 'stick' ? '--peaks-file exp.dat' : '--spectrum-file exp.dat',
    `--delta-assign ${opts.deltaAssign}`,
    `--delta-match ${opts.deltaMatch}`,
    `--kmax ${opts.kmax}`,
    `--max-combos ${opts.maxCombos}`,
    `--J-max ${opts.jmax}`,
    `--intensity-cut ${opts.intensityCut}`,
    `--span ${opts.span}`,
    `--final-span ${opts.finalSpan}`,
    `--mctrl ${opts.mctrl}`,
    `--lambda-rms ${opts.lambdaRms}`,
    `--lambda-bias ${opts.lambdaBias}`,
    `--kappa-col ${opts.kappaCol}`,
    ...(opts.usePowell ? [] : ['--skip-powell']),
    `--sample-temp ${opts.sampleTemp}`,
    `--beam-size ${opts.beamSize}`,
    `--neighbor-k ${opts.neighborK}`,
    ...buildModelParamCliArgs(normalizedModelParams),
    ...buildQuarticScaleCliArgs(quarticScaleState),
    ...(selectedFitParams.length ? [`--fit-params ${selectedFitParams.join(' ')}`] : []),
    ...buildParamBoundCliArgs(customParamBounds),
    `--ls-loss ${opts.lsLoss}`,
    `--ls-f-scale ${opts.lsFScaleMode === 'custom' ? opts.lsFScale : 'auto'}`,
    `--ls-max-nfev ${opts.lsMaxNfev}`,
  ];
  return args.join(' ');
}

function createDownloadLink(name, text) {
  const blob = new Blob([text], { type: 'text/plain' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  link.textContent = `Download ${name}`;
  link.addEventListener('click', () => {
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  return link;
}

function createDownloadLinkFromBlob(name, blob, label = null) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  link.textContent = label || `Download ${name}`;
  link.addEventListener('click', () => {
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  return link;
}

function resolveDirectAssignmentRunState(theorySpectrum, expSpectrum) {
  const useInputAssignments = isInputAssignmentModeEnabled();
  const manualMode = isManualAssignmentMode();
  if (manualMode) {
    syncFitLabelsFromManualAssignments();
  }
  if (useInputAssignments) {
    const assignments = buildInputAssignmentsFromSpectrum(theorySpectrum, expSpectrum);
    return {
      manualMode: true,
      useInputAssignments: true,
      directAssignments: assignments,
      selectedLabels: assignments.map((item) => item.label),
      modeTag: 'input-file',
    };
  }
  if (manualMode) {
    return {
      manualMode: true,
      useInputAssignments: false,
      directAssignments: manualAssignments.slice(),
      selectedLabels: manualAssignments.map((item) => item.label),
      modeTag: 'manual',
    };
  }
  return {
    manualMode: false,
    useInputAssignments: false,
    directAssignments: [],
    selectedLabels: selectedFitLabels.slice(),
    modeTag: 'auto',
  };
}

async function runAutofit() {
  lastAutofitSummary = null;
  void window.WMSPwaEnhancements?.prepareNotificationsFromGesture?.();
  if (!autofitInputText) {
    setAutofitStatus('Please upload the theoretical input file.');
    return;
  }
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  const expSpectrum = getSpectrumById(autofitExpSelect?.value);
  if (!theorySpectrum || !expSpectrum) {
    setAutofitStatus('Select both theoretical and experimental spectra.');
    return;
  }
  if (theorySpectrum.id === expSpectrum.id) {
    setAutofitStatus('Select an experimental spectrum different from the theoretical one.');
    return;
  }
  if (!theorySpectrum.labels || !theorySpectrum.labels.some((l) => l)) {
    setAutofitStatus('Theoretical spectrum has no labels. Use a labeled stick spectrum.');
    return;
  }
  let runState;
  try {
    runState = resolveDirectAssignmentRunState(theorySpectrum, expSpectrum);
  } catch (error) {
    setAutofitStatus(String(error?.message || error));
    return;
  }
  const { manualMode, useInputAssignments, directAssignments, selectedLabels, modeTag } = runState;
  if (manualMode && directAssignments.length < 3) {
    setAutofitStatus(useInputAssignments
      ? 'Use input file assignation requires at least 3 labeled stick lines.'
      : 'Manual mode requires at least 3 theory-experimental pairs.');
    return;
  }
  if (!manualMode && !selectedLabels.length) {
    setAutofitStatus('Select at least one fit label.');
    return;
  }
  const resolved = getAutofitOptionValues(manualMode, directAssignments);
  const runOptions = { ...resolved.values };
  const currentModelParams = getSpectrumModelParams(theorySpectrum);
  const quarticScaleState = getQuarticScaleState(theorySpectrum);
  const fitParams = getSelectedOptionalFitParams();
  let paramBounds = null;
  try {
    paramBounds = getCustomParamBounds();
  } catch (error) {
    setAutofitStatus(String(error?.message || error));
    return;
  }
  if (manualMode && resolved.notes.length) {
    showOverlayToast(resolved.notes.join(' · '));
  }
  lastAutofitContext = {
    manualMode,
    useInputAssignments,
    assignmentSource: modeTag,
    expMode: autofitExpMode?.value || 'profile',
    theorySpectrumId: theorySpectrum.id,
    theorySpectrumName: theorySpectrum.name,
    experimentalSpectrumId: expSpectrum.id,
    experimentalSpectrumName: expSpectrum.name,
    selectedFitLabels: selectedLabels.slice(),
    selectedControlLabels: selectedControlLabels.slice(),
    manualAssignments: cloneJsonValue(directAssignments),
    optionValues: cloneJsonValue(runOptions),
    currentModelParams: cloneJsonValue(currentModelParams),
    quarticScaleState: cloneJsonValue(quarticScaleState),
    fitParams: fitParams.slice(),
    paramBounds: cloneJsonValue(paramBounds),
  };

  setAutofitStatus(currentModelParams
    ? (useInputAssignments
      ? `Starting worker… using ${directAssignments.length} assignments read from the experimental input labels and the latest Hamiltonian parameters as initial guess.`
      : 'Starting worker… using the latest Hamiltonian parameters as initial guess.')
    : (useInputAssignments
      ? `Starting worker… using ${directAssignments.length} assignments read from the experimental input labels.`
      : 'Starting worker…'));
  if (autofitOutput) autofitOutput.textContent = '';
  resetAutofitProgress();
  if (autofitRunButton) autofitRunButton.disabled = true;
  if (autofitCancelButton) autofitCancelButton.disabled = false;
  updateRestoreInitialButtonState();
  setAutofitWarningVisible(true);

  try {
    let expText = expSpectrum.rawText || '';
    let expMode = autofitExpMode?.value || 'profile';
    const labelsText = formatAutofitLabelText(selectedLabels, theorySpectrum);
    const controlLabelsText = manualMode ? '' : formatAutofitLabelText(selectedControlLabels, theorySpectrum);
    if (expMode === 'stick') {
      if (expSpectrum.mode === 'profile' && expSpectrum.stick) {
        expText = formatStickDat(expSpectrum.stick.freqs || [], expSpectrum.stick.baseInts || []);
      }
    }

    const worker = initWorker();
    worker.postMessage({
      type: 'run',
      payload: {
        inputText: autofitInputText,
        labelsText,
        controlLabelsText,
        manualAssignmentsText: manualMode && !useInputAssignments ? formatDirectAssignmentsText(directAssignments) : '',
        useInputAssignments,
        expText,
        expMode,
        options: runOptions,
        currentModelParams,
        quarticScaleState,
        fitParams,
        paramBounds,
      }
    });
    autofitRunInProgress = true;
    updateRestoreInitialButtonState();
    void window.WMSPwaEnhancements?.notifyBackground?.({
      title: 'WMS-FitRot',
      body: 'Autofit started.',
      tag: 'wms-fitrot-autofit'
    });
  } catch (error) {
    autofitRunInProgress = false;
    console.error(error);
    setAutofitStatus('Autofit failed to start.');
    setAutofitWarningVisible(false);
    void window.WMSPwaEnhancements?.notifyBackground?.({
      title: 'WMS-FitRot',
      body: 'Autofit failed to start.',
      tag: 'wms-fitrot-autofit'
    });
    if (autofitOutput) autofitOutput.textContent = String(error);
    if (autofitRunButton) autofitRunButton.disabled = false;
    if (autofitCancelButton) autofitCancelButton.disabled = true;
    updateRestoreInitialButtonState();
  }
}

async function runAutofitLocal() {
  if (!autofitInputText) {
    setAutofitStatus('Please upload the theoretical input file.');
    return;
  }
  const expSpectrum = getSpectrumById(autofitExpSelect?.value);
  if (!expSpectrum) {
    setAutofitStatus('Select an experimental spectrum.');
    return;
  }
  const theorySpectrum = getSpectrumById(autofitTheorySelect?.value);
  if (theorySpectrum && theorySpectrum.id === expSpectrum.id) {
    setAutofitStatus('Select an experimental spectrum different from the theoretical one.');
    return;
  }
  let runState;
  try {
    runState = resolveDirectAssignmentRunState(theorySpectrum, expSpectrum);
  } catch (error) {
    setAutofitStatus(String(error?.message || error));
    return;
  }
  const { manualMode, useInputAssignments, directAssignments, selectedLabels } = runState;
  if (manualMode && directAssignments.length < 3) {
    setAutofitStatus(useInputAssignments
      ? 'Use input file assignation requires at least 3 labeled stick lines.'
      : 'Manual mode requires at least 3 theory-experimental pairs.');
    return;
  }
  if (!manualMode && !selectedLabels.length) {
    setAutofitStatus('Select at least one fit label.');
    return;
  }

  const expMode = autofitExpMode?.value || 'profile';
  let expText = expSpectrum.rawText || '';
  const labelsText = formatAutofitLabelText(selectedLabels, theorySpectrum);
  const controlLabelsText = manualMode ? '' : formatAutofitLabelText(selectedControlLabels, theorySpectrum);
  if (expMode === 'stick') {
    if (expSpectrum.mode === 'profile' && expSpectrum.stick) {
      expText = formatStickDat(expSpectrum.stick.freqs || [], expSpectrum.stick.baseInts || []);
    }
  }
  const resolved = getAutofitOptionValues(manualMode, directAssignments);
  const runOptions = { ...resolved.values };
  const currentModelParams = getSpectrumModelParams(theorySpectrum);
  const quarticScaleState = getQuarticScaleState(theorySpectrum);
  const fitParams = getSelectedOptionalFitParams();
  let paramBounds = null;
  try {
    paramBounds = getCustomParamBounds();
  } catch (error) {
    setAutofitStatus(String(error?.message || error));
    return;
  }
  if (manualMode && resolved.notes.length) {
    showOverlayToast(resolved.notes.join(' · '));
  }

  resetAutofitProgress();
  setAutofitWarningVisible(false);
  if (autofitDownloads) {
    const zip = new JSZip();
    zip.file('input.txt', autofitInputText);
    zip.file('labels.txt', labelsText);
    if (manualMode && !useInputAssignments) {
      zip.file('manual_assignments.txt', formatDirectAssignmentsText(directAssignments));
    } else if (selectedControlLabels.length) {
      zip.file('control_labels.txt', controlLabelsText);
    }
    zip.file('exp.dat', expText);

    const cmd = buildCliArgs(expMode, {
      hasControlLabels: selectedControlLabels.length > 0,
      manualMode,
      useInputAssignments,
      optionValues: runOptions,
      currentModelParams,
      quarticScaleState,
      fitParams,
      customParamBounds: paramBounds,
    });
    zip.file('README.txt', `WMS-FitRot local run\n\n1) Unzip this archive\n2) Run:\n${cmd}\n`);

    const scriptFiles = [
      { path: './fitting_script.txt', name: 'fitting_script.py' },
      { path: './fit_abc.txt', name: 'fit_abc.py' },
      { path: '../vmsrot.txt', name: 'vmsrot.py' },
    ];

    await Promise.all(scriptFiles.map(async (file) => {
      try {
        const response = await fetch(file.path);
        if (!response.ok) throw new Error(`Failed to fetch ${file.path}`);
        const text = await response.text();
        zip.file(file.name, text);
      } catch (error) {
        console.error(error);
      }
    }));

    const zipBlob = await zip.generateAsync({ type: 'blob' });
    autofitDownloads.appendChild(createDownloadLinkFromBlob('fitspectra_local_run.zip', zipBlob, 'Download local run ZIP'));

    const cmdEl = document.createElement('div');
    cmdEl.textContent = `Run locally: ${cmd}`;
    autofitDownloads.appendChild(cmdEl);
  }
  setAutofitStatus('Prepared local run ZIP.');
}

autofitAddLabel?.addEventListener('click', () => {
  if (!autofitLabelInput) return;
  addLabel(autofitLabelInput.value, getClickModeValue());
  autofitLabelInput.value = '';
});

autofitLabelInput?.addEventListener('keydown', (event) => {
  if (event.key === 'Enter') {
    event.preventDefault();
    addLabel(autofitLabelInput.value, getClickModeValue());
    autofitLabelInput.value = '';
  }
});

autofitClearLabels?.addEventListener('click', () => {
  if (autofitUseInputAssignment?.checked) {
    autofitUseInputAssignment.checked = false;
  }
  selectedFitLabels.length = 0;
  selectedControlLabels.length = 0;
  manualAssignments.length = 0;
  pendingManualSelection = null;
  overlayHighlightsFit.clear();
  overlayHighlightsControl.clear();
  renderLabelLists();
  updateOverlayPlot();
});

autofitAssignmentMode?.addEventListener('change', () => {
  if (isManualAssignmentMode()) {
    syncFitLabelsFromManualAssignments();
  } else {
    pendingManualSelection = null;
  }
  renderLabelLists();
  updateOverlayPlot();
});

autofitUseInputAssignment?.addEventListener('change', () => {
  if (autofitUseInputAssignment.checked && autofitExpMode) {
    autofitExpMode.value = 'stick';
  }
  pendingManualSelection = null;
  renderLabelLists();
  updateOverlayPlot();
});

autofitExpMode?.addEventListener('change', () => {
  renderLabelLists();
  updateOverlayPlot();
});

overlayTargetFit?.addEventListener('click', () => {
  setClickModeValue('fit');
  updateOverlayInteractionUi();
});

overlayTargetControl?.addEventListener('click', () => {
  setClickModeValue('control');
  updateOverlayInteractionUi();
});

autofitClickMode?.addEventListener('change', () => {
  updateOverlayInteractionUi();
});

autofitOptions.lsFScaleMode?.addEventListener('change', () => {
  syncLeastSquaresUiState();
});
autofitOptions.span?.addEventListener('input', () => {
  updateAutofitParamBoundPlaceholders();
});
autofitOptions.finalSpan?.addEventListener('input', () => {
  updateAutofitParamBoundPlaceholders();
});

const onAutofitSpectrumSelectionChanged = () => {
  syncQuarticScaleUiFromSpectrum();
  updateAutofitParamBoundPlaceholders();
  updateRestoreInitialButtonState();
  if (isInputAssignmentModeEnabled()) {
    renderLabelLists();
    updateOverlayPlot();
    return;
  }
  if (!isManualAssignmentMode()) return;
  if (!manualAssignments.length && !pendingManualSelection) return;
  clearManualAssignments({ skipRender: true, clearFitLabels: true });
  renderLabelLists();
  updateOverlayPlot();
  showOverlayToast('Manual pairs reset after changing selected spectra');
};

autofitTheorySelect?.addEventListener('change', onAutofitSpectrumSelectionChanged);
autofitExpSelect?.addEventListener('change', onAutofitSpectrumSelectionChanged);
autofitFitQuarticScale?.addEventListener('change', () => {
  syncQuarticFitSelectionUi();
});
QUARTIC_PARAM_KEYS.forEach((key) => {
  autofitFitParamInputs[key]?.addEventListener('change', () => {
    syncQuarticFitSelectionUi(key);
  });
});
autofitQuarticScaleSlider?.addEventListener('input', () => {
  setQuarticScaleUiValue(autofitQuarticScaleSlider.value);
  if (autofitQuarticScaleValue) {
    autofitQuarticScaleValue.value = autofitQuarticScaleSlider.value;
  }
  commitQuarticScaleFromUi();
});
autofitQuarticScaleValue?.addEventListener('input', () => {
  setQuarticScaleUiValue(autofitQuarticScaleValue.value);
});
autofitQuarticScaleValue?.addEventListener('change', () => {
  commitQuarticScaleFromUi();
});
sensitivityHideWeakToggle?.addEventListener('change', () => {
  autofitSensitivityState.hideWeak = Boolean(sensitivityHideWeakToggle.checked);
  ensureSensitivitySelection(lastAutofitFitMetrics);
  renderParameterTable();
  renderSensitivityChart();
  renderCorrelationHeatmap();
  renderParameterDetail();
  syncSensitivityHighlights();
});
sensitivityPanel?.addEventListener('toggle', () => {
  if (!sensitivityPanel.open) return;
  renderSensitivityDashboard(lastAutofitFitMetrics);
});

autofitRestoreInitialButton?.addEventListener('click', () => {
  restoreInitialTheoryParameters();
});
autofitRunButton?.addEventListener('click', () => runAutofit());
autofitRunLocalButton?.addEventListener('click', () => runAutofitLocal());
autofitDownloadLogButton?.addEventListener('click', () => {
  downloadFitStatusLog();
});
autofitUploadLogButton?.addEventListener('click', () => {
  autofitUploadLogInput?.click();
});
autofitUploadLogInput?.addEventListener('change', async (event) => {
  const file = event.target.files?.[0];
  if (!file) return;
  try {
    await handleFitStatusLogUpload(file);
  } catch (error) {
    console.error('Unable to restore Fit Status LOG:', error);
    alert(String(error?.message || error));
  } finally {
    autofitUploadLogInput.value = '';
  }
});
autofitImportSpfitButton?.addEventListener('click', () => autofitImportInput?.click());
autofitImportInput?.addEventListener('change', async (event) => {
  const files = Array.from(event.target.files || []);
  if (!files.length) return;
  try {
    await handleSpfitImportFiles(files);
  } catch (error) {
    console.error('Unable to import SPFIT files:', error);
    pendingSpfitImportState = null;
    setAutofitProgress(1);
    setAutofitStatus(`SPFIT import failed: ${String(error?.message || error)}`);
  } finally {
    autofitImportInput.value = '';
  }
});
autofitCancelButton?.addEventListener('click', () => {
  const wasRunning = autofitRunInProgress;
  if (autofitWorker) {
    autofitWorker.terminate();
    autofitWorker = null;
  }
  autofitRunInProgress = false;
  pendingSpfitImportState = null;
  renderAutofitDiagnostics(null);
  renderSpfitImportReport(null);
  setAutofitStatus('Cancelled.');
  setAutofitWarningVisible(false);
  if (wasRunning) {
    void window.WMSPwaEnhancements?.notifyBackground?.({
      title: 'WMS-FitRot',
      body: 'Autofit cancelled.',
      tag: 'wms-fitrot-autofit'
    });
  }
  if (autofitRunButton) autofitRunButton.disabled = false;
  if (autofitCancelButton) autofitCancelButton.disabled = true;
  updateRestoreInitialButtonState();
});

renderAutofitParamBounds();
syncQuarticScaleUiFromSpectrum();
syncQuarticFitSelectionUi();
renderLabelLists();
syncLeastSquaresUiState();
updateRestoreInitialButtonState();

function addSpectrumFromText({
  id = null,
  name,
  text,
  mode = 'profile',
  setActive = false,
  role = null,
  previewView = null,
  overlayState = null,
  modelParams = null,
  fittedConstants = null,
  quarticScale = null,
  quarticReferenceParams = null,
  initialState = null,
  skipRender = false,
} = {}) {
  const parsed = parseDatSpectrum(text || '');
  if (!parsed.freqs.length) {
    alert(`The file \"${name || 'selected'}\" does not contain valid data.`);
    return null;
  }
  const baseAbs = parsed.intensities.map((value) => Math.abs(Number(value)));
  const spectrumId = id || `spectrum-${++spectrumCounter}`;
  syncSpectrumCounterWithId(spectrumId);
  const spectrum = {
    id: spectrumId,
    name: name || `Spectrum ${spectrumCounter}`,
    mode: mode === 'stick' ? 'stick' : 'profile',
    freqs: parsed.freqs,
    baseAbs,
    labels: parsed.labels,
    rawText: text || '',
    isTheoretical: role === 'theory',
    stick: null,
    previewView: mode === 'profile'
      ? (previewView === 'stick' ? 'stick' : 'profile')
      : 'stick',
    modelParams: normalizeModelParams(modelParams),
    fittedConstants: fittedConstants ? cloneJsonValue(fittedConstants) : null,
    quarticScale: normalizeQuarticScale(quarticScale, 1),
    quarticReferenceParams: normalizeQuarticReferenceParams(quarticReferenceParams),
    initialState: null,
  };
  if (spectrum.mode === 'profile') {
    const stickData = computeDatStickSpectrum(parsed.freqs, parsed.intensities, null, 'absorbance');
    spectrum.stick = {
      freqs: stickData.freqs,
      baseInts: stickData.baseInts,
      labels: [],
    };
  } else {
    spectrum.stick = { freqs: parsed.freqs, baseInts: baseAbs, labels: parsed.labels };
  }
  spectrum.initialState = normalizeSerializedSpectrumState(initialState)
    || (spectrum.isTheoretical ? buildSpectrumInitialStateSnapshot(spectrum) : null);
  spectra.set(spectrum.id, spectrum);
  if (overlayState) {
    applyStoredOverlaySettings(spectrum, overlayState);
  }
  if (setActive) {
    activeSpectrumId = spectrum.id;
  }
  if (!skipRender) {
    renderUploadedTabs();
    renderOverlayList();
    updateOverlayPlot();
    syncInitialTheoryStateFromInputText();
  }
  return spectrum;
}

function ingestQueuedSpectra() {
  const queueRaw = localStorage.getItem(queueKey);
  if (!queueRaw) return;
  let queue = [];
  try {
    queue = JSON.parse(queueRaw) || [];
  } catch (error) {
    console.warn('Unable to parse queued spectra:', error);
  }
  if (!Array.isArray(queue) || !queue.length) {
    localStorage.removeItem(queueKey);
    return;
  }
  queue.forEach((item) => {
    if (item?.mode === 'input') {
      if (typeof item?.text === 'string' && item.text.trim()) {
        autofitInputText = item.text;
        autofitInputName = item?.name || null;
        setAutofitStatus(`Loaded input from WMS-Rot: ${item?.name || 'input file'}`);
      }
      return;
    }
    addSpectrumFromText({
      name: item?.name || 'WMS-Rot spectrum',
      text: item?.text || '',
      mode: item?.mode || 'stick',
      role: item?.role || null,
      setActive: true,
    });
  });
  localStorage.removeItem(queueKey);
  syncAutofitJMaxFromInputText();
  syncInitialTheoryStateFromInputText();
  updateAutofitParamBoundPlaceholders();
}

ingestQueuedSpectra();
