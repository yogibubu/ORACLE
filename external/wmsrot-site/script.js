// script_tmp.js (full-F mixing only for quadrupole; bugfix pyodide out-of-bounds)
// Note: the non-HFS physics/flow follows the original code (same intensity style).
// Note: Qrs is intentionally kept fixed at 22.8174 as requested.

const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.23.4/full/";
const PYODIDE_PACKAGES = ["numpy", "pandas", "matplotlib", "sympy"];

function createPyodideReadyPromise() {
  return loadPyodide({ indexURL: PYODIDE_INDEX_URL });
}

let pyodideReadyPromise = createPyodideReadyPromise();

async function clearPyodideCaches() {
  if (typeof caches === "undefined" || typeof caches.keys !== "function") return;
  try {
    const cacheNames = await caches.keys();
    const toDelete = cacheNames.filter(name => name.startsWith("wms-rot-pyodide-"));
    await Promise.all(toDelete.map(name => caches.delete(name)));
  } catch {
    // no-op
  }
}

const PLOTLY_DOWNLOAD_OPTIONS = { format: "png", scale: 5 };

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
      scale: Number.isFinite(existingScale) ? Math.max(existingScale, PLOTLY_DOWNLOAD_OPTIONS.scale) : PLOTLY_DOWNLOAD_OPTIONS.scale
    }
  };
}

function sanitizePlotFilename(value, fallback = "plotly_figure") {
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
    filename: sanitizePlotFilename(filenameBase, "wms_rot_plot"),
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

const globalLoadingOverlay = (() => {
  const overlay = document.getElementById("loading-overlay");
  if (!overlay) {
    return {
      show() {},
      hide() {},
      setMessage() {},
      showError() {},
    };
  }
  const messageEl = overlay.querySelector('[data-role="loading-message"]');
  const spinnerEl = overlay.querySelector('[data-role="loading-spinner"]');
  return {
    show(message) {
      overlay.classList.remove("loading-overlay--error");
      overlay.classList.add("is-visible");
      overlay.removeAttribute("aria-hidden");
      if (spinnerEl) spinnerEl.hidden = false;
      if (message && messageEl) messageEl.textContent = message;
    },
    hide() {
      overlay.classList.remove("is-visible");
      overlay.setAttribute("aria-hidden", "true");
    },
    setMessage(message) {
      if (message && messageEl) messageEl.textContent = message;
    },
    showError(message) {
      overlay.classList.add("loading-overlay--error");
      overlay.classList.add("is-visible");
      overlay.removeAttribute("aria-hidden");
      if (spinnerEl) spinnerEl.hidden = true;
      if (messageEl) messageEl.textContent = message || "Error during loading.";
    },
  };
})();

async function main() {
  let pyodide;
  try {
    globalLoadingOverlay.show("Loading Pyodide…");
    pyodide = await pyodideReadyPromise;
    globalLoadingOverlay.setMessage("Installing scientific packages…");
    await pyodide.loadPackage(PYODIDE_PACKAGES);

    globalLoadingOverlay.setMessage("Preparing calculation environment…");
    const setupCode = await fetch("vmsrot.txt").then(response => response.text());
    pyodide.runPython(setupCode);
  } catch (err) {
    console.warn("Initial Pyodide bootstrap failed, retrying once after cache cleanup.", err);
    try {
      globalLoadingOverlay.setMessage("Recovering Pyodide cache and retrying…");
      await clearPyodideCaches();
      pyodideReadyPromise = createPyodideReadyPromise();
      pyodide = await pyodideReadyPromise;
      globalLoadingOverlay.setMessage("Installing scientific packages…");
      await pyodide.loadPackage(PYODIDE_PACKAGES);
      globalLoadingOverlay.setMessage("Preparing calculation environment…");
      const setupCode = await fetch("vmsrot.txt").then(response => response.text());
      pyodide.runPython(setupCode);
    } catch (retryErr) {
      console.error(retryErr);
      globalLoadingOverlay.showError("Failed to load Pyodide/packages. Please reload.");
      return;
    }
  }

  const get = id => parseFloat(document.getElementById(id).value);
  const getText = id => document.getElementById(id).value;
  const getChecked = id => document.getElementById(id).checked;
  const DISTORTION_FIELD_IDS = [
    "DJ", "DJK", "DK", "dJ", "dK",
    "HJ", "HJK", "HKJ", "HK", "h1", "h2", "h3"
  ];
  const QUARTIC_DISTORTION_FIELD_IDS = ["DJ", "DJK", "DK", "dJ", "dK"];
  const SEXTIC_DISTORTION_FIELD_IDS = ["HJ", "HJK", "HKJ", "HK", "h1", "h2", "h3"];
  const REDUCTION_KEYS = ["A", "S"];
  const REPRESENTATION_UI_INFO = {
    Ir: {
      label: "Ir (almost-prolate)",
      axes: { z: "a", x: "b", y: "c" }
    },
    IIIl: {
      label: "IIIl (almost-oblate)",
      axes: { z: "c", x: "b", y: "a" }
    }
  };
  let lastPyvars = null;
  let lastSpectrum = null;
  let loadedFileData = null;
  let representationChangedWarningActive = false;
  // Distortion constants are stored independently for Watson A and Watson S.
  // Switching reduction swaps the active set, but does not transform values.
  let activeDistortionReduction = String(getText("reduction") || "A").trim().toUpperCase().startsWith("S") ? "S" : "A";
  let distortionValuesByReductionMHz = { A: null, S: null };
  const simulatedSpectrumSection = document.getElementById("simulated-spectrum-section");
  const reductionSelect = document.getElementById("reduction");
  const representationSelect = document.getElementById("representation");
  const representationStatus = document.getElementById("representation-status");
  const representationHelp = document.getElementById("representation-help");
  const reductionHelp = document.getElementById("reduction-help");
  const distortionRepresentationWarning = document.getElementById("representation-distortion-warning");
  const spectrumDownloadButton = document.getElementById("download-spectrum-png");
  const downloadResultsButton = document.getElementById("download-results");
  const openCompareButton = document.getElementById("open-compare");
  const openWavefunctionViewerButton = document.getElementById("open-wavefunction-viewer");
  const postProcessingStatus = document.getElementById("post-processing-status");

  function setSimulatedSpectrumVisible(visible) {
    if (!simulatedSpectrumSection) return;
    simulatedSpectrumSection.hidden = !visible;
  }

  function setPostProcessingReady(ready) {
    [downloadResultsButton, openWavefunctionViewerButton].forEach((button) => {
      if (button) button.disabled = !ready;
    });
    if (openCompareButton) {
      openCompareButton.disabled = false;
    }
    if (!postProcessingStatus) return;
    postProcessingStatus.textContent = ready
      ? "Current calculated spectrum is ready for fitting, wavefunction inspection, and export."
      : "WMS-FitRot can be opened directly or restored from a Fit Status LOG. Calculate a spectrum in Core Calculation to enable wavefunction inspection and spectrum export.";
    postProcessingStatus.classList.toggle("is-ready", ready);
  }

  setSimulatedSpectrumVisible(false);
  setPostProcessingReady(Boolean(lastPyvars));

  // Keep JS plot ranges consistent with the Python-side (vmsrot.py) "auto" unit
  // heuristic: if FREQ_MAX <= 1000 we interpret the window as GHz, else MHz.
  function normalizeFreqWindowMHz(minRaw, maxRaw) {
    let min = Number(minRaw);
    let max = Number(maxRaw);
    if (Number.isNaN(min) || Number.isNaN(max)) return { minMHz: minRaw, maxMHz: maxRaw };
    if (max > 0 && max <= 1000) {
      max *= 1000;
      if (min > 0) min *= 1000;
    }
    // Mirror vmsrot.py clamping to avoid inverted / degenerate ranges.
    if (min < 5e-5) min = 5e-5;
    if (max < min) max = min;
    return { minMHz: min, maxMHz: max };
  }

  function estimateFirstLinearLineMHz() {
    const B = get("B");
    const DJ = readDistortionFieldMHz("DJ");
    const HJ = readDistortionFieldMHz("HJ");
    if (Number.isNaN(B)) return null;
    // R(0): nu ~= 2B - 4DJ + 8HJ (constants in MHz).
    return 2 * B - 4 * DJ + 8 * HJ;
  }

  function validateRotationalConstantsOrThrow(rotorType, A, B, C) {
    // Small tolerance to avoid false positives from rounding / copy-paste.
    const eps = 1e-9;
    const finite = [A, B, C].every(v => Number.isFinite(v));
    if (!finite) return; // Let existing input validation handle NaNs.

    if (rotorType === "asymmetric") {
      if (!(A > B + eps && B > C + eps)) {
        throw new Error(
          `Asymmetric rotor requires A > B > C (MHz). You entered A=${A}, B=${B}, C=${C}.`
        );
      }
      return;
    }

    // For non-asymmetric rotors, allow equalities (symmetric/linear/spherical),
    // but keep the user convention "largest first": A >= B >= C.
    if (!(A + eps >= B && B + eps >= C)) {
      throw new Error(
        `Please enter rotational constants ordered as A >= B >= C (MHz). You entered A=${A}, B=${B}, C=${C}.`
      );
    }
  }

  const buildHoverTemplate = (label, includeHyperfine = false, includeI12 = false) => {
    const lower = includeHyperfine
      ? (includeI12
        ? "Lower (J,F,I12,Ka,Kc): %{customdata[2]:.0f}, %{customdata[8]}, %{customdata[10]}, %{customdata[3]:.0f}, %{customdata[4]:.0f}<br>"
        : "Lower (J,F,Ka,Kc): %{customdata[2]:.0f}, %{customdata[8]}, %{customdata[3]:.0f}, %{customdata[4]:.0f}<br>")
      : "Lower (J,Ka,Kc): %{customdata[2]:.0f}, %{customdata[3]:.0f}, %{customdata[4]:.0f}<br>";
    const upper = includeHyperfine
      ? (includeI12
        ? "Upper (J,F,I12,Ka,Kc): %{customdata[5]:.0f}, %{customdata[9]}, %{customdata[11]}, %{customdata[6]:.0f}, %{customdata[7]:.0f}<br>"
        : "Upper (J,F,Ka,Kc): %{customdata[5]:.0f}, %{customdata[9]}, %{customdata[6]:.0f}, %{customdata[7]:.0f}<br>")
      : "Upper (J,Ka,Kc): %{customdata[5]:.0f}, %{customdata[6]:.0f}, %{customdata[7]:.0f}<br>";
    return (
      "Frequency: %{customdata[0]:.6f} MHz<br>" +
      `${label}: %{customdata[1]:.3e}<br>` +
      lower +
      upper +
      `Branch: %{customdata[${includeI12 ? 12 : 10}]}<extra></extra>`
    );
  };

  function buildStickTrace(points, color, name, yKey, hoverTemplate) {
    if (!points.length) return null;
    const stickX = [];
    const stickY = [];
    const customdata = [];
    for (const p of points) {
      const yVal = p[yKey];
      const Fl = Number.isFinite(Number(p.Fl)) ? Number(p.Fl) : "";
      const Fu = Number.isFinite(Number(p.Fu)) ? Number(p.Fu) : "";
      const I12l = Number.isFinite(Number(p.I12_l)) ? Number(p.I12_l) : "";
      const I12u = Number.isFinite(Number(p.I12_u)) ? Number(p.I12_u) : "";
      const cd = [p.x, yVal, p.Jl, p.Ka_l, p.Kc_l, p.Ju, p.Ka_u, p.Kc_u, Fl, Fu, I12l, I12u, p.branch];
      stickX.push(p.x, p.x, null);
      stickY.push(0, yVal, null);
      customdata.push(cd, cd, null);
    }
    return {
      x: stickX,
      y: stickY,
      customdata,
      mode: "lines",
      line: { color, width: 1 },
      hovertemplate: hoverTemplate,
      name
    };
  }

  function renderSpectrum(spectrum) {
    const data = [];
    const showP = getChecked("showP");
    const showQ = getChecked("showQ");
    const showR = getChecked("showR");
    const colorP = getText("colorP");
    const colorQ = getText("colorQ");
    const colorR = getText("colorR");
    const scale = getText("intensityScale");
    const yKey = scale === "absolute" ? "yAbs" : "yRel";
    const hoverLabel = scale === "absolute"
      ? "Intensity (nm²·MHz)"
      : "Rel. intensity";
    const hasHyperfineLabels = points =>
      points.some(p => Number.isFinite(Number(p.Fl)) && Number.isFinite(Number(p.Fu)));
    const hasI12Labels = points =>
      points.some(p => Number.isFinite(Number(p.I12_l)) && Number.isFinite(Number(p.I12_u)));

    if (showP) {
      const traceP = buildStickTrace(
        spectrum.branches.P,
        colorP,
        "P",
        yKey,
        buildHoverTemplate(hoverLabel, hasHyperfineLabels(spectrum.branches.P), hasI12Labels(spectrum.branches.P))
      );
      if (traceP) data.push(traceP);
    }
    if (showQ) {
      const traceQ = buildStickTrace(
        spectrum.branches.Q,
        colorQ,
        "Q",
        yKey,
        buildHoverTemplate(hoverLabel, hasHyperfineLabels(spectrum.branches.Q), hasI12Labels(spectrum.branches.Q))
      );
      if (traceQ) data.push(traceQ);
    }
    if (showR) {
      const traceR = buildStickTrace(
        spectrum.branches.R,
        colorR,
        "R",
        yKey,
        buildHoverTemplate(hoverLabel, hasHyperfineLabels(spectrum.branches.R), hasI12Labels(spectrum.branches.R))
      );
      if (traceR) data.push(traceR);
    }

    const xMinRaw = get("freq_min_spcat");
    const xMaxRaw = get("freq_max_spcat");
    const { minMHz: xMin, maxMHz: xMax } = normalizeFreqWindowMHz(xMinRaw, xMaxRaw);
    const layout = {
      title: "",
      xaxis: { title: "Frequency (MHz)" },
      yaxis: {
        title: scale === "absolute"
          ? "Integrated intensity (nm²·MHz)"
          : "Rel. intensity"
      },
      margin: { t: 20, r: 20, b: 50, l: 60 },
      annotations: [],
      showlegend: data.length > 1
    };
    if (!Number.isNaN(xMin) && !Number.isNaN(xMax)) {
      layout.xaxis.range = [xMin, xMax];
    }

    if (!data.length) {
      let msg = spectrum.totalLines
        ? "No lines for selected branches"
        : "No lines above threshold";

      // Extra hint for the most common "linear rotor looks empty" case:
      // frequencies outside the selected window.
      if (!spectrum.totalLines) {
        const rotorType = document.getElementById("rotorType")?.value;
        if (rotorType === "linear" && !Number.isNaN(xMax)) {
          const first = estimateFirstLinearLineMHz();
          if (first !== null && first > xMax) {
            msg = `No lines in frequency window (first line ~${first.toFixed(1)} MHz). Increase FREQ_MAX.`;
          }
        }
      }
      layout.annotations.push({
        text: msg,
        x: 0.5,
        y: 0.5,
        xref: "paper",
        yref: "paper",
        showarrow: false
      });
    }

    setSpectrumDownloadReady(Boolean(spectrum));
    return Plotly.react(
      document.getElementById("spectrum"),
      data,
      layout,
      withHighResDownloadConfig({ responsive: true, displaylogo: false })
    );
  }

  const setValueIfPresent = (id, value) => {
    const el = document.getElementById(id);
    if (!el) return;
    if (value === null || value === undefined || Number.isNaN(value)) return;
    el.value = value;
  };

  const finiteOr = (value, fallback) => (Number.isFinite(value) ? value : fallback);
  const readIntegerField = (id) => {
    const el = document.getElementById(id);
    if (!el) return Number.NaN;
    const raw = String(el.value ?? "").trim();
    if (!raw) return Number.NaN;
    const value = Number(raw);
    return Number.isInteger(value) ? value : Number.NaN;
  };

  function setJRangeStatus(message = "", isError = false) {
    const target = document.getElementById("j-range-status");
    if (!target) return;
    target.textContent = message;
    target.style.color = isError ? "#b91c1c" : "#555";
  }

  function updateJRangeUiState() {
    const autoEstimate = getChecked("auto_estimate_j_range");
    const jMinField = document.getElementById("J_min");
    const jMaxField = document.getElementById("J_max");
    if (jMinField) jMinField.disabled = autoEstimate;
    if (jMaxField) jMaxField.disabled = autoEstimate;
  }

  function collectSimulationBase() {
    const rep = normalizeRepresentationKey();
    const A = get("A");
    const B = get("B");
    const C = get("C");
    const dvibA = Number.isNaN(get("DVibA")) ? 0 : get("DVibA");
    const dvibB = Number.isNaN(get("DVibB")) ? 0 : get("DVibB");
    const dvibC = Number.isNaN(get("DVibC")) ? 0 : get("DVibC");
    const muA = get("mu_a");
    const muB = get("mu_b");
    const muC = get("mu_c");
    const qrotOverride = get("qrot_override");
    const sD2Id = loadedFileData && Number.isFinite(loadedFileData.sD2Id) ? loadedFileData.sD2Id : NaN;
    const Fmax = get("F_max");
    const I_nuc2 = get("I_nuc_2");
    const chi_aa_2 = get("chi_aa_2");
    const chi_bb_2 = get("chi_bb_2");
    const chi_cc_2 = get("chi_cc_2");
    const chi_ab_2 = get("chi_ab_2");
    const chi_ac_2 = get("chi_ac_2");
    const chi_bc_2 = get("chi_bc_2");
    const rotorType = document.getElementById("rotorType").value;

    const base = {
      A: A + dvibA,
      B: B + dvibB,
      C: C + dvibC,
      DJ: readDistortionFieldMHz("DJ"),
      DJK: readDistortionFieldMHz("DJK"),
      DK: readDistortionFieldMHz("DK"),
      dJ: readDistortionFieldMHz("dJ"),
      dK: readDistortionFieldMHz("dK"),
      HJ: readDistortionFieldMHz("HJ"),
      HJK: readDistortionFieldMHz("HJK"),
      HKJ: readDistortionFieldMHz("HKJ"),
      HK: readDistortionFieldMHz("HK"),
      h1: readDistortionFieldMHz("h1"),
      h2: readDistortionFieldMHz("h2"),
      h3: readDistortionFieldMHz("h3"),
      mu_a: muA,
      mu_b: muB,
      mu_c: muC,
      T: get("T"),
      intensity_cut: 1e-20,
      groupSymmetry: getText("groupSymmetry"),
      a_checked: document.getElementById("a_axis").checked,
      b_checked: document.getElementById("b_axis").checked,
      c_checked: document.getElementById("c_axis").checked,
      REDUCTION: getText("reduction"),
      REPRESENTATION: rep,
      QROT_OVERRIDE: Number.isNaN(qrotOverride) ? null : qrotOverride,
      S_D2_ID: Number.isNaN(sD2Id) ? null : Math.round(sD2Id),
      STR0: get("str0"),
      STR1: get("str1"),
      FREQ_MIN: get("freq_min_spcat"),
      FREQ_MAX: get("freq_max_spcat"),
      eeWt: get("eeWt"),
      eoWt: get("eoWt"),
      oeWt: get("oeWt"),
      ooWt: get("ooWt"),
      I_NUC: finiteOr(get("I_nuc"), 0.0),
      I_NUC_1: finiteOr(get("I_nuc"), 0.0),
      I_NUC_2: finiteOr(I_nuc2, 0.0),
      F_MAX: Number.isNaN(Fmax) ? null : Fmax,
      chi_aa: finiteOr(get("chi_aa"), 0.0),
      chi_bb: finiteOr(get("chi_bb"), 0.0),
      chi_cc: finiteOr(get("chi_cc"), 0.0),
      chi_ab: finiteOr(get("chi_ab"), 0.0),
      chi_ac: finiteOr(get("chi_ac"), 0.0),
      chi_bc: finiteOr(get("chi_bc"), 0.0),
      chi_aa_1: finiteOr(get("chi_aa"), 0.0),
      chi_bb_1: finiteOr(get("chi_bb"), 0.0),
      chi_cc_1: finiteOr(get("chi_cc"), 0.0),
      chi_ab_1: finiteOr(get("chi_ab"), 0.0),
      chi_ac_1: finiteOr(get("chi_ac"), 0.0),
      chi_bc_1: finiteOr(get("chi_bc"), 0.0),
      chi_aa_2: finiteOr(chi_aa_2, 0.0),
      chi_bb_2: finiteOr(chi_bb_2, 0.0),
      chi_cc_2: finiteOr(chi_cc_2, 0.0),
      chi_ab_2: finiteOr(chi_ab_2, 0.0),
      chi_ac_2: finiteOr(chi_ac_2, 0.0),
      chi_bc_2: finiteOr(chi_bc_2, 0.0),
      rotorType,
    };

    if (rotorType === "linear") {
      base.A = base.B;
      base.C = base.B;
      base.mu_b = 0;
      base.mu_c = 0;
    } else if (rotorType === "spherical") {
      base.A = base.B;
      base.C = base.B;
      base.mu_a = 0;
      base.mu_b = 0;
      base.mu_c = 0;
    } else if (rotorType === "prolate") {
      base.C = base.B;
      base.mu_b = 0;
      base.mu_c = 0;
    } else if (rotorType === "oblate") {
      base.A = base.B;
      base.mu_a = 0;
      base.mu_b = 0;
    }

    if (!base.a_checked) base.mu_a = 0;
    if (!base.b_checked) base.mu_b = 0;
    if (!base.c_checked) base.mu_c = 0;

    validateRotationalConstantsOrThrow(rotorType, base.A, base.B, base.C);
    return base;
  }

  function buildSimulationPyvars(base, overrides = null) {
    const opts = (overrides && typeof overrides === "object") ? overrides : {};
    const J_min = Number.isFinite(opts.J_min) ? opts.J_min : readIntegerField("J_min");
    const J_max = Number.isFinite(opts.J_max) ? opts.J_max : readIntegerField("J_max");

    if (!Number.isInteger(J_min) || J_min < 0) {
      throw new Error(`J_min must be a non-negative integer. You entered ${document.getElementById("J_min")?.value ?? ""}.`);
    }
    if (!Number.isInteger(J_max) || J_max < J_min) {
      throw new Error(`J_max must be an integer greater than or equal to J_min. You entered J_min=${J_min}, J_max=${document.getElementById("J_max")?.value ?? J_max}.`);
    }

    return {
      ...base,
      J_min,
      J_max,
    };
  }

  async function estimateJRangeFromWindow(base = null, { updateFields = true } = {}) {
    const simBase = base || collectSimulationBase();
    const currentJmax = readIntegerField("J_max");
    const scanStart = Number.isInteger(currentJmax) && currentJmax > 0 ? currentJmax : 12;
    const scanCap = Math.max(scanStart, 1200);
    const estimatePyvars = buildSimulationPyvars(simBase, { J_min: 0, J_max: scanStart });

    Object.entries(estimatePyvars).forEach(([k, v]) => pyodide.globals.set(k, v));
    pyodide.globals.set("ESTIMATE_SCAN_J_START", scanStart);
    pyodide.globals.set("ESTIMATE_SCAN_J_CAP", scanCap);

    try {
      const estimateJson = await pyodide.runPythonAsync(`
import json
_estimate = estimate_j_range_from_frequency_window(
    T, A, B, C,
    DJ, DJK, DK, dJ, dK,
    HJ, HJK, HKJ, HK, h1, h2, h3,
    mu_a, mu_b, mu_c,
    rotorType=rotorType,
    groupSymmetry=groupSymmetry,
    a_checked=a_checked, b_checked=b_checked, c_checked=c_checked,
    scan_J_start=ESTIMATE_SCAN_J_START,
    scan_J_cap=ESTIMATE_SCAN_J_CAP,
)
json.dumps(_estimate)
`);
      const estimate = JSON.parse(estimateJson);
      if (updateFields) {
        document.getElementById("J_min").value = String(estimate.J_min);
        document.getElementById("J_max").value = String(estimate.J_max);
      }
      const sourceLabel = Array.isArray(estimate.rotor_types) && estimate.rotor_types.length
        ? estimate.rotor_types.join(" + ")
        : "selected approximation";
      const suffix = estimate.truncated ? ` Scan capped at J=${estimate.scan_J_max}.` : "";
      setJRangeStatus(`Estimated J range: ${estimate.J_min} - ${estimate.J_max} (${sourceLabel}).${suffix}`);
      return estimate;
    } catch (err) {
      setJRangeStatus(err.message || "Unable to estimate J range.", true);
      throw err;
    }
  }

  function normalizeDistortionUnit(unitRaw) {
    const raw = String(unitRaw || "").trim().toLowerCase();
    if (raw === "mhz") return "MHz";
    return "kHz";
  }

  function normalizeReductionKey(reductionRaw = getText("reduction")) {
    return String(reductionRaw || "A").trim().toUpperCase().startsWith("S") ? "S" : "A";
  }

  function normalizeRepresentationKey(representationRaw = getText("representation")) {
    const raw = String(representationRaw || "").trim().replace(/\s+/g, "").toLowerCase();
    if (!raw) return "Ir";
    if (raw === "il" || raw === "iil" || raw === "iiil") return "IIIl";
    if (raw === "ir" || raw === "iir" || raw === "iiir") return "Ir";
    return raw.endsWith("l") ? "IIIl" : "Ir";
  }

  function getRepresentationUiInfo(representationRaw = getText("representation")) {
    return REPRESENTATION_UI_INFO[normalizeRepresentationKey(representationRaw)] || REPRESENTATION_UI_INFO.Ir;
  }

  function formatAxisRoleLetter(letter) {
    const safe = String(letter || "").trim().toLowerCase();
    if (!safe) return "";
    if (safe === "z") return '<span class="axis-role-badge__mapped-axis axis-role-badge__mapped-axis--z">z</span>';
    return `<span class="axis-role-badge__mapped-axis">${safe}</span>`;
  }

  function formatAxisRoleBadge(source, target) {
    const safeSource = String(source || "").trim().toLowerCase();
    const safeTarget = String(target || "").trim().toLowerCase();
    if (!safeSource || !safeTarget) return "";
    return `${safeSource} &rarr; ${formatAxisRoleLetter(safeTarget)}`;
  }

  function formatAxisPairBadge(sourcePair, targetPair) {
    const safeSource = String(sourcePair || "").trim().toLowerCase();
    const safeTarget = String(targetPair || "").trim().toLowerCase();
    if (!safeSource || !safeTarget || safeSource.length !== safeTarget.length) return "";
    const mapped = safeTarget
      .split("")
      .map((letter) => formatAxisRoleLetter(letter))
      .join("");
    return `${safeSource} &rarr; ${mapped}`;
  }

  function updateAxisRoleBadges() {
    const rotorType = document.getElementById("rotorType")?.value || "asymmetric";
    const isAsymmetric = rotorType === "asymmetric";
    const repInfo = getRepresentationUiInfo();
    const axisToRole = {
      [repInfo.axes.z]: "z",
      [repInfo.axes.x]: "x",
      [repInfo.axes.y]: "y"
    };
    document.querySelectorAll(".axis-role-badge[data-axis-role]").forEach((node) => {
      if (!(node instanceof HTMLElement)) return;
      const axis = String(node.dataset.axisRole || "").toLowerCase();
      const mappedRole = axisToRole[axis] || "";
      const html = isAsymmetric ? formatAxisRoleBadge(axis, mappedRole) : "";
      node.innerHTML = html;
      node.style.display = html ? "inline-flex" : "none";
    });
    document.querySelectorAll(".axis-role-badge[data-axis-pair]").forEach((node) => {
      if (!(node instanceof HTMLElement)) return;
      const pair = String(node.dataset.axisPair || "").toLowerCase();
      const mapped = pair
        .split("")
        .map((axis) => axisToRole[axis] || "")
        .join("");
      const html = isAsymmetric && mapped.length === pair.length ? formatAxisPairBadge(pair, mapped) : "";
      node.innerHTML = html;
      node.style.display = html ? "inline-flex" : "none";
    });
  }

  function setRepresentationCompatibilityWarningsVisible(visible) {
    const rotorType = document.getElementById("rotorType")?.value || "asymmetric";
    const shouldShow = Boolean(visible) && rotorType === "asymmetric";
    if (distortionRepresentationWarning) {
      distortionRepresentationWarning.hidden = !shouldShow;
    }
  }

  function cloneDistortionValues(values) {
    if (!values || typeof values !== "object") return null;
    const out = {};
    DISTORTION_FIELD_IDS.forEach((id) => {
      const value = values[id];
      out[id] = Number.isFinite(value) ? value : null;
    });
    return out;
  }

  function hasAnyDistortionValue(values) {
    return DISTORTION_FIELD_IDS.some((id) => Number.isFinite(values?.[id]));
  }

  function persistActiveReductionState() {
    // Keep exactly what the user typed for the currently visible reduction.
    const reduction = activeDistortionReduction || normalizeReductionKey();
    distortionValuesByReductionMHz[reduction] = captureDistortionValuesForStateMHz();
  }

  function getCurrentDistortionUnit() {
    const el = document.getElementById("distortion_unit");
    return normalizeDistortionUnit(el ? el.value : "kHz");
  }

  function distortionToMHzFactor(unitRaw = getCurrentDistortionUnit()) {
    return normalizeDistortionUnit(unitRaw) === "MHz" ? 1.0 : 0.001;
  }

  function distortionFromMHzFactor(unitRaw = getCurrentDistortionUnit()) {
    return 1.0 / distortionToMHzFactor(unitRaw);
  }

  function readDistortionFieldMHz(id) {
    const raw = get(id);
    if (Number.isNaN(raw)) return 0.0;
    return raw * distortionToMHzFactor();
  }

  function readDistortionFieldStateMHz(id) {
    const el = document.getElementById(id);
    if (!el) return null;
    const rawText = String(el.value ?? "").trim();
    if (!rawText) return null;
    const rawNumber = parseFloat(rawText);
    if (Number.isNaN(rawNumber)) return null;
    return rawNumber * distortionToMHzFactor();
  }

  function distortionValueFromMHz(mhzValue) {
    return Number.isFinite(mhzValue) ? mhzValue * distortionFromMHzFactor() : mhzValue;
  }

  function captureDistortionValuesForStateMHz() {
    const out = {};
    DISTORTION_FIELD_IDS.forEach((id) => {
      out[id] = readDistortionFieldStateMHz(id);
    });
    return out;
  }

  function captureDisplayedDistortionValues() {
    const out = {};
    DISTORTION_FIELD_IDS.forEach((id) => {
      const el = document.getElementById(id);
      if (!el) return;
      out[id] = String(el.value ?? "");
    });
    return out;
  }

  function restoreDisplayedDistortionValues(values) {
    DISTORTION_FIELD_IDS.forEach((id) => {
      const el = document.getElementById(id);
      if (!el || !values || !Object.prototype.hasOwnProperty.call(values, id)) return;
      el.value = values[id];
    });
  }

  function restoreDistortionValuesFromMHz(valuesMHz, { clearMissing = false } = {}) {
    DISTORTION_FIELD_IDS.forEach((id) => {
      const el = document.getElementById(id);
      if (!el) return;
      const value = valuesMHz ? valuesMHz[id] : null;
      if (!Number.isFinite(value)) {
        if (clearMissing) el.value = "";
        return;
      }
      el.value = String(distortionValueFromMHz(value));
    });
  }

  function updateRepresentationChangeControls() {
    const rotorType = document.getElementById("rotorType")?.value || "asymmetric";
    const isAsymmetric = rotorType === "asymmetric";
    const selectedReduction = normalizeReductionKey();
    const repInfo = getRepresentationUiInfo();

    updateAxisRoleBadges();

    if (representationStatus) {
      representationStatus.style.display = isAsymmetric ? "" : "none";
      if (isAsymmetric) {
        representationStatus.textContent =
          `${repInfo.label}: Jz follows axis ${repInfo.axes.z}, ` +
          `Jx axis ${repInfo.axes.x}, and Jy axis ${repInfo.axes.y}. ` +
          `Use this same convention for rotational, dipole, quadrupole, and distortion constants ` +
          `with Watson ${selectedReduction}.`;
      } else {
        representationStatus.textContent = "";
      }
    }
    if (representationHelp) representationHelp.style.display = isAsymmetric ? "" : "none";
    if (reductionHelp) reductionHelp.style.display = isAsymmetric ? "" : "none";
    setRepresentationCompatibilityWarningsVisible(representationChangedWarningActive);
  }

  function handleDistortionFieldInput() {
    distortionValuesByReductionMHz[activeDistortionReduction] = captureDistortionValuesForStateMHz();
    updateRepresentationChangeControls();
  }

  function bindDistortionFieldListeners() {
    DISTORTION_FIELD_IDS.forEach((id) => {
      const input = document.getElementById(id);
      if (!input || input.dataset.pendingRepresentationBound === "true") return;
      input.addEventListener("input", handleDistortionFieldInput);
      input.dataset.pendingRepresentationBound = "true";
    });
  }

  function setSpectrumDownloadReady(ready) {
    if (!spectrumDownloadButton) return;
    spectrumDownloadButton.disabled = !ready;
  }

  function convertDistortionInputs(fromUnitRaw, toUnitRaw) {
    const fromFactor = distortionToMHzFactor(fromUnitRaw);
    const toFactor = distortionFromMHzFactor(toUnitRaw);
    DISTORTION_FIELD_IDS.forEach((id) => {
      const el = document.getElementById(id);
      if (!el) return;
      const raw = parseFloat(el.value);
      if (!Number.isFinite(raw)) return;
      const mhz = raw * fromFactor;
      el.value = String(mhz * toFactor);
    });
  }

  function updateDistortionUnitLabel() {
    const unit = getCurrentDistortionUnit();
    const unitLabel = document.getElementById("distortion-unit-label");
    if (unitLabel) unitLabel.textContent = unit;
    const unitLabelSextic = document.getElementById("distortion-unit-label-sextic");
    if (unitLabelSextic) unitLabelSextic.textContent = unit;
  }

  function setQrotAutoInfo(message = "") {
    const target = document.getElementById("qrot-auto-info");
    if (!target) return;
    target.innerHTML = `<small>${message || ""}</small>`;
  }

  function downloadFile(contents, filename, type) {
    const blob = new Blob([contents], { type });
    const url = URL.createObjectURL(blob);
    const link = Object.assign(document.createElement("a"), {
      href: url,
      download: filename
    });
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  function parseRotationalSection(text) {
    const lines = text.split(/\r?\n/);
    let inSection = false;
    const data = {};
    for (const line of lines) {
      if (!inSection) {
        if (line.startsWith("#ROTATIONAL")) {
          inSection = true;
        }
        continue;
      }
      if (line.startsWith("#")) break;
      if (!line.trim()) continue;
      const eqIndex = line.indexOf("=");
      if (eqIndex === -1) continue;
      const key = line.slice(0, eqIndex).trim();
      const value = line.slice(eqIndex + 1).trim();
      data[key] = value;
    }
    if (!inSection) {
      throw new Error("Sezione #ROTATIONAL non trovata nel file.");
    }
    return data;
  }

  function parseSciNumber(value) {
    if (value === null || value === undefined) return null;
    const cleaned = String(value).replace(/D/gi, "E");
    const num = parseFloat(cleaned);
    return Number.isNaN(num) ? null : num;
  }

  function inferRotorType(A, B, C) {
    if (![A, B, C].every(Number.isFinite)) return null;
    const tolRel = 1e-9;
    const tolAbs = 1e-9;
    const close = (x, y) => Math.abs(x - y) <= Math.max(tolAbs, tolRel * Math.max(Math.abs(x), Math.abs(y)));
    const minVal = Math.min(A, B, C);
    if (minVal <= 1e-3) return "linear";
    if (close(A, B) && close(B, C)) return "spherical";
    if (close(B, C) && A > B) return "prolate";
    if (close(A, B) && B > C) return "oblate";
    return "asymmetric";
  }

  function gaussianLogToInputText(text) {
    const pickettBlocks = text.match(/-{5,}[\s\S]*?-{5,}/g) || [];
    let pickettBlock = null;
    for (const block of pickettBlocks) {
      if (block.match(/^\s*\d+\s+[-\d+.DE]+\s+[-\d+.DE]+\s+\/\s*[^/]+?\s*\//m)) {
        pickettBlock = block;
      }
    }
    if (!pickettBlock) {
      const globalMatch = text.match(/^\s*\d+\s+[-\d+.DE]+\s+[-\d+.DE]+\s+\/\s*[^/]+?\s*\/.*$/gm);
      if (globalMatch) {
        pickettBlock = globalMatch.join("\n");
      }
    }

    const pickettParams = {};
    const codeParams = {};
    let reduction = null;
    if (pickettBlock) {
      const lines = pickettBlock.split(/\r?\n/);
      for (const line of lines) {
        const match = line.match(/^\s*(\d+)\s+([-\d+.DE]+)\s+([-\d+.DE]+)\s+\/\s*([^/]+?)\s*\//);
        if (!match) continue;
        const code = parseInt(match[1], 10);
        const value = parseSciNumber(match[2]);
        if (value === null) continue;
        const rawLabel = match[4].trim();
        const label = rawLabel.replace(/\s+/g, " ").trim();
        const labelNoSign = label.replace(/^[-+]\s*/, "");
        let canonicalLabel = labelNoSign;
        if (/^DELTA\s+N$/i.test(canonicalLabel)) canonicalLabel = "DELTA J";
        if (/^DELTA\s+NK$/i.test(canonicalLabel)) canonicalLabel = "DELTA JK";
        if (/^delta\s+N$/i.test(canonicalLabel)) canonicalLabel = "delta J";
        if (/^delta\s+NK$/i.test(canonicalLabel)) canonicalLabel = "delta JK";
        pickettParams[canonicalLabel] = value;
        if (Number.isFinite(code)) codeParams[code] = value;
        const labelUpper = canonicalLabel.toUpperCase();
        if (labelUpper.includes("DELTA") || labelUpper.startsWith("PHI")) reduction = "A";
        if (labelUpper === "D1" || labelUpper === "D2" || labelUpper === "DJ" || labelUpper === "DJK" || labelUpper === "DK") {
          if (!reduction) reduction = "S";
        }
        if (/^H\s*(N|K|NK|KN|1|2|3)$/.test(labelUpper) || labelUpper === "H1" || labelUpper === "H2" || labelUpper === "H3") {
          if (!reduction) reduction = "S";
        }
        if (code === 41000) reduction = "A";
        if (code === 50000) reduction = "S";
        if (code === 50100 || code === 51000 || code === 60000) reduction = "S";
        if (code === 40200 || code === 41100 || code === 42000) reduction = "A";
      }
    }

    const constantsRegex = /Rotational constants\s+\((MHZ|GHZ)\):\s*([-\d+.DE]+)\s+([-\d+.DE]+)\s+([-\d+.DE]+)/gi;
    let constantsMatch = null;
    let m;
    while ((m = constantsRegex.exec(text)) !== null) constantsMatch = m;
    let A = null;
    let B = null;
    let C = null;
    if (constantsMatch) {
      const unit = constantsMatch[1].toUpperCase();
      const conv = unit === "GHZ" ? 1000 : 1;
      A = parseSciNumber(constantsMatch[2]);
      B = parseSciNumber(constantsMatch[3]);
      C = parseSciNumber(constantsMatch[4]);
      if (A !== null) A *= conv;
      if (B !== null) B *= conv;
      if (C !== null) C *= conv;
    }

    if (codeParams[10000] !== undefined) A = codeParams[10000];
    if (codeParams[20000] !== undefined) B = codeParams[20000];
    if (codeParams[30000] !== undefined) C = codeParams[30000];
    if (pickettParams.A !== undefined && !Number.isFinite(A)) A = pickettParams.A;
    if (pickettParams.B !== undefined && !Number.isFinite(B)) B = pickettParams.B;
    if (pickettParams.C !== undefined && !Number.isFinite(C)) C = pickettParams.C;

    const dipoleRegex = /Dipole moment \(Debye\):\s*[\r\n]+\s*([-\d+.DE]+)\s+([-\d+.DE]+)\s+([-\d+.DE]+)/gi;
    let dipoleMatch = null;
    while ((m = dipoleRegex.exec(text)) !== null) dipoleMatch = m;
    let mu_a = null;
    let mu_b = null;
    let mu_c = null;
    if (dipoleMatch) {
      mu_a = parseSciNumber(dipoleMatch[1]);
      mu_b = parseSciNumber(dipoleMatch[2]);
      mu_c = parseSciNumber(dipoleMatch[3]);
    }

    const pointGroupRegex = /Full point group\s+([A-Za-z0-9]+)/gi;
    let pointGroupMatch = null;
    while ((m = pointGroupRegex.exec(text)) !== null) pointGroupMatch = m;
    const pointGroup = pointGroupMatch ? pointGroupMatch[1] : null;

    const rotorType = inferRotorType(A, B, C);

    const dvibRegex = /Rotational Constants\s*\(in MHz\)[\s\S]*?Ae=\s*([-\d+.DE]+)\s*A00=\s*([-\d+.DE]+)[\s\S]*?Be=\s*([-\d+.DE]+)\s*B00=\s*([-\d+.DE]+)[\s\S]*?Ce=\s*([-\d+.DE]+)\s*C00=\s*([-\d+.DE]+)/gi;
    let dvibMatch = null;
    while ((m = dvibRegex.exec(text)) !== null) dvibMatch = m;
    let DVibA = null;
    let DVibB = null;
    let DVibC = null;
    if (dvibMatch) {
      const Ae = parseSciNumber(dvibMatch[1]);
      const A00 = parseSciNumber(dvibMatch[2]);
      const Be = parseSciNumber(dvibMatch[3]);
      const B00 = parseSciNumber(dvibMatch[4]);
      const Ce = parseSciNumber(dvibMatch[5]);
      const C00 = parseSciNumber(dvibMatch[6]);
      if (Number.isFinite(Ae) && Number.isFinite(A00)) DVibA = A00 - Ae;
      if (Number.isFinite(Be) && Number.isFinite(B00)) DVibB = B00 - Be;
      if (Number.isFinite(Ce) && Number.isFinite(C00)) DVibC = C00 - Ce;
    }

    const lines = ["#ROTATIONAL"];
    if (rotorType) lines.push(`rotor_type = ${rotorType}`);
    if (pointGroup) lines.push(`Point Group = ${pointGroup}`);
    if (reduction) lines.push(`Watson Reduction = ${reduction}`);
    if (Number.isFinite(A)) lines.push(`A_MHz = ${A}`);
    if (Number.isFinite(B)) lines.push(`B_MHz = ${B}`);
    if (Number.isFinite(C)) lines.push(`C_MHz = ${C}`);
    if (Number.isFinite(mu_a)) lines.push(`Dipole_a_D = ${mu_a}`);
    if (Number.isFinite(mu_b)) lines.push(`Dipole_b_D = ${mu_b}`);
    if (Number.isFinite(mu_c)) lines.push(`Dipole_c_D = ${mu_c}`);
    if (Number.isFinite(DVibA)) lines.push(`DVibA_MHz = ${DVibA}`);
    if (Number.isFinite(DVibB)) lines.push(`DVibB_MHz = ${DVibB}`);
    if (Number.isFinite(DVibC)) lines.push(`DVibC_MHz = ${DVibC}`);

    const pickettLabel = (label) => Object.prototype.hasOwnProperty.call(pickettParams, label) ? pickettParams[label] : null;
    const pickettCode = (code) => Object.prototype.hasOwnProperty.call(codeParams, code) ? codeParams[code] : null;
    const addDistortion = (key, value) => {
      if (Number.isFinite(value)) lines.push(`${key} = ${value}`);
    };

    const deltaJ = pickettCode(200) ?? pickettLabel("DELTA J") ?? pickettLabel("DELTAJ");
    const deltaJK = pickettCode(1100) ?? pickettLabel("DELTA JK") ?? pickettLabel("DELTAJK");
    const deltaK = pickettCode(2000) ?? pickettLabel("DELTA K") ?? pickettLabel("DELTAK");
    const deltaj = (reduction === "A" ? pickettCode(40100) : null) ?? pickettLabel("delta J") ?? pickettLabel("deltaJ");
    const deltak = pickettCode(41000) ?? pickettLabel("delta K") ?? pickettLabel("deltaK");
    const DJ = (reduction === "S" ? pickettCode(200) : null) ?? pickettLabel("DJ");
    const DJK = (reduction === "S" ? pickettCode(1100) : null) ?? pickettLabel("DJK");
    const DK = (reduction === "S" ? pickettCode(2000) : null) ?? pickettLabel("DK");
    const d1 = (reduction === "S" ? pickettCode(40100) : null) ?? pickettLabel("d1") ?? pickettLabel("D1");
    const d2 = (reduction === "S" ? pickettCode(50000) : null) ?? pickettLabel("d2") ?? pickettLabel("D2");
    const phiN = (reduction === "A" ? pickettCode(40200) : null) ?? pickettLabel("phi N") ?? pickettLabel("phiN");
    const phiNK = (reduction === "A" ? pickettCode(41100) : null) ?? pickettLabel("phi NK") ?? pickettLabel("phiNK");
    const phiK = (reduction === "A" ? pickettCode(42000) : null) ?? pickettLabel("phi K") ?? pickettLabel("phiK");
    const PHIN = (reduction === "A" ? pickettCode(300) : null) ?? pickettLabel("PHI N") ?? pickettLabel("PHIN");
    const PHINK = (reduction === "A" ? pickettCode(1200) : null) ?? pickettLabel("PHI NK") ?? pickettLabel("PHINK");
    const PHIKN = (reduction === "A" ? pickettCode(2100) : null) ?? pickettLabel("PHI KN") ?? pickettLabel("PHIKN");
    const PHIK = (reduction === "A" ? pickettCode(3000) : null) ?? pickettLabel("PHI K") ?? pickettLabel("PHIK");
    const HN = (reduction === "S" ? pickettCode(300) : null) ?? pickettLabel("H N") ?? pickettLabel("HN");
    const HNK = (reduction === "S" ? pickettCode(1200) : null) ?? pickettLabel("H NK") ?? pickettLabel("HNK");
    const HKN = (reduction === "S" ? pickettCode(2100) : null) ?? pickettLabel("H KN") ?? pickettLabel("HKN");
    const HK = (reduction === "S" ? pickettCode(3000) : null) ?? pickettLabel("H K") ?? pickettLabel("HK");
    const h1 = (reduction === "S" ? pickettCode(50100) : null) ?? pickettLabel("h 1") ?? pickettLabel("h1") ?? pickettLabel("H1");
    const h2 = (reduction === "S" ? pickettCode(51000) : null) ?? pickettLabel("h 2") ?? pickettLabel("h2") ?? pickettLabel("H2");
    const h3 = (reduction === "S" ? pickettCode(60000) : null) ?? pickettLabel("h 3") ?? pickettLabel("h3") ?? pickettLabel("H3");

    if (reduction === "A") {
      addDistortion("DELTA J_MHz", deltaJ);
      addDistortion("DELTA JK_MHz", deltaJK);
      addDistortion("DELTA K_MHz", deltaK);
      addDistortion("delta J_MHz", deltaj);
      addDistortion("delta K_MHz", deltak);
      addDistortion("PHI N_MHz", PHIN);
      addDistortion("PHI NK_MHz", PHINK);
      addDistortion("PHI KN_MHz", PHIKN);
      addDistortion("PHI K_MHz", PHIK);
      addDistortion("phi N_MHz", phiN);
      addDistortion("phi NK_MHz", phiNK);
      addDistortion("phi K_MHz", phiK);
    } else if (reduction === "S") {
      addDistortion("DJ_MHz", DJ);
      addDistortion("DJK_MHz", DJK);
      addDistortion("DK_MHz", DK);
      addDistortion("d1_MHz", d1);
      addDistortion("d2_MHz", d2);
      addDistortion("H N_MHz", HN);
      addDistortion("H NK_MHz", HNK);
      addDistortion("H KN_MHz", HKN);
      addDistortion("H K_MHz", HK);
      addDistortion("h1_MHz", h1);
      addDistortion("h2_MHz", h2);
      addDistortion("h3_MHz", h3);
    } else {
      addDistortion("DJ_MHz", DJ);
      addDistortion("DJK_MHz", DJK);
      addDistortion("DK_MHz", DK);
      addDistortion("d1_MHz", d1);
      addDistortion("d2_MHz", d2);
      addDistortion("H N_MHz", HN);
      addDistortion("H NK_MHz", HNK);
      addDistortion("H KN_MHz", HKN);
      addDistortion("H K_MHz", HK);
      addDistortion("h1_MHz", h1);
      addDistortion("h2_MHz", h2);
      addDistortion("h3_MHz", h3);
      addDistortion("DELTA J_MHz", deltaJ);
      addDistortion("DELTA JK_MHz", deltaJK);
      addDistortion("DELTA K_MHz", deltaK);
      addDistortion("delta J_MHz", deltaj);
      addDistortion("delta K_MHz", deltak);
      addDistortion("PHI N_MHz", PHIN);
      addDistortion("PHI NK_MHz", PHINK);
      addDistortion("PHI KN_MHz", PHIKN);
      addDistortion("PHI K_MHz", PHIK);
      addDistortion("phi N_MHz", phiN);
      addDistortion("phi NK_MHz", phiNK);
      addDistortion("phi K_MHz", phiK);
    }

    if (lines.length === 1) {
      throw new Error("Unable to extract data from the Gaussian file.");
    }
    return lines.join("\n") + "\n";
  }

  function buildFileData(raw) {
    const readRaw = (key, alternates = []) => {
      if (Object.prototype.hasOwnProperty.call(raw, key)) return raw[key];
      for (const alt of alternates) {
        if (Object.prototype.hasOwnProperty.call(raw, alt)) return raw[alt];
      }
      return undefined;
    };
    const readNum = (key, alternates = []) => {
      const value = readRaw(key, alternates);
      if (value === undefined) return null;
      return parseSciNumber(value);
    };
    const readBool = (key, alternates = []) => {
      const value = readRaw(key, alternates);
      if (value === undefined || value === null) return null;
      const normalized = String(value).trim().toLowerCase();
      if (["1", "true", "yes", "on"].includes(normalized)) return true;
      if (["0", "false", "no", "off"].includes(normalized)) return false;
      return null;
    };
    const declaredDistortionUnitRaw = readRaw(
      "Distortion_Unit",
      [
        "Distortion Unit",
        "Distortion Units",
        "Distortion units",
        "distortion_unit",
        "distortion unit",
        "distortion units"
      ]
    );
    const hasDeclaredDistortionUnit = (
      declaredDistortionUnitRaw !== undefined
      && declaredDistortionUnitRaw !== null
      && String(declaredDistortionUnitRaw).trim() !== ""
    );
    const declaredDistortionUnit = hasDeclaredDistortionUnit
      ? normalizeDistortionUnit(declaredDistortionUnitRaw)
      : null;
    const declaredDistortionToMHz = declaredDistortionUnit
      ? distortionToMHzFactor(declaredDistortionUnit)
      : null;
    const explicitDistortionUnits = new Set();
    const readDistortionWithUnitsMHz = ({ mhzKeys = [], khzKeys = [], unitlessKeys = [] } = {}) => {
      for (const key of mhzKeys) {
        const num = readNum(key);
        if (num !== null) {
          explicitDistortionUnits.add("MHz");
          return num;
        }
      }
      for (const key of khzKeys) {
        const num = readNum(key);
        if (num !== null) {
          explicitDistortionUnits.add("kHz");
          return num * 0.001;
        }
      }
      if (declaredDistortionToMHz !== null) {
        for (const key of unitlessKeys) {
          const num = readNum(key);
          if (num !== null) return num * declaredDistortionToMHz;
        }
      }
      return null;
    };

    const rotorRaw = readRaw("rotor_type");
    let rotorType = null;
    if (rotorRaw) {
      const rotorLower = rotorRaw.toLowerCase();
      if (rotorLower.includes("asymmetric")) rotorType = "asymmetric";
      else if (rotorLower.includes("linear")) rotorType = "linear";
      else if (rotorLower.includes("spherical")) rotorType = "spherical";
      else if (rotorLower.includes("oblate")) rotorType = "oblate";
      else if (rotorLower.includes("prolate")) rotorType = "prolate";
    }

    const reductionRaw = readRaw("Watson Reduction");
    let reduction = null;
    if (reductionRaw) {
      const r = reductionRaw.trim().toUpperCase();
      if (r.startsWith("S")) reduction = "S";
      if (r.startsWith("A")) reduction = "A";
    }

    const distortion = {
      A: {
        DJ: readDistortionWithUnitsMHz({
          mhzKeys: ["DELTA J_MHz", "DELTA_J_MHz"],
          khzKeys: ["DELTA J_kHz", "DELTA_J_kHz", "DELTA J_KHz", "DELTA_J_KHz"],
          unitlessKeys: ["DELTA J", "DELTA_J"]
        }),
        DJK: readDistortionWithUnitsMHz({
          mhzKeys: ["DELTA JK_MHz", "DELTA_JK_MHz"],
          khzKeys: ["DELTA JK_kHz", "DELTA_JK_kHz", "DELTA JK_KHz", "DELTA_JK_KHz"],
          unitlessKeys: ["DELTA JK", "DELTA_JK"]
        }),
        DK: readDistortionWithUnitsMHz({
          mhzKeys: ["DELTA K_MHz", "DELTA_K_MHz"],
          khzKeys: ["DELTA K_kHz", "DELTA_K_kHz", "DELTA K_KHz", "DELTA_K_KHz"],
          unitlessKeys: ["DELTA K", "DELTA_K"]
        }),
        dJ: readDistortionWithUnitsMHz({
          mhzKeys: ["delta J_MHz", "delta_J_MHz"],
          khzKeys: ["delta J_kHz", "delta_J_kHz", "delta J_KHz", "delta_J_KHz"],
          unitlessKeys: ["delta J", "delta_J"]
        }),
        dK: readDistortionWithUnitsMHz({
          mhzKeys: ["delta K_MHz", "delta_K_MHz"],
          khzKeys: ["delta K_kHz", "delta_K_kHz", "delta K_KHz", "delta_K_KHz"],
          unitlessKeys: ["delta K", "delta_K"]
        }),
        HJ: readDistortionWithUnitsMHz({
          mhzKeys: ["PHI N_MHz", "PHI_N_MHz"],
          khzKeys: ["PHI N_kHz", "PHI_N_kHz", "PHI N_KHz", "PHI_N_KHz"],
          unitlessKeys: ["PHI N", "PHI_N"]
        }),
        HJK: readDistortionWithUnitsMHz({
          mhzKeys: ["PHI NK_MHz", "PHI_NK_MHz"],
          khzKeys: ["PHI NK_kHz", "PHI_NK_kHz", "PHI NK_KHz", "PHI_NK_KHz"],
          unitlessKeys: ["PHI NK", "PHI_NK"]
        }),
        HKJ: readDistortionWithUnitsMHz({
          mhzKeys: ["PHI KN_MHz", "PHI_KN_MHz"],
          khzKeys: ["PHI KN_kHz", "PHI_KN_kHz", "PHI KN_KHz", "PHI_KN_KHz"],
          unitlessKeys: ["PHI KN", "PHI_KN"]
        }),
        HK: readDistortionWithUnitsMHz({
          mhzKeys: ["PHI K_MHz", "PHI_K_MHz"],
          khzKeys: ["PHI K_kHz", "PHI_K_kHz", "PHI K_KHz", "PHI_K_KHz"],
          unitlessKeys: ["PHI K", "PHI_K"]
        }),
        h1: readDistortionWithUnitsMHz({
          mhzKeys: ["phi N_MHz", "phi_N_MHz"],
          khzKeys: ["phi N_kHz", "phi_N_kHz", "phi N_KHz", "phi_N_KHz"],
          unitlessKeys: ["phi N", "phi_N"]
        }),
        h2: readDistortionWithUnitsMHz({
          mhzKeys: ["phi NK_MHz", "phi_NK_MHz"],
          khzKeys: ["phi NK_kHz", "phi_NK_kHz", "phi NK_KHz", "phi_NK_KHz"],
          unitlessKeys: ["phi NK", "phi_NK"]
        }),
        h3: readDistortionWithUnitsMHz({
          mhzKeys: ["phi K_MHz", "phi_K_MHz"],
          khzKeys: ["phi K_kHz", "phi_K_kHz", "phi K_KHz", "phi_K_KHz"],
          unitlessKeys: ["phi K", "phi_K"]
        })
      },
      S: {
        DJ: readDistortionWithUnitsMHz({
          mhzKeys: ["DJ_MHz"],
          khzKeys: ["DJ_kHz", "DJ_KHz"],
          unitlessKeys: ["DJ"]
        }),
        DJK: readDistortionWithUnitsMHz({
          mhzKeys: ["DJK_MHz"],
          khzKeys: ["DJK_kHz", "DJK_KHz"],
          unitlessKeys: ["DJK"]
        }),
        DK: readDistortionWithUnitsMHz({
          mhzKeys: ["DK_MHz"],
          khzKeys: ["DK_kHz", "DK_KHz"],
          unitlessKeys: ["DK"]
        }),
        dJ: readDistortionWithUnitsMHz({
          mhzKeys: ["d1_MHz"],
          khzKeys: ["d1_kHz", "d1_KHz"],
          unitlessKeys: ["d1", "D1"]
        }),
        dK: readDistortionWithUnitsMHz({
          mhzKeys: ["d2_MHz"],
          khzKeys: ["d2_kHz", "d2_KHz"],
          unitlessKeys: ["d2", "D2"]
        }),
        HJ: readDistortionWithUnitsMHz({
          mhzKeys: ["H N_MHz", "H_N_MHz"],
          khzKeys: ["H N_kHz", "H_N_kHz", "H N_KHz", "H_N_KHz"],
          unitlessKeys: ["H N", "H_N"]
        }),
        HJK: readDistortionWithUnitsMHz({
          mhzKeys: ["H NK_MHz", "H_NK_MHz"],
          khzKeys: ["H NK_kHz", "H_NK_kHz", "H NK_KHz", "H_NK_KHz"],
          unitlessKeys: ["H NK", "H_NK"]
        }),
        HKJ: readDistortionWithUnitsMHz({
          mhzKeys: ["H KN_MHz", "H_KN_MHz"],
          khzKeys: ["H KN_kHz", "H_KN_kHz", "H KN_KHz", "H_KN_KHz"],
          unitlessKeys: ["H KN", "H_KN"]
        }),
        HK: readDistortionWithUnitsMHz({
          mhzKeys: ["H K_MHz", "H_K_MHz"],
          khzKeys: ["H K_kHz", "H_K_kHz", "H K_KHz", "H_K_KHz"],
          unitlessKeys: ["H K", "H_K"]
        }),
        h1: readDistortionWithUnitsMHz({
          mhzKeys: ["h1_MHz"],
          khzKeys: ["h1_kHz", "h1_KHz"],
          unitlessKeys: ["h1", "H1", "h 1"]
        }),
        h2: readDistortionWithUnitsMHz({
          mhzKeys: ["h2_MHz"],
          khzKeys: ["h2_kHz", "h2_KHz"],
          unitlessKeys: ["h2", "H2", "h 2"]
        }),
        h3: readDistortionWithUnitsMHz({
          mhzKeys: ["h3_MHz"],
          khzKeys: ["h3_kHz", "h3_KHz"],
          unitlessKeys: ["h3", "H3", "h 3"]
        })
      }
    };
    const inferredDistortionUnit = declaredDistortionUnit || (
      explicitDistortionUnits.size === 1
        ? Array.from(explicitDistortionUnits)[0]
        : null
    );

    return {
      raw,
      rotorType,
      representation: normalizeRepresentationKey(readRaw("representation") || null),
      reduction,
      pointGroup: readRaw("Point Group") || null,
      T: readNum("T_K"),
      A: readNum("A_MHz"),
      B: readNum("B_MHz"),
      C: readNum("C_MHz"),
      DVibA: readNum("DVibA_MHz"),
      DVibB: readNum("DVibB_MHz"),
      DVibC: readNum("DVibC_MHz"),
      mu_a: readNum("Dipole_a_D"),
      mu_b: readNum("Dipole_b_D"),
      mu_c: readNum("Dipole_c_D"),
      qrot: readNum("Q_rot"),
      sD2Id: readNum("S_D2_ID", ["s_d2_id", "S_d2_ID"]),
      J_min: readNum("J_MIN", ["J_min"]),
      J_max: readNum("J_MAX", ["J_max"]),
      autoEstimateJRange: readBool("AUTO_ESTIMATE_J_RANGE", ["auto_estimate_j_range"]),
      distortionUnit: inferredDistortionUnit,
      distortion,
      hfs: {
        I_nuc: readNum("I_NUC_1", ["I_nuc_1", "I_NUC", "I_nuc"]),
        I_nuc_1: readNum("I_NUC_1", ["I_nuc_1", "I_NUC", "I_nuc"]),
        I_nuc_2: readNum("I_NUC_2", ["I_nuc_2"]),
        F_max: readNum("F_MAX", ["F_max"]),
        chi_aa: readNum("chi_aa_1", ["chi_aa_1_MHz", "chi_aa", "chi_aa_MHz"]),
        chi_bb: readNum("chi_bb_1", ["chi_bb_1_MHz", "chi_bb", "chi_bb_MHz"]),
        chi_cc: readNum("chi_cc_1", ["chi_cc_1_MHz", "chi_cc", "chi_cc_MHz"]),
        chi_ab: readNum("chi_ab_1", ["chi_ab_1_MHz", "chi_ab", "chi_ab_MHz"]),
        chi_ac: readNum("chi_ac_1", ["chi_ac_1_MHz", "chi_ac", "chi_ac_MHz"]),
        chi_bc: readNum("chi_bc_1", ["chi_bc_1_MHz", "chi_bc", "chi_bc_MHz"]),
        chi_aa_1: readNum("chi_aa_1", ["chi_aa_1_MHz", "chi_aa", "chi_aa_MHz"]),
        chi_bb_1: readNum("chi_bb_1", ["chi_bb_1_MHz", "chi_bb", "chi_bb_MHz"]),
        chi_cc_1: readNum("chi_cc_1", ["chi_cc_1_MHz", "chi_cc", "chi_cc_MHz"]),
        chi_ab_1: readNum("chi_ab_1", ["chi_ab_1_MHz", "chi_ab", "chi_ab_MHz"]),
        chi_ac_1: readNum("chi_ac_1", ["chi_ac_1_MHz", "chi_ac", "chi_ac_MHz"]),
        chi_bc_1: readNum("chi_bc_1", ["chi_bc_1_MHz", "chi_bc", "chi_bc_MHz"]),
        chi_aa_2: readNum("chi_aa_2", ["chi_aa_2_MHz"]),
        chi_bb_2: readNum("chi_bb_2", ["chi_bb_2_MHz"]),
        chi_cc_2: readNum("chi_cc_2", ["chi_cc_2_MHz"]),
        chi_ab_2: readNum("chi_ab_2", ["chi_ab_2_MHz"]),
        chi_ac_2: readNum("chi_ac_2", ["chi_ac_2_MHz"]),
        chi_bc_2: readNum("chi_bc_2", ["chi_bc_2_MHz"])
      }
    };
  }

  function applyDistortionForReduction() {
    // Reduction switching is now a pure state swap between stored A and S sets.
    const reduction = normalizeReductionKey();
    let set = distortionValuesByReductionMHz[reduction];
    if (!set && loadedFileData?.distortion?.[reduction]) {
      set = cloneDistortionValues(loadedFileData.distortion[reduction]);
      distortionValuesByReductionMHz[reduction] = set;
    }
    activeDistortionReduction = reduction;
    restoreDistortionValuesFromMHz(set, { clearMissing: true });
    updateRepresentationChangeControls();
  }

  function applyFileDataToForm(data) {
    representationChangedWarningActive = false;
    if (data.rotorType) {
      document.getElementById("rotorType").value = data.rotorType;
    }
    if (data.representation) {
      document.getElementById("representation").value = normalizeRepresentationKey(data.representation);
    }
    if (data.reduction) {
      document.getElementById("reduction").value = data.reduction;
    }
    if (data.distortionUnit && distortionUnitSelect) {
      distortionUnitSelect.value = normalizeDistortionUnit(data.distortionUnit);
      previousDistortionUnit = getCurrentDistortionUnit();
      updateDistortionUnitLabel();
    }
    REDUCTION_KEYS.forEach((reduction) => {
      const set = cloneDistortionValues(data.distortion?.[reduction]);
      distortionValuesByReductionMHz[reduction] = hasAnyDistortionValue(set) ? set : null;
    });
    activeDistortionReduction = normalizeReductionKey(data.reduction || getText("reduction"));
    updateEquations();

    if (data.pointGroup) {
      const group = document.getElementById("groupSymmetry");
      if (group) {
        const hasOption = Array.from(group.options).some(opt => opt.value === data.pointGroup);
        if (!hasOption) {
          const opt = document.createElement("option");
          opt.value = data.pointGroup;
          opt.textContent = data.pointGroup;
          group.appendChild(opt);
        }
        group.value = data.pointGroup;
      }
    }

    setValueIfPresent("T", data.T);
    setValueIfPresent("qrot_override", data.qrot);
    setValueIfPresent("A", data.A);
    setValueIfPresent("B", data.B);
    setValueIfPresent("C", data.C);
    setValueIfPresent("J_min", data.J_min);
    setValueIfPresent("J_max", data.J_max);
    if (typeof data.autoEstimateJRange === "boolean") {
      document.getElementById("auto_estimate_j_range").checked = data.autoEstimateJRange;
    }
    setValueIfPresent("DVibA", data.DVibA);
    setValueIfPresent("DVibB", data.DVibB);
    setValueIfPresent("DVibC", data.DVibC);
    setValueIfPresent("mu_a", data.mu_a);
    setValueIfPresent("mu_b", data.mu_b);
    setValueIfPresent("mu_c", data.mu_c);
    if (data.hfs) {
      setValueIfPresent("I_nuc", data.hfs.I_nuc_1 ?? data.hfs.I_nuc);
      setValueIfPresent("I_nuc_2", data.hfs.I_nuc_2);
      setValueIfPresent("F_max", data.hfs.F_max);
      setValueIfPresent("chi_aa", data.hfs.chi_aa_1 ?? data.hfs.chi_aa);
      setValueIfPresent("chi_bb", data.hfs.chi_bb_1 ?? data.hfs.chi_bb);
      setValueIfPresent("chi_cc", data.hfs.chi_cc_1 ?? data.hfs.chi_cc);
      setValueIfPresent("chi_ab", data.hfs.chi_ab_1 ?? data.hfs.chi_ab);
      setValueIfPresent("chi_ac", data.hfs.chi_ac_1 ?? data.hfs.chi_ac);
      setValueIfPresent("chi_bc", data.hfs.chi_bc_1 ?? data.hfs.chi_bc);
      setValueIfPresent("chi_aa_2", data.hfs.chi_aa_2);
      setValueIfPresent("chi_bb_2", data.hfs.chi_bb_2);
      setValueIfPresent("chi_cc_2", data.hfs.chi_cc_2);
      setValueIfPresent("chi_ab_2", data.hfs.chi_ab_2);
      setValueIfPresent("chi_ac_2", data.hfs.chi_ac_2);
      setValueIfPresent("chi_bc_2", data.hfs.chi_bc_2);
    }
    applyDistortionForReduction();
    updateRepresentationChangeControls();
    updateJRangeUiState();
    setJRangeStatus("");

    // If the default frequency window is too low for a linear rotor, auto-expand
    // it to something usable (otherwise CO-like cases show "No lines" by default).
    const type = document.getElementById("rotorType")?.value;
    if (type === "linear") {
      const xMaxRaw = get("freq_max_spcat");
      const { maxMHz: xMaxMHz } = normalizeFreqWindowMHz(get("freq_min_spcat"), xMaxRaw);
      const first = estimateFirstLinearLineMHz();
      // Trigger only if the user didn't tweak the default window (or it's clearly too small).
      if (!Number.isNaN(xMaxRaw) && (Math.abs(xMaxRaw - 20000.0) < 1e-9 || (first !== null && first > xMaxMHz))) {
        const Jmax = readIntegerField("J_max");
        const B = get("B");
        const safeJ = Number.isFinite(Jmax) ? Jmax : 30;
        if (!Number.isNaN(B)) {
          const estMax = 2 * B * (safeJ + 1);
          document.getElementById("freq_min_spcat").value = "0.0";
          document.getElementById("freq_max_spcat").value = String(Math.ceil(estMax));
        }
      }
    }
  }

  function buildInputTextFromForm() {
    const rotorType = document.getElementById("rotorType").value;
    const representation = normalizeRepresentationKey();
    const reduction = getText("reduction");
    const pointGroup = getText("groupSymmetry");
    const T = get("T");
    const A = get("A");
    const B = get("B");
    const C = get("C");
    const DVibA = get("DVibA");
    const DVibB = get("DVibB");
    const DVibC = get("DVibC");
    const mu_a = get("mu_a");
    const mu_b = get("mu_b");
    const mu_c = get("mu_c");
    const qrotOverride = get("qrot_override");
    const sD2Id = loadedFileData && Number.isFinite(loadedFileData.sD2Id) ? loadedFileData.sD2Id : NaN;
    const J_min = readIntegerField("J_min");
    const J_max = readIntegerField("J_max");
    const autoEstimateJRange = getChecked("auto_estimate_j_range");
    const I_nuc = get("I_nuc");
    const I_nuc_2 = get("I_nuc_2");
    const F_max = get("F_max");
    const chi_aa = get("chi_aa");
    const chi_bb = get("chi_bb");
    const chi_cc = get("chi_cc");
    const chi_ab = get("chi_ab");
    const chi_ac = get("chi_ac");
    const chi_bc = get("chi_bc");
    const chi_aa_2 = get("chi_aa_2");
    const chi_bb_2 = get("chi_bb_2");
    const chi_cc_2 = get("chi_cc_2");
    const chi_ab_2 = get("chi_ab_2");
    const chi_ac_2 = get("chi_ac_2");
    const chi_bc_2 = get("chi_bc_2");

    // Always export distortion in MHz for compatibility with WMS-FitRot parser.
    const DJ = readDistortionFieldMHz("DJ");
    const DJK = readDistortionFieldMHz("DJK");
    const DK = readDistortionFieldMHz("DK");
    const dJ = readDistortionFieldMHz("dJ");
    const dK = readDistortionFieldMHz("dK");
    const HJ = readDistortionFieldMHz("HJ");
    const HJK = readDistortionFieldMHz("HJK");
    const HKJ = readDistortionFieldMHz("HKJ");
    const HK = readDistortionFieldMHz("HK");
    const h1 = readDistortionFieldMHz("h1");
    const h2 = readDistortionFieldMHz("h2");
    const h3 = readDistortionFieldMHz("h3");

    const lines = [];
    lines.push("#ROTATIONAL");
    lines.push(`rotor_type = ${rotorType}`);
    lines.push(`representation = ${representation}`);
    lines.push(`Watson Reduction = ${reduction}`);
    lines.push(`Point Group = ${pointGroup}`);
    lines.push(`T_K = ${Number.isFinite(T) ? T : ""}`);
    lines.push(`A_MHz = ${Number.isFinite(A) ? A : ""}`);
    lines.push(`B_MHz = ${Number.isFinite(B) ? B : ""}`);
    lines.push(`C_MHz = ${Number.isFinite(C) ? C : ""}`);
    lines.push(`DVibA_MHz= ${Number.isFinite(DVibA) ? DVibA : 0}`);
    lines.push(`DVibB_MHz= ${Number.isFinite(DVibB) ? DVibB : 0}`);
    lines.push(`DVibC_MHz= ${Number.isFinite(DVibC) ? DVibC : 0}`);
    lines.push(`Dipole_a_D = ${Number.isFinite(mu_a) ? mu_a : ""}`);
    lines.push(`Dipole_b_D = ${Number.isFinite(mu_b) ? mu_b : ""}`);
    lines.push(`Dipole_c_D = ${Number.isFinite(mu_c) ? mu_c : ""}`);
    lines.push(`I_NUC = ${Number.isFinite(I_nuc) ? I_nuc : 0}`);
    lines.push(`I_NUC_1 = ${Number.isFinite(I_nuc) ? I_nuc : 0}`);
    lines.push(`I_NUC_2 = ${Number.isFinite(I_nuc_2) ? I_nuc_2 : 0}`);
    if (Number.isFinite(F_max)) lines.push(`F_MAX = ${F_max}`);
    lines.push(`chi_aa = ${Number.isFinite(chi_aa) ? chi_aa : 0}`);
    lines.push(`chi_bb = ${Number.isFinite(chi_bb) ? chi_bb : 0}`);
    lines.push(`chi_cc = ${Number.isFinite(chi_cc) ? chi_cc : 0}`);
    lines.push(`chi_ab = ${Number.isFinite(chi_ab) ? chi_ab : 0}`);
    lines.push(`chi_ac = ${Number.isFinite(chi_ac) ? chi_ac : 0}`);
    lines.push(`chi_bc = ${Number.isFinite(chi_bc) ? chi_bc : 0}`);
    lines.push(`chi_aa_1 = ${Number.isFinite(chi_aa) ? chi_aa : 0}`);
    lines.push(`chi_bb_1 = ${Number.isFinite(chi_bb) ? chi_bb : 0}`);
    lines.push(`chi_cc_1 = ${Number.isFinite(chi_cc) ? chi_cc : 0}`);
    lines.push(`chi_ab_1 = ${Number.isFinite(chi_ab) ? chi_ab : 0}`);
    lines.push(`chi_ac_1 = ${Number.isFinite(chi_ac) ? chi_ac : 0}`);
    lines.push(`chi_bc_1 = ${Number.isFinite(chi_bc) ? chi_bc : 0}`);
    lines.push(`chi_aa_2 = ${Number.isFinite(chi_aa_2) ? chi_aa_2 : 0}`);
    lines.push(`chi_bb_2 = ${Number.isFinite(chi_bb_2) ? chi_bb_2 : 0}`);
    lines.push(`chi_cc_2 = ${Number.isFinite(chi_cc_2) ? chi_cc_2 : 0}`);
    lines.push(`chi_ab_2 = ${Number.isFinite(chi_ab_2) ? chi_ab_2 : 0}`);
    lines.push(`chi_ac_2 = ${Number.isFinite(chi_ac_2) ? chi_ac_2 : 0}`);
    lines.push(`chi_bc_2 = ${Number.isFinite(chi_bc_2) ? chi_bc_2 : 0}`);
    if (Number.isInteger(J_min)) lines.push(`J_MIN = ${J_min}`);
    if (Number.isInteger(J_max)) lines.push(`J_MAX = ${J_max}`);
    lines.push("FREQ_UNIT = MHz");
    lines.push(`AUTO_ESTIMATE_J_RANGE = ${autoEstimateJRange ? "true" : "false"}`);

    if (String(reduction).toUpperCase().startsWith("A")) {
      lines.push(`DELTA J_MHz =  ${Number.isFinite(DJ) ? DJ : ""}`);
      lines.push(`DELTA JK_MHz = ${Number.isFinite(DJK) ? DJK : ""}`);
      lines.push(`DELTA K_MHz = ${Number.isFinite(DK) ? DK : ""}`);
      lines.push(`delta J_MHz = ${Number.isFinite(dJ) ? dJ : ""}`);
      lines.push(`delta K_MHz = ${Number.isFinite(dK) ? dK : ""}`);
      lines.push(`PHI N_MHz = ${Number.isFinite(HJ) ? HJ : ""}`);
      lines.push(`PHI NK_MHz = ${Number.isFinite(HJK) ? HJK : ""}`);
      lines.push(`PHI KN_MHz = ${Number.isFinite(HKJ) ? HKJ : ""}`);
      lines.push(`PHI K_MHz = ${Number.isFinite(HK) ? HK : ""}`);
      lines.push(`phi N_MHz = ${Number.isFinite(h1) ? h1 : ""}`);
      lines.push(`phi NK_MHz = ${Number.isFinite(h2) ? h2 : ""}`);
      lines.push(`phi K_MHz = ${Number.isFinite(h3) ? h3 : ""}`);
    } else {
      lines.push(`DJ_MHz = ${Number.isFinite(DJ) ? DJ : ""}`);
      lines.push(`DJK_MHz = ${Number.isFinite(DJK) ? DJK : ""}`);
      lines.push(`DK_MHz = ${Number.isFinite(DK) ? DK : ""}`);
      lines.push(`d1_MHz = ${Number.isFinite(dJ) ? dJ : ""}`);
      lines.push(`d2_MHz = ${Number.isFinite(dK) ? dK : ""}`);
      lines.push(`H N_MHz = ${Number.isFinite(HJ) ? HJ : ""}`);
      lines.push(`H NK_MHz = ${Number.isFinite(HJK) ? HJK : ""}`);
      lines.push(`H KN_MHz = ${Number.isFinite(HKJ) ? HKJ : ""}`);
      lines.push(`H K_MHz = ${Number.isFinite(HK) ? HK : ""}`);
      lines.push(`h1_MHz = ${Number.isFinite(h1) ? h1 : ""}`);
      lines.push(`h2_MHz = ${Number.isFinite(h2) ? h2 : ""}`);
      lines.push(`h3_MHz = ${Number.isFinite(h3) ? h3 : ""}`);
    }
    if (Number.isFinite(qrotOverride)) {
      lines.push(`Q_rot=${qrotOverride}`);
    }
    if (Number.isFinite(sD2Id)) {
      lines.push(`S_D2_ID=${Math.round(sD2Id)}`);
    }
    return lines.join("\n") + "\n";
  }

  function applyRotorUiState() {
    const type = document.getElementById("rotorType").value;
    const AField = document.getElementById("A").parentElement;
    const BField = document.getElementById("B").parentElement;
    const CField = document.getElementById("C").parentElement;

    AField.style.display = '';
    BField.style.display = '';
    CField.style.display = '';
    if (type === 'linear' || type === 'spherical') {
      AField.style.display = 'none';
      CField.style.display = 'none';
    } else if (type === 'prolate') {
      CField.style.display = 'none';
    } else if (type === 'oblate') {
      AField.style.display = 'none';
    }

    const repLabel = document.getElementById("representation")?.parentElement;
    const redLabel = document.getElementById("reduction")?.parentElement;
    if (repLabel) repLabel.style.display = type === "asymmetric" ? "" : "none";
    if (redLabel) redLabel.style.display = type === "asymmetric" ? "" : "none";

    const aAxis = document.getElementById("a_axis");
    const bAxis = document.getElementById("b_axis");
    const cAxis = document.getElementById("c_axis");
    const setAxis = (el, visible, checked, disabled) => {
      if (!el) return;
      if (el.parentElement) el.parentElement.style.display = visible ? "" : "none";
      if (typeof checked === "boolean") el.checked = checked;
      if (typeof disabled === "boolean") el.disabled = disabled;
    };
    if (type === "linear") {
      setAxis(aAxis, true, true, false);
      setAxis(bAxis, false, false, true);
      setAxis(cAxis, false, false, true);
    } else if (type === "prolate") {
      setAxis(aAxis, true, true, false);
      setAxis(bAxis, false, false, true);
      setAxis(cAxis, false, false, true);
    } else if (type === "oblate") {
      setAxis(aAxis, false, false, true);
      setAxis(bAxis, false, false, true);
      setAxis(cAxis, true, true, false);
    } else if (type === "spherical") {
      setAxis(aAxis, false, false, true);
      setAxis(bAxis, false, false, true);
      setAxis(cAxis, false, false, true);
    } else {
      setAxis(aAxis, true, aAxis ? aAxis.checked : undefined, false);
      setAxis(bAxis, true, bAxis ? bAxis.checked : undefined, false);
      setAxis(cAxis, true, cAxis ? cAxis.checked : undefined, false);
    }

    const showQ = document.getElementById("showQ")?.parentElement;
    const colorQ = document.getElementById("colorQ")?.parentElement;
    if (showQ) showQ.style.display = type === "linear" ? "none" : "";
    if (colorQ) colorQ.style.display = type === "linear" ? "none" : "";

    const rotorNote = document.getElementById("rotor-constants-note");
    if (rotorNote) {
      let note = "";
      if (type === "prolate") {
        note = "Prolate symmetric top: A > B = C (symmetry axis along A).";
      } else if (type === "oblate") {
        note = "Oblate symmetric top: A = B > C (symmetry axis along C).";
      } else if (type === "spherical") {
        note = "Spherical top: A = B = C.";
      } else if (type === "linear") {
        note = "Linear rotor: A and C are not used; B is the rotational constant; dipole is along the molecular axis (μ∥).";
      } else {
        note = "Asymmetric top: A, B, C are all distinct (no symmetry axis).";
      }
      rotorNote.innerHTML = `<small>${note}</small>`;
    }

    const muALabelText = document.getElementById("mu_a_label_text");
    const muBLabel = document.getElementById("mu_b_label");
    const muCLabel = document.getElementById("mu_c_label");
    const muAField = document.getElementById("mu_a")?.parentElement;
    const muBField = document.getElementById("mu_b")?.parentElement;
    const muCField = document.getElementById("mu_c")?.parentElement;
    const muSphereWrap = document.getElementById("mu_spherical_wrap");
    if (muAField && muBField && muCField) {
      if (type === "linear") {
        if (muALabelText) muALabelText.innerHTML = "μ<sub>∥</sub>";
        muAField.style.display = "";
        muBField.style.display = "none";
        muCField.style.display = "none";
        if (muSphereWrap) muSphereWrap.style.display = "none";
      } else if (type === "prolate") {
        if (muALabelText) muALabelText.innerHTML = "μ<sub>a</sub>";
        muAField.style.display = "";
        muBField.style.display = "none";
        muCField.style.display = "none";
        if (muSphereWrap) muSphereWrap.style.display = "none";
      } else if (type === "oblate") {
        if (muALabelText) muALabelText.innerHTML = "μ<sub>a</sub>";
        muAField.style.display = "none";
        muBField.style.display = "none";
        muCField.style.display = "";
        if (muSphereWrap) muSphereWrap.style.display = "none";
      } else if (type === "spherical") {
        if (muALabelText) muALabelText.innerHTML = "μ<sub>a</sub>";
        muAField.style.display = "none";
        muBField.style.display = "none";
        muCField.style.display = "none";
        if (muSphereWrap) muSphereWrap.style.display = "";
      } else {
        if (muALabelText) muALabelText.innerHTML = "μ<sub>a</sub>";
        muAField.style.display = "";
        muBField.style.display = "";
        muCField.style.display = "";
        if (muSphereWrap) muSphereWrap.style.display = "none";
      }
    }

    const setRowDisplay = (id, show) => {
      const row = document.getElementById(id);
      if (row) row.style.display = show ? "" : "none";
    };
    const disableInput = (id, flag, reset) => {
      const el = document.getElementById(id);
      if (!el) return;
      el.disabled = flag;
      if (flag && reset) el.value = "0";
    };
    if (type === "linear") {
      setRowDisplay("DJ/DeltaJ", true);
      setRowDisplay("DJK/DeltaJK", false);
      setRowDisplay("DK/DeltaK", false);
      setRowDisplay("d1/dJ", false);
      setRowDisplay("d2/dK", false);
      setRowDisplay("HJ/PHIN", true);
      setRowDisplay("HJK/PHINK", false);
      setRowDisplay("HKJ/PHIKN", false);
      setRowDisplay("HK/PHIK", false);
      setRowDisplay("h1/phiN", false);
      setRowDisplay("h2/phiNK", false);
      setRowDisplay("h3/phiK", false);
      disableInput("DJK", true, true);
      disableInput("DK", true, true);
      disableInput("dJ", true, true);
      disableInput("dK", true, true);
      disableInput("DJ", false);
      disableInput("HJK", true, true);
      disableInput("HKJ", true, true);
      disableInput("HK", true, true);
      disableInput("h1", true, true);
      disableInput("h2", true, true);
      disableInput("h3", true, true);
      disableInput("HJ", false);
    } else if (type === "spherical") {
      setRowDisplay("DJ/DeltaJ", true);
      setRowDisplay("DJK/DeltaJK", false);
      setRowDisplay("DK/DeltaK", false);
      setRowDisplay("d1/dJ", false);
      setRowDisplay("d2/dK", false);
      setRowDisplay("HJ/PHIN", true);
      setRowDisplay("HJK/PHINK", false);
      setRowDisplay("HKJ/PHIKN", false);
      setRowDisplay("HK/PHIK", false);
      setRowDisplay("h1/phiN", false);
      setRowDisplay("h2/phiNK", false);
      setRowDisplay("h3/phiK", false);
      disableInput("DJK", true, true);
      disableInput("DK", true, true);
      disableInput("dJ", true, true);
      disableInput("dK", true, true);
      disableInput("DJ", false);
      disableInput("HJK", true, true);
      disableInput("HKJ", true, true);
      disableInput("HK", true, true);
      disableInput("h1", true, true);
      disableInput("h2", true, true);
      disableInput("h3", true, true);
      disableInput("HJ", false);
    } else {
      setRowDisplay("DJ/DeltaJ", true);
      setRowDisplay("DJK/DeltaJK", true);
      setRowDisplay("DK/DeltaK", true);
      setRowDisplay("d1/dJ", type === "asymmetric");
      setRowDisplay("d2/dK", type === "asymmetric");
      setRowDisplay("HJ/PHIN", true);
      setRowDisplay("HJK/PHINK", true);
      setRowDisplay("HKJ/PHIKN", true);
      setRowDisplay("HK/PHIK", true);
      setRowDisplay("h1/phiN", type === "asymmetric");
      setRowDisplay("h2/phiNK", type === "asymmetric");
      setRowDisplay("h3/phiK", type === "asymmetric");
      disableInput("DJ", false);
      disableInput("DJK", false);
      disableInput("DK", false);
      disableInput("dJ", false);
      disableInput("dK", false);
      disableInput("HJ", false);
      disableInput("HJK", false);
      disableInput("HKJ", false);
      disableInput("HK", false);
      disableInput("h1", false);
      disableInput("h2", false);
      disableInput("h3", false);
    }
  }
  function updateEquations() {
    const displayedDistortionValues = captureDisplayedDistortionValues();
    const currentGroupSymmetry = getText("groupSymmetry");
    const reduction = getText("reduction");
    const rep = normalizeRepresentationKey();
    const rotorType = document.getElementById("rotorType").value;
    const redLabel = reduction === "S" ? "S-reduction" : "A-reduction";
    if (rotorType === "asymmetric") {
      reduction === "S" ? document.getElementById("d1/dJ").innerHTML = '<label>d<sub>1</sub>:  <input type="number" id="dJ"  step="any" value="0"></label>' : document.getElementById("d1/dJ").innerHTML = '<label>δ<sub>J</sub>:  <input type="number" id="dJ"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("d2/dK").innerHTML = '<label>d<sub>2</sub>:  <input type="number" id="dK"  step="any" value="0"></label>' : document.getElementById("d2/dK").innerHTML = '<label>δ<sub>K</sub>:  <input type="number" id="dK"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("DJ/DeltaJ").innerHTML = '<label>D<sub>J</sub>:  <input type="number" id="DJ"  step="any" value="0"></label>' : document.getElementById("DJ/DeltaJ").innerHTML = '<label>Δ<sub>J</sub>:  <input type="number" id="DJ"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("DJK/DeltaJK").innerHTML = '<label>D<sub>JK</sub>: <input type="number" id="DJK" step="any" value="0"></label>' : document.getElementById("DJK/DeltaJK").innerHTML = '<label>Δ<sub>JK</sub>: <input type="number" id="DJK" step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("DK/DeltaK").innerHTML = '<label>D<sub>K</sub>:  <input type="number" id="DK"  step="any" value="0"></label>' : document.getElementById("DK/DeltaK").innerHTML = '<label>Δ<sub>K</sub>:  <input type="number" id="DK"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("HJ/PHIN").innerHTML = '<label>H<sub>N</sub>:  <input type="number" id="HJ"  step="any" value="0"></label>' : document.getElementById("HJ/PHIN").innerHTML = '<label>Φ<sub>N</sub>:  <input type="number" id="HJ"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("HJK/PHINK").innerHTML = '<label>H<sub>NK</sub>: <input type="number" id="HJK" step="any" value="0"></label>' : document.getElementById("HJK/PHINK").innerHTML = '<label>Φ<sub>NK</sub>: <input type="number" id="HJK" step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("HKJ/PHIKN").innerHTML = '<label>H<sub>KN</sub>: <input type="number" id="HKJ" step="any" value="0"></label>' : document.getElementById("HKJ/PHIKN").innerHTML = '<label>Φ<sub>KN</sub>: <input type="number" id="HKJ" step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("HK/PHIK").innerHTML = '<label>H<sub>K</sub>:  <input type="number" id="HK"  step="any" value="0"></label>' : document.getElementById("HK/PHIK").innerHTML = '<label>Φ<sub>K</sub>:  <input type="number" id="HK"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("h1/phiN").innerHTML = '<label>h<sub>1</sub>:  <input type="number" id="h1"  step="any" value="0"></label>' : document.getElementById("h1/phiN").innerHTML = '<label>φ<sub>N</sub>:  <input type="number" id="h1"  step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("h2/phiNK").innerHTML = '<label>h<sub>2</sub>: <input type="number" id="h2" step="any" value="0"></label>' : document.getElementById("h2/phiNK").innerHTML = '<label>φ<sub>NK</sub>: <input type="number" id="h2" step="any" value="0"></label>';
      reduction === "S" ? document.getElementById("h3/phiK").innerHTML = '<label>h<sub>3</sub>:  <input type="number" id="h3"  step="any" value="0"></label>' : document.getElementById("h3/phiK").innerHTML = '<label>φ<sub>K</sub>:  <input type="number" id="h3"  step="any" value="0"></label>';
    } else {
      document.getElementById("d1/dJ").innerHTML = '<label>d<sub>J</sub>:  <input type="number" id="dJ"  step="any" value="0"></label>';
      document.getElementById("d2/dK").innerHTML = '<label>d<sub>K</sub>:  <input type="number" id="dK"  step="any" value="0"></label>';
      document.getElementById("DJ/DeltaJ").innerHTML = '<label>D<sub>J</sub>:  <input type="number" id="DJ"  step="any" value="0"></label>';
      document.getElementById("DJK/DeltaJK").innerHTML = '<label>D<sub>JK</sub>: <input type="number" id="DJK" step="any" value="0"></label>';
      document.getElementById("DK/DeltaK").innerHTML = '<label>D<sub>K</sub>:  <input type="number" id="DK"  step="any" value="0"></label>';
      document.getElementById("HJ/PHIN").innerHTML = '<label>H<sub>N</sub>:  <input type="number" id="HJ"  step="any" value="0"></label>';
      document.getElementById("HJK/PHINK").innerHTML = '<label>H<sub>NK</sub>: <input type="number" id="HJK" step="any" value="0"></label>';
      document.getElementById("HKJ/PHIKN").innerHTML = '<label>H<sub>KN</sub>: <input type="number" id="HKJ" step="any" value="0"></label>';
      document.getElementById("HK/PHIK").innerHTML = '<label>H<sub>K</sub>:  <input type="number" id="HK"  step="any" value="0"></label>';
      document.getElementById("h1/phiN").innerHTML = '<label>h<sub>1</sub>:  <input type="number" id="h1"  step="any" value="0"></label>';
      document.getElementById("h2/phiNK").innerHTML = '<label>h<sub>2</sub>: <input type="number" id="h2" step="any" value="0"></label>';
      document.getElementById("h3/phiK").innerHTML = '<label>h<sub>3</sub>:  <input type="number" id="h3"  step="any" value="0"></label>';
    }
      
      switch (rotorType) {
              case "asymmetric":
                document.getElementById("groupSymmetry").innerHTML = '<option value="C1" selected>C1</option><option value="Cs">Cs</option><option value="Ci">Ci</option><option value="C2">C2</option><option value="C2v">C2v</option><option value="C2h">C2h</option><option value="D2">D2</option><option value="D2h">D2h</option><option value="D2d">D2d</option>'
                break;
              case "linear":
                document.getElementById("groupSymmetry").innerHTML = '<option value="C1" selected>C1</option><option value="Cinfv">Cinfv</option><option value="Dinfh">(Dinfh)</option>'
                break;
              case "spherical":
                document.getElementById("groupSymmetry").innerHTML = '<option value="Td">Td</option><option value="Oh">Oh</option><option value="Ih">Ih</option>'
                break;
              case "prolate":
                document.getElementById("groupSymmetry").innerHTML = '<option value="C3">C3</option><option value="C4">C4</option><option value="C5">C5</option><option value="C6">C6</option><option value="C3v">C3v</option><option value="C4v">C4v</option><option value="C6v">C6v</option><option value="C3h">C3h</option><option value="C4h">C4h</option><option value="C6h">C6h</option><option value="D3">D3</option><option value="D4">D4</option><option value="D6">D6</option><option value="D3h">D3h</option><option value="D4h">D4h</option><option value="D6h">D6h</option><option value="D3d">D3d</option><option value="D4d">D4d</option><option value="S4">S4</option><option value="S6">S6</option><option value="S8">S8</option>'
                break;
              case "oblate":
                document.getElementById("groupSymmetry").innerHTML = '<option value="C3">C3</option><option value="C4">C4</option><option value="C5">C5</option><option value="C6">C6</option><option value="C3v">C3v</option><option value="C4v">C4v</option><option value="C6v">C6v</option><option value="C3h">C3h</option><option value="C4h">C4h</option><option value="C6h">C6h</option><option value="D3">D3</option><option value="D4">D4</option><option value="D6">D6</option><option value="D3h">D3h</option><option value="D4h">D4h</option><option value="D6h">D6h</option><option value="D3d">D3d</option><option value="D4d">D4d</option><option value="S4">S4</option><option value="S6">S6</option><option value="S8">S8</option>'
                break;
      }

    const groupSymmetrySelect = document.getElementById("groupSymmetry");
    if (groupSymmetrySelect) {
      const hasCurrentGroup = Array.from(groupSymmetrySelect.options).some((option) => option.value === currentGroupSymmetry);
      if (hasCurrentGroup) {
        groupSymmetrySelect.value = currentGroupSymmetry;
      }
    }
      
    const repInfo = getRepresentationUiInfo(rep);
    const repAxis = repInfo.axes.z;
    const repX = repInfo.axes.x;
    const repY = repInfo.axes.y;
    const kDef = "K^2 = \\hat J_" + repAxis + "^2";
    const jpmDef = "\\hat J_\\pm = \\hat J_" + repX + " \\pm i\\hat J_" + repY;
    const basisNote = "representation: " + repInfo.label + ", \\hat J_z\\equiv\\hat J_" + repAxis +
      ",\\;\\hat J_x\\equiv\\hat J_" + repX + ",\\;\\hat J_y\\equiv\\hat J_" + repY;

    const quarticTerms = reduction === "S"
      ? `$$\\hat H_{\\mathrm{cd}} = D_J\\hat J^4 + D_{JK}\\hat J^2\\hat K^2 + D_K\\hat K^4
          + d_1\\hat J^2(\\hat J_+^2 + \\hat J_-^2)
          + d_2(\\hat J_+^4 + \\hat J_-^4)$$`
      : `$$\\hat H_{\\mathrm{cd}} = Δ_J\\hat J^4 + Δ_{JK}\\hat J^2\\hat K^2 + Δ_K\\hat K^4
          + δ_J\\hat J^2(\\hat J_+^2 + \\hat J_-^2)
          + \\frac{δ_K}{2}\\left[\\hat K^2(\\hat J_+^2 + \\hat J_-^2) + (\\hat J_+^2 + \\hat J_-^2)\\hat K^2\\right]$$`;
    const sexticTerms = reduction === "S"
      ? `$$\\hat H_{\\mathrm{cd}}^{(6)} = H_N\\hat J^6 + H_{NK}\\hat J^4\\hat K^2 + H_{KN}\\hat J^2\\hat K^4 + H_K\\hat K^6
          + h_1\\hat J^2(\\hat J_+^4 + \\hat J_-^4)
          + \\frac{h_2}{2}\\left[\\hat K^2(\\hat J_+^4 + \\hat J_-^4) + (\\hat J_+^4 + \\hat J_-^4)\\hat K^2\\right]
          + h_3(\\hat J_+^6 + \\hat J_-^6)$$`
      : `$$\\hat H_{\\mathrm{cd}}^{(6)} = \\Phi_N\\hat J^6 + \\Phi_{NK}\\hat J^4\\hat K^2 + \\Phi_{KN}\\hat J^2\\hat K^4 + \\Phi_K\\hat K^6
          + \\phi_N\\hat J^4(\\hat J_+^2 + \\hat J_-^2)
          + \\frac{\\phi_{NK}}{2}\\left[\\hat J^2\\hat K^2(\\hat J_+^2 + \\hat J_-^2) + (\\hat J_+^2 + \\hat J_-^2)\\hat K^2\\hat J^2\\right]
          + \\frac{\\phi_K}{2}\\left[\\hat K^4(\\hat J_+^2 + \\hat J_-^2) + (\\hat J_+^2 + \\hat J_-^2)\\hat K^4\\right]$$`;

    const isAsym = rotorType === "asymmetric";
    const DJLabel = isAsym ? (reduction === "S" ? "D_J" : "\\Delta_J") : "D_J";
    const DJKLabel = isAsym ? (reduction === "S" ? "D_{JK}" : "\\Delta_{JK}") : "D_{JK}";
    const DKLabel = isAsym ? (reduction === "S" ? "D_K" : "\\Delta_K") : "D_K";
    const HJLabel = isAsym ? (reduction === "S" ? "H_N" : "\\Phi_N") : "H_J";
    const HJKLabel = isAsym ? (reduction === "S" ? "H_{NK}" : "\\Phi_{NK}") : "H_{JK}";
    const HKJLabel = isAsym ? (reduction === "S" ? "H_{KN}" : "\\Phi_{KN}") : "H_{KJ}";
    const HKLabel = isAsym ? (reduction === "S" ? "H_K" : "\\Phi_K") : "H_K";

    const ecdSym = `$$E_{\\mathrm{cd}} = +${DJLabel} J^2(J+1)^2 + ${DJKLabel} J(J+1)K^2 + ${DKLabel}K^4
      + ${HJLabel}[J(J+1)]^3 + ${HJKLabel}[J(J+1)]^2K^2 + ${HKJLabel}J(J+1)K^4 + ${HKLabel}K^6$$`;
    const lineStrength = `$$S_{if} = (2J_u+1)\\left|\\begin{pmatrix} J_u & 1 & J_l \\\\ -K_u & q & K_l \\end{pmatrix}\\right|^2,\\quad q=\\Delta K$$`;
    const perpStrength = `$$S_{if}^{(\\perp)} = \\frac{1}{2}S_{if}\\quad(\\Delta K=\\pm1)$$`;
    const qrotSym = `$$Q_{\\mathrm{rot}}(T)=\\frac{1}{\\sigma}\\sum_{J,K}(2J+1)(2-\\delta_{K,0})\\,g_{\\mathrm{ns}}(sp_{J,K})\\,e^{-E_{J,K}/kT}$$`;
    const intensitySym = `$$I_{if}\\propto \\frac{g_{\\mathrm{ns}}(sp_l)}{Q_{\\mathrm{rot}}}\\,S_{if}\\,\\mu_{\\mathrm{axis}}^2\\,\\nu\\,(1-e^{-h\\nu/kT})\\,e^{-E_l/kT}$$`;

    let eqHtml = "";
    if (rotorType === "asymmetric") {
      eqHtml = `
        <p><strong>Hamiltonian (summary)</strong></p>
        <p>$$\\hat H = \\hat H_{\\mathrm{rot}} + \\hat H_{\\mathrm{cd}} + \\hat H_Q$$</p>
        <p>$$\\hat H_{\\mathrm{rot}} = A\\hat J_a^2 + B\\hat J_b^2 + C\\hat J_c^2$$</p>
        <p>$$\\hat J^2 = \\hat J_a^2 + \\hat J_b^2 + \\hat J_c^2,\\quad ${kDef},\\quad ${jpmDef}$$</p>
        <p>$$\\text{Axis convention: }${basisNote}$$</p>
        <p>$$\\text{Watson ${redLabel} quartics:}$$</p>
        <p>${quarticTerms}</p>
        <p>$$\\text{Watson ${redLabel} sextics:}$$</p>
        <p>${sexticTerms}</p>
        <p><strong>Hyperfine quadrupole (implemented)</strong></p>
        <p>Enabled when one or two nuclei have \\(I_n>1/2\\) and at least one \\(\\chi^{(n)}\\) component is non-zero.</p>
        <p>$$\\chi_{ij}^{(n)}\\rightarrow \\chi_{ij}^{\\prime(n)}=\\chi_{ij}^{(n)}-\\tfrac{1}{3}\\mathrm{Tr}(\\chi^{(n)})\\,\\delta_{ij},\\qquad
            \\tilde\\chi_q^{(n)}=\\chi_{q,\\mathrm{traceless}}^{(n)}/f_{\\mathrm{red}}(I_n)$$</p>
        <p>$$\\hat H_Q=\\sum_{n=1}^{N_Q}\\sum_{q=-2}^{2}\\tilde\\chi_q^{(n)}\\,\\hat C_q^{(2)}\\otimes\\hat T_q^{(2)}(I_n),\\qquad N_Q\\le 2$$</p>
        <p>$$\\langle J,K|\\hat C_q^{(2)}|J',K'\\rangle
            =(-1)^{J-K}\\sqrt{(2J+1)(2J'+1)}
            \\begin{pmatrix} J&2&J'\\\\-K&q&K'\\end{pmatrix},\\quad K'=K-q$$</p>
        <p>$$\\text{Diagonalization HFS: for each }F\\text{ build the direct Wang-spin basis }|J,w,I_{12};F\\rangle$$</p>
        <p>$$\\text{with }I_{12}=I_1+I_2\\text{ if two nuclei are active, then diagonalize the full }H_F\\text{ block.}$$</p>
        <p>Rotational eigenstates are still used for labeling and wavefunction analysis, but the hyperfine Hamiltonian is assembled directly in the coupled Wang-spin basis.</p>
        <p><strong>Intensity (summary)</strong></p>
        <p>$$M_{(J_l,w_l,I_{12}),(J_u,w_u,I_{12})}=
            \\sqrt{(2F_l+1)(2F_u+1)}
            \\begin{Bmatrix}J_l&F_l&I_{12}\\\\F_u&J_u&1\\end{Bmatrix}
            \\sum_{g\\in\\{a,b,c\\}}\\mu_g\\,\\langle J_lw_l|T_g^{(1)}|J_uw_u\\rangle$$</p>
        <p>$$A=C_{F_l}^\\dagger M C_{F_u},\\quad S_{\\mathrm{rot}}=|A|^2,\\quad
            S_{\\mathrm{tot}}=S_{\\mathrm{rot}}\\times f_{\\mathrm{HFS}}$$</p>
        <p>$$\\Delta F=0,\\pm1\\quad(0\\nleftrightarrow0)$$</p>
        <p>$$S_{if} \\propto \\left|\\mu_a\\langle i|\\hat T_a^{(1)}|f\\rangle + \\mu_b\\langle i|\\hat T_b^{(1)}|f\\rangle + \\mu_c\\langle i|\\hat T_c^{(1)}|f\\rangle\\right|^2$$</p>
        <p>$$P_i \\propto (2J_i+1)\\,g_{\\mathrm{ns}}(\\Gamma_{i})\\,e^{-E_i/kT}/Q_{\\mathrm{rot}}$$</p>
        <p>$$Q_{\\mathrm{rot}}(T)=\\frac{1}{\\sigma}\\sum_{J,\\alpha}(2J+1)\\,g_{\\mathrm{ns}}\\,e^{-E_{J\\alpha}/kT}$$</p>
      `;
    } else {
      let energyHtml = "";
      let selectionHtml = "";
      if (rotorType === "linear") {
        energyHtml = `
          <p><strong>Hamiltonian (linear)</strong></p>
          <p>$$K=0$$</p>
          <p>$$E_J = B J(J+1) - D_J J^2(J+1)^2 + H_J J^3(J+1)^3$$</p>
          <p>$$\\text{Simple mode: quadrupole off.}$$</p>
          <p>$$\\text{With quadrupole: exact direct-}F\\text{ matrix projected onto the physical }K=0\\text{ subspace, then diagonalized.}$$</p>
        `;
        selectionHtml = "<p><strong>Selection rules</strong></p><p>Linear (a-type only): \\(\\Delta J=\\pm1\\), \\(\\Delta K=0\\). With quadrupole active, the same \\(\\mu_{\\parallel}\\) operator is evaluated between the mixed direct-\\(F\\) states.</p>";
      } else if (rotorType === "spherical") {
        energyHtml = `
          <p><strong>Hamiltonian (spherical)</strong></p>
          <p>$$A=B=C$$</p>
          <p>$$E_J = B J(J+1) + D_J J^2(J+1)^2 + H_J J^3(J+1)^3$$</p>
          <p>$$\\text{Simple mode: }d_J=d_K=h_1=h_2=h_3=0\\text{ and quadrupole off.}$$</p>
          <p>$$\\text{Otherwise: full matrix diagonalization (Watson quartic/sextic, optional }\\hat H_Q\\text{).}$$</p>
        `;
        selectionHtml = "<p><strong>Selection rules</strong></p><p>As symmetric-top in the \\(|J,K\\rangle\\) basis: a-type \\(\\Delta K=0\\), b/c-type \\(\\Delta K=\\pm1\\), with \\(\\Delta J=0,\\pm1\\).</p>";
      } else {
        const kAxis = rotorType === "oblate" ? "K \\equiv K_c" : "K \\equiv K_a";
        const e0 = rotorType === "prolate"
          ? "B J(J+1) + (A-B)K^2"
          : rotorType === "oblate"
            ? "B J(J+1) + (C-B)K^2"
            : "B J(J+1)\\; (A=B=C)";
        const modeNote = rotorType === "oblate"
          ? `
          <p>$$\\text{Simple mode: }d_J=d_K=h_1=h_2=h_3=0\\text{ and quadrupole off.}$$</p>
          <p>$$\\text{Otherwise: full matrix diagonalization in the selected oblate-axis convention.}$$</p>
          `
          : `
          <p>$$\\text{Simple mode is evaluated on the canonical symmetry axis }(K\\equiv K_a).$$</p>
          <p>$$\\text{Simple mode: }d_J=d_K=h_1=h_2=h_3=0\\text{ and quadrupole off.}$$</p>
          <p>$$\\text{Otherwise: full matrix diagonalization (Watson quartic/sextic, optional }\\hat H_Q\\text{).}$$</p>
          `;
        energyHtml = `
          <p><strong>Hamiltonian (symmetric, hybrid)</strong></p>
          <p>$$${kAxis},\\; 0\\le K\\le J$$</p>
          <p>$$E_{J,K} = ${e0} + E_{\\mathrm{cd}}$$</p>
          <p>${ecdSym}</p>
          ${modeNote}
        `;
        if (rotorType === "oblate") {
          selectionHtml = "<p><strong>Selection rules</strong></p><p>Simple limit: c-type (\\(\\mu_c\\)) \\(\\Delta K=0\\), \\(\\Delta J=\\pm1\\); a/b-type (\\(\\mu_a,\\mu_b\\)) \\(\\Delta K=\\pm1\\), \\(\\Delta J=0,\\pm1\\). With off-diagonal/quadrupole terms, transitions are computed from mixed eigenvectors after full diagonalization.</p>";
        } else {
          selectionHtml = "<p><strong>Selection rules</strong></p><p>a-type (\\(\\mu_a\\)): \\(\\Delta K=0\\), \\(\\Delta J=\\pm1\\); b/c-type (\\(\\mu_b,\\mu_c\\)): \\(\\Delta K=\\pm1\\), \\(\\Delta J=0,\\pm1\\).</p>";
        }
      }

      eqHtml = `
        ${energyHtml}
        ${selectionHtml}
        <p><strong>Intensity (used here)</strong></p>
        <p>${lineStrength}</p>
        <p>${perpStrength}</p>
        <p>${intensitySym}</p>
        <p>${qrotSym}</p>
      `;
    }

    const eqBox = document.getElementById("equations");
    eqBox.innerHTML = eqHtml;
    restoreDisplayedDistortionValues(displayedDistortionValues);
    applyRotorUiState();
    bindDistortionFieldListeners();
    updateRepresentationChangeControls();
    if (window.MathJax && window.MathJax.typesetPromise) {
      window.MathJax.typesetPromise([eqBox]).catch(() => {});
    }
  }
  document.getElementById("rotorType").addEventListener("change", updateEquations);
  if (reductionSelect) {
    reductionSelect.addEventListener("change", () => {
      persistActiveReductionState();
      activeDistortionReduction = normalizeReductionKey();
      updateEquations();
      applyDistortionForReduction();
      updateRepresentationChangeControls();
    });
  }
  if (representationSelect) {
    representationSelect.addEventListener("change", (event) => {
      void event;
      representationChangedWarningActive = true;
      distortionValuesByReductionMHz[activeDistortionReduction] = captureDistortionValuesForStateMHz();
      updateEquations();
      updateRepresentationChangeControls();
    });
  }
  let previousDistortionUnit = getCurrentDistortionUnit();
  const distortionUnitSelect = document.getElementById("distortion_unit");
  if (distortionUnitSelect) {
    updateDistortionUnitLabel();
    distortionUnitSelect.addEventListener("change", () => {
      const nextUnit = getCurrentDistortionUnit();
      convertDistortionInputs(previousDistortionUnit, nextUnit);
      distortionValuesByReductionMHz[activeDistortionReduction] = captureDistortionValuesForStateMHz();
      previousDistortionUnit = nextUnit;
      updateDistortionUnitLabel();
    });
  }
  bindDistortionFieldListeners();
  updateEquations();
  updateRepresentationChangeControls();
  setSpectrumDownloadReady(Boolean(lastSpectrum));
  updateJRangeUiState();

  const autoEstimateCheckbox = document.getElementById("auto_estimate_j_range");
  if (autoEstimateCheckbox) {
    autoEstimateCheckbox.addEventListener("change", () => {
      updateJRangeUiState();
      setJRangeStatus("");
    });
  }

  const estimateJRangeBtn = document.getElementById("estimate-j-range");
  if (estimateJRangeBtn) {
    estimateJRangeBtn.addEventListener("click", async () => {
      try {
        const base = collectSimulationBase();
        setJRangeStatus("Estimating J range...");
        await estimateJRangeFromWindow(base);
      } catch (err) {
        console.error("Unable to estimate J range:", err);
        alert(err.message || "Unable to estimate J range.");
      }
    });
  }

  ["showP", "showQ", "showR", "intensityScale"].forEach(id => {
    document.getElementById(id).addEventListener("change", () => {
      if (lastSpectrum) void renderSpectrum(lastSpectrum);
    });
  });
  ["colorP", "colorQ", "colorR"].forEach(id => {
    document.getElementById(id).addEventListener("input", () => {
      if (lastSpectrum) void renderSpectrum(lastSpectrum);
    });
  });
  if (spectrumDownloadButton) {
    spectrumDownloadButton.addEventListener("click", async () => {
      spectrumDownloadButton.disabled = true;
      try {
        await downloadPlotlyFigureAsPng(
          document.getElementById("spectrum"),
          `wms_rot_spectrum_${normalizeRepresentationKey()}_${getText("reduction")}`
        );
      } catch (err) {
        console.error("Unable to download spectrum PNG:", err);
        alert(err.message || "Unable to download the high-resolution PNG.");
      } finally {
        setSpectrumDownloadReady(Boolean(lastSpectrum));
      }
    });
  }

  document.getElementById("calculate").addEventListener("click", async e => {
    e.preventDefault();
    void window.WMSPwaEnhancements?.prepareNotificationsFromGesture?.();

    const calculateBtn = document.getElementById("calculate");
    setSimulatedSpectrumVisible(true);
    globalLoadingOverlay.show("Calculating spectrum…");
    setQrotAutoInfo("");
    if (calculateBtn) calculateBtn.disabled = true;
    void window.WMSPwaEnhancements?.notifyBackground?.({
      title: "WMS-Rot",
      body: "Spectrum calculation started.",
      tag: "wms-rot-calculate"
    });

    try {
      const base = collectSimulationBase();
      let jRange = null;
      if (getChecked("auto_estimate_j_range")) {
        globalLoadingOverlay.setMessage("Estimating J range from frequency window...");
        const estimate = await estimateJRangeFromWindow(base);
        jRange = { J_min: estimate.J_min, J_max: estimate.J_max };
      }
      globalLoadingOverlay.setMessage("Calculating spectrum...");
      const pyvars = buildSimulationPyvars(base, jRange);

      Object.entries(pyvars).forEach(([k, v]) => pyodide.globals.set(k, v));
      lastPyvars = { ...pyvars };

      // per vedere errori veri di Pyodide:
      const plotJson = await pyodide.runPythonAsync(`
clear_wigner_cache()
df = simulate_rigid_spectrum(
    T, A, B, C,
    DJ, DJK, DK, dJ, dK,
    HJ, HJK, HKJ, HK, h1, h2, h3,
    mu_a, mu_b, mu_c,
    J_max, intensity_cut,
    groupSymmetry,
    rotorType,
    a_checked, b_checked, c_checked,
    True, J_min
)
import json
x = df['Frequency (MHz)'].tolist() if not df.empty else []
if 'Intensity' in df.columns:
    y_abs = df['Intensity'].tolist()
else:
    y_abs = df['Relative intensity'].tolist()
if 'Relative intensity' in df.columns:
    y_rel = df['Relative intensity'].tolist()
else:
    y_rel = y_abs
Jl = df['Jl'].tolist() if not df.empty else []
Ka_l = df['Ka_l'].tolist() if not df.empty else []
Kc_l = df['Kc_l'].tolist() if not df.empty else []
Ju = df['Ju'].tolist() if not df.empty else []
Ka_u = df['Ka_u'].tolist() if not df.empty else []
Kc_u = df['Kc_u'].tolist() if not df.empty else []
Fl = df['Fl'].tolist() if ('Fl' in df.columns and not df.empty) else [None] * len(x)
Fu = df['Fu'].tolist() if ('Fu' in df.columns and not df.empty) else [None] * len(x)
I12_l = df['I12_l'].tolist() if ('I12_l' in df.columns and not df.empty) else [None] * len(x)
I12_u = df['I12_u'].tolist() if ('I12_u' in df.columns and not df.empty) else [None] * len(x)
Branch = df['Branch'].tolist() if not df.empty else []
json.dumps({
  "x": x, "y_abs": y_abs, "y_rel": y_rel,
  "Jl": Jl, "Ka_l": Ka_l, "Kc_l": Kc_l,
  "Ju": Ju, "Ka_u": Ka_u, "Kc_u": Kc_u,
  "Fl": Fl, "Fu": Fu, "I12_l": I12_l, "I12_u": I12_u,
  "Branch": Branch,
  "qrot_used": globals().get("LAST_QROT_USED", None),
  "qrot_source": globals().get("LAST_QROT_SOURCE", None)
})
`);
      const {
        x,
        y_abs,
        y_rel,
        Jl,
        Ka_l,
        Kc_l,
        Ju,
        Ka_u,
        Kc_u,
        Fl,
        Fu,
        I12_l,
        I12_u,
        Branch,
        qrot_used,
        qrot_source
      } = JSON.parse(plotJson);
      const branches = { P: [], Q: [], R: [] };
      for (let i = 0; i < x.length; i += 1) {
        const branch = Branch[i] || "";
        if (!branches[branch]) continue;
        branches[branch].push({
          x: x[i],
          yAbs: y_abs[i],
          yRel: y_rel[i],
          Jl: Jl[i],
          Ka_l: Ka_l[i],
          Kc_l: Kc_l[i],
          Ju: Ju[i],
          Ka_u: Ka_u[i],
          Kc_u: Kc_u[i],
          Fl: Fl[i],
          Fu: Fu[i],
          I12_l: I12_l[i],
          I12_u: I12_u[i],
          branch: branch
        });
      }
      lastSpectrum = { branches, totalLines: x.length };
      await Promise.resolve(renderSpectrum(lastSpectrum));
      setPostProcessingReady(true);
      if (pyvars.QROT_OVERRIDE === null && Number.isFinite(qrot_used)) {
        const qrotStr = Number(qrot_used).toPrecision(8);
        const autoTag = String(qrot_source || "").toLowerCase() === "auto" ? "auto" : "calculated";
        setQrotAutoInfo(`Qrot (${autoTag}) = ${qrotStr}`);
      } else {
        setQrotAutoInfo("");
      }
      void window.WMSPwaEnhancements?.notifyBackground?.({
        title: "WMS-Rot",
        body: `Spectrum calculation completed (${x.length} lines).`,
        tag: "wms-rot-calculate"
      });
      // Give Plotly a tick to paint before removing the overlay.
      await new Promise(r => requestAnimationFrame(() => r()));
    } catch (err) {
      console.error("CALC ERROR:", err);
      void window.WMSPwaEnhancements?.notifyBackground?.({
        title: "WMS-Rot",
        body: "Spectrum calculation failed.",
        tag: "wms-rot-calculate"
      });
      alert(String(err));
    } finally {
      globalLoadingOverlay.hide();
      if (calculateBtn) calculateBtn.disabled = false;
    }

  });

  const downloadInputFileBtn = document.getElementById("download-input-file");
  if (downloadInputFileBtn) {
    downloadInputFileBtn.addEventListener("click", () => {
      try {
        const inputText = buildInputTextFromForm();
        const timestamp = new Date().toISOString().replace(/:/g, "-").replace("T", "_").slice(0, 19);
        downloadFile(inputText, `wms-rot_input_${timestamp}.txt`, "text/plain;charset=utf-8");
      } catch (err) {
        console.error("Unable to export input file:", err);
        alert("Unable to export input file.");
      }
    });
  }

  document.getElementById("download-results").addEventListener("click", async () => {
    if (!lastPyvars) {
      alert("Please calculate a spectrum first.");
      return;
    }

    Object.entries(lastPyvars).forEach(([k, v]) => pyodide.globals.set(k, v));
      
      
    const sleep = (ms) => new Promise(r => setTimeout(r, ms));


    const csvString = await pyodide.runPythonAsync(`df.to_csv(index=False)`);
    downloadFile(csvString, "spectrum.csv", "text/csv");
      
    await sleep(250);

    const levelsCsv = await pyodide.runPythonAsync(`
import pandas as pd
import numpy as np
CMC_SPCAT = 29979.2458
J_min = int(globals().get("J_min", 0))
if LAST_DF_HF is not None and not LAST_DF_HF.empty:
    out = LAST_DF_HF.copy()
    out['Energy (MHz)'] = out['Energy(cm^-1)'] * CMC_SPCAT
    out = out.rename(columns={
        'J_dom': 'J',
        'Energy(cm^-1)': 'Energy (cm-1, egy)',
        'Parity': 'sp_l'
    })
    keep = ['J', 'F', 'alpha', 'Ka', 'Kc', 'Energy (MHz)', 'Energy (cm-1, egy)', 'sp_l', 'tau_dom', 'purity']
    if 'I12_dom' in out.columns:
        keep.append('I12_dom')
    csv_str = out[keep].sort_values(['J', 'F', 'Energy (MHz)']).to_csv(index=False)
else:
    levels_list = []

    # For linear rotors SPCAT writes a single K=0 level per J. Using the full
    # asymmetric-top Hamiltonian would generate (2J+1) degenerate levels, which
    # does not match .egy; so we handle linear explicitly here.
    if rotorType == "linear":
        for J in range(J_min, J_max + 1):
            E_val = (B * J * (J + 1)) - (DJ * (J**2) * ((J + 1)**2)) + (HJ * (J**3) * ((J + 1)**3))
            E_cm = float(E_val) / CMC_SPCAT
            levels_list.append((J, 0, 0, float(E_val), f"{E_cm:.6f}", "ee"))
    else:
        rep_user = _normalize_rep_key(REPRESENTATION)
        A_int, B_int, C_int = float(A), float(B), float(C)
        rep_int = rep_user
        use_wang_labels = _pickett_rep_axes(rep_user)[0] in ("a", "c")
        for J in range(J_min, J_max + 1):
            H = H_matrix(A_int, B_int, C_int, DJ, DJK, DK, dJ, dK, J, reduction=REDUCTION, rep=rep_int,
                         HJ=HJ, HJK=HJK, HKJ=HKJ, HK=HK, h1=h1, h2=h2, h3=h3)
            W = wang_transform(J)
            H_w = W.T @ H @ W
            E, U_w = np.linalg.eigh(H_w)
            idx = np.argsort(E.real)
            E = E.real[idx]
            U_w = U_w[:, idx]
            if not use_wang_labels:
                U_k = W @ U_w
                Ja, Jb, Jc, Jz = _rot_ops_rep(J, rep=rep_int)
                Ja2 = Ja @ Ja
                Jc2 = Jc @ Jc
            for i, E_val in enumerate(E):
                if use_wang_labels:
                    Ka, Kc, tau, species = assign_KaKc_tau_species_from_wang(U_w[:, i], J, rep_user)
                else:
                    Ka, Kc, tau, species, _, _ = assign_KaKc_tau_species_from_expectation(U_k[:, i], J, Ja2, Jc2)
                E_cm = float(E_val) / CMC_SPCAT
                levels_list.append((J, Ka, Kc, float(E_val), f"{E_cm:.6f}", species))
    csv_str = pd.DataFrame(levels_list, columns=['J','Ka','Kc','Energy (MHz)','Energy (cm-1, egy)','Parity'])\
               .sort_values(['J','Ka','Kc','Energy (MHz)']).to_csv(index=False)
csv_str
`);
    downloadFile(levelsCsv, "energies.csv", "text/csv");
      
    await sleep(250);
    const datString = await pyodide.runPythonAsync(`
import pandas as pd
def _label_row(r):
    if 'Fl' in df.columns and 'Fu' in df.columns and 'I12_l' in df.columns and 'I12_u' in df.columns:
        return f"J={int(r.Jl)} F={float(r.Fl):g} I12={float(r.I12_l):g} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} F={float(r.Fu):g} I12={float(r.I12_u):g} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
    if 'Fl' in df.columns and 'Fu' in df.columns:
        return f"J={int(r.Jl)} F={float(r.Fl):g} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} F={float(r.Fu):g} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
    return f"J={int(r.Jl)} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
df_copy = df.copy()
df_copy["Label"] = df_copy.apply(_label_row, axis=1)
df_copy[['Frequency (MHz)','Relative intensity','Label']].to_csv(index=False, header=False, sep=' ')
`);
    downloadFile(datString, "spectrum.dat", "text/plain");
  });

  openCompareButton.addEventListener("click", async () => {
    if (!lastPyvars) {
      window.open("./WMS-FitRot/index.html", "_blank");
      return;
    }
    Object.entries(lastPyvars).forEach(([k, v]) => pyodide.globals.set(k, v));
    try {
      const datString = await pyodide.runPythonAsync(`
import pandas as pd
def _label_row(r):
    if 'Fl' in df.columns and 'Fu' in df.columns and 'I12_l' in df.columns and 'I12_u' in df.columns:
        return f"J={int(r.Jl)} F={float(r.Fl):g} I12={float(r.I12_l):g} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} F={float(r.Fu):g} I12={float(r.I12_u):g} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
    if 'Fl' in df.columns and 'Fu' in df.columns:
        return f"J={int(r.Jl)} F={float(r.Fl):g} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} F={float(r.Fu):g} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
    return f"J={int(r.Jl)} Ka={int(r.Ka_l)} Kc={int(r.Kc_l)} -> J={int(r.Ju)} Ka={int(r.Ka_u)} Kc={int(r.Kc_u)} ({r.Branch})"
df_copy = df.copy()
df_copy["Label"] = df_copy.apply(_label_row, axis=1)
df_copy[['Frequency (MHz)','Relative intensity','Label']].to_csv(index=False, header=False, sep=' ')
`);
      const timestamp = new Date();
      const label = `WMS-Rot spectrum ${timestamp.toISOString().slice(0, 19).replace('T', ' ')}`;
      const queueRaw = localStorage.getItem("fitSpectraQueue");
      const queue = queueRaw ? JSON.parse(queueRaw) : [];
      queue.push({
        name: label,
        mode: "stick",
        role: "theory",
        text: datString
      });
      queue.push({
        name: `${label} input`,
        mode: "input",
        text: buildInputTextFromForm()
      });
      localStorage.setItem("fitSpectraQueue", JSON.stringify(queue));
      window.open("./WMS-FitRot/index.html", "_blank");
    } catch (err) {
      console.error("Unable to export spectrum for WMS-FitRot:", err);
      alert("Unable to send the spectrum to WMS-FitRot.");
    }
  });

  if (openWavefunctionViewerButton) {
    openWavefunctionViewerButton.addEventListener("click", async () => {
      if (!lastPyvars) {
        alert("Please calculate a spectrum first.");
        return;
      }
      Object.entries(lastPyvars).forEach(([k, v]) => pyodide.globals.set(k, v));
      try {
        const payloadJson = await pyodide.runPythonAsync(`
import json, datetime
cache = globals().get("LAST_WAVEFUNC_CACHE", None)
if cache is None:
    raise RuntimeError("Wavefunction cache not available. Recalculate first.")
meta = cache.setdefault("meta", {})
meta["generated_at_utc"] = datetime.datetime.utcnow().isoformat() + "Z"
json.dumps(cache)
`);
        localStorage.setItem("wmsRotWavefunctionPayload", payloadJson);
        window.open("./WMS-Wavefunctions/index.html", "_blank");
      } catch (err) {
        console.error("Unable to open wavefunction viewer:", err);
        alert("Unable to open wavefunction viewer. Try lowering Jmax and recalculating.");
      }
    });
  }

  const fileButton = document.getElementById("load-input-file");
  const fileInput = document.getElementById("input-file");
  const fileName = document.getElementById("input-file-name");
  const fileStatus = document.getElementById("input-file-status");
  if (fileButton && fileInput) {
    fileButton.addEventListener("click", () => fileInput.click());
    fileInput.addEventListener("change", async () => {
      const file = fileInput.files && fileInput.files[0];
      if (!file) return;
      if (fileName) fileName.textContent = file.name;
      if (fileStatus) fileStatus.textContent = "Loading file...";
      try {
        const text = await file.text();
        const raw = parseRotationalSection(text);
        const data = buildFileData(raw);
        loadedFileData = data;
        applyFileDataToForm(data);
        if (fileStatus) fileStatus.textContent = "Loaded data from #ROTATIONAL section.";
      } catch (err) {
        console.error(err);
        if (fileStatus) fileStatus.textContent = "Error while parsing file.";
        alert(String(err));
      }
    });
  }

  // Initial page ready (equations + Pyodide environment).
  try {
    if (window.MathJax && window.MathJax.typesetPromise) {
      globalLoadingOverlay.setMessage("Rendering equations…");
      await window.MathJax.typesetPromise();
    }
    globalLoadingOverlay.hide();
  } catch (err) {
    console.error(err);
    globalLoadingOverlay.showError("Failed to initialize the page. Please reload.");
  }
}

main();
