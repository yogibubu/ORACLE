"use strict";

const REP_MAP = {
  Ir: ["a", "b", "c"],
  IIr: ["b", "c", "a"],
  IIIr: ["c", "a", "b"],
  Il: ["a", "c", "b"],
  IIl: ["b", "a", "c"],
  IIIl: ["c", "b", "a"],
};

const AXIS_INDEX = { a: 0, b: 1, c: 2 };
const GJF_TEMPLATE_PATH = "../G16_Inputs/pipeline.txt";
const GJF_XYZ_PLACEHOLDER = "<INSERT XYZ COORDINATES HERE>";
const GJF_NAME_PLACEHOLDER = "<INSERT NAME HERE>";
const GJF_REPRESENTATION_PLACEHOLDER = "<INSERT REPRESENTATION HERE>";
const STEP1_SMILES_SOURCE_NAME = "smiles_input.xyz";
const ASYM_REDUCTION_HEADER = "Constants in the Asymmetrically reduced Hamiltonian";
const SYM_REDUCTION_HEADER = "Constants in the Symmetrically Reduced Hamiltonian";
const A_SEXTIC_REDUCTION_HEADER = "Constants in the A reduced Hamiltonian";
const S_SEXTIC_REDUCTION_HEADER = "Constants in the S reduced Hamiltonian";
const DIPOLE_HEADER = "Electric dipole moment (input orientation):";
const BDPCS3_PRELOAD_STORAGE_KEY = "bdpcs3_preload_xyz";
const NUMBER_TOKEN_PATTERN =
  "[-+]?(?:\\d+\\.\\d*|\\d*\\.\\d+|\\d+)(?:[DEde][-+]?\\d+)?";
const FLOW_TARGETS = ["geometry", "dvib", "quartic", "sextic", "quadrupole"];
const FLOW_SOURCE_LANE_ORDER = ["geometry", "dvib", "quartic", "sextic", "quadrupole"];
const FLOW_TARGET_ROW_ORDER = ["geometry", "dvib", "quartic", "sextic", "quadrupole"];
const FLOW_CONNECTABLE_TARGET_SET = new Set(FLOW_TARGETS);
const FLOW_TARGET_LABELS = {
  geometry: "Geometry",
  hessian: "Hessian",
  quartic: "Quartic",
  sextic: "Sextic",
  dvib: "DVib",
  quadrupole: "Quadrupole",
};
const CARTESIAN_AXIS_LABELS = ["x", "y", "z"];
const PRINCIPAL_AXIS_LABELS = ["a", "b", "c"];
const INPUT_REP_ALIGNMENT_MIN = 0.92;
const INPUT_REP_MOMENT_DEGENERACY_REL = 1e-5;
const BDPCS3_PAIR_GAP_PX = 78;
const BDPCS3_PAIR_MIN_LANE_WIDTH = 540;
const XYZ_INPUT_DEBOUNCE_MS = 180;
const AMU_TO_AU = 1822.888486209;
const ANGSTROM_TO_BOHR = 1.8897261254535;
const BOHR_TO_ANGSTROM = 1 / ANGSTROM_TO_BOHR;
const AU_TO_DEBYE = 2.541746473;
const CM1_TO_MHZ = 29979.2458;
const PLANCK_CONSTANT_SI = 6.62607015e-34;
const BOLTZMANN_CONSTANT_SI = 1.380649e-23;
const C_AU = 137.035999084;
const TWO_PI = Math.PI * 2;
const HARMONIC_CORIOLIS_ENABLED = true;
const HARMONIC_CORIOLIS_SCALE = 0.5;
const HARMONIC_CORIOLIS_MIN_DENOM = 1e-12;
const CATREF_REDUCTIONS = ["A", "S"];
const DEFAULT_CATREF_REDUCTION = "A";
const DEFAULT_CATREF_TEMPERATURE_K = 4.0;
const DEFAULT_CATREF_STR0 = -12.0;
const DEFAULT_CATREF_STR1 = -10.0;
const DEFAULT_CATREF_FREQ_MAX_GHZ = 20.0;
const DEFAULT_CATREF_FLAGS = 1;
const DEFAULT_CATREF_TAG = 1;
const DEFAULT_CATREF_MIN_JMAX = 25;
const DEFAULT_CATREF_MAX_JMAX = 250;
const PYODIDE_INDEX_URL = "../vendor/pyodide/v0.23.4/full/";
const PYODIDE_CDN_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v0.23.4/full/";
const PYODIDE_CDN_SCRIPT_URL = `${PYODIDE_CDN_INDEX_URL}pyodide.js`;
const PYODIDE_DISTORTION_MODULES = [
  {
    fileName: "wms_prep_representation_utils.py",
    url: "./pyodide_distortion/representation_utils.txt",
  },
  {
    fileName: "wms_prep_distortion_models.py",
    url: "./pyodide_distortion/distortion_models.txt",
  },
  {
    fileName: "wms_prep_bridge.py",
    url: "./pyodide_distortion/wms_prep_bridge.txt",
  },
];
const TAU_COMPONENT_KEYS = [
  "taaaa",
  "tbbbb",
  "tcccc",
  "taabb",
  "tbbcc",
  "tccaa",
  "tabab",
  "tbcbc",
  "tcaca",
];
const SEXTIC_INERTIAL_KEYS = [
  "eaaa",
  "eaab",
  "eaac",
  "eabb",
  "eabc",
  "eacc",
  "ebbb",
  "ebbc",
  "ebcc",
  "eccc",
];
const QUAD_CART_KEYS = ["chi_aa", "chi_bb", "chi_cc", "chi_ab", "chi_ac", "chi_bc"];
const TAU_REP_TO_RIGHT_HANDED = {
  Ir: "Ir",
  Il: "Ir",
  IIr: "IIr",
  IIl: "IIr",
  IIIr: "IIIr",
  IIIl: "IIIr",
};
const TAU_RIGHT_HANDED_AXIS_MAP = {
  Ir: { a: "a", b: "b", c: "c" },
  IIr: { a: "b", b: "c", c: "a" },
  IIIr: { a: "c", b: "a", c: "b" },
};
const TAU_CANONICAL_INDEX_LOOKUP = buildTauCanonicalIndexLookup();
const SEXTIC_CANONICAL_INDEX_LOOKUP = buildSexticCanonicalIndexLookup();
const state = {
  fileName: null,
  sourceComment: "",
  loadedAtoms: null,
  principalAtoms: null,
  analysisData: null,
  outputAtoms: null,
  gjfTemplate: null,
  logFileName: null,
  logFiles: [],
  logEntries: [],
  nextSourceBatchId: 0,
  logSources: [],
  sourceById: {},
  connections: {
    geometry: null,
    dvib: null,
    quartic: null,
    sextic: null,
    quadrupole: null,
  },
  nodePositions: {},
  laneBounds: {},
  activeLinkDrag: null,
  combinedParsedData: null,
  flowSummaryText: "No source loaded.",
  wmsInputText: "",
  catrefInputText: "",
  catrefBundles: [],
  catrefVarText: "",
  catrefIntText: "",
  dpcs3XYZText: "",
  xyzBundles: {},
  fchkBundles: {},
  hessianBundles: {},
  tauBundles: {},
  sexticBundles: {},
  outputRefreshToken: 0,
};

const ui = {};
let xyzInputDebounceTimer = null;
let rdkitModule = null;
let rdkitReadyPromise = null;
let distortionPyodide = null;
let distortionPyodideReady = null;
const loadingOverlay = (() => {
  let visibleDepth = 0;

  function getElements() {
    if (typeof document === "undefined") return null;
    const overlay = document.getElementById("loading-overlay");
    if (!overlay) return null;
    return {
      overlay,
      messageEl: overlay.querySelector('[data-role="loading-message"]'),
      spinnerEl: overlay.querySelector('[data-role="loading-spinner"]'),
    };
  }

  return {
    show(message) {
      const parts = getElements();
      if (!parts) return;
      visibleDepth += 1;
      parts.overlay.classList.remove("loading-overlay--error");
      parts.overlay.classList.add("is-visible");
      parts.overlay.removeAttribute("aria-hidden");
      if (parts.spinnerEl) parts.spinnerEl.hidden = false;
      if (message && parts.messageEl) parts.messageEl.textContent = message;
    },
    hide() {
      const parts = getElements();
      if (!parts) return;
      visibleDepth = Math.max(0, visibleDepth - 1);
      if (visibleDepth > 0) return;
      parts.overlay.classList.remove("is-visible");
      parts.overlay.setAttribute("aria-hidden", "true");
    },
    setMessage(message) {
      const parts = getElements();
      if (!parts || !message || !parts.messageEl) return;
      parts.messageEl.textContent = message;
    },
    showError(message) {
      const parts = getElements();
      if (!parts) return;
      visibleDepth = Math.max(1, visibleDepth);
      parts.overlay.classList.add("loading-overlay--error");
      parts.overlay.classList.add("is-visible");
      parts.overlay.removeAttribute("aria-hidden");
      if (parts.spinnerEl) parts.spinnerEl.hidden = true;
      if (parts.messageEl) {
        parts.messageEl.textContent = message || "Error during loading.";
      }
    },
  };
})();

if (typeof window !== "undefined" && typeof document !== "undefined") {
  window.addEventListener("DOMContentLoaded", init);
}

const STEP1_REPRESENTATION_CHANGE_DISABLED_MESSAGE =
  "Representation changes are disabled in WMS-Prep. XYZ and GJF exports keep the detected input representation.";

function syncStep1OutputRepresentationControl() {
  if (!ui.step1OutputRepresentation) return;
  const lockedRep = normalizeFullWatsonRepresentationToken(getStep1InputRepresentation(), "Ir");
  ui.step1OutputRepresentation.disabled = true;
  ui.step1OutputRepresentation.title = STEP1_REPRESENTATION_CHANGE_DISABLED_MESSAGE;
  if (ui.step1OutputRepresentation.value !== lockedRep) {
    ui.step1OutputRepresentation.value = lockedRep;
  }
}

function syncWMSOutputRepresentationControl() {
  if (!ui.wmsOutputRepresentation) return;
  ui.wmsOutputRepresentation.disabled = false;
  ui.wmsOutputRepresentation.title =
    "Select the final representation used for the WMS-Rot export.";
  if (ui.wmsOutputRepresentation.value !== "Ir" && ui.wmsOutputRepresentation.value !== "IIIl") {
    ui.wmsOutputRepresentation.value = "Ir";
  }
}

function syncRepresentationChangeControls() {
  syncStep1OutputRepresentationControl();
  syncWMSOutputRepresentationControl();
}

function buildInitialWMSOutputWarnings(parsed, targetRepresentation) {
  const lines = [];
  if (!parsed || typeof parsed !== "object") {
    return lines;
  }
  const reps = collectConnectedRepresentations(parsed);
  if (reps.length > 1) {
    lines.push(
      `Connected sources use multiple input representations (${reps.join(", ")}); WMS-Prep will transform them to ${targetRepresentation}.`
    );
  }
  return lines;
}

function init() {
  ui.xyzFile = document.getElementById("xyz-file");
  ui.smilesInput = document.getElementById("smiles-input");
  ui.convertSmilesBtn = document.getElementById("convert-smiles-btn");
  ui.xyzInput = document.getElementById("xyz-input");
  ui.gaussianLogFile = document.getElementById("gaussian-log-file");
  ui.loadLogBtn = document.getElementById("load-log-btn");
  ui.step1DetectedRepresentation = document.getElementById("step1-detected-representation");
  ui.step1OutputRepresentation = document.getElementById("step1-output-representation");
  ui.step2Representation = document.getElementById("step2-representation");
  ui.downloadBtn = document.getElementById("download-btn");
  ui.downloadGjfBtn = document.getElementById("download-gjf-btn");
  ui.downloadBasisSetBtn = document.getElementById("download-basis-set-btn");
  ui.downloadWmsInputBtn = document.getElementById("download-wms-input-btn");
  ui.downloadCatrefInputBtn = document.getElementById("download-catref-input-btn");
  ui.openBdpcs3Btn = document.getElementById("open-bdpcs3-btn");
  ui.clearLogsBtn = document.getElementById("clear-logs-btn");
  ui.wmsOutputRepresentation = document.getElementById("wms-output-representation");
  ui.wmsOutputRepresentationWarning = document.getElementById("wms-output-representation-warning");
  ui.loadStep2XYZBtn = document.getElementById("load-step2-xyz-btn");
  ui.useStep1XYZBtn = document.getElementById("use-step1-xyz-btn");
  ui.step2XYZFile = document.getElementById("step2-xyz-file");
  ui.loadFchkBtn = document.getElementById("load-fchk-btn");
  ui.fchkFile = document.getElementById("fchk-file");
  ui.loadHessianXYZBtn = document.getElementById("load-hessian-xyz-btn");
  ui.hessianXYZFile = document.getElementById("hessian-xyz-file");
  ui.status = document.getElementById("status");
  ui.analysis = document.getElementById("analysis");
  ui.xyzOutput = document.getElementById("xyz-output");
  ui.wmsOutput = document.getElementById("wms-output");
  ui.catrefOutput = document.getElementById("catref-output");
  ui.dpcs3XYZOutput = document.getElementById("dpcs3-xyz-output");
  ui.logSummary = document.getElementById("log-summary");
  ui.loadedLogList = document.getElementById("loaded-log-list");
  ui.flowBoard = document.getElementById("flow-board");
  ui.flowNodes = document.getElementById("flow-nodes");
  ui.flowLanes = document.getElementById("flow-lanes");
  ui.flowLines = document.getElementById("flow-lines");
  ui.flowBoardWrap = ui.flowBoard ? ui.flowBoard.parentElement : null;
  ui.flowNodeModal = document.getElementById("flow-node-modal");
  ui.flowNodeModalTitle = document.getElementById("flow-node-modal-title");
  ui.flowNodeModalBody = document.getElementById("flow-node-modal-body");
  ui.flowNodeModalClose = document.getElementById("flow-node-modal-close");

  if (typeof Molecule === "undefined" || typeof atomicMass === "undefined") {
    setStatus("chem.js did not load correctly.", true);
    ui.downloadBtn.disabled = true;
    ui.downloadGjfBtn.disabled = true;
    ui.downloadBasisSetBtn.disabled = true;
    ui.downloadWmsInputBtn.disabled = true;
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
    if (ui.openBdpcs3Btn) ui.openBdpcs3Btn.disabled = true;
    return;
  }

  ui.xyzFile.addEventListener("change", onFileSelected);
  if (ui.convertSmilesBtn) {
    ui.convertSmilesBtn.addEventListener("click", onConvertSmilesClicked);
  }
  if (ui.xyzInput) {
    ui.xyzInput.addEventListener("input", onXYZInputChanged);
  }
  if (ui.loadLogBtn) {
    ui.loadLogBtn.addEventListener("click", onLoadLogsClicked);
  }
  if (ui.wmsOutputRepresentation) {
    ui.wmsOutputRepresentation.addEventListener("change", () => {
      void handleWMSOutputRepresentationChange();
    });
  }
  ui.gaussianLogFile.addEventListener("change", onGaussianLogSelected);
  if (ui.loadStep2XYZBtn) {
    ui.loadStep2XYZBtn.addEventListener("click", onLoadStep2XYZClicked);
  }
  if (ui.useStep1XYZBtn) {
    ui.useStep1XYZBtn.addEventListener("click", onUseStep1XYZClicked);
  }
  if (ui.step2XYZFile) {
    ui.step2XYZFile.addEventListener("change", onStep2XYZSelected);
  }
  if (ui.loadFchkBtn) {
    ui.loadFchkBtn.addEventListener("click", onLoadFchkClicked);
  }
  if (ui.fchkFile) {
    ui.fchkFile.addEventListener("change", onFchkSelected);
  }
  if (ui.loadHessianXYZBtn) {
    ui.loadHessianXYZBtn.addEventListener("click", onLoadHessianXYZClicked);
  }
  if (ui.hessianXYZFile) {
    ui.hessianXYZFile.addEventListener("change", onHessianXYZSelected);
  }
  if (ui.step1OutputRepresentation) {
    ui.step1OutputRepresentation.addEventListener("change", () => {
      if (state.principalAtoms) {
        applyRepresentationAndRefresh();
      }
    });
  }
  ui.downloadBtn.addEventListener("click", downloadOutputXYZ);
  ui.downloadGjfBtn.addEventListener("click", onDownloadGJFClicked);
  ui.downloadBasisSetBtn.addEventListener("click",
      onDownloadBasisSetClicked);
  ui.downloadWmsInputBtn.addEventListener("click", downloadWMSInput);
  if (ui.downloadCatrefInputBtn) {
    ui.downloadCatrefInputBtn.addEventListener("click", downloadCATREFInput);
  }
  if (ui.openBdpcs3Btn) {
    ui.openBdpcs3Btn.addEventListener("click", onOpenInBDPCS3Clicked);
  }
  ui.clearLogsBtn.addEventListener("click", () => {
    if (ui.gaussianLogFile) ui.gaussianLogFile.value = "";
    if (ui.step2XYZFile) ui.step2XYZFile.value = "";
    if (ui.fchkFile) ui.fchkFile.value = "";
    if (ui.hessianXYZFile) ui.hessianXYZFile.value = "";
    clearGaussianLogState();
    setStatus("Loaded sources cleared.");
  });
  window.addEventListener("resize", scheduleFlowLineRedraw);
  if (ui.flowBoardWrap) {
    ui.flowBoardWrap.addEventListener("scroll", scheduleFlowLineRedraw);
  }
  if (ui.flowNodeModal) {
    ui.flowNodeModal.addEventListener("click", (event) => {
      if (event.target === ui.flowNodeModal) {
        closeFlowNodeModal();
      }
    });
  }
  if (ui.flowNodeModalClose) {
    ui.flowNodeModalClose.addEventListener("click", closeFlowNodeModal);
  }
  document.addEventListener("keydown", onGlobalKeyDown);

  renderFlowBoard();
  renderLoadedLogList();
  updateFlowSummary();
  syncRepresentationChangeControls();
  setStatus("Ready. Load XYZ or SMILES for geometry inspection, or load Gaussian logs and sources to generate WMS-Rot input in Ir or IIIl.");
}

function waitForNextPaint() {
  if (typeof window === "undefined" || typeof window.requestAnimationFrame !== "function") {
    return Promise.resolve();
  }
  return new Promise((resolve) => {
    window.requestAnimationFrame(() => {
      window.requestAnimationFrame(resolve);
    });
  });
}

async function handleWMSOutputRepresentationChange() {
  renderFlowBoard();
  loadingOverlay.show("Refreshing outputs…");
  try {
    await waitForNextPaint();
    await refreshWMSInputOutput();
  } finally {
    loadingOverlay.hide();
  }
}

function getStep1InputRepresentation() {
  const detection = state.analysisData && state.analysisData.inputRepresentationDetection
    ? state.analysisData.inputRepresentationDetection
    : null;
  return detection && detection.representation ? detection.representation : null;
}

function getStep1OutputRepresentation() {
  syncStep1OutputRepresentationControl();
  return normalizeFullWatsonRepresentationToken(getStep1InputRepresentation(), "Ir");
}

function getStep2Representation() {
  return ui.step2Representation ? ui.step2Representation.value || "Ir" : "Ir";
}

function onFileSelected(event) {
  const file = event.target.files && event.target.files[0];
  if (!file) return;

  const reader = new FileReader();
  reader.onload = () => {
    try {
      const text = String(reader.result || "");
      if (ui.xyzInput) {
        ui.xyzInput.value = text;
      }
      const parsed = loadXYZFromText(text, file.name);
      const detectedRep = getStep1InputRepresentation();
      const repSuffix = detectedRep
        ? ` Detected input representation: ${detectedRep}.`
        : " Input representation could not be determined.";
      setStatus(`Loaded ${parsed.atoms.length} atoms from ${file.name}.${repSuffix}`);
    } catch (error) {
      clearState();
      setStatus(error instanceof Error ? error.message : String(error), true);
    }
  };
  reader.onerror = () => {
    clearState();
    setStatus("Error while reading the file.", true);
  };
  reader.readAsText(file);
}

function onXYZInputChanged() {
  if (typeof window === "undefined") return;
  if (xyzInputDebounceTimer !== null) {
    window.clearTimeout(xyzInputDebounceTimer);
  }
  xyzInputDebounceTimer = window.setTimeout(() => {
    xyzInputDebounceTimer = null;
    const text = ui.xyzInput ? String(ui.xyzInput.value || "") : "";
    if (!text.trim()) {
      clearState();
      setStatus("XYZ input cleared.");
      return;
    }
    try {
      const parsed = loadXYZFromText(text, "Input XYZ");
      const detectedRep = getStep1InputRepresentation();
      const repSuffix = detectedRep
        ? ` Detected input representation: ${detectedRep}.`
        : " Input representation could not be determined.";
      setStatus(`Parsed ${parsed.atoms.length} atoms from Input XYZ.${repSuffix}`);
    } catch (error) {
      clearState();
      setStatus(error instanceof Error ? error.message : String(error), true);
    }
  }, XYZ_INPUT_DEBOUNCE_MS);
}

async function onConvertSmilesClicked() {
  const smiles = ui.smilesInput ? String(ui.smilesInput.value || "").trim() : "";
  if (!smiles) {
    setStatus("SMILES input is empty.", true);
    return;
  }

  if (ui.convertSmilesBtn) {
    ui.convertSmilesBtn.disabled = true;
  }

  try {
    const converted = await convertSmilesToXYZ(smiles);
    if (ui.xyzInput) {
      ui.xyzInput.value = converted.xyzText;
    }
    const parsed = loadXYZFromText(converted.xyzText, STEP1_SMILES_SOURCE_NAME);
    const detectedRep = getStep1InputRepresentation();
    const repSuffix = detectedRep
      ? ` Detected input representation: ${detectedRep}.`
      : " Input representation could not be determined.";
    const uffSuffix = converted.uffApplied
      ? " UFF optimization applied."
      : " UFF optimization unavailable in this RDKit build; using RDKit-generated coordinates.";
    setStatus(`Converted SMILES to XYZ (${parsed.atoms.length} atoms).${repSuffix}${uffSuffix}`);
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  } finally {
    if (ui.convertSmilesBtn) {
      ui.convertSmilesBtn.disabled = false;
    }
  }
}

function ensureRDKitReady() {
  if (rdkitModule) {
    return Promise.resolve(rdkitModule);
  }
  if (rdkitReadyPromise) {
    return rdkitReadyPromise;
  }
  if (
    typeof window === "undefined" ||
    typeof window.initRDKitModule !== "function"
  ) {
    return Promise.reject(
      new Error("RDKit is not available on this page. Unable to process SMILES.")
    );
  }
  rdkitReadyPromise = window
    .initRDKitModule()
    .then((module) => {
      if (!module || typeof module.get_mol !== "function") {
        throw new Error("RDKit module loaded, but get_mol is unavailable.");
      }
      rdkitModule = module;
      return module;
    })
    .catch((error) => {
      rdkitReadyPromise = null;
      throw new Error(
        `Unable to initialize RDKit: ${error instanceof Error ? error.message : String(error)}`
      );
    });
  return rdkitReadyPromise;
}

async function convertSmilesToXYZ(smiles) {
  const rdkit = await ensureRDKitReady();
  const molsToDelete = [];
  const registerMol = (mol) => {
    if (mol && typeof mol.delete === "function") {
      molsToDelete.push(mol);
    }
    return mol;
  };

  try {
    let workingMol = registerMol(rdkit.get_mol(smiles));
    if (!workingMol) {
      throw new Error("Invalid SMILES string. RDKit could not build a molecule.");
    }

    const withHsBlock =
      typeof workingMol.add_hs === "function" ? String(workingMol.add_hs() || "") : "";
    if (withHsBlock.trim()) {
      const withHsMol = registerMol(rdkit.get_mol(withHsBlock));
      if (withHsMol) {
        workingMol = withHsMol;
      }
    }

    const recoordMol = buildMolWithGeneratedCoords(rdkit, workingMol);
    if (recoordMol) {
      workingMol = registerMol(recoordMol);
    }

    const optimized = tryApplyUFFOptimization(rdkit, workingMol);
    if (optimized.mol && optimized.mol !== workingMol) {
      workingMol = registerMol(optimized.mol);
    }

    const molBlock =
      typeof workingMol.get_molblock === "function" ? String(workingMol.get_molblock() || "") : "";
    if (!molBlock.trim()) {
      throw new Error("RDKit returned an empty MolBlock while converting SMILES to XYZ.");
    }

    const atoms = parseMolBlockAtoms(molBlock);
    if (!atoms.length) {
      throw new Error("RDKit did not provide atom coordinates for the supplied SMILES.");
    }

    const note = optimized.uffApplied ? "RDKit SMILES conversion + UFF" : "RDKit SMILES conversion";
    const xyzText = buildXYZTextFromAtomList(atoms, note);
    return {
      xyzText,
      uffApplied: optimized.uffApplied,
    };
  } finally {
    for (let i = molsToDelete.length - 1; i >= 0; i -= 1) {
      try {
        molsToDelete[i].delete();
      } catch (_) {
        // no-op
      }
    }
  }
}

function buildMolWithGeneratedCoords(rdkit, mol) {
  if (!mol) return null;
  if (typeof mol.get_new_coords === "function") {
    try {
      const block = String(mol.get_new_coords(true) || "");
      if (block.trim()) {
        const newMol = rdkit.get_mol(block);
        if (newMol) return newMol;
      }
    } catch (_) {
      // fallback below
    }
  }
  if (typeof mol.set_new_coords === "function") {
    try {
      mol.set_new_coords(true);
    } catch (_) {
      // keep original geometry if coord generation fails
    }
  }
  return null;
}

function tryApplyUFFOptimization(rdkit, mol) {
  if (!mol) return { mol, uffApplied: false };

  const methodCandidates = [
    { target: mol, method: "uff_optimize", args: [] },
    { target: mol, method: "UFFOptimizeMolecule", args: [] },
    { target: mol, method: "optimize_uff", args: [] },
    { target: rdkit, method: "uff_optimize", args: [mol] },
    { target: rdkit, method: "UFFOptimizeMolecule", args: [mol] },
    { target: rdkit, method: "optimize_uff", args: [mol] },
  ];

  for (let i = 0; i < methodCandidates.length; i += 1) {
    const candidate = methodCandidates[i];
    const fn = candidate.target && candidate.target[candidate.method];
    if (typeof fn !== "function") continue;
    try {
      const result = fn.apply(candidate.target, candidate.args);
      if (typeof result === "string" && result.trim()) {
        const optimizedMol = rdkit.get_mol(result);
        if (optimizedMol) {
          return { mol: optimizedMol, uffApplied: true };
        }
      }
      return { mol, uffApplied: true };
    } catch (_) {
      // keep searching for another available signature
    }
  }
  return { mol, uffApplied: false };
}

function parseMolBlockAtoms(molBlockText) {
  const lines = String(molBlockText || "").split(/\r?\n/);
  if (lines.length < 4) {
    throw new Error("Invalid MolBlock received from RDKit.");
  }
  const countsLine = String(lines[3] || "");
  if (/V3000/i.test(countsLine)) {
    return parseV3000MolBlockAtoms(lines);
  }
  return parseV2000MolBlockAtoms(lines);
}

function parseV2000MolBlockAtoms(lines) {
  const countsLine = String(lines[3] || "");
  const atomCount = Number.parseInt(countsLine.slice(0, 3).trim(), 10);
  if (!Number.isFinite(atomCount) || atomCount <= 0) {
    throw new Error("RDKit MolBlock has an invalid V2000 atom count.");
  }
  const atoms = [];
  for (let i = 0; i < atomCount; i += 1) {
    const line = String(lines[4 + i] || "");
    const x = parseRDKitCoordinateToken(line.slice(0, 10));
    const y = parseRDKitCoordinateToken(line.slice(10, 20));
    const z = parseRDKitCoordinateToken(line.slice(20, 30));
    const label = normalizeElementLabel(line.slice(31, 34));
    atoms.push({ label, coords: [x, y, z] });
  }
  return atoms;
}

function parseV3000MolBlockAtoms(lines) {
  const atoms = [];
  let inAtomSection = false;
  for (let i = 0; i < lines.length; i += 1) {
    const raw = String(lines[i] || "").trim();
    if (/^M\s+V30\s+BEGIN\s+ATOM$/i.test(raw)) {
      inAtomSection = true;
      continue;
    }
    if (/^M\s+V30\s+END\s+ATOM$/i.test(raw)) {
      inAtomSection = false;
      continue;
    }
    if (!inAtomSection) continue;
    const match = raw.match(
      /^M\s+V30\s+\d+\s+([A-Za-z][A-Za-z]?)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\b/
    );
    if (!match) continue;
    atoms.push({
      label: normalizeElementLabel(match[1]),
      coords: [
        parseRDKitCoordinateToken(match[2]),
        parseRDKitCoordinateToken(match[3]),
        parseRDKitCoordinateToken(match[4]),
      ],
    });
  }
  if (!atoms.length) {
    throw new Error("RDKit MolBlock has an empty V3000 atom section.");
  }
  return atoms;
}

function parseRDKitCoordinateToken(value) {
  const parsed = Number.parseFloat(String(value || "").replace(/D/i, "E").trim());
  if (!Number.isFinite(parsed)) {
    throw new Error(`Invalid coordinate token from RDKit MolBlock: "${value}"`);
  }
  return parsed;
}

function normalizeElementLabel(value) {
  const token = String(value || "").trim();
  if (!token) return "X";
  if (token.length === 1) return token.toUpperCase();
  return `${token[0].toUpperCase()}${token.slice(1).toLowerCase()}`;
}

function buildXYZTextFromAtomList(atoms, comment) {
  const lines = [];
  lines.push(String(atoms.length));
  lines.push(String(comment || "RDKit SMILES conversion").trim());
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    lines.push(
      `${atom.label} ${formatNumber(atom.coords[0])} ${formatNumber(atom.coords[1])} ${formatNumber(atom.coords[2])}`
    );
  }
  return `${lines.join("\n")}\n`;
}

function onLoadStep2XYZClicked() {
  if (!ui.step2XYZFile) return;
  ui.step2XYZFile.click();
}

function onLoadLogsClicked() {
  if (!ui.gaussianLogFile) return;
  ui.gaussianLogFile.click();
}

function onUseStep1XYZClicked() {
  const text = ui.xyzOutput ? String(ui.xyzOutput.value || "") : "";
  if (!text.trim()) {
    setStatus("Pre-processing output XYZ is empty. Generate it first, then retry.", true);
    return;
  }
  try {
    const parsed = parseXYZ(text);
    const created = addXYZGeometrySource({
      fileName: "PreProcessing_Output_XYZ",
      atoms: parsed.atoms,
    });
    setStatus(
      `Added Geometry source from pre-processing output (detected representation ${created.representation}; dipole defaults 1/1/1).`
    );
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  }
}

function onStep2XYZSelected(event) {
  const file = event?.target?.files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const text = String(reader.result || "");
      const parsed = parseXYZ(text);
      const created = addXYZGeometrySource({
        fileName: file.name,
        atoms: parsed.atoms,
      });
      setStatus(
        `Added XYZ Geometry source "${file.name}" (detected representation ${created.representation}; dipole defaults 1/1/1).`
      );
    } catch (error) {
      setStatus(error instanceof Error ? error.message : String(error), true);
    } finally {
      if (ui.step2XYZFile) {
        ui.step2XYZFile.value = "";
      }
    }
  };
  reader.onerror = () => {
    setStatus(`Error while reading "${file.name}".`, true);
    if (ui.step2XYZFile) {
      ui.step2XYZFile.value = "";
    }
  };
  reader.readAsText(file);
}

function onLoadFchkClicked() {
  if (!ui.fchkFile) return;
  ui.fchkFile.click();
}

function onLoadHessianXYZClicked() {
  if (!ui.hessianXYZFile) return;
  ui.hessianXYZFile.click();
}

function onFchkSelected(event) {
  const file = event?.target?.files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const text = String(reader.result || "");
      const geometry = parseFchkGeometry(text);
      const hessian = parseFchkCartesianForceConstants(text);
      const dipoleMomentAu = parseFchkDipoleMomentOptional(text);
      const created = addFchkSources({
        fileName: file.name,
        atoms: geometry.atoms,
        hessianMatrix: hessian.matrix,
        lowerTriangularCount: hessian.lowerTriangularCount,
        dipoleMomentAu,
      });
      setStatus(
        `Added FCHK source "${file.name}" with Geometry + Quartic (rotated ${created.representation}).`
      );
    } catch (error) {
      setStatus(error instanceof Error ? error.message : String(error), true);
    } finally {
      if (ui.fchkFile) {
        ui.fchkFile.value = "";
      }
    }
  };
  reader.onerror = () => {
    setStatus(`Error while reading "${file.name}".`, true);
    if (ui.fchkFile) {
      ui.fchkFile.value = "";
    }
  };
  reader.readAsText(file);
}

function onHessianXYZSelected(event) {
  const file = event?.target?.files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const text = String(reader.result || "");
      const parsed = parseHessianXYZBundleText(text);
      const created = addHessianXYZSources({
        fileName: file.name,
        atoms: parsed.atoms,
        hessianMatrix: parsed.hessianMatrix,
        pointGroup: parsed.pointGroup,
        dipoleVectorDebye: parsed.dipoleVectorDebye,
      });
      const metadataTags = [];
      if (parsed.pointGroup) metadataTags.push(`point group ${parsed.pointGroup}`);
      if (Array.isArray(parsed.dipoleVectorDebye)) metadataTags.push("dipole");
      const metadataSuffix = metadataTags.length ? `; metadata: ${metadataTags.join(", ")}` : "";
      setStatus(
        `Added Hessian+XYZ source "${file.name}" with Geometry + Hessian + Quartic (rotated ${created.representation}${metadataSuffix}).`
      );
    } catch (error) {
      setStatus(error instanceof Error ? error.message : String(error), true);
    } finally {
      if (ui.hessianXYZFile) {
        ui.hessianXYZFile.value = "";
      }
    }
  };
  reader.onerror = () => {
    setStatus(`Error while reading "${file.name}".`, true);
    if (ui.hessianXYZFile) {
      ui.hessianXYZFile.value = "";
    }
  };
  reader.readAsText(file);
}

function addXYZGeometrySource({ fileName, atoms }) {
  const batchId = state.nextSourceBatchId;
  state.nextSourceBatchId += 1;
  const bundleId = `xyz-${batchId}`;

  state.xyzBundles[bundleId] = {
    id: bundleId,
    fileName,
    atoms: cloneAtoms(atoms),
  };

  const geometryPayload = buildXYZDerivedGeometryPayload(state.xyzBundles[bundleId]);
  const rep = normalizeWatsonRepresentationToken(geometryPayload.inputRepresentation);
  const geometrySource = {
    id: `${batchId}-geometry-1`,
    type: "geometry",
    fileName,
    title: `${fileName} - Geometry [XYZ${geometryPayload.inputRepresentation ? ` ${geometryPayload.inputRepresentation}` : ""}]`,
    details: summarizeGeometryComponent(geometryPayload),
    payload: geometryPayload,
    defaultForType: true,
    contextLabel: geometryPayload.contextLabel || "",
    sourceIndex: geometryPayload.sourceIndex || Date.now(),
    xyzBundleId: bundleId,
    xyzRole: "geometry",
    originType: "B",
  };

  const entry = {
    id: `xyz-${batchId}`,
    fileName,
    hasCompletePipeline: false,
    componentsLabel: "Geometry x1 (XYZ input)",
  };
  state.logEntries.push(entry);
  state.logSources.push(geometrySource);
  state.sourceById[geometrySource.id] = geometrySource;
  if (!state.connections.geometry) {
    state.connections.geometry = geometrySource.id;
  }
  state.logFileName = fileName;
  state.flowSummaryText = buildFlowSummaryFromEntries();
  renderLoadedLogList();
  renderFlowBoard();
  void refreshWMSInputOutput();
  return {
    source: geometrySource,
    representation: rep || "undetermined",
  };
}

function buildXYZDerivedGeometryPayload(bundle) {
  const rawAtoms = normalizeAtomsForQuartic(bundle.atoms);
  if (!rawAtoms.length) {
    throw new Error(`Unable to build XYZ geometry payload for "${bundle.fileName}".`);
  }
  const detectedRepresentation = detectGeometryRepresentationFromAtoms(rawAtoms);
  const rotationalConstantsMHz = estimateRotationalConstantsMHz(rawAtoms);
  const repLabel = detectedRepresentation || "undetermined";
  const sourceIndex = Date.now();
  return {
    contextLabel: `XYZ geometry (loaded as-is; detected rep ${repLabel}; dipole default 1/1/1)`,
    sourceIndex,
    inputRepresentation: detectedRepresentation,
    pointGroup: "C1",
    rotationalConstantsMHz,
    dipole: {
      a: toFortranScientific(1, 10),
      b: toFortranScientific(1, 10),
      c: toFortranScientific(1, 10),
    },
    dpcs3Atoms: rawAtoms.map((atom) => ({
      z: atom.z,
      coords: [atom.coords[0], atom.coords[1], atom.coords[2]],
    })),
    dpcs3XYZText: buildAtomicNumberXYZText(rawAtoms),
  };
}

function refreshXYZSourcesForRepresentation() {
  const bundleIds = Object.keys(state.xyzBundles || {});
  if (!bundleIds.length) return;
  for (let i = 0; i < bundleIds.length; i += 1) {
    const bundleId = bundleIds[i];
    const bundle = state.xyzBundles[bundleId];
    if (!bundle) continue;
    let geometryPayload;
    try {
      geometryPayload = buildXYZDerivedGeometryPayload(bundle);
    } catch (error) {
      continue;
    }
    const sourceIds = Object.keys(state.sourceById);
    for (let j = 0; j < sourceIds.length; j += 1) {
      const source = state.sourceById[sourceIds[j]];
      if (!source || source.xyzBundleId !== bundleId) continue;
      if (source.xyzRole === "geometry") {
        source.payload = geometryPayload;
        source.details = summarizeGeometryComponent(geometryPayload);
        source.title = `${bundle.fileName} - Geometry [XYZ${geometryPayload.inputRepresentation ? ` ${geometryPayload.inputRepresentation}` : ""}]`;
        source.contextLabel = geometryPayload.contextLabel || "";
      }
    }
  }
}

function addHessianXYZSources({
  fileName,
  atoms,
  hessianMatrix,
  pointGroup = null,
  dipoleVectorDebye = null,
}) {
  const rep = getStep2Representation();
  const batchId = state.nextSourceBatchId;
  state.nextSourceBatchId += 1;
  const bundleId = `hessian-${batchId}`;

  state.hessianBundles[bundleId] = {
    id: bundleId,
    fileName,
    atoms: cloneAtoms(atoms),
    hessianMatrix: cloneMatrix(hessianMatrix),
    pointGroup: pointGroup ? String(pointGroup) : null,
    dipoleVectorDebye: Array.isArray(dipoleVectorDebye)
      ? dipoleVectorDebye.slice(0, 3).map((value) => Number(value))
      : null,
  };

  const derived = buildHessianDerivedPayloads(state.hessianBundles[bundleId], rep);
  const geometrySource = {
    id: `${batchId}-geometry-1`,
    type: "geometry",
    fileName,
    title: `${fileName} - Geometry [Hessian+XYZ rotated ${rep}]`,
    details: summarizeGeometryComponent(derived.geometry),
    payload: derived.geometry,
    defaultForType: true,
    contextLabel: derived.geometry.contextLabel || "",
    sourceIndex: derived.geometry.sourceIndex || Date.now(),
    hessianBundleId: bundleId,
    hessianRole: "geometry",
    originType: "D",
  };
  const hessianSource = {
    id: `${batchId}-hessian-1`,
    type: "hessian",
    fileName,
    title: `${fileName} - Hessian [matrix]`,
    details: summarizeHessianComponent(derived.hessian),
    payload: derived.hessian,
    defaultForType: false,
    contextLabel: derived.hessian.contextLabel || "",
    sourceIndex: derived.hessian.sourceIndex || Date.now(),
    hessianBundleId: bundleId,
    hessianRole: "hessian",
    originType: "D",
  };
  const quarticSource = {
    id: `${batchId}-quartic-1`,
    type: "quartic",
    fileName,
    title: `${fileName} - Quartic [harmonic rotated ${rep}]`,
    details: summarizeQuarticComponent(derived.quartic),
    payload: derived.quartic,
    defaultForType: true,
    contextLabel: derived.quartic.contextLabel || "",
    sourceIndex: derived.quartic.sourceIndex || Date.now(),
    hessianBundleId: bundleId,
    hessianRole: "quartic",
    originType: "D",
  };

  const createdSources = [geometrySource, hessianSource, quarticSource];
  const entry = {
    id: `hessian-${batchId}`,
    fileName,
    hasCompletePipeline: false,
    componentsLabel: "Geometry x1, Hessian x1, Quartic x1 (Hessian+XYZ)",
  };

  state.logEntries.push(entry);
  state.logSources.push(...createdSources);
  for (let i = 0; i < createdSources.length; i += 1) {
    state.sourceById[createdSources[i].id] = createdSources[i];
  }
  if (!state.connections.geometry) {
    state.connections.geometry = geometrySource.id;
  }
  if (!state.connections.quartic) {
    state.connections.quartic = quarticSource.id;
  }
  state.logFileName = fileName;
  state.flowSummaryText = buildFlowSummaryFromEntries();
  renderLoadedLogList();
  renderFlowBoard();
  void refreshWMSInputOutput();

  return {
    sources: createdSources,
    representation: rep,
  };
}

function buildHessianDerivedPayloads(bundle, representation) {
  const rawAtoms = normalizeAtomsForQuartic(bundle.atoms);
  if (!rawAtoms.length) {
    throw new Error(`Unable to build Hessian+XYZ payload for "${bundle.fileName}".`);
  }
  const oriented = orientGeometryForRepresentation(rawAtoms, representation);
  const rotatedAtoms = oriented.atoms;
  const rotationalConstantsMHz = estimateRotationalConstantsMHz(rotatedAtoms);
  const pointGroup = normalizePointGroupToken(bundle.pointGroup) || "C1";
  const dipoleFromFile = formatInputDipoleForRepresentation(
    bundle.dipoleVectorDebye,
    oriented.rotation3
  );
  const dipole =
    dipoleFromFile ||
    {
      a: toFortranScientific(1, 10),
      b: toFortranScientific(1, 10),
      c: toFortranScientific(1, 10),
    };
  const dipoleTag = dipoleFromFile ? "dipole from file" : "dipole default 1/1/1";
  const sourceIndex = Date.now();
  const geometryVariant = {
    contextLabel: `Hessian+XYZ geometry (rotated ${representation}; ${dipoleTag})`,
    sourceIndex,
    inputRepresentation: normalizeWatsonRepresentationToken(representation) || "Ir",
    pointGroup,
    rotationalConstantsMHz,
    dipole,
    dpcs3Atoms: rotatedAtoms.map((atom) => ({
      z: atom.z,
      coords: [atom.coords[0], atom.coords[1], atom.coords[2]],
    })),
    dpcs3XYZText: buildAtomicNumberXYZText(rotatedAtoms),
  };

  const harmonic = computeHarmonicQuarticFromHessian({
    hessianCartesian: bundle.hessianMatrix,
    geometryAtoms: rawAtoms,
    representation,
  });
  const coriolisLabel =
    harmonic && harmonic.coriolisModel && harmonic.coriolisModel.enabled
      ? " + Coriolis(harmonic)"
      : "";
  const quarticVariant = {
    contextLabel: `Hessian+XYZ quartic harmonic${coriolisLabel} (rotated ${representation})`,
    sourceIndex: sourceIndex + 1,
    asym: { ...harmonic.asymFormatted },
    sym: { ...harmonic.symFormatted },
    inputRepresentation: normalizeWatsonRepresentationToken(representation) || "Ir",
  };

  const hessianVariant = {
    contextLabel: "Input Hessian matrix",
    sourceIndex: sourceIndex + 2,
    atomCount: rawAtoms.length,
    matrixSize: Array.isArray(bundle.hessianMatrix) ? bundle.hessianMatrix.length : 0,
    symmetric: isApproximatelySymmetricMatrix(bundle.hessianMatrix),
    trace: matrixTrace(bundle.hessianMatrix),
    frobeniusNorm: matrixFrobeniusNorm(bundle.hessianMatrix),
  };

  return {
    geometry: geometryVariant,
    hessian: hessianVariant,
    quartic: quarticVariant,
  };
}

function refreshHessianSourcesForRepresentation() {
  const rep = getStep2Representation();
  const bundleIds = Object.keys(state.hessianBundles || {});
  if (!bundleIds.length) return;
  for (let i = 0; i < bundleIds.length; i += 1) {
    const bundleId = bundleIds[i];
    const bundle = state.hessianBundles[bundleId];
    if (!bundle) continue;
    let derived;
    try {
      derived = buildHessianDerivedPayloads(bundle, rep);
    } catch (error) {
      continue;
    }
    const sourceIds = Object.keys(state.sourceById);
    for (let j = 0; j < sourceIds.length; j += 1) {
      const source = state.sourceById[sourceIds[j]];
      if (!source || source.hessianBundleId !== bundleId) continue;
      if (source.hessianRole === "geometry") {
        source.payload = derived.geometry;
        source.details = summarizeGeometryComponent(derived.geometry);
        source.title = `${bundle.fileName} - Geometry [Hessian+XYZ rotated ${rep}]`;
        source.contextLabel = derived.geometry.contextLabel || "";
      } else if (source.hessianRole === "hessian") {
        source.payload = derived.hessian;
        source.details = summarizeHessianComponent(derived.hessian);
        source.title = `${bundle.fileName} - Hessian [matrix]`;
        source.contextLabel = derived.hessian.contextLabel || "";
      } else if (source.hessianRole === "quartic") {
        source.payload = derived.quartic;
        source.details = summarizeQuarticComponent(derived.quartic);
        source.title = `${bundle.fileName} - Quartic [harmonic rotated ${rep}]`;
        source.contextLabel = derived.quartic.contextLabel || "";
      }
    }
  }
}

function addFchkSources({ fileName, atoms, hessianMatrix, lowerTriangularCount, dipoleMomentAu = null }) {
  const rep = getStep2Representation();
  const batchId = state.nextSourceBatchId;
  state.nextSourceBatchId += 1;
  const bundleId = `fchk-${batchId}`;

  state.fchkBundles[bundleId] = {
    id: bundleId,
    fileName,
    atoms: cloneAtoms(atoms),
    hessianMatrix: cloneMatrix(hessianMatrix),
    lowerTriangularCount,
    dipoleMomentAu: Array.isArray(dipoleMomentAu) ? dipoleMomentAu.slice(0, 3) : null,
  };

  const derived = buildFchkDerivedPayloads(state.fchkBundles[bundleId], rep);
  const geometrySource = {
    id: `${batchId}-geometry-1`,
    type: "geometry",
    fileName,
    title: `${fileName} - Geometry [FCHK rotated ${rep}]`,
    details: summarizeGeometryComponent(derived.geometry),
    payload: derived.geometry,
    defaultForType: true,
    contextLabel: derived.geometry.contextLabel || "",
    sourceIndex: derived.geometry.sourceIndex || Date.now(),
    fchkBundleId: bundleId,
    fchkRole: "geometry",
    originType: "C",
  };
  const quarticSource = {
    id: `${batchId}-quartic-1`,
    type: "quartic",
    fileName,
    title: `${fileName} - Quartic [FCHK rotated ${rep}]`,
    details: summarizeQuarticComponent(derived.quartic),
    payload: derived.quartic,
    defaultForType: true,
    contextLabel: derived.quartic.contextLabel || "",
    sourceIndex: derived.quartic.sourceIndex || Date.now(),
    fchkBundleId: bundleId,
    fchkRole: "quartic",
    originType: "C",
  };

  state.logFileName = fileName;
  state.logFiles.push(fileName);
  const createdSources = [geometrySource, quarticSource];
  const geometryCount = createdSources.filter((item) => item.type === "geometry").length;
  const quarticCount = createdSources.filter((item) => item.type === "quartic").length;
  state.logEntries.push({
    id: `log-${batchId}`,
    fileName,
    hasCompletePipeline: false,
    componentsLabel: `Geometry x${geometryCount}, Quartic x${quarticCount} (FCHK rotated)`,
  });
  state.logSources.push(...createdSources);
  for (let i = 0; i < createdSources.length; i += 1) {
    state.sourceById[createdSources[i].id] = createdSources[i];
  }
  state.flowSummaryText = buildFlowSummaryFromEntries();

  if (!state.connections.geometry) {
    state.connections.geometry = geometrySource.id;
  }
  if (!state.connections.quartic) {
    state.connections.quartic = quarticSource.id;
  }

  renderLoadedLogList();
  renderFlowBoard();
  void refreshWMSInputOutput();
  return { representation: rep };
}

function addTauSources({ fileName, tauText }) {
  const parsedTau = parseTauPayloadText(tauText);
  const batchId = state.nextSourceBatchId;
  state.nextSourceBatchId += 1;
  const bundleId = `tau-${batchId}`;

  state.tauBundles[bundleId] = {
    id: bundleId,
    fileName,
    inputRepresentation: parsedTau.inputRepresentation,
    tau: { ...parsedTau.tau },
    rotationalConstantsMHz: parsedTau.rotationalConstantsMHz
      ? { ...parsedTau.rotationalConstantsMHz }
      : null,
    moments: parsedTau.moments ? { ...parsedTau.moments } : null,
    sigma: Number.isFinite(parsedTau.sigma) ? parsedTau.sigma : null,
    sigmaSource: parsedTau.sigmaSource || "none",
  };

  const sourceRep = normalizeWatsonRepresentationToken(parsedTau.inputRepresentation) || "Ir";
  const derived = buildTauDerivedPayloads(state.tauBundles[bundleId], sourceRep);
  const quarticSource = {
    id: `${batchId}-quartic-1`,
    type: "quartic",
    fileName,
    title: `${fileName} - Quartic [Tau input ${sourceRep}]`,
    details: summarizeQuarticComponent(derived.quartic),
    payload: derived.quartic,
    defaultForType: true,
    contextLabel: derived.quartic.contextLabel || "",
    sourceIndex: derived.quartic.sourceIndex || Date.now(),
    tauBundleId: bundleId,
    tauRole: "quartic",
  };

  state.logFileName = fileName;
  state.logFiles.push(fileName);
  state.logEntries.push({
    id: `log-${batchId}`,
    fileName,
    hasCompletePipeline: false,
    componentsLabel: `Quartic x1 (tau input ${sourceRep})`,
  });
  state.logSources.push(quarticSource);
  state.sourceById[quarticSource.id] = quarticSource;
  state.flowSummaryText = buildFlowSummaryFromEntries();

  if (!state.connections.quartic) {
    state.connections.quartic = quarticSource.id;
  }

  renderLoadedLogList();
  renderFlowBoard();
  void refreshWMSInputOutput();
  return {
    representation: sourceRep,
    inputRepresentation: parsedTau.inputRepresentation,
  };
}

function addSexticSources({ fileName, sexticText }) {
  const parsedSextic = parseSexticPayloadText(sexticText);
  const batchId = state.nextSourceBatchId;
  state.nextSourceBatchId += 1;
  const bundleId = `sextic-${batchId}`;

  state.sexticBundles[bundleId] = {
    id: bundleId,
    fileName,
    inputRepresentation: parsedSextic.inputRepresentation,
    inertial: { ...parsedSextic.inertial },
    rotationalConstantsMHz: parsedSextic.rotationalConstantsMHz
      ? { ...parsedSextic.rotationalConstantsMHz }
      : null,
    moments: parsedSextic.moments ? { ...parsedSextic.moments } : null,
    sigma: Number.isFinite(parsedSextic.sigma) ? parsedSextic.sigma : null,
    sigmaSource: parsedSextic.sigmaSource || "none",
  };

  const sourceRep = normalizeWatsonRepresentationToken(parsedSextic.inputRepresentation) || "Ir";
  const derived = buildSexticDerivedPayloads(state.sexticBundles[bundleId], sourceRep);
  const sexticSource = {
    id: `${batchId}-sextic-1`,
    type: "sextic",
    fileName,
    title: `${fileName} - Sextic [input ${sourceRep}]`,
    details: summarizeSexticComponent(derived.sextic),
    payload: derived.sextic,
    defaultForType: true,
    contextLabel: derived.sextic.contextLabel || "",
    sourceIndex: derived.sextic.sourceIndex || Date.now(),
    sexticBundleId: bundleId,
    sexticRole: "sextic",
  };

  state.logFileName = fileName;
  state.logFiles.push(fileName);
  state.logEntries.push({
    id: `log-${batchId}`,
    fileName,
    hasCompletePipeline: false,
    componentsLabel: `Sextic x1 (input ${sourceRep})`,
  });
  state.logSources.push(sexticSource);
  state.sourceById[sexticSource.id] = sexticSource;
  state.flowSummaryText = buildFlowSummaryFromEntries();

  if (!state.connections.sextic) {
    state.connections.sextic = sexticSource.id;
  }

  renderLoadedLogList();
  renderFlowBoard();
  void refreshWMSInputOutput();
  return {
    representation: sourceRep,
    inputRepresentation: parsedSextic.inputRepresentation,
  };
}

function refreshFchkSourcesForRepresentation() {
  const rep = getStep2Representation();
  const bundleIds = Object.keys(state.fchkBundles || {});
  if (!bundleIds.length) return;

  for (let i = 0; i < bundleIds.length; i += 1) {
    const bundleId = bundleIds[i];
    const bundle = state.fchkBundles[bundleId];
    if (!bundle) continue;
    let derived;
    try {
      derived = buildFchkDerivedPayloads(bundle, rep);
    } catch (error) {
      continue;
    }
    const sourceIds = Object.keys(state.sourceById);
    const hasBDPCS3GeometrySource = sourceIds.some((sourceId) => {
      const source = state.sourceById[sourceId];
      return source && source.fchkBundleId === bundleId && source.fchkRole === "bdpcs3-geometry";
    });
    const bdpcs3Geometry = hasBDPCS3GeometrySource
      ? buildBDPCS3DerivedGeometryPayload(derived.geometry)
      : null;
    for (let j = 0; j < sourceIds.length; j += 1) {
      const source = state.sourceById[sourceIds[j]];
      if (!source || source.fchkBundleId !== bundleId) continue;
      if (source.fchkRole === "geometry") {
        source.payload = derived.geometry;
        source.details = summarizeGeometryComponent(derived.geometry);
        source.title = `${bundle.fileName} - Geometry [FCHK rotated ${rep}]`;
        source.contextLabel = derived.geometry.contextLabel || "";
      } else if (source.fchkRole === "quartic") {
        source.payload = derived.quartic;
        source.details = summarizeQuarticComponent(derived.quartic);
        source.title = `${bundle.fileName} - Quartic [FCHK rotated ${rep}]`;
        source.contextLabel = derived.quartic.contextLabel || "";
      } else if (source.fchkRole === "bdpcs3-geometry" && bdpcs3Geometry) {
        source.payload = bdpcs3Geometry;
        source.details = summarizeGeometryComponent(bdpcs3Geometry);
        source.title = `${bundle.fileName} - Geometry [BDPCS3 corrected] - FCHK rotated ${rep}`;
        source.contextLabel = bdpcs3Geometry.contextLabel || "";
      }
    }
  }
}

function refreshTauSourcesForRepresentation() {
  const rep = getStep2Representation();
  const bundleIds = Object.keys(state.tauBundles || {});
  if (!bundleIds.length) return;
  for (let i = 0; i < bundleIds.length; i += 1) {
    const bundleId = bundleIds[i];
    const bundle = state.tauBundles[bundleId];
    if (!bundle) continue;
    let derived;
    try {
      derived = buildTauDerivedPayloads(bundle, rep);
    } catch (error) {
      continue;
    }
    const sourceIds = Object.keys(state.sourceById);
    for (let j = 0; j < sourceIds.length; j += 1) {
      const source = state.sourceById[sourceIds[j]];
      if (!source || source.tauBundleId !== bundleId) continue;
      if (source.tauRole === "quartic") {
        source.payload = derived.quartic;
        source.details = summarizeQuarticComponent(derived.quartic);
        source.title = `${bundle.fileName} - Quartic [Tau rotated ${rep}]`;
        source.contextLabel = derived.quartic.contextLabel || "";
      }
    }
  }
}

function refreshSexticSourcesForRepresentation() {
  const rep = getStep2Representation();
  const bundleIds = Object.keys(state.sexticBundles || {});
  if (!bundleIds.length) return;
  for (let i = 0; i < bundleIds.length; i += 1) {
    const bundleId = bundleIds[i];
    const bundle = state.sexticBundles[bundleId];
    if (!bundle) continue;
    let derived;
    try {
      derived = buildSexticDerivedPayloads(bundle, rep);
    } catch (error) {
      continue;
    }
    const sourceIds = Object.keys(state.sourceById);
    for (let j = 0; j < sourceIds.length; j += 1) {
      const source = state.sourceById[sourceIds[j]];
      if (!source || source.sexticBundleId !== bundleId) continue;
      if (source.sexticRole === "sextic") {
        source.payload = derived.sextic;
        source.details = summarizeSexticComponent(derived.sextic);
        source.title = `${bundle.fileName} - Sextic [rotated ${rep}]`;
        source.contextLabel = derived.sextic.contextLabel || "";
      }
    }
  }
}

function rotateTauPrimeComponentsForRepresentation(
  tauPrime,
  inputRepresentation,
  outputRepresentation
) {
  const remap = buildTauAxisRemap("Ir", outputRepresentation);
  const outKeys = ["taaaa", "tbbaa", "tbbbb", "tccaa", "tccbb", "tcccc"];
  const rotated = {};
  for (let i = 0; i < outKeys.length; i += 1) {
    const key = outKeys[i];
    const outputIndex = key.slice(1);
    const mappedInputIndex =
      remap.axisOutToIn[outputIndex[0]] +
      remap.axisOutToIn[outputIndex[1]] +
      remap.axisOutToIn[outputIndex[2]] +
      remap.axisOutToIn[outputIndex[3]];
    const value = getTauPrimeValueByIndices(tauPrime, mappedInputIndex);
    if (!Number.isFinite(value)) {
      throw new Error(
        `Unable to rotate TauP component "${key}" to representation ${outputRepresentation}.`
      );
    }
    rotated[key] = value;
  }
  return {
    tauPrime: rotated,
    remap,
  };
}

function getTauPrimeValueByIndices(tauPrime, indexToken) {
  const key = canonicalTauPrimeKeyFromToken(indexToken);
  if (!key) return null;
  const value = tauPrime ? tauPrime[key] : null;
  return Number.isFinite(value) ? value : null;
}

function canonicalTauPrimeKeyFromToken(value) {
  const token = String(value || "")
    .toLowerCase()
    .replace(/[^abc]/g, "");
  if (token.length !== 4) return null;
  if (token[0] === token[1] && token[1] === token[2] && token[2] === token[3]) {
    return `t${token}`;
  }
  if (token[0] !== token[1] || token[2] !== token[3]) return null;
  const left = token[0];
  const right = token[2];
  const sorted = [left, right].sort();
  const low = sorted[0];
  const high = sorted[1];
  return `t${high}${high}${low}${low}`;
}

function loadXYZFromText(text, sourceName) {
  const parsed = parseXYZ(text);
  state.fileName = sourceName || "Input XYZ";
  state.sourceComment = parsed.comment || "";
  state.loadedAtoms = parsed.atoms;
  computePrincipalFrame();
  applyRepresentationAndRefresh();
  return parsed;
}

function onGaussianLogSelected(event) {
  const files = Array.from((event.target && event.target.files) || []);
  if (!files.length) {
    return;
  }

  Promise.all(
    files.map((file) =>
      readFileAsText(file).then((text) => ({
        name: file.name,
        parsed: parseGaussianLogComponents(text),
      }))
    )
  )
    .then((parsedLogs) => {
      const createdSources = [];
      const createdEntries = [];
      for (let i = 0; i < parsedLogs.length; i += 1) {
        const entry = parsedLogs[i];
        const batchId = state.nextSourceBatchId;
        state.nextSourceBatchId += 1;
        const created = createSourcesFromParsedLog(entry.parsed, entry.name, batchId);
        createdSources.push(...created);
        const typeCounts = {
          geometry: created.filter((source) => source.type === "geometry").length,
          dvib: created.filter((source) => source.type === "dvib").length,
          quartic: created.filter((source) => source.type === "quartic").length,
          sextic: created.filter((source) => source.type === "sextic").length,
          quadrupole: created.filter((source) => source.type === "quadrupole").length,
        };
        const presentTypes = FLOW_SOURCE_LANE_ORDER
          .filter((type) => typeCounts[type] > 0)
          .map((type) => `${FLOW_TARGET_LABELS[type]} x${typeCounts[type]}`);
        createdEntries.push({
          id: `log-${batchId}`,
          fileName: entry.name,
          hasCompletePipeline: Boolean(entry.parsed.summary && entry.parsed.summary.hasCompletePipeline),
          componentsLabel: presentTypes.length ? presentTypes.join(", ") : "no usable WMS components",
        });
      }

      state.logFileName = files[files.length - 1].name;
      state.logFiles.push(...files.map((f) => f.name));
      state.logEntries.push(...createdEntries);
      state.logSources.push(...createdSources);
      for (let i = 0; i < createdSources.length; i += 1) {
        state.sourceById[createdSources[i].id] = createdSources[i];
      }
      refreshGaussianLogDistortionSourcePresentation();
      state.flowSummaryText = buildFlowSummaryFromEntries();

      const suggestInput = state.logEntries.map((item) => ({
        name: item.fileName,
        parsed: { summary: { hasCompletePipeline: item.hasCompletePipeline } },
      }));
      const suggestedConnections = autoConnectFlowSources(suggestInput, state.logSources);
      for (let i = 0; i < FLOW_TARGETS.length; i += 1) {
        const targetType = FLOW_TARGETS[i];
        const current = state.connections[targetType];
        if (current && state.sourceById[current]) continue;
        if (suggestedConnections[targetType]) {
          state.connections[targetType] = suggestedConnections[targetType];
        }
      }

      renderLoadedLogList();
      renderFlowBoard();
      void refreshWMSInputOutput();

      if (!createdSources.length) {
        setStatus("No compatible WMS-Rot data blocks found in the selected Gaussian log files.", true);
        return;
      }

      const completePipeline = parsedLogs.find((entry) => entry.parsed.summary.hasCompletePipeline);
      if (completePipeline) {
        setStatus(
          `Added ${files.length} log file(s). Auto-connected missing targets from ${completePipeline.name}.`
        );
      } else {
        setStatus(
          `Added ${files.length} log file(s). Drag hooks to connect or reconnect blocks.`
        );
      }
    })
    .catch((error) => {
      clearGaussianLogState();
      setStatus(error instanceof Error ? error.message : String(error), true);
    })
    .finally(() => {
      if (ui.gaussianLogFile) {
        ui.gaussianLogFile.value = "";
      }
    });
}

function readFileAsText(file) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(String(reader.result || ""));
    reader.onerror = () => reject(new Error(`Error while reading "${file.name}".`));
    reader.readAsText(file);
  });
}

function clearState() {
  state.fileName = null;
  state.sourceComment = "";
  state.loadedAtoms = null;
  state.principalAtoms = null;
  state.analysisData = null;
  state.outputAtoms = null;
  ui.xyzOutput.value = "";
  ui.analysis.textContent = "No file loaded.";
  updateStep1DetectedRepresentationBadge(null);
  syncStep1OutputRepresentationControl();
  ui.downloadBtn.disabled = true;
  ui.downloadGjfBtn.disabled = true;
  ui.downloadBasisSetBtn.disabled = true;
}

function updateStep1DetectedRepresentationBadge(detection) {
  if (!ui.step1DetectedRepresentation) return;
  if (!detection || !detection.representation) {
    const reason = detection && detection.reason ? ` (${detection.reason})` : "";
    ui.step1DetectedRepresentation.textContent = `Undetermined${reason}`;
    return;
  }
  const reasonSuffix = detection.reason ? ` (${detection.reason})` : "";
  const alignmentSuffix = Number.isFinite(detection.alignment)
    ? ` (${(detection.alignment * 100).toFixed(1)}% axis alignment)`
    : "";
  ui.step1DetectedRepresentation.textContent =
    `${detection.representation}${reasonSuffix}${alignmentSuffix}`;
}

function clearGaussianLogState() {
  clearActiveLinkDrag();
  state.logFileName = null;
  state.logFiles = [];
  state.logEntries = [];
  state.nextSourceBatchId = 0;
  state.logSources = [];
  state.sourceById = {};
  state.connections = {
    geometry: null,
    dvib: null,
    quartic: null,
    sextic: null,
    quadrupole: null,
  };
  state.nodePositions = {};
  state.laneBounds = {};
  state.activeLinkDrag = null;
  state.combinedParsedData = null;
  state.flowSummaryText = "No source loaded.";
  state.wmsInputText = "";
  state.catrefInputText = "";
  state.catrefBundles = [];
  state.catrefVarText = "";
  state.catrefIntText = "";
  state.dpcs3XYZText = "";
  state.xyzBundles = {};
  state.fchkBundles = {};
  state.hessianBundles = {};
  state.tauBundles = {};
  state.sexticBundles = {};
  state.outputRefreshToken += 1;
  ui.wmsOutput.value = "";
  if (ui.catrefOutput) ui.catrefOutput.value = "";
  if (ui.dpcs3XYZOutput) ui.dpcs3XYZOutput.value = "";
  renderWMSOutputRepresentationWarning([]);
  syncWMSOutputRepresentationControl();
  ui.downloadWmsInputBtn.disabled = true;
  if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
  if (ui.openBdpcs3Btn) ui.openBdpcs3Btn.disabled = true;
  renderLoadedLogList();
  updateFlowSummary();
  renderFlowBoard();
}

function createSourcesFromParsedLog(parsedLog, fileName, fileIndex) {
  const sources = [];
  appendSourcesForType({
    sources,
    fileName,
    fileIndex,
    type: "geometry",
    variants: parsedLog.geometryVariants || [],
    defaultVariant: parsedLog.geometry || null,
    summarize: summarizeGeometryComponent,
  });
  appendSourcesForType({
    sources,
    fileName,
    fileIndex,
    type: "quartic",
    variants: parsedLog.quarticVariants || [],
    defaultVariant: parsedLog.quartic || null,
    summarize: summarizeQuarticComponent,
  });
  appendSourcesForType({
    sources,
    fileName,
    fileIndex,
    type: "sextic",
    variants: parsedLog.sexticVariants || [],
    defaultVariant: parsedLog.sextic || null,
    summarize: summarizeSexticComponent,
  });
  appendSourcesForType({
    sources,
    fileName,
    fileIndex,
    type: "dvib",
    variants: parsedLog.dvibVariants || [],
    defaultVariant: parsedLog.dvib || null,
    summarize: summarizeDVibComponent,
  });
  appendSourcesForType({
    sources,
    fileName,
    fileIndex,
    type: "quadrupole",
    variants: parsedLog.quadrupoleVariants || [],
    defaultVariant: parsedLog.quadrupole || null,
    summarize: summarizeQuadrupoleComponent,
  });
  return sources;
}

function appendSourcesForType({
  sources,
  fileName,
  fileIndex,
  type,
  variants,
  defaultVariant,
  summarize,
}) {
  if (!Array.isArray(variants) || !variants.length) return;
  for (let i = 0; i < variants.length; i += 1) {
    const variant = variants[i];
    const isDefault = variant === defaultVariant;
    const context = variant.contextLabel ? ` - ${variant.contextLabel}` : "";
    const source = {
      id: `${fileIndex}-${type}-${i + 1}`,
      type,
      fileName,
      title: `${fileName} - ${FLOW_TARGET_LABELS[type]}${context}${isDefault ? " [default]" : ""}`,
      details: summarize(variant),
      payload: variant,
      defaultForType: isDefault,
      contextLabel: variant.contextLabel || "",
      sourceIndex: Number.isFinite(variant.sourceIndex) ? variant.sourceIndex : i,
      originType: "A",
    };
    applyGaussianLogDistortionSourcePresentation(source);
    sources.push(source);
  }
}

function createBDPCS3DerivedSourceForGeometry(source) {
  if (!source || source.type !== "geometry" || source.bdpcs3Derived) return null;
  const payload = buildBDPCS3DerivedGeometryPayload(source.payload);
  if (!payload) return null;
  const contextText = source.contextLabel ? ` - ${source.contextLabel}` : "";
  const derivedId = buildUniqueSourceId(`${source.id}-bdpcs3`);
  const pairGroupId = source.bdpcs3PairGroupId || `bdpcs3-pair-${source.id}`;
  const derived = {
    id: derivedId,
    type: "geometry",
    fileName: source.fileName,
    title: `${source.fileName} - Geometry [BDPCS3 corrected]${contextText}`,
    details: summarizeGeometryComponent(payload),
    payload,
    defaultForType: false,
    contextLabel: payload.contextLabel || "",
    sourceIndex: Number.isFinite(source.sourceIndex) ? source.sourceIndex + 0.0001 : Date.now(),
    bdpcs3Derived: true,
    bdpcs3DerivedFromSourceId: source.id,
    bdpcs3PairSourceId: source.id,
    bdpcs3PairGroupId: pairGroupId,
    fchkBundleId: source.fchkBundleId || null,
    fchkRole: source.fchkBundleId ? "bdpcs3-geometry" : null,
    originType: source.originType || "A",
  };
  source.bdpcs3PairSourceId = derived.id;
  source.bdpcs3PairGroupId = pairGroupId;
  return derived;
}

function createBDPCS3PairFromSourceId(sourceId) {
  const source = state.sourceById[sourceId];
  if (!source || source.type !== "geometry") {
    setStatus("Select a valid Geometry source to create a BDPCS3 Δ pair.", true);
    return;
  }
  if (source.bdpcs3Derived) {
    setStatus("This geometry is already a BDPCS3 corrected source.");
    return;
  }
  if (source.bdpcs3PairSourceId && state.sourceById[source.bdpcs3PairSourceId]) {
    setStatus("BDPCS3 Δ source already created for this geometry.");
    return;
  }
  const derived = createBDPCS3DerivedSourceForGeometry(source);
  if (!derived) {
    setStatus(
      "Unable to compute BDPCS3 corrected geometry from this source. Check geometry content and retry.",
      true
    );
    return;
  }
  insertSourceAfter(state.logSources, source.id, derived);
  state.sourceById[derived.id] = derived;
  renderFlowBoard();
  void refreshWMSInputOutput();
  setStatus(`Created BDPCS3 Δ geometry for "${source.fileName}".`);
}

function buildUniqueSourceId(baseId) {
  let out = String(baseId || "source");
  let suffix = 2;
  while (state.sourceById[out]) {
    out = `${baseId}-${suffix}`;
    suffix += 1;
  }
  return out;
}

function insertSourceAfter(list, anchorId, item) {
  if (!Array.isArray(list)) return;
  const idx = list.findIndex((entry) => entry && entry.id === anchorId);
  if (idx < 0) {
    list.push(item);
    return;
  }
  list.splice(idx + 1, 0, item);
}

function buildBDPCS3DerivedGeometryPayload(baseGeometry) {
  if (!baseGeometry || !Array.isArray(baseGeometry.dpcs3Atoms) || !baseGeometry.dpcs3Atoms.length) {
    return null;
  }
  const baseAtoms = cloneAtoms(baseGeometry.dpcs3Atoms);
  if (!baseAtoms.length) return null;
  const corrected = runBDPCS3CorrectionOnAtoms(baseAtoms);
  if (!corrected || !Array.isArray(corrected.atoms) || !corrected.atoms.length) return null;

  const refinedAtoms = corrected.atoms;
  const refinedRotational =
    corrected.rotationalConstantsMHz || estimateRotationalConstantsMHz(refinedAtoms);
  const dipole = normalizeGeometryDipole(baseGeometry.dipole);
  const contextPrefix = String(baseGeometry.contextLabel || "").trim();
  const contextLabel = contextPrefix
    ? `${contextPrefix} | BDPCS3 corrected geometry (Δ)`
    : "BDPCS3 corrected geometry (Δ)";

  return {
    contextLabel,
    sourceIndex: Number.isFinite(baseGeometry.sourceIndex)
      ? baseGeometry.sourceIndex + 0.0001
      : Date.now(),
    pointGroup: baseGeometry.pointGroup || "CS",
    rotationalConstantsMHz: refinedRotational,
    dipole,
    dpcs3Atoms: refinedAtoms,
    dpcs3XYZText: buildAtomicNumberXYZText(refinedAtoms),
    bdpcs3Info: corrected.info || null,
  };
}

function runBDPCS3CorrectionOnAtoms(atoms) {
  if (!canRunBDPCS3CorrectionLocally()) return null;
  // BDPCS3's Molecule.readFromString() expects standard XYZ with first two header lines.
  const xyzText = buildBDPCS3XYZFileText(atoms);
  if (!xyzText.trim()) return null;
  try {
    const molecule = new Molecule();
    molecule.readFromString(xyzText);
    if (molecule.nAtmsOriginal !== atoms.length) {
      throw new Error(
        `BDPCS3 input atom mismatch: expected ${atoms.length}, parsed ${molecule.nAtmsOriginal}.`
      );
    }
    if (!Number.isFinite(molecule.nAtmsOriginal) || molecule.nAtmsOriginal < 1) {
      return null;
    }
    molecule.eckartOrientation(false);
    const beforeRotationalMHz = extractBDPCS3RotationalConstantsMHz(
      molecule.computeRotationalConstant()
    );
    molecule.addLinearAnglesLonePairs();
    const optimizer = new BDPCS3Optimizer(molecule, true);
    const optimizerStatus = String(optimizer.gradient_descent(0.4, 10000, 0) || "");
    const afterRotationalMHz = extractBDPCS3RotationalConstantsMHz(
      molecule.computeRotationalConstant()
    );
    const refinedAtoms = [];
    for (let i = 0; i < molecule.nAtmsOriginal; i += 1) {
      refinedAtoms.push({
        label: String(molecule.atomicNumbers[i]),
        z: Number(molecule.atomicNumbers[i]),
        coords: [
          Number(molecule.coordinates[i][0]),
          Number(molecule.coordinates[i][1]),
          Number(molecule.coordinates[i][2]),
        ],
      });
    }
    return {
      atoms: cloneAtoms(refinedAtoms),
      rotationalConstantsMHz: afterRotationalMHz,
      info: {
        status: optimizerStatus,
        inputAtomCount: atoms.length,
        outputAtomCount: molecule.nAtmsOriginal,
        beforeRotationalMHz,
        afterRotationalMHz,
      },
    };
  } catch (error) {
    return null;
  }
}

function extractBDPCS3RotationalConstantsMHz(raw) {
  if (!Array.isArray(raw)) return null;
  if (raw.length >= 4) {
    const A = Number(raw[1]);
    const B = Number(raw[2]);
    const C = Number(raw[3]);
    if (Number.isFinite(A) && Number.isFinite(B) && Number.isFinite(C)) {
      return { A, B, C };
    }
  }
  if (raw.length >= 2) {
    const value = Number(raw[1]);
    if (Number.isFinite(value)) {
      return { A: value, B: value, C: Number.NaN };
    }
  }
  return null;
}

function canRunBDPCS3CorrectionLocally() {
  return (
    typeof Molecule !== "undefined" &&
    typeof BDPCS3Optimizer !== "undefined"
  );
}

function normalizeGeometryDipole(dipole) {
  if (!dipole || typeof dipole !== "object") {
    return { a: "", b: "", c: "" };
  }
  return {
    a: String(dipole.a || "").trim(),
    b: String(dipole.b || "").trim(),
    c: String(dipole.c || "").trim(),
  };
}

function autoConnectFlowSources(parsedLogs, sources) {
  const connections = {
    geometry: null,
    dvib: null,
    quartic: null,
    sextic: null,
    quadrupole: null,
  };

  let preferredFileName = null;
  for (let i = parsedLogs.length - 1; i >= 0; i -= 1) {
    const item = parsedLogs[i];
    if (item.parsed.summary.hasCompletePipeline) {
      preferredFileName = item.name;
      break;
    }
  }

  for (let i = 0; i < FLOW_TARGETS.length; i += 1) {
    const type = FLOW_TARGETS[i];
    const preferredPool = preferredFileName
      ? sources.filter((source) => source.type === type && source.fileName === preferredFileName)
      : sources.filter((source) => source.type === type);
    const fallbackPool = sources.filter((source) => source.type === type);

    let chosen = pickBestSourceForType(preferredPool);
    if (!chosen) {
      chosen = pickBestSourceForType(fallbackPool);
    }
    if (chosen) {
      connections[type] = chosen.id;
    }
  }
  return connections;
}

function pickBestSourceForType(pool) {
  if (!Array.isArray(pool) || !pool.length) return null;
  const defaults = pool.filter((source) => source.defaultForType);
  const list = defaults.length ? defaults : pool;
  return list[list.length - 1] || null;
}

function summarizeGeometryComponent(component) {
  const lines = [];
  if (component.contextLabel) lines.push(component.contextLabel);
  if (component.pointGroup) lines.push(`Point Group: ${component.pointGroup}`);
  if (component.rotationalConstantsMHz) {
    lines.push(
      `A/B/C MHz: ${formatCompactNumber(component.rotationalConstantsMHz.A)} / ${formatCompactNumber(component.rotationalConstantsMHz.B)} / ${formatCompactNumber(component.rotationalConstantsMHz.C)}`
    );
  }
  if (component.dipole) {
    lines.push(`Dipole a/b/c D: ${component.dipole.a} / ${component.dipole.b} / ${component.dipole.c}`);
  }
  if (Array.isArray(component.dpcs3Atoms)) {
    lines.push(`Atoms: ${component.dpcs3Atoms.length}`);
  }
  if (component.bdpcs3Info) {
    if (component.bdpcs3Info.status) {
      lines.push(`BDPCS3 status: ${component.bdpcs3Info.status}`);
    }
    const before = component.bdpcs3Info.beforeRotationalMHz;
    if (before) {
      lines.push(
        `BDPCS3 before A/B/C MHz: ${formatCompactNumber(before.A)} / ${formatCompactNumber(before.B)} / ${formatCompactNumber(before.C)}`
      );
    }
    const after = component.bdpcs3Info.afterRotationalMHz;
    if (after) {
      lines.push(
        `BDPCS3 after A/B/C MHz: ${formatCompactNumber(after.A)} / ${formatCompactNumber(after.B)} / ${formatCompactNumber(after.C)}`
      );
    }
  }
  return lines.join("\n");
}

function summarizeQuarticComponent(component, options = {}) {
  const lines = [];
  const prefixLines = Array.isArray(options.prefixLines) ? options.prefixLines : [];
  const includeContext = options.includeContext !== false;
  if (includeContext && component.contextLabel) lines.push(component.contextLabel);
  for (let i = 0; i < prefixLines.length; i += 1) {
    const line = String(prefixLines[i] || "").trim();
    if (line) lines.push(line);
  }
  const asym = component.asym || {};
  const sym = component.sym || {};
  lines.push(
    `A: DELTA J=${formatPreviewValue(asym.deltaJ)} | DELTA JK=${formatPreviewValue(asym.deltaJK)} | DELTA K=${formatPreviewValue(asym.deltaK)} | delta J=${formatPreviewValue(asym.smallDeltaJ)} | delta K=${formatPreviewValue(asym.smallDeltaK)}`
  );
  lines.push(
    `S: DJ=${formatPreviewValue(sym.DJ)} | DJK=${formatPreviewValue(sym.DJK)} | DK=${formatPreviewValue(sym.DK)} | d1=${formatPreviewValue(sym.d1)} | d2=${formatPreviewValue(sym.d2)}`
  );
  return lines.join("\n");
}

function summarizeSexticComponent(component, options = {}) {
  const lines = [];
  const prefixLines = Array.isArray(options.prefixLines) ? options.prefixLines : [];
  const includeContext = options.includeContext !== false;
  if (includeContext && component.contextLabel) lines.push(component.contextLabel);
  for (let i = 0; i < prefixLines.length; i += 1) {
    const line = String(prefixLines[i] || "").trim();
    if (line) lines.push(line);
  }
  const asym = component.asymSextic || {};
  const sym = component.symSextic || {};
  lines.push(
    `A: PHI N=${formatPreviewValue(asym.PHIN)} | PHI NK=${formatPreviewValue(asym.PHINK)} | PHI KN=${formatPreviewValue(asym.PHIKN)} | PHI K=${formatPreviewValue(asym.PHIK)} | phi N=${formatPreviewValue(asym.phiN)} | phi NK=${formatPreviewValue(asym.phiNK)} | phi K=${formatPreviewValue(asym.phiK)}`
  );
  lines.push(
    `S: H N=${formatPreviewValue(sym.HN)} | H NK=${formatPreviewValue(sym.HNK)} | H KN=${formatPreviewValue(sym.HKN)} | H K=${formatPreviewValue(sym.HK)} | h1=${formatPreviewValue(sym.h1)} | h2=${formatPreviewValue(sym.h2)} | h3=${formatPreviewValue(sym.h3)}`
  );
  return lines.join("\n");
}

function summarizeHessianComponent(component) {
  if (!component || typeof component !== "object") return "";
  const lines = [];
  if (component.contextLabel) lines.push(component.contextLabel);
  if (Number.isFinite(component.atomCount)) {
    lines.push(`Atoms: ${component.atomCount}`);
  }
  if (Number.isFinite(component.matrixSize)) {
    lines.push(`Matrix: ${component.matrixSize} x ${component.matrixSize}`);
  }
  lines.push(`Symmetric: ${component.symmetric ? "yes" : "no"}`);
  if (Number.isFinite(component.trace)) {
    lines.push(`Trace: ${formatCompactNumber(component.trace)}`);
  }
  if (Number.isFinite(component.frobeniusNorm)) {
    lines.push(`Frobenius norm: ${formatCompactNumber(component.frobeniusNorm)}`);
  }
  return lines.join("\n");
}

function formatPreviewValue(value) {
  if (value === null || value === undefined) return "—";
  const text = String(value).trim();
  return text || "—";
}

function stripLegacyLogDistortionContextLabel(text) {
  return String(text || "")
    .replace(/\s*\|\s*TauP-derived quartic\s*$/i, "")
    .replace(/\s*\|\s*Phi\/eta-derived sextic\s*$/i, "")
    .trim();
}

function buildGaussianLogDistortionRepresentationLine(component) {
  if (!component || typeof component !== "object") return "";
  const sourceRep = normalizeWatsonRepresentationToken(component.inputRepresentation);
  if (!sourceRep) return "";
  return `Representation from log: ${sourceRep}`;
}

function buildGaussianLogDistortionSourceTitle(source) {
  const payload = source && source.payload ? source.payload : {};
  const context = stripLegacyLogDistortionContextLabel(payload.contextLabel || source.contextLabel || "");
  const sourceRep = normalizeWatsonRepresentationToken(payload.inputRepresentation);
  const parts = [`${source.fileName} - ${FLOW_TARGET_LABELS[source.type]}`];
  if (context) parts.push(`- ${context}`);
  if (sourceRep) {
    parts.push(`[${sourceRep}]`);
  }
  if (source.defaultForType) parts.push("[default]");
  return parts.join(" ");
}

function applyGaussianLogDistortionSourcePresentation(source) {
  if (!source || String(source.originType || "").toUpperCase() !== "A") return;
  if (source.type !== "quartic" && source.type !== "sextic") return;
  const displayPayload = {
    ...(source.payload || {}),
    asym: { ...((source.payload && source.payload.asym) || {}) },
    sym: { ...((source.payload && source.payload.sym) || {}) },
    asymSextic: { ...((source.payload && source.payload.asymSextic) || {}) },
    symSextic: { ...((source.payload && source.payload.symSextic) || {}) },
  };
  const representationLine = buildGaussianLogDistortionRepresentationLine(displayPayload);
  source.title = buildGaussianLogDistortionSourceTitle(source);
  source.details =
    source.type === "quartic"
      ? summarizeQuarticComponent(displayPayload, { prefixLines: [representationLine] })
      : summarizeSexticComponent(displayPayload, { prefixLines: [representationLine] });
  source.contextLabel = displayPayload.contextLabel || "";
}

function refreshGaussianLogDistortionSourcePresentation() {
  const sourceIds = Object.keys(state.sourceById || {});
  for (let i = 0; i < sourceIds.length; i += 1) {
    const source = state.sourceById[sourceIds[i]];
    applyGaussianLogDistortionSourcePresentation(source);
  }
}

function summarizeDVibComponent(component) {
  const prefix = component.contextLabel ? `${component.contextLabel}\n` : "";
  return `${prefix}DVib MHz: ${formatCompactNumber(component.A)} / ${formatCompactNumber(component.B)} / ${formatCompactNumber(component.C)}`;
}

function summarizeQuadrupoleComponent(component) {
  const lines = [];
  if (component.contextLabel) lines.push(component.contextLabel);
  const inputRep = normalizeWatsonRepresentationToken(component.inputRepresentation);
  if (inputRep) {
    lines.push(`Representation from log: ${inputRep}`);
  }
  const nuclei = Array.isArray(component.nuclei) ? component.nuclei : [];
  if (!nuclei.length) {
    lines.push("No quadrupole tensors detected.");
    return lines.join("\n");
  }
  for (let i = 0; i < nuclei.length; i += 1) {
    const nucleus = nuclei[i];
    const chi = nucleus && nucleus.chi ? nucleus.chi : {};
    const label = nucleus && nucleus.label ? nucleus.label : `Nucleus ${i + 1}`;
    lines.push(
      `${label}: aa=${formatCompactNumber(chi.chi_aa)} | bb=${formatCompactNumber(chi.chi_bb)} | cc=${formatCompactNumber(chi.chi_cc)} | ab=${formatCompactNumber(chi.chi_ab)} | ac=${formatCompactNumber(chi.chi_ac)} | bc=${formatCompactNumber(chi.chi_bc)}`
    );
  }
  return lines.join("\n");
}

function hasAnyNonEmptyValue(obj) {
  if (!obj || typeof obj !== "object") return false;
  const values = Object.values(obj);
  for (let i = 0; i < values.length; i += 1) {
    const value = values[i];
    if (value === null || value === undefined) continue;
    if (typeof value === "string" && !value.trim()) continue;
    if (Number.isFinite(value)) return true;
    if (typeof value === "string") return true;
  }
  return false;
}

function updateFlowSummary() {
  if (!ui.logSummary) return;
  ui.logSummary.textContent = state.flowSummaryText || "No source loaded.";
}

function renderLoadedLogList() {
  if (!ui.loadedLogList) return;
  ui.loadedLogList.innerHTML = "";
  if (!state.logEntries.length) {
    const li = document.createElement("li");
    li.textContent = "No source loaded.";
    ui.loadedLogList.appendChild(li);
    return;
  }
  for (let i = 0; i < state.logEntries.length; i += 1) {
    const entry = state.logEntries[i];
    const li = document.createElement("li");
    li.textContent = `${entry.fileName} - ${entry.componentsLabel}`;
    ui.loadedLogList.appendChild(li);
  }
}

function buildFlowSummaryFromEntries() {
  if (!state.logEntries.length) return "No source loaded.";
  const lines = [];
  for (let i = 0; i < state.logEntries.length; i += 1) {
    const entry = state.logEntries[i];
    lines.push(`${entry.fileName}: ${entry.componentsLabel}`);
  }
  return lines.join("\n");
}

function renderFlowBoard() {
  if (!ui.flowBoard || !ui.flowNodes || !ui.flowLines || !ui.flowLanes) return;
  clearActiveLinkDrag();
  ui.flowNodes.innerHTML = "";
  ui.flowLines.innerHTML = "";
  ui.flowLanes.innerHTML = "";
  updateFlowSummary();
  const validNodeIds = new Set();
  state.laneBounds = {};
  ui.flowBoard.style.width = "1180px";
  ui.flowBoard.style.height = "420px";

  if (!state.logSources.length) {
    const empty = document.createElement("div");
    empty.className = "flow-empty";
    empty.textContent = "Load Gaussian log or XYZ data. Detected components appear as source boxes that can be linked to the WMS-Rot targets.";
    ui.flowNodes.appendChild(empty);
    pruneNodePositions(validNodeIds);
    ensureFlowBoardSize();
    scheduleFlowLineRedraw();
    return;
  }

  const laneLeftPadding = 16;
  const laneGap = 14;
  const sourceVerticalGap = 16;
  const sourceStartY = 56;
  const laneCount = FLOW_SOURCE_LANE_ORDER.length;
  const groupedSources = {
    geometry: [],
    hessian: [],
    dvib: [],
    quartic: [],
    sextic: [],
    quadrupole: [],
  };
  for (let i = 0; i < state.logSources.length; i += 1) {
    const source = state.logSources[i];
    if (!groupedSources[source.type]) continue;
    groupedSources[source.type].push(source);
  }

  const hasGeometryPairs = groupedSources.geometry.some(
    (source) => source.bdpcs3PairSourceId || source.bdpcs3DerivedFromSourceId
  );
  const baseLaneWidths = FLOW_SOURCE_LANE_ORDER.map((type) => {
    if (type === "geometry" && hasGeometryPairs) {
      return BDPCS3_PAIR_MIN_LANE_WIDTH;
    }
    return 250;
  });
  const baseLaneWidthSum = baseLaneWidths.reduce((sum, width) => sum + width, 0);
  const minBoardWidth =
    laneLeftPadding * 2 + baseLaneWidthSum + (laneCount - 1) * laneGap;
  const wrapWidth = ui.flowBoardWrap ? Math.floor(ui.flowBoardWrap.clientWidth) : 0;
  const boardWidth = Math.max(minBoardWidth, wrapWidth || 0);
  const extraWidth = Math.max(0, boardWidth - minBoardWidth);
  const laneExtra = Math.floor(extraWidth / laneCount);
  let laneExtraRemainder = extraWidth - laneExtra * laneCount;

  const laneMetrics = [];
  let laneCursorX = laneLeftPadding;
  for (let laneIndex = 0; laneIndex < FLOW_SOURCE_LANE_ORDER.length; laneIndex += 1) {
    const type = FLOW_SOURCE_LANE_ORDER[laneIndex];
    const width =
      baseLaneWidths[laneIndex] +
      laneExtra +
      (laneExtraRemainder > 0 ? 1 : 0);
    if (laneExtraRemainder > 0) {
      laneExtraRemainder -= 1;
    }
    laneMetrics.push({
      type,
      x: laneCursorX,
      width,
    });
    laneCursorX += width + laneGap;
  }

  const sourceNodeWidth = 188;
  let maxSourceBottom = sourceStartY;
  for (let laneIndex = 0; laneIndex < laneMetrics.length; laneIndex += 1) {
    const lane = laneMetrics[laneIndex];
    state.laneBounds[lane.type] = {
      x: lane.x,
      width: lane.width,
      minY: sourceStartY,
      maxY: null,
    };
    const defaultSourceX =
      lane.x + Math.max(10, Math.floor((lane.width - sourceNodeWidth) / 2));
    const sources = groupedSources[lane.type] || [];
    const renderGroups = buildLaneRenderGroups(sources);
    let laneCursorY = sourceStartY;
    for (let groupIndex = 0; groupIndex < renderGroups.length; groupIndex += 1) {
      const group = renderGroups[groupIndex];
      const primarySource = group.primary;
      const secondarySource = group.secondary;
      if (!primarySource) continue;

      if (secondarySource) {
        const pairWidth = sourceNodeWidth * 2 + BDPCS3_PAIR_GAP_PX;
        const primaryX =
          lane.x + Math.max(10, Math.floor((lane.width - pairWidth) / 2));
        const secondaryX = primaryX + sourceNodeWidth + BDPCS3_PAIR_GAP_PX;
        const primaryNode = buildFlowNode({
          className: getSourceNodeClassName(primarySource),
          x: primaryX,
          y: laneCursorY,
          title: primarySource.title,
          details: primarySource.details,
          id: primarySource.id,
          extraElement: buildGeometryNodeActions(primarySource),
          outputHandle: buildOutputHandleForSource(primarySource),
          laneType: primarySource.type,
        });
        const secondaryNode = buildFlowNode({
          className: getSourceNodeClassName(secondarySource),
          x: secondaryX,
          y: laneCursorY,
          title: secondarySource.title,
          details: secondarySource.details,
          id: secondarySource.id,
          extraElement: buildGeometryNodeActions(secondarySource),
          outputHandle: buildOutputHandleForSource(secondarySource),
          laneType: secondarySource.type,
        });
        ui.flowNodes.appendChild(primaryNode);
        ui.flowNodes.appendChild(secondaryNode);
        markRigidPairNodeData(
          primaryNode,
          secondaryNode,
          primarySource.bdpcs3PairGroupId || `bdpcs3-pair-${primarySource.id}`
        );
        constrainRigidPairNodes(primaryNode, secondaryNode, lane.type, true);
        validNodeIds.add(primarySource.id);
        validNodeIds.add(secondarySource.id);
        const nodeTop = parseFloat(primaryNode.style.top) || laneCursorY;
        const rowHeight = Math.max(primaryNode.offsetHeight || 120, secondaryNode.offsetHeight || 120);
        maxSourceBottom = Math.max(maxSourceBottom, nodeTop + rowHeight);
        laneCursorY = nodeTop + rowHeight + sourceVerticalGap;
        continue;
      }

      const node = buildFlowNode({
        className: getSourceNodeClassName(primarySource),
        x: defaultSourceX,
        y: laneCursorY,
        title: primarySource.title,
        details: primarySource.details,
        id: primarySource.id,
        extraElement: buildGeometryNodeActions(primarySource),
        outputHandle: buildOutputHandleForSource(primarySource),
        laneType: primarySource.type,
      });
      ui.flowNodes.appendChild(node);
      constrainNodeToLane(node, primarySource.type, true);
      validNodeIds.add(primarySource.id);
      const nodeTop = parseFloat(node.style.top) || laneCursorY;
      const nodeHeight = node.offsetHeight || 120;
      maxSourceBottom = Math.max(maxSourceBottom, nodeTop + nodeHeight);
      laneCursorY = nodeTop + nodeHeight + sourceVerticalGap;
    }
  }

  const laneMinimumHeight = 260;
  const laneGuideBottom = Math.max(maxSourceBottom + 12, sourceStartY + laneMinimumHeight);
  const generatorTopGap = 24;
  const generatorHorizontalMargin = 24;
  const laneAreaWidth = boardWidth - laneLeftPadding * 2;
  const generatorWidth = Math.max(700, laneAreaWidth - generatorHorizontalMargin * 2);
  const generatorX =
    laneLeftPadding + Math.max(0, (laneAreaWidth - generatorWidth) / 2);
  const generatorY = laneGuideBottom + generatorTopGap;
  const generatorHeight = 180;
  const boardHeight = Math.max(560, generatorY + generatorHeight + 56);
  ui.flowBoard.style.width = `${boardWidth}px`;
  ui.flowBoard.style.height = `${boardHeight}px`;

  for (let laneIndex = 0; laneIndex < laneMetrics.length; laneIndex += 1) {
    const lane = laneMetrics[laneIndex];
    const type = lane.type;
    state.laneBounds[type] = {
      x: lane.x,
      width: lane.width,
      minY: sourceStartY,
      maxY: laneGuideBottom,
    };
    renderFlowLaneGuide(type, lane.x, lane.width, 8, laneGuideBottom);
  }

  const laneNodes = ui.flowNodes.querySelectorAll(".flow-node.source[data-lane-type]");
  const constrainedRigidGroups = new Set();
  laneNodes.forEach((node) => {
    const rigidGroupId = node.dataset.rigidGroupId || "";
    const rigidRole = node.dataset.rigidRole || "";
    if (rigidGroupId && rigidRole === "primary" && !constrainedRigidGroups.has(rigidGroupId)) {
      const mateId = node.dataset.rigidMateId;
      const mate = mateId
        ? ui.flowNodes.querySelector(`[data-node-id="${cssEscapeAttr(mateId)}"]`)
        : null;
      if (mate) {
        const laneType = node.dataset.laneType;
        constrainRigidPairNodes(node, mate, laneType, true);
        constrainedRigidGroups.add(rigidGroupId);
        return;
      }
    }
    if (rigidGroupId) return;
    const laneType = node.dataset.laneType;
    constrainNodeToLane(node, laneType, true);
  });

  state.nodePositions["target-generator"] = {
    x: generatorX,
    y: generatorY,
  };
  const generatorNode = buildGeneratorNode({
    id: "target-generator",
    x: generatorX,
    y: generatorY,
    width: generatorWidth,
  });
  ui.flowNodes.appendChild(generatorNode);
  validNodeIds.add("target-generator");

  pruneNodePositions(validNodeIds);
  ensureFlowBoardSize();
  scheduleFlowLineRedraw();
}

function buildLaneRenderGroups(sources) {
  if (!Array.isArray(sources) || !sources.length) return [];
  const sourceById = {};
  for (let i = 0; i < sources.length; i += 1) {
    const source = sources[i];
    if (source && source.id) {
      sourceById[source.id] = source;
    }
  }
  const visited = new Set();
  const groups = [];
  for (let i = 0; i < sources.length; i += 1) {
    const source = sources[i];
    if (!source || visited.has(source.id)) continue;
    if (source.bdpcs3DerivedFromSourceId) {
      const parent = sourceById[source.bdpcs3DerivedFromSourceId];
      if (parent && !visited.has(parent.id)) {
        continue;
      }
    }
    const pairId = source.bdpcs3PairSourceId;
    const pair = pairId ? sourceById[pairId] : null;
    if (pair && pair.bdpcs3DerivedFromSourceId === source.id && !visited.has(pair.id)) {
      visited.add(source.id);
      visited.add(pair.id);
      groups.push({ primary: source, secondary: pair });
      continue;
    }
    visited.add(source.id);
    groups.push({ primary: source, secondary: null });
  }
  return groups;
}

function getSourceNodeClassName(source) {
  let className = "flow-node source";
  const origin = source && source.originType ? String(source.originType).toUpperCase() : "A";
  if (origin === "A") className += " origin-a";
  else if (origin === "B") className += " origin-b";
  else if (origin === "C") className += " origin-c";
  else if (origin === "D") className += " origin-d";
  if (source && source.type === "hessian") {
    className += " hessian-source";
  }
  if (source && source.bdpcs3Derived) {
    className += " bdpcs3-derived-source";
  } else if (source && source.bdpcs3PairSourceId) {
    className += " bdpcs3-parent-source";
  }
  return className;
}

function buildOutputHandleForSource(source) {
  if (!source || !FLOW_CONNECTABLE_TARGET_SET.has(source.type)) {
    return null;
  }
  return {
    sourceId: source.id,
    sourceType: source.type,
  };
}

function buildGeometryNodeActions(source) {
  if (!source || source.type !== "geometry") return null;
  const container = document.createElement("div");
  container.className = "flow-node-actions";
  const button = document.createElement("button");
  button.type = "button";
  button.className = "delta-create-btn";
  button.textContent = "+Δ";
  if (source.bdpcs3Derived) {
    button.disabled = true;
    button.title = "This is already a BDPCS3 corrected geometry.";
  } else if (source.bdpcs3PairSourceId && state.sourceById[source.bdpcs3PairSourceId]) {
    button.disabled = true;
    button.title = "BDPCS3 Δ geometry already created for this source.";
  } else {
    button.title = "Create BDPCS3 corrected geometry (Δ)";
  }
  button.addEventListener("click", (event) => {
    event.preventDefault();
    event.stopPropagation();
    createBDPCS3PairFromSourceId(source.id);
  });
  container.appendChild(button);
  return container;
}

function markRigidPairNodeData(primaryNode, secondaryNode, groupId) {
  if (!primaryNode || !secondaryNode || !groupId) return;
  const primaryId = primaryNode.dataset.nodeId;
  const secondaryId = secondaryNode.dataset.nodeId;
  if (!primaryId || !secondaryId) return;
  primaryNode.dataset.rigidGroupId = groupId;
  primaryNode.dataset.rigidMateId = secondaryId;
  primaryNode.dataset.rigidRole = "primary";
  secondaryNode.dataset.rigidGroupId = groupId;
  secondaryNode.dataset.rigidMateId = primaryId;
  secondaryNode.dataset.rigidRole = "secondary";
}

function constrainRigidPairNodes(primaryNode, secondaryNode, laneType, persistState = false) {
  if (!primaryNode || !secondaryNode || !laneType) return;
  const bounds = state.laneBounds[laneType];
  if (!bounds) return;
  const primaryId = primaryNode.dataset.nodeId;
  const secondaryId = secondaryNode.dataset.nodeId;
  const primaryWidth = primaryNode.offsetWidth || 188;
  const secondaryWidth = secondaryNode.offsetWidth || 188;
  const primaryHeight = primaryNode.offsetHeight || 120;
  const secondaryHeight = secondaryNode.offsetHeight || 120;
  const pairOffsetX = primaryWidth + BDPCS3_PAIR_GAP_PX;
  const pairWidth = pairOffsetX + secondaryWidth;
  const minX = bounds.x + 10;
  const maxX = Math.max(minX, bounds.x + bounds.width - pairWidth - 10);
  const primaryRawX = parseFloat(primaryNode.style.left);
  const secondaryRawX = parseFloat(secondaryNode.style.left);
  let primaryX = Number.isFinite(primaryRawX)
    ? primaryRawX
    : Number.isFinite(secondaryRawX)
      ? secondaryRawX - pairOffsetX
      : minX;
  primaryX = Math.min(maxX, Math.max(minX, primaryX));

  const primaryRawY = parseFloat(primaryNode.style.top);
  const secondaryRawY = parseFloat(secondaryNode.style.top);
  let top = Number.isFinite(primaryRawY)
    ? primaryRawY
    : Number.isFinite(secondaryRawY)
      ? secondaryRawY
      : bounds.minY || 0;
  const minY = bounds.minY || 0;
  const pairHeight = Math.max(primaryHeight, secondaryHeight);
  const maxY = Number.isFinite(bounds.maxY)
    ? Math.max(minY, bounds.maxY - pairHeight - 8)
    : null;
  const clampedTop =
    typeof maxY === "number" ? Math.min(maxY, Math.max(minY, top)) : Math.max(minY, top);

  primaryNode.style.left = `${primaryX}px`;
  primaryNode.style.top = `${clampedTop}px`;
  secondaryNode.style.left = `${primaryX + pairOffsetX}px`;
  secondaryNode.style.top = `${clampedTop}px`;

  if (persistState && primaryId) {
    state.nodePositions[primaryId] = { x: primaryX, y: clampedTop };
  }
  if (persistState && secondaryId) {
    state.nodePositions[secondaryId] = { x: primaryX + pairOffsetX, y: clampedTop };
  }
}

function renderFlowLaneGuide(type, laneX, laneWidth, top = 8, bottom = null) {
  if (!ui.flowLanes) return;
  const lane = document.createElement("div");
  lane.className = "flow-lane";
  lane.style.left = `${laneX}px`;
  lane.style.width = `${laneWidth}px`;
  lane.style.top = `${Math.max(0, top)}px`;
  if (typeof bottom === "number" && Number.isFinite(bottom)) {
    lane.style.bottom = "auto";
    lane.style.height = `${Math.max(120, bottom - top)}px`;
  }
  const label = document.createElement("span");
  label.className = "flow-lane-label";
  label.textContent = FLOW_TARGET_LABELS[type];
  lane.appendChild(label);
  ui.flowLanes.appendChild(lane);
}

function buildGeneratorNode({ id, x, y, width }) {
  const rows = document.createElement("div");
  rows.className = "flow-generator-rows";
  for (let i = 0; i < FLOW_TARGET_ROW_ORDER.length; i += 1) {
    const type = FLOW_TARGET_ROW_ORDER[i];
    const connected = Boolean(getConnectedSource(type));
    const row = document.createElement("div");
    row.className = `flow-generator-row${connected ? " connected" : ""}`;
    row.dataset.targetType = type;
    const handle = document.createElement("span");
    handle.className = `flow-handle input${connected ? " connected" : ""}`;
    handle.title = `${FLOW_TARGET_LABELS[type]} input`;
    handle.dataset.handleKind = "input";
    handle.dataset.targetType = type;
    attachFlowHandleEvents(handle);
    const text = document.createElement("span");
    text.textContent = FLOW_TARGET_LABELS[type];
    row.appendChild(handle);
    row.appendChild(text);
    rows.appendChild(row);
  }
  const note = buildFlowGeneratorRepresentationNote();
  return buildFlowNode({
    className: "flow-node generator flow-generator",
    x,
    y,
    width,
    title: "Input Generator",
    details: note,
    id,
    extraElement: rows,
    draggable: false,
  });
}

function buildFlowNode({
  className,
  x,
  y,
  width = null,
  title,
  details,
  id,
  extraElement = null,
  inputHandle = null,
  outputHandle = null,
  laneType = null,
  draggable = true,
}) {
  const node = document.createElement("article");
  node.className = className;
  node.dataset.nodeId = id;
  if (laneType) {
    node.dataset.laneType = laneType;
  }
  const position = resolveNodePosition(id, x, y);
  node.style.left = `${position.x}px`;
  node.style.top = `${position.y}px`;
  if (typeof width === "number" && Number.isFinite(width) && width > 0) {
    node.style.width = `${width}px`;
  }

  const detailText = String(details || "").trim();
  if (detailText) {
    const expandButton = document.createElement("button");
    expandButton.type = "button";
    expandButton.className = "flow-node-expand-btn";
    expandButton.textContent = "+";
    expandButton.title = "View details";
    expandButton.setAttribute("aria-label", `View details for ${String(title || "panel")}`);
    expandButton.addEventListener("click", (event) => {
      event.preventDefault();
      event.stopPropagation();
      openFlowNodeModal(title, detailText);
    });
    node.appendChild(expandButton);
  }

  const heading = document.createElement("h4");
  heading.textContent = title;
  heading.title = String(title || "");
  node.appendChild(heading);

  if (extraElement) {
    node.appendChild(extraElement);
  }
  if (inputHandle) {
    const handle = document.createElement("span");
    handle.className = `flow-handle input${inputHandle.connected ? " connected" : ""}`;
    handle.title = `Input ${FLOW_TARGET_LABELS[inputHandle.targetType]}`;
    handle.dataset.handleKind = "input";
    handle.dataset.targetType = inputHandle.targetType;
    attachFlowHandleEvents(handle);
    node.appendChild(handle);
  }
  if (outputHandle) {
    const handle = document.createElement("span");
    handle.className = "flow-handle output";
    handle.title = `Output ${FLOW_TARGET_LABELS[outputHandle.sourceType]}`;
    handle.dataset.handleKind = "output";
    handle.dataset.sourceId = outputHandle.sourceId;
    handle.dataset.sourceType = outputHandle.sourceType;
    attachFlowHandleEvents(handle);
    node.appendChild(handle);
  }
  if (draggable) {
    makeNodeDraggable(node);
  }
  return node;
}

function resolveNodePosition(nodeId, fallbackX, fallbackY) {
  if (!Object.prototype.hasOwnProperty.call(state.nodePositions, nodeId)) {
    state.nodePositions[nodeId] = { x: fallbackX, y: fallbackY };
  }
  return state.nodePositions[nodeId];
}

function openFlowNodeModal(title, details) {
  if (!ui.flowNodeModal || !ui.flowNodeModalTitle || !ui.flowNodeModalBody) return;
  const safeTitle = String(title || "Panel details");
  const safeDetails = String(details || "No additional details available.");
  ui.flowNodeModalTitle.textContent = safeTitle;
  ui.flowNodeModalBody.textContent = `${safeTitle}\n\n${safeDetails}`;
  ui.flowNodeModal.hidden = false;
  ui.flowNodeModal.setAttribute("aria-hidden", "false");
  if (typeof document !== "undefined" && document.body) {
    document.body.classList.add("modal-open");
  }
}

function closeFlowNodeModal() {
  if (!ui.flowNodeModal) return;
  ui.flowNodeModal.hidden = true;
  ui.flowNodeModal.setAttribute("aria-hidden", "true");
  if (typeof document !== "undefined" && document.body) {
    document.body.classList.remove("modal-open");
  }
}

function onGlobalKeyDown(event) {
  if (!event || event.key !== "Escape") return;
  if (ui.flowNodeModal && !ui.flowNodeModal.hidden) {
    closeFlowNodeModal();
  }
}

function pruneNodePositions(validNodeIds) {
  const ids = Object.keys(state.nodePositions);
  for (let i = 0; i < ids.length; i += 1) {
    if (!validNodeIds.has(ids[i])) {
      delete state.nodePositions[ids[i]];
    }
  }
}

function constrainNodeToLane(node, laneType, persistState = false) {
  if (!node || !laneType) return;
  const bounds = state.laneBounds[laneType];
  if (!bounds) return;
  const nodeId = node.dataset.nodeId;
  const nodeWidth = node.offsetWidth || 188;
  const nodeHeight = node.offsetHeight || 120;
  const minX = bounds.x + 10;
  const maxX = Math.max(minX, bounds.x + bounds.width - nodeWidth - 10);
  const currentX = parseFloat(node.style.left) || minX;
  const clampedX = Math.min(maxX, Math.max(minX, currentX));
  node.style.left = `${clampedX}px`;

  const currentY = parseFloat(node.style.top) || bounds.minY || 0;
  const minY = bounds.minY || 0;
  const maxY = Number.isFinite(bounds.maxY)
    ? Math.max(minY, bounds.maxY - nodeHeight - 8)
    : null;
  const clampedY =
    typeof maxY === "number" ? Math.min(maxY, Math.max(minY, currentY)) : Math.max(minY, currentY);
  node.style.top = `${clampedY}px`;

  if (persistState && nodeId) {
    state.nodePositions[nodeId] = { x: clampedX, y: clampedY };
  }
}

function ensureFlowBoardSize() {
  if (!ui.flowBoard || !ui.flowNodes) return;
  let maxRight = 920;
  let maxBottom = 420;
  const currentWidth = parseFloat(ui.flowBoard.style.width) || 0;
  const currentHeight = parseFloat(ui.flowBoard.style.height) || 0;
  const wrapWidth = ui.flowBoardWrap ? Math.floor(ui.flowBoardWrap.clientWidth) : 0;
  if (currentWidth > 0) maxRight = Math.max(maxRight, currentWidth);
  if (currentHeight > 0) maxBottom = Math.max(maxBottom, currentHeight);
  if (wrapWidth > 0) maxRight = Math.max(maxRight, wrapWidth - 2);
  const nodes = ui.flowNodes.querySelectorAll(".flow-node");
  nodes.forEach((node) => {
    const left = parseFloat(node.style.left) || 0;
    const top = parseFloat(node.style.top) || 0;
    const width = node.offsetWidth || 188;
    const height = node.offsetHeight || 120;
    maxRight = Math.max(maxRight, left + width + 40);
    maxBottom = Math.max(maxBottom, top + height + 40);
  });
  ui.flowBoard.style.width = `${Math.ceil(maxRight)}px`;
  ui.flowBoard.style.height = `${Math.ceil(maxBottom)}px`;
}

function makeNodeDraggable(node) {
  if (!node || !ui.flowBoard) return;

  node.addEventListener("pointerdown", (event) => {
    if (event.pointerType === "mouse" && event.button !== 0) return;
    if (!event.isPrimary) return;
    if (
      event.target &&
      typeof event.target.closest === "function" &&
      event.target.closest("select, option, input, textarea, button, a, .flow-handle")
    ) {
      return;
    }

    const nodeId = node.dataset.nodeId;
    const laneType = node.dataset.laneType || null;
    if (!nodeId) return;

    const rigidGroupId = node.dataset.rigidGroupId || "";
    const rigidMateId = node.dataset.rigidMateId || "";
    const rigidMateNode =
      rigidGroupId && rigidMateId && ui.flowNodes
        ? ui.flowNodes.querySelector(`[data-node-id="${cssEscapeAttr(rigidMateId)}"]`)
        : null;
    const nodesForDrag = [node];
    if (
      rigidMateNode &&
      rigidMateNode !== node &&
      rigidMateNode.dataset.rigidGroupId === rigidGroupId
    ) {
      nodesForDrag.push(rigidMateNode);
    }

    const dragEntries = [];
    const seenNodeIds = new Set();
    for (let i = 0; i < nodesForDrag.length; i += 1) {
      const dragNode = nodesForDrag[i];
      if (!dragNode) continue;
      const dragNodeId = dragNode.dataset.nodeId;
      if (!dragNodeId || seenNodeIds.has(dragNodeId)) continue;
      seenNodeIds.add(dragNodeId);
      dragEntries.push({
        id: dragNodeId,
        node: dragNode,
        startLeft: parseFloat(dragNode.style.left) || 0,
        startTop: parseFloat(dragNode.style.top) || 0,
        width: dragNode.offsetWidth || 188,
        height: dragNode.offsetHeight || 120,
        offsetLeft: 0,
        offsetTop: 0,
      });
    }
    if (!dragEntries.length) return;
    const minStartLeft = Math.min(...dragEntries.map((entry) => entry.startLeft));
    const minStartTop = Math.min(...dragEntries.map((entry) => entry.startTop));
    for (let i = 0; i < dragEntries.length; i += 1) {
      dragEntries[i].offsetLeft = dragEntries[i].startLeft - minStartLeft;
      dragEntries[i].offsetTop = dragEntries[i].startTop - minStartTop;
    }
    const groupWidth = Math.max(...dragEntries.map((entry) => entry.offsetLeft + entry.width));
    const groupHeight = Math.max(...dragEntries.map((entry) => entry.offsetTop + entry.height));

    const startX = event.clientX;
    const startY = event.clientY;
    for (let i = 0; i < dragEntries.length; i += 1) {
      dragEntries[i].node.classList.add("dragging");
    }
    event.preventDefault();

    if (typeof node.setPointerCapture === "function") {
      try {
        node.setPointerCapture(event.pointerId);
      } catch (error) {
        // no-op: keep dragging even if capture fails.
      }
    }

    const onPointerMove = (moveEvent) => {
      const dx = moveEvent.clientX - startX;
      const dy = moveEvent.clientY - startY;

      let groupLeft = Math.max(8, minStartLeft + dx);
      let groupTop = Math.max(8, minStartTop + dy);
      if (laneType) {
        const bounds = state.laneBounds[laneType];
        if (bounds) {
          const minX = bounds.x + 10;
          const maxX = Math.max(minX, bounds.x + bounds.width - groupWidth - 10);
          groupLeft = Math.min(maxX, Math.max(minX, groupLeft));
          const minY = bounds.minY || 0;
          const maxY = Number.isFinite(bounds.maxY)
            ? Math.max(minY, bounds.maxY - groupHeight - 8)
            : null;
          groupTop =
            typeof maxY === "number" ? Math.min(maxY, Math.max(minY, groupTop)) : Math.max(minY, groupTop);
        }
      }

      for (let i = 0; i < dragEntries.length; i += 1) {
        const entry = dragEntries[i];
        const left = groupLeft + entry.offsetLeft;
        const top = groupTop + entry.offsetTop;
        entry.node.style.left = `${left}px`;
        entry.node.style.top = `${top}px`;
        state.nodePositions[entry.id] = { x: left, y: top };
      }

      ensureFlowBoardSize();
      scheduleFlowLineRedraw();
    };

    const onPointerEnd = (endEvent) => {
      for (let i = 0; i < dragEntries.length; i += 1) {
        dragEntries[i].node.classList.remove("dragging");
      }
      node.removeEventListener("pointermove", onPointerMove);
      node.removeEventListener("pointerup", onPointerEnd);
      node.removeEventListener("pointercancel", onPointerEnd);
      if (typeof node.releasePointerCapture === "function") {
        try {
          node.releasePointerCapture(endEvent.pointerId);
        } catch (error) {
          // no-op
        }
      }
    };

    node.addEventListener("pointermove", onPointerMove);
    node.addEventListener("pointerup", onPointerEnd);
    node.addEventListener("pointercancel", onPointerEnd);
  });
}

function attachFlowHandleEvents(handle) {
  if (!handle) return;
  handle.addEventListener("pointerdown", (event) => {
    if (event.pointerType === "mouse" && event.button !== 0) return;
    if (!event.isPrimary) return;
    const handleKind = handle.dataset.handleKind;
    if (handleKind === "output") {
      startActiveLinkFromOutput(handle, event);
      return;
    }
    if (handleKind === "input") {
      startActiveLinkFromInput(handle, event);
    }
  });
}

function startActiveLinkFromOutput(handle, event) {
  const sourceId = handle.dataset.sourceId;
  const sourceType = handle.dataset.sourceType;
  if (!sourceId || !sourceType) return;
  beginActiveLinkDrag({
    pointerId: event.pointerId,
    mode: "from-output",
    sourceId,
    sourceType,
    targetType: null,
    startHandle: handle,
    startPoint: getHandleBoardPoint(handle),
  });
  event.preventDefault();
}

function startActiveLinkFromInput(handle, event) {
  const targetType = handle.dataset.targetType;
  if (!targetType) return;
  beginActiveLinkDrag({
    pointerId: event.pointerId,
    mode: "from-input",
    sourceId: state.connections[targetType] || null,
    sourceType: targetType,
    targetType,
    startHandle: handle,
    startPoint: getHandleBoardPoint(handle),
  });
  event.preventDefault();
}

function beginActiveLinkDrag({
  pointerId,
  mode,
  sourceId,
  sourceType,
  targetType,
  startHandle,
  startPoint,
}) {
  if (!ui.flowLines || !ui.flowBoard || !startPoint) return;
  clearActiveLinkDrag();
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("fill", "none");
  path.setAttribute("stroke", "#1f7dd6");
  path.setAttribute("stroke-width", "2");
  path.setAttribute("stroke-linecap", "round");
  path.setAttribute("stroke-dasharray", "4 4");
  path.setAttribute("opacity", "0.85");
  ui.flowLines.appendChild(path);

  state.activeLinkDrag = {
    pointerId,
    mode,
    sourceId,
    sourceType,
    targetType,
    startHandle,
    path,
    startPoint,
    startClientX: null,
    startClientY: null,
    moved: false,
    lastClientX: null,
    lastClientY: null,
    moveHandler: null,
    upHandler: null,
  };

  const onPointerMove = (moveEvent) => {
    if (!state.activeLinkDrag || moveEvent.pointerId !== pointerId) return;
    if (state.activeLinkDrag.startClientX === null || state.activeLinkDrag.startClientY === null) {
      state.activeLinkDrag.startClientX = moveEvent.clientX;
      state.activeLinkDrag.startClientY = moveEvent.clientY;
    }
    state.activeLinkDrag.lastClientX = moveEvent.clientX;
    state.activeLinkDrag.lastClientY = moveEvent.clientY;
    if (
      Math.abs(moveEvent.clientX - state.activeLinkDrag.startClientX) > 3 ||
      Math.abs(moveEvent.clientY - state.activeLinkDrag.startClientY) > 3
    ) {
      state.activeLinkDrag.moved = true;
    }
    const point = clientToBoardPoint(moveEvent.clientX, moveEvent.clientY);
    if (!point) return;
    drawLinkPath(state.activeLinkDrag.path, state.activeLinkDrag.startPoint, point);
  };

  const onPointerEnd = (endEvent) => {
    if (!state.activeLinkDrag || endEvent.pointerId !== pointerId) return;
    completeActiveLinkDrag(endEvent);
  };

  state.activeLinkDrag.moveHandler = onPointerMove;
  state.activeLinkDrag.upHandler = onPointerEnd;
  window.addEventListener("pointermove", onPointerMove);
  window.addEventListener("pointerup", onPointerEnd);
  window.addEventListener("pointercancel", onPointerEnd);

  if (typeof startHandle.setPointerCapture === "function") {
    try {
      startHandle.setPointerCapture(pointerId);
    } catch (error) {
      // no-op
    }
  }
}

function completeActiveLinkDrag(event) {
  const active = state.activeLinkDrag;
  if (!active) return;
  const targetEl = document.elementFromPoint(event.clientX, event.clientY);
  let didConnect = false;

  if (active.mode === "from-output") {
    if (!active.moved) {
      clearActiveLinkDrag();
      scheduleFlowLineRedraw();
      return;
    }
    const inputHandle = targetEl && targetEl.closest
      ? targetEl.closest('.flow-handle.input[data-handle-kind="input"]')
      : null;
    if (inputHandle) {
      const targetType = inputHandle.dataset.targetType;
      if (targetType && targetType === active.sourceType) {
        connectSourceToTarget(active.sourceId, targetType);
        didConnect = true;
      }
    }
  } else if (active.mode === "from-input") {
    const outputHandle = targetEl && targetEl.closest
      ? targetEl.closest('.flow-handle.output[data-handle-kind="output"]')
      : null;
    if (outputHandle) {
      const sourceId = outputHandle.dataset.sourceId;
      const sourceType = outputHandle.dataset.sourceType;
      if (sourceId && sourceType && sourceType === active.targetType) {
        connectSourceToTarget(sourceId, active.targetType);
        didConnect = true;
      }
    }
    if (!didConnect && active.targetType && active.moved) {
      disconnectTarget(active.targetType);
      didConnect = true;
    }
  }

  clearActiveLinkDrag();
  if (didConnect) {
    renderFlowBoard();
    void refreshWMSInputOutput();
  } else {
    scheduleFlowLineRedraw();
  }
}

function clearActiveLinkDrag() {
  const active = state.activeLinkDrag;
  if (!active) return;
  if (active.moveHandler) {
    window.removeEventListener("pointermove", active.moveHandler);
    window.removeEventListener("pointerup", active.upHandler);
    window.removeEventListener("pointercancel", active.upHandler);
  }
  if (active.path && active.path.parentNode) {
    active.path.parentNode.removeChild(active.path);
  }
  if (active.startHandle && typeof active.startHandle.releasePointerCapture === "function") {
    try {
      active.startHandle.releasePointerCapture(active.pointerId);
    } catch (error) {
      // no-op
    }
  }
  state.activeLinkDrag = null;
}

function connectSourceToTarget(sourceId, targetType) {
  if (!sourceId || !targetType) return;
  const source = state.sourceById[sourceId];
  if (!source || source.type !== targetType) return;
  state.connections[targetType] = sourceId;
}

function disconnectTarget(targetType) {
  if (!targetType) return;
  state.connections[targetType] = null;
}

function getHandleBoardPoint(handle) {
  if (!handle || !ui.flowBoard) return null;
  const rect = handle.getBoundingClientRect();
  return clientToBoardPoint(rect.left + rect.width / 2, rect.top + rect.height / 2);
}

function clientToBoardPoint(clientX, clientY) {
  if (!ui.flowBoard) return null;
  const boardRect = ui.flowBoard.getBoundingClientRect();
  const scrollX = ui.flowBoard.scrollLeft || 0;
  const scrollY = ui.flowBoard.scrollTop || 0;
  return {
    x: clientX - boardRect.left + scrollX,
    y: clientY - boardRect.top + scrollY,
  };
}

function drawLinkPath(path, fromPoint, toPoint) {
  if (!path || !fromPoint || !toPoint) return;
  const dx = toPoint.x - fromPoint.x;
  const dy = toPoint.y - fromPoint.y;
  if (Math.abs(dy) > Math.abs(dx)) {
    const ctrlY = Math.max(42, Math.abs(dy) * 0.45);
    const dirY = dy >= 0 ? 1 : -1;
    path.setAttribute(
      "d",
      `M ${fromPoint.x} ${fromPoint.y} C ${fromPoint.x} ${fromPoint.y + dirY * ctrlY}, ${toPoint.x} ${toPoint.y - dirY * ctrlY}, ${toPoint.x} ${toPoint.y}`
    );
    return;
  }
  const ctrlX = Math.max(36, Math.abs(dx) * 0.45);
  const dirX = dx >= 0 ? 1 : -1;
  path.setAttribute(
    "d",
    `M ${fromPoint.x} ${fromPoint.y} C ${fromPoint.x + dirX * ctrlX} ${fromPoint.y}, ${toPoint.x - dirX * ctrlX} ${toPoint.y}, ${toPoint.x} ${toPoint.y}`
  );
}

function scheduleFlowLineRedraw() {
  if (typeof window === "undefined") return;
  window.requestAnimationFrame(drawFlowLines);
}

function drawFlowLines() {
  if (!ui.flowBoard || !ui.flowNodes || !ui.flowLines) return;
  const activeTempPath = state.activeLinkDrag ? state.activeLinkDrag.path : null;
  ui.flowLines.innerHTML = "";
  appendFlowLineDefs();
  if (activeTempPath) {
    ui.flowLines.appendChild(activeTempPath);
  }
  drawRigidPairArrows();

  const generator = ui.flowNodes.querySelector('[data-node-id="target-generator"]');
  if (!generator) return;

  for (let i = 0; i < FLOW_TARGET_ROW_ORDER.length; i += 1) {
    const type = FLOW_TARGET_ROW_ORDER[i];
    const targetInput = ui.flowNodes.querySelector(
      `.flow-handle.input[data-handle-kind="input"][data-target-type="${type}"]`
    );
    if (!targetInput) continue;
    const sourceId = state.connections[type];
    if (!sourceId) continue;
    const sourceHandle = ui.flowNodes.querySelector(
      `.flow-handle.output[data-handle-kind="output"][data-source-id="${cssEscapeAttr(sourceId)}"]`
    );
    if (sourceHandle && targetInput) {
      drawLineBetweenHandles(sourceHandle, targetInput, "#0f67b8");
      continue;
    }
    const sourceNode = ui.flowNodes.querySelector(`[data-node-id="${cssEscapeAttr(sourceId)}"]`);
    if (!sourceNode) continue;
    drawLineBetweenNodes(sourceNode, generator, "#0f67b8");
  }
}

function appendFlowLineDefs() {
  if (!ui.flowLines) return;
  const defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
  const marker = document.createElementNS("http://www.w3.org/2000/svg", "marker");
  marker.setAttribute("id", "flow-rigid-arrowhead");
  marker.setAttribute("viewBox", "0 0 10 10");
  marker.setAttribute("refX", "8");
  marker.setAttribute("refY", "5");
  marker.setAttribute("markerWidth", "7");
  marker.setAttribute("markerHeight", "7");
  marker.setAttribute("orient", "auto");
  marker.setAttribute("markerUnits", "strokeWidth");
  const arrowPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
  arrowPath.setAttribute("d", "M 0 0 L 10 5 L 0 10 z");
  arrowPath.setAttribute("fill", "#2a7a61");
  marker.appendChild(arrowPath);
  defs.appendChild(marker);
  ui.flowLines.appendChild(defs);
}

function drawRigidPairArrows() {
  if (!ui.flowNodes || !ui.flowLines) return;
  const primaryNodes = ui.flowNodes.querySelectorAll(
    '.flow-node.source[data-rigid-role="primary"][data-rigid-mate-id]'
  );
  primaryNodes.forEach((primaryNode) => {
    const mateId = primaryNode.dataset.rigidMateId;
    if (!mateId) return;
    const secondaryNode = ui.flowNodes.querySelector(
      `[data-node-id="${cssEscapeAttr(mateId)}"]`
    );
    if (!secondaryNode) return;
    drawRigidPairArrowBetweenNodes(primaryNode, secondaryNode);
  });
}

function drawRigidPairArrowBetweenNodes(fromNode, toNode) {
  if (!ui.flowBoard || !ui.flowLines || !fromNode || !toNode) return;
  const boardRect = ui.flowBoard.getBoundingClientRect();
  const fromRect = fromNode.getBoundingClientRect();
  const toRect = toNode.getBoundingClientRect();
  const scrollX = ui.flowBoard ? ui.flowBoard.scrollLeft : 0;
  const scrollY = ui.flowBoard ? ui.flowBoard.scrollTop : 0;
  const x1 = fromRect.right - boardRect.left + scrollX + 3;
  const x2 = toRect.left - boardRect.left + scrollX - 3;
  const y = fromRect.top + fromRect.height / 2 - boardRect.top + scrollY;
  if (!(x2 > x1 + 10)) return;

  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("fill", "none");
  path.setAttribute("stroke", "#2a7a61");
  path.setAttribute("stroke-width", "2.4");
  path.setAttribute("stroke-linecap", "round");
  path.setAttribute("opacity", "0.94");
  path.setAttribute("marker-end", "url(#flow-rigid-arrowhead)");
  path.setAttribute("d", `M ${x1} ${y} L ${x2} ${y}`);
  ui.flowLines.appendChild(path);

  const delta = document.createElementNS("http://www.w3.org/2000/svg", "text");
  delta.setAttribute("x", `${(x1 + x2) / 2}`);
  delta.setAttribute("y", `${y - 8}`);
  delta.setAttribute("fill", "#174f3d");
  delta.setAttribute("font-size", "15");
  delta.setAttribute("font-weight", "700");
  delta.setAttribute("text-anchor", "middle");
  delta.textContent = "Δ";
  ui.flowLines.appendChild(delta);
}

function drawLineBetweenHandles(fromHandle, toHandle, strokeColor) {
  const fromPoint = getHandleBoardPoint(fromHandle);
  const toPoint = getHandleBoardPoint(toHandle);
  if (!fromPoint || !toPoint) return;
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("fill", "none");
  path.setAttribute("stroke", strokeColor);
  path.setAttribute("stroke-width", "2");
  path.setAttribute("stroke-linecap", "round");
  path.setAttribute("stroke-linejoin", "round");
  path.setAttribute("opacity", "0.92");
  drawLinkPath(path, fromPoint, toPoint);
  ui.flowLines.appendChild(path);
}

function drawLineBetweenNodes(fromNode, toNode, strokeColor) {
  if (!ui.flowLines) return;
  const boardRect = ui.flowBoard.getBoundingClientRect();
  const fromRect = fromNode.getBoundingClientRect();
  const toRect = toNode.getBoundingClientRect();
  const scrollX = ui.flowBoard ? ui.flowBoard.scrollLeft : 0;
  const scrollY = ui.flowBoard ? ui.flowBoard.scrollTop : 0;
  const x1 = fromRect.right - boardRect.left + scrollX;
  const y1 = fromRect.top + fromRect.height / 2 - boardRect.top + scrollY;
  const x2 = toRect.left - boardRect.left + scrollX;
  const y2 = toRect.top + toRect.height / 2 - boardRect.top + scrollY;
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("fill", "none");
  path.setAttribute("stroke", strokeColor);
  path.setAttribute("stroke-width", "2");
  path.setAttribute("stroke-linecap", "round");
  path.setAttribute("stroke-linejoin", "round");
  path.setAttribute("opacity", "0.92");
  drawLinkPath(path, { x: x1, y: y1 }, { x: x2, y: y2 });
  ui.flowLines.appendChild(path);
}

function cssEscapeAttr(value) {
  const text = String(value || "");
  return text.replace(/["\\]/g, "\\$&");
}

function parseXYZ(text) {
  const clean = text.replace(/^\uFEFF/, "");
  const lines = clean.split(/\r?\n/);
  const atoms = [];

  const first = (lines[0] || "").trim();
  const nAtoms = Number.parseInt(first, 10);
  if (/^[+-]?\d+$/.test(first) && Number.isFinite(nAtoms) && nAtoms > 0) {
    const comment = lines[1] || "";
    for (let i = 2; i < lines.length && atoms.length < nAtoms; i += 1) {
      const parsed = parseAtomLine(lines[i], i + 1, true);
      if (parsed) atoms.push(parsed);
    }
    if (atoms.length !== nAtoms) {
      throw new Error(
        `Invalid XYZ file: expected ${nAtoms} atoms, found ${atoms.length}.`
      );
    }
    return { atoms, comment };
  }

  // Fallback: parse all lines that look like atom rows.
  for (let i = 0; i < lines.length; i += 1) {
    const parsed = parseAtomLine(lines[i], i + 1, false);
    if (parsed) atoms.push(parsed);
  }
  if (!atoms.length) {
    throw new Error("No valid atomic coordinates found in the file.");
  }
  return { atoms, comment: "" };
}

function parseAtomLine(line, lineNo, strict) {
  const trimmed = String(line || "").trim();
  if (!trimmed) return null;

  const parts = trimmed.split(/\s+/);
  if (parts.length < 4) {
    if (!strict) return null;
    throw new Error(`Line ${lineNo}: invalid atomic row format.`);
  }

  const rawLabel = parts[0];
  let atomicNumber;
  if (/^[+-]?\d+$/.test(rawLabel)) {
    atomicNumber = Number.parseInt(rawLabel, 10);
  } else {
    atomicNumber = lookupAtomicNumberByLabel(rawLabel);
  }

  if (!Number.isFinite(atomicNumber) || atomicNumber <= 0) {
    if (!strict) return null;
    throw new Error(`Line ${lineNo}: unrecognized element "${rawLabel}".`);
  }

  const x = Number.parseFloat(parts[1]);
  const y = Number.parseFloat(parts[2]);
  const z = Number.parseFloat(parts[3]);
  if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
    if (!strict) return null;
    throw new Error(`Line ${lineNo}: non-numeric coordinates.`);
  }

  return {
    label: rawLabel,
    z: atomicNumber,
    coords: [x, y, z],
  };
}

function lookupAtomicNumberByLabel(label) {
  const token = String(label || "").trim().toUpperCase();
  if (!token) return null;
  if (typeof atomicNumbersTable !== "undefined" && atomicNumbersTable) {
    const value = atomicNumbersTable[token];
    if (Number.isFinite(value) && value > 0) return Number(value);
  }
  const fallback = {
    H: 1, HE: 2, LI: 3, BE: 4, B: 5, C: 6, N: 7, O: 8, F: 9, NE: 10,
    NA: 11, MG: 12, AL: 13, SI: 14, P: 15, S: 16, CL: 17, AR: 18, K: 19, CA: 20,
    SC: 21, TI: 22, V: 23, CR: 24, MN: 25, FE: 26, CO: 27, NI: 28, CU: 29, ZN: 30,
    GA: 31, GE: 32, AS: 33, SE: 34, BR: 35, KR: 36, RB: 37, SR: 38, Y: 39, ZR: 40,
    NB: 41, MO: 42, TC: 43, RU: 44, RH: 45, PD: 46, AG: 47, CD: 48, IN: 49, SN: 50,
    SB: 51, TE: 52, I: 53, XE: 54, CS: 55, BA: 56, LA: 57, CE: 58, PR: 59, ND: 60,
    PM: 61, SM: 62, EU: 63, GD: 64, TB: 65, DY: 66, HO: 67, ER: 68, TM: 69, YB: 70,
    LU: 71, HF: 72, TA: 73, W: 74, RE: 75, OS: 76, IR: 77, PT: 78, AU: 79, HG: 80,
    TL: 81, PB: 82, BI: 83, PO: 84, AT: 85, RN: 86, FR: 87, RA: 88, AC: 89, TH: 90,
    PA: 91, U: 92, NP: 93, PU: 94, AM: 95, CM: 96, BK: 97, CF: 98, ES: 99, FM: 100,
    MD: 101, NO: 102, LR: 103, RF: 104, DB: 105, SG: 106, BH: 107, HS: 108, MT: 109,
    DS: 110, RG: 111, CN: 112, NH: 113, FL: 114, MC: 115, LV: 116, TS: 117, OG: 118,
    UUT: 113, UUP: 115, UUS: 117, UUO: 118,
  };
  return Object.prototype.hasOwnProperty.call(fallback, token) ? fallback[token] : null;
}

function computePrincipalFrame() {
  if (!state.loadedAtoms || !state.loadedAtoms.length) return;

  const center = computeCenterOfMass(state.loadedAtoms);
  const inertiaTensor = computeInertiaTensor(state.loadedAtoms, center);
  const { moments, basis } = diagonalizeSymmetricTensor(inertiaTensor);

  const principalAtoms = state.loadedAtoms.map((atom) => {
    const shifted = subtractVec(atom.coords, center);
    return {
      label: atom.label,
      z: atom.z,
      coords: rotateToBasis(shifted, basis),
    };
  });

  state.principalAtoms = principalAtoms;
  const inputRepresentationDetection = detectInputRepresentationFromPrincipalBasis({
    basis,
    moments,
  });
  state.analysisData = {
    center,
    inertiaTensor,
    moments,
    basis,
    inputRepresentationDetection,
  };
}

function detectInputRepresentationFromPrincipalBasis({ basis, moments }) {
  if (!Array.isArray(basis) || basis.length !== 3) {
    return { representation: null, reason: "missing basis", alignment: null, axisMap: null };
  }
  const relDegeneracy = computeMomentRelativeDegeneracy(moments);
  if (Number.isFinite(relDegeneracy) && relDegeneracy < INPUT_REP_MOMENT_DEGENERACY_REL) {
    return {
      representation: null,
      reason: "near-degenerate inertia moments",
      alignment: null,
      axisMap: null,
    };
  }

  const axisMap = { x: null, y: null, z: null };
  const usedCartesian = new Set();
  let minAlignment = Number.POSITIVE_INFINITY;

  for (let i = 0; i < PRINCIPAL_AXIS_LABELS.length; i += 1) {
    const principalAxis = PRINCIPAL_AXIS_LABELS[i];
    const vec = Array.isArray(basis[i]) ? basis[i] : null;
    if (!vec || vec.length < 3) {
      return { representation: null, reason: "invalid basis vector", alignment: null, axisMap: null };
    }
    let bestCartesianIndex = 0;
    let bestAbs = Math.abs(vec[0]);
    for (let j = 1; j < 3; j += 1) {
      const absVal = Math.abs(vec[j]);
      if (absVal > bestAbs) {
        bestAbs = absVal;
        bestCartesianIndex = j;
      }
    }
    if (usedCartesian.has(bestCartesianIndex)) {
      return {
        representation: null,
        reason: "axes are not uniquely aligned",
        alignment: bestAbs,
        axisMap: null,
      };
    }
    usedCartesian.add(bestCartesianIndex);
    const cartesianAxis = CARTESIAN_AXIS_LABELS[bestCartesianIndex];
    axisMap[cartesianAxis] = principalAxis;
    minAlignment = Math.min(minAlignment, bestAbs);
  }

  const representation = inferRepresentationFromAxisMap(axisMap);
  if (!representation) {
    return {
      representation: null,
      reason: "axis permutation not recognized",
      alignment: minAlignment,
      axisMap,
    };
  }
  const lowAlignment = !Number.isFinite(minAlignment) || minAlignment < INPUT_REP_ALIGNMENT_MIN;
  return {
    representation,
    reason: lowAlignment ? "low axis alignment" : null,
    alignment: minAlignment,
    axisMap,
  };
}

function detectGeometryRepresentationFromAtoms(atoms) {
  const normalizedAtoms = normalizeAtomsForQuartic(atoms);
  if (!normalizedAtoms.length) return null;
  try {
    const center = computeCenterOfMass(normalizedAtoms);
    const inertiaTensor = computeInertiaTensor(normalizedAtoms, center);
    const principal = diagonalizeSymmetricTensor(inertiaTensor);
    const detection = detectInputRepresentationFromPrincipalBasis({
      basis: principal.basis,
      moments: principal.moments,
    });
    return detection && detection.representation ? detection.representation : null;
  } catch (_) {
    return null;
  }
}

function inferRepresentationFromAxisMap(axisMap) {
  if (!axisMap || !axisMap.x || !axisMap.y || !axisMap.z) return null;
  const reps = Object.keys(REP_MAP);
  for (let i = 0; i < reps.length; i += 1) {
    const rep = reps[i];
    const mapping = REP_MAP[rep];
    if (!mapping) continue;
    const expectedAxisMap = {
      x: mapping[1],
      y: mapping[2],
      z: mapping[0],
    };
    if (
      expectedAxisMap.x === axisMap.x &&
      expectedAxisMap.y === axisMap.y &&
      expectedAxisMap.z === axisMap.z
    ) {
      return rep;
    }
  }
  return null;
}

function computeMomentRelativeDegeneracy(moments) {
  if (!Array.isArray(moments) || moments.length < 3) return Number.NaN;
  const m0 = Number(moments[0]);
  const m1 = Number(moments[1]);
  const m2 = Number(moments[2]);
  if (!Number.isFinite(m0) || !Number.isFinite(m1) || !Number.isFinite(m2)) return Number.NaN;
  const scale = Math.max(Math.abs(m0), Math.abs(m1), Math.abs(m2), 1);
  const d01 = Math.abs(m1 - m0);
  const d12 = Math.abs(m2 - m1);
  return Math.min(d01, d12) / scale;
}

function applyRepresentationAndRefresh() {
  if (!state.principalAtoms || !state.analysisData) return;

  const inputRep = getStep1InputRepresentation();
  const outputRep = getStep1OutputRepresentation();
  state.outputAtoms = applyRepresentation(state.principalAtoms, outputRep);

  const xyzText = buildXYZText(state.outputAtoms, outputRep, inputRep || "undetermined");
  ui.xyzOutput.value = xyzText;
  ui.downloadBtn.disabled = false;
  ui.downloadGjfBtn.disabled = false;
  ui.downloadBasisSetBtn.disabled = false;

  updateStep1DetectedRepresentationBadge(state.analysisData.inputRepresentationDetection || null);
  ui.analysis.textContent = formatAnalysis(state.analysisData, outputRep);
}

function applyRepresentation(principalAtoms, rep) {
  const mapping = REP_MAP[rep] || REP_MAP.Ir;
  const zAxis = mapping[0];
  const xAxis = mapping[1];
  const yAxis = mapping[2];

  return principalAtoms.map((atom) => {
    const p = atom.coords;
    return {
      label: atom.label,
      z: atom.z,
      coords: [p[AXIS_INDEX[xAxis]], p[AXIS_INDEX[yAxis]], p[AXIS_INDEX[zAxis]]],
    };
  });
}

function computeCenterOfMass(atoms) {
  let totalMass = 0;
  const center = [0, 0, 0];

  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    const mass = atomicMass[atom.z];
    if (!Number.isFinite(mass)) {
      throw new Error(`Atomic mass not available for Z=${atom.z}.`);
    }
    totalMass += mass;
    center[0] += atom.coords[0] * mass;
    center[1] += atom.coords[1] * mass;
    center[2] += atom.coords[2] * mass;
  }

  if (!Number.isFinite(totalMass) || totalMass <= 0) {
    throw new Error("Unable to compute the center of mass.");
  }

  center[0] /= totalMass;
  center[1] /= totalMass;
  center[2] /= totalMass;
  return center;
}

function computeInertiaTensor(atoms, center) {
  let ixx = 0;
  let iyy = 0;
  let izz = 0;
  let ixy = 0;
  let ixz = 0;
  let iyz = 0;

  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    const mass = atomicMass[atom.z];
    const x = atom.coords[0] - center[0];
    const y = atom.coords[1] - center[1];
    const z = atom.coords[2] - center[2];

    ixx += mass * (y * y + z * z);
    iyy += mass * (x * x + z * z);
    izz += mass * (x * x + y * y);
    ixy -= mass * x * y;
    ixz -= mass * x * z;
    iyz -= mass * y * z;
  }

  return [
    [ixx, ixy, ixz],
    [ixy, iyy, iyz],
    [ixz, iyz, izz],
  ];
}

function diagonalizeSymmetricTensor(tensor) {
  const eig = math.eigs(math.matrix(tensor));
  const values = toArray(eig.values);
  const vectors = toArray(eig.vectors);

  if (!Array.isArray(values) || !Array.isArray(vectors) || vectors.length !== 3) {
    throw new Error("Inertia tensor diagonalization failed.");
  }

  const pairs = [];
  for (let i = 0; i < 3; i += 1) {
    const vec = [
      asReal(vectors[0][i]),
      asReal(vectors[1][i]),
      asReal(vectors[2][i]),
    ];
    pairs.push({
      value: asReal(values[i]),
      vector: normalizeVec(vec),
    });
  }

  pairs.sort((a, b) => a.value - b.value);

  let eA = normalizeVec(pairs[0].vector);
  let eB = orthonormalize(pairs[1].vector, eA);
  let eC = normalizeVec(crossVec(eA, eB));
  if (dotVec(eC, pairs[2].vector) < 0) {
    eC = scaleVec(eC, -1);
  }

  eA = stabilizeSign(eA);
  eB = stabilizeSign(eB);
  eC = stabilizeSign(eC);

  if (determinantFromColumns(eA, eB, eC) < 0) {
    eC = scaleVec(eC, -1);
  }

  return {
    moments: [pairs[0].value, pairs[1].value, pairs[2].value],
    basis: [eA, eB, eC],
  };
}

function toArray(value) {
  if (Array.isArray(value)) return value;
  if (value && typeof value.toArray === "function") return value.toArray();
  return [];
}

function asReal(value) {
  if (typeof value === "number") return value;
  if (value && typeof value.re === "number") return value.re;
  const coerced = Number(value);
  return Number.isFinite(coerced) ? coerced : 0;
}

function rotateToBasis(vector, basis) {
  return [
    dotVec(vector, basis[0]),
    dotVec(vector, basis[1]),
    dotVec(vector, basis[2]),
  ];
}

function dotVec(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function crossVec(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function normVec(v) {
  return Math.sqrt(dotVec(v, v));
}

function normalizeVec(v) {
  const n = normVec(v);
  if (n < 1e-14) return [1, 0, 0];
  return [v[0] / n, v[1] / n, v[2] / n];
}

function scaleVec(v, factor) {
  return [v[0] * factor, v[1] * factor, v[2] * factor];
}

function subtractVec(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function orthonormalize(vector, reference) {
  const projection = scaleVec(reference, dotVec(vector, reference));
  let result = subtractVec(vector, projection);
  if (normVec(result) < 1e-12) {
    const trial = Math.abs(reference[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    result = crossVec(reference, trial);
  }
  return normalizeVec(result);
}

function stabilizeSign(vector) {
  const absVals = vector.map((x) => Math.abs(x));
  let idx = 0;
  if (absVals[1] > absVals[idx]) idx = 1;
  if (absVals[2] > absVals[idx]) idx = 2;
  if (vector[idx] < 0) return scaleVec(vector, -1);
  return vector;
}

function determinantFromColumns(c1, c2, c3) {
  return (
    c1[0] * (c2[1] * c3[2] - c2[2] * c3[1]) -
    c2[0] * (c1[1] * c3[2] - c1[2] * c3[1]) +
    c3[0] * (c1[1] * c2[2] - c1[2] * c2[1])
  );
}

function buildXYZText(atoms, outputRep, inputRep = "undetermined") {
  const commentParts = [
    `WMS-Prep input_representation=${inputRep} output_representation=${outputRep}`,
  ];
  if (state.sourceComment.trim()) {
    commentParts.push(`source="${state.sourceComment.trim()}"`);
  }

  const out = [];
  out.push(String(atoms.length));
  out.push(commentParts.join(" | "));

  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    out.push(
      `${atom.label} ${formatNumber(atom.coords[0])} ${formatNumber(atom.coords[1])} ${formatNumber(atom.coords[2])}`
    );
  }

  return `${out.join("\n")}\n`;
}

function formatNumber(value) {
  const v = Math.abs(value) < 5e-13 ? 0 : value;
  return v.toFixed(8);
}

function formatAnalysis(data, outputRep) {
  const zxy = REP_MAP[outputRep] || REP_MAP.Ir;
  const zAxis = zxy[0];
  const xAxis = zxy[1];
  const yAxis = zxy[2];
  const detection = data && data.inputRepresentationDetection
    ? data.inputRepresentationDetection
    : null;
  const inputRepLine = detection && detection.representation
    ? detection.representation
    : "Undetermined";
  const alignmentLine = detection && Number.isFinite(detection.alignment)
    ? `${(detection.alignment * 100).toFixed(1)}%`
    : "n/a";
  const detectionReason = detection && detection.reason ? detection.reason : "n/a";
  const axisMap = detection && detection.axisMap ? detection.axisMap : null;

  const lines = [];
  lines.push(`Detected input representation: ${inputRepLine}`);
  lines.push(`Input detection alignment: ${alignmentLine}`);
  lines.push(`Input detection note: ${detectionReason}`);
  lines.push(`Output representation: ${outputRep}`);
  lines.push("");
  lines.push("Center of mass (Angstrom):");
  lines.push(`  X = ${formatNumber(data.center[0])}`);
  lines.push(`  Y = ${formatNumber(data.center[1])}`);
  lines.push(`  Z = ${formatNumber(data.center[2])}`);
  lines.push("");
  lines.push("Inertia tensor (amu*Angstrom^2) with respect to the center of mass:");
  lines.push(matrixRow(data.inertiaTensor[0]));
  lines.push(matrixRow(data.inertiaTensor[1]));
  lines.push(matrixRow(data.inertiaTensor[2]));
  lines.push("");
  lines.push("Principal moments (Ia, Ib, Ic):");
  lines.push(`  Ia = ${formatNumber(data.moments[0])}`);
  lines.push(`  Ib = ${formatNumber(data.moments[1])}`);
  lines.push(`  Ic = ${formatNumber(data.moments[2])}`);
  lines.push("");
  lines.push("Principal axes in original coordinates (unit vectors):");
  lines.push(`  a = [${vecRow(data.basis[0])}]`);
  lines.push(`  b = [${vecRow(data.basis[1])}]`);
  lines.push(`  c = [${vecRow(data.basis[2])}]`);
  if (axisMap && axisMap.x && axisMap.y && axisMap.z) {
    lines.push("");
    lines.push("Detected principal-axis map in input XYZ:");
    lines.push(`  X ~= ${axisMap.x}`);
    lines.push(`  Y ~= ${axisMap.y}`);
    lines.push(`  Z ~= ${axisMap.z}`);
  }
  lines.push("");
  lines.push("Axis permutation by representation (from principal frame a,b,c):");
  lines.push(`  X <- ${xAxis}`);
  lines.push(`  Y <- ${yAxis}`);
  lines.push(`  Z <- ${zAxis}`);

  return lines.join("\n");
}

function matrixRow(row) {
  return `  [ ${formatNumber(row[0])}  ${formatNumber(row[1])}  ${formatNumber(row[2])} ]`;
}

function vecRow(row) {
  return `${formatNumber(row[0])}, ${formatNumber(row[1])}, ${formatNumber(row[2])}`;
}

function downloadOutputXYZ() {
  if (!state.outputAtoms || !state.outputAtoms.length) return;
  const rep = getStep1OutputRepresentation();
  const fileName = buildOutputFileName(rep, "xyz");
  downloadTextFile(ui.xyzOutput.value, fileName);
}

async function onDownloadGJFClicked() {
  try {
    await downloadOutputGJF();
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  }
}

async function onDownloadBasisSetClicked() {
  try {
    await downloadBasisSet();
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  }
}

async function downloadOutputGJF() {
  if (!state.outputAtoms || !state.outputAtoms.length) return;

  const template = await getGJFTemplate();
  const rep = getStep1OutputRepresentation();
  const fileName = buildOutputFileName(rep, "gjf");
  const gjfBaseName = fileName.replace(/\.gjf$/i, "");
  const gjfText = buildGJFTextFromTemplate(template, state.outputAtoms, gjfBaseName, rep);
  downloadTextFile(gjfText, fileName);
}

async function downloadBasisSet() {
  const a = document.createElement("a");
  a.href = "../G16_Inputs/3F12red.gbs.txt";
  a.download = "3F12red.gbs" || "";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}

function getGJFTemplate() {
  if (state.gjfTemplate) {
    return Promise.resolve(state.gjfTemplate);
  }

  return fetch(GJF_TEMPLATE_PATH, { cache: "no-store" })
    .then((response) => {
      if (!response.ok) {
        throw new Error(
          `Unable to load the GJF template (${response.status} ${response.statusText}).`
        );
      }
      return response.text();
    })
    .then((templateText) => {
      if (!templateText.includes(GJF_XYZ_PLACEHOLDER)) {
        throw new Error(
          `Template missing placeholder: ${GJF_XYZ_PLACEHOLDER}`
        );
      }
      if (!templateText.includes(GJF_NAME_PLACEHOLDER)) {
        throw new Error(
          `Template missing placeholder: ${GJF_NAME_PLACEHOLDER}`
        );
      }
      state.gjfTemplate = templateText;
      return templateText;
    });
}

function buildGJFTextFromTemplate(template, atoms, gjfBaseName, representation = "Ir") {
  const xyzBlock = buildGJFXYZBlock(atoms);
  let out = template
    .replace(GJF_XYZ_PLACEHOLDER, xyzBlock)
    .split(GJF_NAME_PLACEHOLDER)
    .join(gjfBaseName);
  if (out.includes(GJF_REPRESENTATION_PLACEHOLDER)) {
    out = out.split(GJF_REPRESENTATION_PLACEHOLDER).join(representation);
  }
  return out;
}

function buildGJFXYZBlock(atoms) {
  const lines = [];
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    lines.push(
      `${atom.label} ${formatNumber(atom.coords[0])} ${formatNumber(atom.coords[1])} ${formatNumber(atom.coords[2])}`
    );
  }
  return lines.join("\n");
}

const EMPTY_QUARTIC_ASYM = Object.freeze({
  deltaJ: "",
  deltaJK: "",
  deltaK: "",
  smallDeltaJ: "",
  smallDeltaK: "",
});
const EMPTY_QUARTIC_SYM = Object.freeze({
  DJ: "",
  DJK: "",
  DK: "",
  d1: "",
  d2: "",
});
const EMPTY_SEXTIC_ASYM = Object.freeze({
  PHIN: "",
  PHINK: "",
  PHIKN: "",
  PHIK: "",
  phiN: "",
  phiNK: "",
  phiK: "",
});
const EMPTY_SEXTIC_SYM = Object.freeze({
  HN: "",
  HNK: "",
  HKN: "",
  HK: "",
  h1: "",
  h2: "",
  h3: "",
});
const EMPTY_QUADRUPOLE_WMS_FIELDS = Object.freeze({
  I_NUC: "0",
  I_NUC_1: "0",
  I_NUC_2: "0",
  chi_aa: "0",
  chi_bb: "0",
  chi_cc: "0",
  chi_ab: "0",
  chi_ac: "0",
  chi_bc: "0",
  chi_aa_1: "0",
  chi_bb_1: "0",
  chi_cc_1: "0",
  chi_ab_1: "0",
  chi_ac_1: "0",
  chi_bc_1: "0",
  chi_aa_2: "0",
  chi_bb_2: "0",
  chi_cc_2: "0",
  chi_ab_2: "0",
  chi_ac_2: "0",
  chi_bc_2: "0",
});

function parseGaussianLogComponents(logText) {
  const pointGroup = extractPointGroup(logText) || "CS";
  const contextResolver = createContextResolver(logText);
  const geometryVariants = extractGeometryVariantsFromLog(logText, pointGroup, contextResolver);
  const quarticVariants = extractQuarticVariantsFromLog(logText, contextResolver);
  const sexticVariants = extractSexticVariantsFromLog(logText, contextResolver);
  const dvibVariants = extractDVibVariantsFromLog(logText, contextResolver);
  const quadrupoleVariants = extractQuadrupoleVariantsFromLog(logText, contextResolver);

  const geometry = pickDefaultVariantForType("geometry", geometryVariants);
  const quartic = pickDefaultVariantForType("quartic", quarticVariants);
  const sextic = pickDefaultVariantForType("sextic", sexticVariants);
  const dvib = pickDefaultVariantForType("dvib", dvibVariants);
  const quadrupole = pickDefaultVariantForType("quadrupole", quadrupoleVariants);

  return {
    pointGroup,
    geometryVariants,
    quarticVariants,
    sexticVariants,
    dvibVariants,
    quadrupoleVariants,
    geometry,
    quartic,
    sextic,
    dvib,
    quadrupole,
    summary: {
      hasGeometry: Boolean(geometry),
      hasQuartic: Boolean(quartic),
      hasSextic: Boolean(sextic),
      hasDVib: Boolean(dvib),
      hasQuadrupole: Boolean(quadrupole),
      hasCompletePipeline: Boolean(geometry && quartic && sextic && dvib),
    },
  };
}

function parseGaussianLogForWMSInput(logText) {
  const parsed = parseGaussianLogComponents(logText);
  const combined = combineParsedComponents({
    geometry: parsed.geometry,
    quartic: parsed.quartic,
    sextic: parsed.sextic,
    dvib: parsed.dvib,
    quadrupole: parsed.quadrupole,
  });
  if (!combined) {
    throw new Error("Unable to parse enough Gaussian data to build a WMS-Rot input.");
  }
  return combined;
}

function extractGeometryVariantsFromLog(logText, fallbackPointGroup, contextResolver) {
  const rotationalEvents = extractAllRotationalConstantsOccurrences(logText);
  if (!rotationalEvents.length) {
    const single = extractGeometryComponentFromLog(logText, fallbackPointGroup);
    return single ? [single] : [];
  }

  const orientationEvents = extractAllInputOrientationOccurrences(logText, contextResolver);
  const pointGroupEvents = extractAllPointGroupOccurrences(logText, contextResolver);
  const dipoleEvents = extractAllDipoleOccurrences(logText, contextResolver);
  const grouped = new Map();

  for (let i = 0; i < rotationalEvents.length; i += 1) {
    const event = rotationalEvents[i];
    const context = extractContextInfo(logText, event.index, contextResolver);
    const key = `geometry|${context.groupKey || "global"}`;
    const inputRepresentation = resolveLogSectionInputRepresentation(
      logText,
      event.index,
      logText.slice(event.index, Math.min(logText.length, event.index + 2400))
    );
    const atomsEvent =
      findLastBeforeWithGroup(
        orientationEvents,
        event.index,
        context.groupKey,
        context.segmentId
      ) ||
      findLastBefore(orientationEvents, event.index);
    const pointGroupEvent =
      findNearestByGroup(pointGroupEvents, event.index, context.groupKey, context.segmentId) ||
      findLastBefore(pointGroupEvents, event.index);
    const dipoleEvent =
      findNearestByGroup(dipoleEvents, event.index, context.groupKey, context.segmentId) ||
      findNearestOccurrence(dipoleEvents, event.index);
    const atoms = atomsEvent && Array.isArray(atomsEvent.atoms) ? atomsEvent.atoms : [];

    grouped.set(key, {
      contextLabel: context.contextLabel,
      sourceIndex: event.index,
      inputRepresentation,
      pointGroup:
        (pointGroupEvent && pointGroupEvent.pointGroup) ||
        fallbackPointGroup ||
        "CS",
      rotationalConstantsMHz: event.rotationalConstantsMHz,
      dipole: dipoleEvent && dipoleEvent.dipole
        ? dipoleEvent.dipole
        : { a: "", b: "", c: "" },
      dpcs3Atoms: atoms,
      dpcs3XYZText: atoms.length ? buildAtomicNumberXYZText(atoms) : "",
    });
  }

  const variants = Array.from(grouped.values());
  if (variants.length) return variants;
  const single = extractGeometryComponentFromLog(logText, fallbackPointGroup);
  return single ? [single] : [];
}

function extractQuarticVariantsFromLog(logText, contextResolver) {
  const tauPrimeVariants = extractQuarticVariantsFromTauPrimeSections(logText, contextResolver);
  const reducedVariants = extractQuarticVariantsFromReducedSections(logText, contextResolver);
  return tauPrimeVariants.length ? tauPrimeVariants : reducedVariants;
}

function extractQuarticVariantsFromTauPrimeSections(logText, contextResolver) {
  const rotationalEvents = extractAllRotationalConstantContextEvents(logText, contextResolver);
  const tauSections = extractAllTauPrimeSections(logText);
  const variants = [];
  const contextCount = {};

  for (let i = 0; i < tauSections.length; i += 1) {
    const section = tauSections[i];
    const tauPrime = parseQuarticTauPrimeFromSection(section.text);
    if (!tauPrime) continue;
    const sigma = parseSigmaFromSection(section.text);
    const context = extractContextInfo(logText, section.index, contextResolver);
    const inputRepresentation = resolveLogSectionInputRepresentation(
      logText,
      section.index,
      section.text
    );
    const rotatedTauPrime = rotateTauPrimeComponentsForRepresentation(
      tauPrime,
      "Ir",
      inputRepresentation
    ).tauPrime;
    const preview = buildQuarticVariantFromTauPrimeSection({
      tauPrime: rotatedTauPrime,
      sigma,
      nielsen: parseNielsenQuarticFromSection(section.text),
    });
    if (!preview) continue;
    variants.push({
      contextLabel: buildContextLabelWithOccurrence(context, contextCount, "TauP"),
      sourceIndex: section.index,
      inputRepresentation,
      rotationalConstantsMHz: findRotationalConstantsNearContext(
        rotationalEvents,
        section.index,
        context.groupKey,
        context.segmentId
      ),
      sigma: Number.isFinite(sigma) ? sigma : null,
      sigma1: Number.isFinite(preview.sigma1) ? preview.sigma1 : null,
      tauPrime: { ...tauPrime },
      T: preview.T ? { ...preview.T } : null,
      nielsen: preview.nielsen ? { ...preview.nielsen } : null,
      asym: { ...(preview.asym || EMPTY_QUARTIC_ASYM) },
      sym: { ...(preview.sym || EMPTY_QUARTIC_SYM) },
    });
  }

  return variants;
}

function extractQuarticVariantsFromReducedSections(logText, contextResolver) {
  const rotationalEvents = extractAllRotationalConstantContextEvents(logText, contextResolver);
  const asymSections = extractAllSections(logText, ASYM_REDUCTION_HEADER, 2600);
  const symSections = extractAllSections(logText, SYM_REDUCTION_HEADER, 2600);
  const grouped = new Map();

  for (let i = 0; i < asymSections.length; i += 1) {
    const section = asymSections[i];
    const context = extractContextInfo(logText, section.index, contextResolver);
    const key = `quartic|${context.groupKey || "global"}`;
    const current = grouped.get(key) || {
      contextLabel: context.contextLabel,
      sourceIndex: section.index,
      context,
      rotationalConstantsMHz: null,
      inputRepresentation: resolveLogSectionInputRepresentation(
        logText,
        section.index,
        section.text
      ),
      asym: { ...EMPTY_QUARTIC_ASYM },
      sym: { ...EMPTY_QUARTIC_SYM },
    };
    current.sourceIndex = Math.max(current.sourceIndex, section.index);
    current.contextLabel = context.contextLabel || current.contextLabel;
    current.context = context;
    current.inputRepresentation =
      current.inputRepresentation ||
      resolveLogSectionInputRepresentation(logText, section.index, section.text);
    current.rotationalConstantsMHz =
      current.rotationalConstantsMHz ||
      findRotationalConstantsNearContext(
        rotationalEvents,
        section.index,
        context.groupKey,
        context.segmentId
      );
    current.asym = {
      deltaJ: extractMHzFromSectionOptional(section.text, "DELTA\\s+J"),
      deltaJK: extractMHzFromSectionOptional(section.text, "DELTA\\s+JK"),
      deltaK: extractMHzFromSectionOptional(section.text, "DELTA\\s+K"),
      smallDeltaJ: extractMHzFromSectionOptional(section.text, "delta\\s+J"),
      smallDeltaK: extractMHzFromSectionOptional(section.text, "delta\\s+K"),
    };
    grouped.set(key, current);
  }

  for (let i = 0; i < symSections.length; i += 1) {
    const section = symSections[i];
    const context = extractContextInfo(logText, section.index, contextResolver);
    const key = `quartic|${context.groupKey || "global"}`;
    const current = grouped.get(key) || {
      contextLabel: context.contextLabel,
      sourceIndex: section.index,
      context,
      rotationalConstantsMHz: null,
      inputRepresentation: resolveLogSectionInputRepresentation(
        logText,
        section.index,
        section.text
      ),
      asym: { ...EMPTY_QUARTIC_ASYM },
      sym: { ...EMPTY_QUARTIC_SYM },
    };
    current.sourceIndex = Math.max(current.sourceIndex, section.index);
    current.contextLabel = context.contextLabel || current.contextLabel;
    current.context = context;
    current.inputRepresentation =
      current.inputRepresentation ||
      resolveLogSectionInputRepresentation(logText, section.index, section.text);
    current.rotationalConstantsMHz =
      current.rotationalConstantsMHz ||
      findRotationalConstantsNearContext(
        rotationalEvents,
        section.index,
        context.groupKey,
        context.segmentId
      );
    current.sym = {
      DJ: extractMHzFromSectionOptional(section.text, "D\\s+J"),
      DJK: extractMHzFromSectionOptional(section.text, "D\\s+JK"),
      DK: extractMHzFromSectionOptional(section.text, "D\\s+K"),
      d1: extractMHzFromSectionOptional(section.text, "d\\s+1"),
      d2: extractMHzFromSectionOptional(section.text, "d\\s+2"),
    };
    grouped.set(key, current);
  }

  const variants = Array.from(grouped.values()).filter(
    (item) => hasAnyNonEmptyValue(item.asym) || hasAnyNonEmptyValue(item.sym)
  );
  const contextCount = {};
  for (let i = 0; i < variants.length; i += 1) {
    variants[i].contextLabel = buildContextLabelWithOccurrence(
      variants[i].context,
      contextCount
    );
    delete variants[i].context;
  }
  return variants;
}

function extractSexticVariantsFromLog(logText, contextResolver) {
  const reducedVariants = extractSexticVariantsFromReducedSections(logText, contextResolver);
  const phiVariants = extractSexticVariantsFromPhiSections(logText, contextResolver);
  return phiVariants.length ? phiVariants : reducedVariants;
}

function extractSexticVariantsFromPhiSections(logText, contextResolver) {
  const rotationalEvents = extractAllRotationalConstantContextEvents(logText, contextResolver);
  const quarticSections = extractAllTauPrimeSections(logText);
  const sexticSections = extractAllSexticDistortionSections(logText);
  const variants = [];
  const contextCount = {};

  for (let i = 0; i < sexticSections.length; i += 1) {
    const section = sexticSections[i];
    const inertial = parseSexticInertialFromPhiSection(section.text);
    if (!inertial) continue;
    const sigma = parseSigmaFromSection(section.text);
    const context = extractContextInfo(logText, section.index, contextResolver);
    const inputRepresentation = resolveLogSectionInputRepresentation(
      logText,
      section.index,
      section.text
    );
    const rotatedInertial = rotateSexticInertialComponents(
      inertial,
      "Ir",
      inputRepresentation
    ).inertial;
    const rotationalConstantsMHz = findRotationalConstantsNearContext(
      rotationalEvents,
      section.index,
      context.groupKey,
      context.segmentId
    );
    const quarticT = findQuarticTensorForSexticSection(
      quarticSections,
      logText,
      section.index,
      contextResolver,
      context.groupKey,
      context.segmentId,
      inputRepresentation
    );
    const preview = buildSexticVariantFromInertialSection({
      inertial: rotatedInertial,
      sigma,
      quarticT,
      rotationalConstantsMHz: buildWatsonABCFromPrincipalConstants(
        rotationalConstantsMHz,
        inverseFullWatsonRepresentationToken(inputRepresentation)
      ),
    });
    if (!preview) continue;
    variants.push({
      contextLabel: buildContextLabelWithOccurrence(context, contextCount, "Phi"),
      sourceIndex: section.index,
      inputRepresentation,
      rotationalConstantsMHz,
      sigma: Number.isFinite(sigma) ? sigma : null,
      inertial: { ...inertial },
      asymSextic: { ...(preview.asymSextic || EMPTY_SEXTIC_ASYM) },
      symSextic: { ...(preview.symSextic || EMPTY_SEXTIC_SYM) },
    });
  }

  return variants;
}

function buildWatsonABCFromPrincipalConstants(values, representation) {
  if (!values || typeof values !== "object") return null;
  const principal = {
    a: getRotationalConstantForSextic(values, "A"),
    b: getRotationalConstantForSextic(values, "B"),
    c: getRotationalConstantForSextic(values, "C"),
  };
  if (!Number.isFinite(principal.a) || !Number.isFinite(principal.b) || !Number.isFinite(principal.c)) {
    return null;
  }
  const mapping = REP_MAP[normalizeFullWatsonRepresentationToken(representation)] || REP_MAP.Ir;
  return {
    A: principal[mapping[0]],
    B: principal[mapping[1]],
    C: principal[mapping[2]],
  };
}

function findQuarticTensorForSexticSection(
  quarticSections,
  logText,
  beforeIndex,
  contextResolver,
  groupKey,
  segmentId,
  targetRepresentation
) {
  if (!Array.isArray(quarticSections) || !quarticSections.length) return null;
  let best = null;
  for (let i = 0; i < quarticSections.length; i += 1) {
    const section = quarticSections[i];
    if (!section || !Number.isFinite(section.index) || section.index > beforeIndex) continue;
    const context = extractContextInfo(logText, section.index, contextResolver);
    if (groupKey && context.groupKey && context.groupKey !== groupKey) continue;
    if (segmentId && context.segmentId && context.segmentId !== segmentId) continue;
    const tauPrime = parseQuarticTauPrimeFromSection(section.text);
    const inputRepresentation = resolveLogSectionInputRepresentation(
      logText,
      section.index,
      section.text
    );
    const rotatedTauPrime = tauPrime
      ? rotateTauPrimeComponentsForRepresentation(
          tauPrime,
          inputRepresentation,
          targetRepresentation
        ).tauPrime
      : null;
    const T = rotatedTauPrime ? deriveQuarticTensorFromTauPrime(rotatedTauPrime) : null;
    if (T) best = T;
  }
  return best;
}

function extractQuadrupoleVariantsFromLog(logText, contextResolver) {
  const sections = extractSectionsByHeaderUntilMarkers(
    logText,
    "Nuclear quadrupole coupling constants [Chi] (MHz):",
    [
      "Nuclear quadrupole coupling constants [Chi] (MHz):",
      "Dipole moment (Debye):",
      "Electric dipole moment (input orientation):",
      "Quartic Centrifugal Distortion Constants",
      "Normal termination of Gaussian",
      "Link1:  Proceeding to internal job step number",
      "Initial command:",
    ],
    3200
  );
  const variants = [];
  const contextCount = {};

  for (let i = 0; i < sections.length; i += 1) {
    const section = sections[i];
    const nuclei = parseQuadrupoleNucleiFromSection(section.text);
    if (!nuclei.length) continue;
    const context = extractContextInfo(logText, section.index, contextResolver);
    const inputRepresentation = resolveLogSectionInputRepresentation(
      logText,
      section.index,
      section.text
    );
    variants.push({
      contextLabel: buildContextLabelWithOccurrence(context, contextCount, "Quadrupole"),
      sourceIndex: section.index,
      inputRepresentation,
      nuclei,
    });
  }

  return variants;
}

function extractSexticVariantsFromReducedSections(logText, contextResolver) {
  const rotationalEvents = extractAllRotationalConstantContextEvents(logText, contextResolver);
  const asymSections = extractAllSections(logText, A_SEXTIC_REDUCTION_HEADER, 2600);
  const symSections = extractAllSections(logText, S_SEXTIC_REDUCTION_HEADER, 2600);
  const grouped = new Map();

  for (let i = 0; i < asymSections.length; i += 1) {
    const section = asymSections[i];
    const context = extractContextInfo(logText, section.index, contextResolver);
    const key = `sextic|${context.groupKey || "global"}`;
    const current = grouped.get(key) || {
      contextLabel: context.contextLabel,
      sourceIndex: section.index,
      context,
      rotationalConstantsMHz: null,
      inputRepresentation: resolveLogSectionInputRepresentation(
        logText,
        section.index,
        section.text
      ),
      asymSextic: { ...EMPTY_SEXTIC_ASYM },
      symSextic: { ...EMPTY_SEXTIC_SYM },
    };
    current.sourceIndex = Math.max(current.sourceIndex, section.index);
    current.contextLabel = context.contextLabel || current.contextLabel;
    current.context = context;
    current.inputRepresentation =
      current.inputRepresentation ||
      resolveLogSectionInputRepresentation(logText, section.index, section.text);
    current.rotationalConstantsMHz =
      current.rotationalConstantsMHz ||
      findRotationalConstantsNearContext(
        rotationalEvents,
        section.index,
        context.groupKey,
        context.segmentId
      );
    current.asymSextic = {
      PHIN: extractHzAsMHzFromSection(section.text, "Phi\\s+J"),
      PHINK: extractHzAsMHzFromSection(section.text, "Phi\\s+JK"),
      PHIKN: extractHzAsMHzFromSection(section.text, "Phi\\s+KJ"),
      PHIK: extractHzAsMHzFromSection(section.text, "Phi\\s+K"),
      phiN: extractHzAsMHzFromSection(section.text, "phi\\s+j"),
      phiNK: extractHzAsMHzFromSection(section.text, "phi\\s+jk"),
      phiK: extractHzAsMHzFromSection(section.text, "phi\\s+k"),
    };
    grouped.set(key, current);
  }

  for (let i = 0; i < symSections.length; i += 1) {
    const section = symSections[i];
    const context = extractContextInfo(logText, section.index, contextResolver);
    const key = `sextic|${context.groupKey || "global"}`;
    const current = grouped.get(key) || {
      contextLabel: context.contextLabel,
      sourceIndex: section.index,
      context,
      rotationalConstantsMHz: null,
      inputRepresentation: resolveLogSectionInputRepresentation(
        logText,
        section.index,
        section.text
      ),
      asymSextic: { ...EMPTY_SEXTIC_ASYM },
      symSextic: { ...EMPTY_SEXTIC_SYM },
    };
    current.sourceIndex = Math.max(current.sourceIndex, section.index);
    current.contextLabel = context.contextLabel || current.contextLabel;
    current.context = context;
    current.inputRepresentation =
      current.inputRepresentation ||
      resolveLogSectionInputRepresentation(logText, section.index, section.text);
    current.rotationalConstantsMHz =
      current.rotationalConstantsMHz ||
      findRotationalConstantsNearContext(
        rotationalEvents,
        section.index,
        context.groupKey,
        context.segmentId
      );
    current.symSextic = {
      HN: extractHzAsMHzFromSection(section.text, "H\\s+J"),
      HNK: extractHzAsMHzFromSection(section.text, "H\\s+JK"),
      HKN: extractHzAsMHzFromSection(section.text, "H\\s+KJ"),
      HK: extractHzAsMHzFromSection(section.text, "H\\s+K"),
      h1: extractHzAsMHzFromSection(section.text, "h\\s+1"),
      h2: extractHzAsMHzFromSection(section.text, "h\\s+2"),
      h3: extractHzAsMHzFromSection(section.text, "h\\s+3"),
    };
    grouped.set(key, current);
  }

  const variants = Array.from(grouped.values()).filter(
    (item) =>
      hasAnyNonEmptyValue(item.asymSextic) || hasAnyNonEmptyValue(item.symSextic)
  );
  const contextCount = {};
  for (let i = 0; i < variants.length; i += 1) {
    variants[i].contextLabel = buildContextLabelWithOccurrence(
      variants[i].context,
      contextCount
    );
    delete variants[i].context;
  }
  return variants;
}

function buildContextLabelWithOccurrence(context, counter, suffix = "") {
  const groupKey = context && context.groupKey ? context.groupKey : "global";
  const nextCount = (counter[groupKey] || 0) + 1;
  counter[groupKey] = nextCount;
  const base = context && context.contextLabel ? context.contextLabel : "Gaussian block";
  const occurrence = nextCount > 1 ? ` [set ${nextCount}]` : "";
  if (!suffix) {
    return `${base}${occurrence}`;
  }
  return `${base}${occurrence} | ${suffix}`;
}

function extractAllRotationalConstantContextEvents(logText, resolver = null) {
  const events = extractAllRotationalConstantsOccurrences(logText);
  const out = [];
  for (let i = 0; i < events.length; i += 1) {
    const event = events[i];
    const context = extractContextInfo(logText, event.index, resolver);
    out.push({
      ...event,
      groupKey: context.groupKey,
      segmentId: context.segmentId,
    });
  }
  return out;
}

function findRotationalConstantsNearContext(events, index, groupKey, segmentId = null) {
  const event =
    findLastBeforeWithGroup(events, index, groupKey, segmentId) ||
    findNearestByGroup(events, index, groupKey, segmentId) ||
    findLastBefore(events, index) ||
    findNearestOccurrence(events, index);
  if (!event || !event.rotationalConstantsMHz) return null;
  return {
    A: Number(event.rotationalConstantsMHz.A),
    B: Number(event.rotationalConstantsMHz.B),
    C: Number(event.rotationalConstantsMHz.C),
  };
}

function extractAllTauPrimeSections(logText) {
  return extractSectionsByHeaderUntilMarkers(
    logText,
    "Quartic Centrifugal Distortion Constants Tau Prime",
    [
      "Quartic Centrifugal Distortion Constants Tau Prime",
      "Sextic Centrifugal Distortion Constants",
      "Average Coordinates and Mean Square Amplitudes",
      "Normal termination of Gaussian",
      "Link1:  Proceeding to internal job step number",
      "Initial command:",
    ],
    7000
  );
}

function extractAllSexticDistortionSections(logText) {
  return extractSectionsByHeaderUntilMarkers(
    logText,
    "Sextic Distortion Constants",
    [
      "Sextic Distortion Constants",
      "Average Coordinates and Mean Square Amplitudes",
      "Normal termination of Gaussian",
      "Link1:  Proceeding to internal job step number",
      "Initial command:",
    ],
    8000
  );
}

function extractSectionsByHeaderUntilMarkers(logText, header, markers, fallbackLength = 5000) {
  const sections = [];
  let from = 0;
  while (from < logText.length) {
    const idx = logText.indexOf(header, from);
    if (idx < 0) break;
    const start = idx;
    let end = Math.min(logText.length, idx + fallbackLength);
    for (let i = 0; i < markers.length; i += 1) {
      const marker = markers[i];
      if (!marker) continue;
      const markerIdx = logText.indexOf(marker, idx + header.length);
      if (markerIdx >= 0 && markerIdx < end) {
        end = markerIdx;
      }
    }
    sections.push({
      index: start,
      text: logText.slice(start, end),
    });
    from = idx + header.length;
  }
  return sections;
}

function parseQuarticTauPrimeFromSection(sectionText) {
  if (!sectionText) return null;
  const out = {
    taaaa: null,
    tbbaa: null,
    tbbbb: null,
    tccaa: null,
    tccbb: null,
    tcccc: null,
  };
  const regex = new RegExp(
    `^\\s*TauP\\s+([abc]{4})\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "gmi"
  );
  let match;
  while ((match = regex.exec(sectionText)) !== null) {
    const token = String(match[1] || "").toLowerCase();
    const valueMHz = parseSciNumber(match[3]);
    if (!Number.isFinite(valueMHz)) continue;
    const key = canonicalTauPrimeKeyFromToken(token);
    if (!key || !Object.prototype.hasOwnProperty.call(out, key)) continue;
    out[key] = valueMHz;
  }
  const keys = Object.keys(out);
  for (let i = 0; i < keys.length; i += 1) {
    if (!Number.isFinite(out[keys[i]])) return null;
  }
  return out;
}

function parseSigmaFromSection(sectionText) {
  if (!sectionText) return null;
  const regex = new RegExp(
    `^\\s*Sigma\\s*:\\s*(${NUMBER_TOKEN_PATTERN})(?:\\s*\\|.*)?$`,
    "mi"
  );
  const match = sectionText.match(regex);
  if (!match) return null;
  const sigma = parseSciNumber(match[1]);
  return Number.isFinite(sigma) ? sigma : null;
}

function parseNielsenQuarticFromSection(sectionText) {
  if (!sectionText) return null;
  const block = extractSubsectionByHeader(sectionText, "Nielsen Centrifugal Distortion Constants", 1600);
  if (!block) return null;
  const DJ = extractLabeledTableValueSecondColumn(block, "DJ");
  const DJK = extractLabeledTableValueSecondColumn(block, "DJK");
  const DK = extractLabeledTableValueSecondColumn(block, "DK");
  const dJ = extractLabeledTableValueSecondColumn(block, "dJ");
  const R5 = extractLabeledTableValueSecondColumn(block, "R5");
  const R6 = extractLabeledTableValueSecondColumn(block, "R6");
  if (
    !Number.isFinite(DJ) ||
    !Number.isFinite(DJK) ||
    !Number.isFinite(DK) ||
    !Number.isFinite(dJ) ||
    !Number.isFinite(R5) ||
    !Number.isFinite(R6)
  ) {
    return null;
  }
  return { DJ, DJK, DK, dJ, R5, R6 };
}

function extractSubsectionByHeader(text, header, fallbackLength = 1000) {
  if (!text || !header) return null;
  const idx = text.indexOf(header);
  if (idx < 0) return null;
  const end = Math.min(text.length, idx + Math.max(120, fallbackLength));
  return text.slice(idx, end);
}

function extractLabeledTableValueSecondColumn(sectionText, label) {
  if (!sectionText || !label) return null;
  const escapedLabel = escapeRegExp(label);
  const regex = new RegExp(
    `^\\s*${escapedLabel}\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "m"
  );
  const match = sectionText.match(regex);
  if (!match) return null;
  const parsed = parseSciNumber(match[2]);
  return Number.isFinite(parsed) ? parsed : null;
}

function deriveQuarticTensorFromTauPrime(tauPrime) {
  const taaaa = Number(tauPrime && tauPrime.taaaa);
  const taabb = Number(tauPrime && tauPrime.tbbaa);
  const tbbbb = Number(tauPrime && tauPrime.tbbbb);
  const taacc = Number(tauPrime && tauPrime.tccaa);
  const tbbcc = Number(tauPrime && tauPrime.tccbb);
  const tcccc = Number(tauPrime && tauPrime.tcccc);

  if (
    !Number.isFinite(taaaa) ||
    !Number.isFinite(taabb) ||
    !Number.isFinite(tbbbb) ||
    !Number.isFinite(taacc) ||
    !Number.isFinite(tbbcc) ||
    !Number.isFinite(tcccc)
  ) {
    return null;
  }

  const DJ = -(3.0 * tbbbb + 3.0 * tcccc + 2.0 * tbbcc) / 32.0;
  const DK = DJ - 0.25 * (taaaa - taabb - taacc);
  const DJK = -DJ - DK - 0.25 * taaaa;
  const dJ = -(tbbbb - tcccc) / 16.0;
  const R5 = -(tbbbb - tcccc - 2.0 * taabb + 2.0 * taacc) / 32.0;
  const R6 = (tbbbb + tcccc - 2.0 * tbbcc) / 64.0;

  const Taa = taaaa / 4.0;
  const Tbb = tbbbb / 4.0;
  const Tcc = tcccc / 4.0;
  const Tab = taabb / 8.0;
  const Tac = taacc / 8.0;
  const Tbc = tbbcc / 8.0;
  const T400 = -DJ;
  const T220 = -DJK;
  const T040 = -DK;
  const T202 = -dJ;
  const T022 = 2.0 * R5;
  const T004 = R6;

  return {
    Taa,
    Tbb,
    Tcc,
    Tab,
    Tac,
    Tbc,
    T400,
    T220,
    T040,
    T202,
    T022,
    T004,
    DJ,
    DJK,
    DK,
    dJ,
    R5,
    R6,
  };
}

function buildQuarticVariantFromTauPrimeSection({ tauPrime, sigma, nielsen }) {
  const T = deriveQuarticTensorFromTauPrime(tauPrime);
  if (!T) return null;

  const s = Number.isFinite(sigma) ? sigma : 0;
  const sigma1 = Math.abs(s) > 1.0e-16 ? 1.0 / s : 0.0;

  const asymMHz = {
    deltaJ: T.DJ - 2.0 * T.R6,
    deltaJK: T.DJK + 12.0 * T.R6,
    deltaK: T.DK - 10.0 * T.R6,
    smallDeltaJ: T.dJ,
    smallDeltaK: -2.0 * T.R5 - 4.0 * s * T.R6,
  };

  const symMHz = {
    DJ: T.DJ + T.R5 * sigma1,
    DJK: T.DJK - 6.0 * T.R5 * sigma1,
    DK: T.DK + 5.0 * T.R5 * sigma1,
    d1: -T.dJ,
    d2: T.R6 + 0.5 * T.R5 * sigma1,
  };

  return {
    sigma: s,
    sigma1,
    tauPrime: { ...tauPrime },
    T: { ...T },
    nielsen: nielsen ? { ...nielsen } : null,
    asym: {
      deltaJ: toFortranScientific(asymMHz.deltaJ, 10),
      deltaJK: toFortranScientific(asymMHz.deltaJK, 10),
      deltaK: toFortranScientific(asymMHz.deltaK, 10),
      smallDeltaJ: toFortranScientific(asymMHz.smallDeltaJ, 10),
      smallDeltaK: toFortranScientific(asymMHz.smallDeltaK, 10),
    },
    sym: {
      DJ: toFortranScientific(symMHz.DJ, 10),
      DJK: toFortranScientific(symMHz.DJK, 10),
      DK: toFortranScientific(symMHz.DK, 10),
      d1: toFortranScientific(symMHz.d1, 10),
      d2: toFortranScientific(symMHz.d2, 10),
    },
  };
}

function parseSexticInertialFromPhiSection(sectionText) {
  if (!sectionText) return null;
  const inertial = {};
  for (let i = 0; i < SEXTIC_INERTIAL_KEYS.length; i += 1) {
    inertial[SEXTIC_INERTIAL_KEYS[i]] = null;
  }
  const regex = new RegExp(
    `^\\s*Phi\\s+([abc]{3})\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "gmi"
  );
  let match;
  while ((match = regex.exec(sectionText)) !== null) {
    const key = canonicalSexticInertialKeyFromToken(match[1]);
    if (!key) continue;
    const hzValue = parseSciNumber(match[3]);
    if (!Number.isFinite(hzValue)) continue;
    inertial[key] = hzValue / 1e6;
  }
  const missing = SEXTIC_INERTIAL_KEYS.filter((key) => !Number.isFinite(inertial[key]));
  if (missing.length) return null;
  return inertial;
}

function parseQuadrupoleNucleiFromSection(sectionText) {
  const lines = String(sectionText || "").split(/\r?\n/);
  const nuclei = [];

  for (let i = 0; i < lines.length; i += 1) {
    const headerMatch = lines[i].match(/^\s*(\d+)\s+([A-Z][A-Za-z]?(?:\(\d+\)|-\d+)?)\s*$/);
    if (!headerMatch) continue;
    const centerIndex = Number.parseInt(headerMatch[1], 10);
    const rawLabel = String(headerMatch[2] || "").trim();
    const tensor = {};
    let consumed = 0;

    for (let j = i + 1; j < lines.length && consumed < 3; j += 1) {
      const line = lines[j];
      if (!line || !line.trim()) continue;
      consumed += 1;
      const pairRegex = new RegExp(`([abc]{2})\\s*=\\s*(${NUMBER_TOKEN_PATTERN})`, "gi");
      let pairMatch;
      while ((pairMatch = pairRegex.exec(line)) !== null) {
        const pair = String(pairMatch[1] || "").toLowerCase();
        const value = parseSciNumber(pairMatch[2]);
        if (!Number.isFinite(value)) continue;
        tensor[pair] = value;
      }
    }

    const chi = {
      chi_aa: tensor.aa,
      chi_bb: tensor.bb,
      chi_cc: tensor.cc,
      chi_ab: Number.isFinite(tensor.ab) ? tensor.ab : tensor.ba,
      chi_ac: Number.isFinite(tensor.ac) ? tensor.ac : tensor.ca,
      chi_bc: Number.isFinite(tensor.bc) ? tensor.bc : tensor.cb,
    };
    if (!QUAD_CART_KEYS.every((key) => Number.isFinite(chi[key]))) {
      continue;
    }

    const isotopeMatch = rawLabel.match(/^([A-Z][A-Za-z]?)(?:\((\d+)\)|-(\d+))?/);
    nuclei.push({
      centerIndex,
      label: rawLabel,
      element: isotopeMatch ? isotopeMatch[1] : "",
      massNumber: isotopeMatch
        ? Number.parseInt(isotopeMatch[2] || isotopeMatch[3] || "", 10) || null
        : null,
      chi,
    });
  }

  return nuclei;
}

function buildSexticVariantFromInertialSection({
  inertial,
  sigma,
  quarticT = null,
  rotationalConstantsMHz = null,
}) {
  if (!inertial) return null;
  const s = Number.isFinite(sigma) ? sigma : 0;
  const sextic = buildSexticFromInertialComponents({
    inertial,
    sigma: s,
    quarticT,
    rotationalConstantsMHz,
  });
  return {
    sigma: s,
    inertial: { ...inertial },
    asymSextic: { ...sextic.asymFormatted },
    symSextic: { ...sextic.symFormatted },
  };
}

function extractDVibVariantsFromLog(logText, contextResolver) {
  const regex =
    /Rotational Constants\s*\(in MHz\)[\s\S]*?Ae=\s*([-\d+.DEde]+)\s*A00=\s*([-\d+.DEde]+)[\s\S]*?Be=\s*([-\d+.DEde]+)\s*B00=\s*([-\d+.DEde]+)[\s\S]*?Ce=\s*([-\d+.DEde]+)\s*C00=\s*([-\d+.DEde]+)/gi;
  const grouped = new Map();
  let match;
  while ((match = regex.exec(logText)) !== null) {
    const Ae = parseSciNumber(match[1]);
    const A00 = parseSciNumber(match[2]);
    const Be = parseSciNumber(match[3]);
    const B00 = parseSciNumber(match[4]);
    const Ce = parseSciNumber(match[5]);
    const C00 = parseSciNumber(match[6]);
    if (
      !Number.isFinite(Ae) ||
      !Number.isFinite(A00) ||
      !Number.isFinite(Be) ||
      !Number.isFinite(B00) ||
      !Number.isFinite(Ce) ||
      !Number.isFinite(C00)
    ) {
      continue;
    }
    const context = extractContextInfo(logText, match.index, contextResolver);
    const key = `dvib|${context.groupKey || "global"}`;
    grouped.set(key, {
      contextLabel: context.contextLabel,
      sourceIndex: match.index,
      inputRepresentation: resolveLogSectionInputRepresentation(
        logText,
        match.index,
        logText.slice(match.index, Math.min(logText.length, match.index + 2200))
      ),
      A: A00 - Ae,
      B: B00 - Be,
      C: C00 - Ce,
    });
  }
  return Array.from(grouped.values());
}

function pickDefaultVariantForType(type, variants) {
  if (!Array.isArray(variants) || !variants.length) return null;
  let best = variants[0];
  let bestScore = computeVariantPriority(type, best);
  for (let i = 1; i < variants.length; i += 1) {
    const candidate = variants[i];
    const score = computeVariantPriority(type, candidate);
    if (score > bestScore) {
      best = candidate;
      bestScore = score;
      continue;
    }
    if (score === bestScore && (candidate.sourceIndex || 0) > (best.sourceIndex || 0)) {
      best = candidate;
      bestScore = score;
    }
  }
  return best;
}

function computeVariantPriority(type, variant) {
  const text = String(variant.contextLabel || "").toLowerCase();
  const sourceScore = Number.isFinite(variant.sourceIndex)
    ? variant.sourceIndex / 100000000
    : 0;
  let score = sourceScore;
  if (type === "geometry" || type === "quartic") {
    if (/(step\s*3)/.test(text)) score += 120;
    if (/(dpcs3|rdsd|dsdpbep86|3f12)/.test(text)) score += 100;
    if (/(step\s*2|anharm)/.test(text)) score -= 25;
  } else if (type === "sextic" || type === "dvib") {
    if (/(step\s*2)/.test(text)) score += 120;
    if (/(anharm|b3lyp)/.test(text)) score += 100;
    if (/(step\s*3|dpcs3|rdsd)/.test(text)) score -= 20;
  }
  return score;
}

function extractGeometryComponentFromLog(logText, pointGroup) {
  const rotationalConstantsMHz = tryExtractPenultimateRotationalConstantsMHz(logText);
  const dipole = tryExtractDipoleFromLog(logText);
  const dpcs3Atoms = tryExtractLastInputOrientationAtoms(logText);

  if (!rotationalConstantsMHz && !dipole && !dpcs3Atoms) {
    return null;
  }

  const safeAtoms = Array.isArray(dpcs3Atoms) ? dpcs3Atoms : [];
  return {
    pointGroup: pointGroup || "CS",
    rotationalConstantsMHz: rotationalConstantsMHz || { A: NaN, B: NaN, C: NaN },
    dipole: dipole || { a: "", b: "", c: "" },
    dpcs3Atoms: safeAtoms,
    dpcs3XYZText: safeAtoms.length ? buildAtomicNumberXYZText(safeAtoms) : "",
  };
}

function extractQuarticComponentFromLog(logText) {
  const asymSection = tryGetLastSection(logText, ASYM_REDUCTION_HEADER, 2600);
  const symSection = tryGetLastSection(logText, SYM_REDUCTION_HEADER, 2600);
  const asym = {
    deltaJ: extractMHzFromSectionOptional(asymSection, "DELTA\\s+J"),
    deltaJK: extractMHzFromSectionOptional(asymSection, "DELTA\\s+JK"),
    deltaK: extractMHzFromSectionOptional(asymSection, "DELTA\\s+K"),
    smallDeltaJ: extractMHzFromSectionOptional(asymSection, "delta\\s+J"),
    smallDeltaK: extractMHzFromSectionOptional(asymSection, "delta\\s+K"),
  };
  const sym = {
    DJ: extractMHzFromSectionOptional(symSection, "D\\s+J"),
    DJK: extractMHzFromSectionOptional(symSection, "D\\s+JK"),
    DK: extractMHzFromSectionOptional(symSection, "D\\s+K"),
    d1: extractMHzFromSectionOptional(symSection, "d\\s+1"),
    d2: extractMHzFromSectionOptional(symSection, "d\\s+2"),
  };
  if (!hasAnyNonEmptyValue(asym) && !hasAnyNonEmptyValue(sym)) {
    return null;
  }
  return { asym, sym };
}

function extractSexticComponentFromLog(logText) {
  const asymSexticSection = tryGetLastSection(logText, A_SEXTIC_REDUCTION_HEADER, 2600);
  const symSexticSection = tryGetLastSection(logText, S_SEXTIC_REDUCTION_HEADER, 2600);
  const asymSextic = {
    PHIN: extractHzAsMHzFromSection(asymSexticSection, "Phi\\s+J"),
    PHINK: extractHzAsMHzFromSection(asymSexticSection, "Phi\\s+JK"),
    PHIKN: extractHzAsMHzFromSection(asymSexticSection, "Phi\\s+KJ"),
    PHIK: extractHzAsMHzFromSection(asymSexticSection, "Phi\\s+K"),
    phiN: extractHzAsMHzFromSection(asymSexticSection, "phi\\s+j"),
    phiNK: extractHzAsMHzFromSection(asymSexticSection, "phi\\s+jk"),
    phiK: extractHzAsMHzFromSection(asymSexticSection, "phi\\s+k"),
  };
  const symSextic = {
    HN: extractHzAsMHzFromSection(symSexticSection, "H\\s+J"),
    HNK: extractHzAsMHzFromSection(symSexticSection, "H\\s+JK"),
    HKN: extractHzAsMHzFromSection(symSexticSection, "H\\s+KJ"),
    HK: extractHzAsMHzFromSection(symSexticSection, "H\\s+K"),
    h1: extractHzAsMHzFromSection(symSexticSection, "h\\s+1"),
    h2: extractHzAsMHzFromSection(symSexticSection, "h\\s+2"),
    h3: extractHzAsMHzFromSection(symSexticSection, "h\\s+3"),
  };
  if (!hasAnyNonEmptyValue(asymSextic) && !hasAnyNonEmptyValue(symSextic)) {
    return null;
  }
  return { asymSextic, symSextic };
}

function getLastSection(text, header, fallbackLength) {
  const section = tryGetLastSection(text, header, fallbackLength);
  if (!section) {
    throw new Error(`Section not found in Gaussian log: ${header}`);
  }
  return section;
}

function tryGetLastSection(text, header, fallbackLength) {
  const idx = text.lastIndexOf(header);
  if (idx < 0) {
    return null;
  }

  const tail = text.slice(idx + header.length);
  const nextHeaderOffset = tail.search(/\n\s*Constants in the [^\n\r]*Hamiltonian/);
  if (nextHeaderOffset >= 0) {
    return text.slice(idx, idx + header.length + nextHeaderOffset);
  }

  const end = Math.min(text.length, idx + fallbackLength);
  return text.slice(idx, end);
}

function extractAllSections(text, header, fallbackLength) {
  const sections = [];
  let from = 0;
  while (from < text.length) {
    const idx = text.indexOf(header, from);
    if (idx < 0) break;
    const tail = text.slice(idx + header.length);
    const nextHeaderOffset = tail.search(/\n\s*Constants in the [^\n\r]*Hamiltonian/);
    const sectionText =
      nextHeaderOffset >= 0
        ? text.slice(idx, idx + header.length + nextHeaderOffset)
        : text.slice(idx, Math.min(text.length, idx + fallbackLength));
    sections.push({ index: idx, text: sectionText });
    from = idx + header.length;
  }
  return sections;
}

function createContextResolver(logText) {
  const segmentMarkers = [{ index: 0, segmentId: 0 }];
  const segmentRegex = /Normal termination of Gaussian|Link1:\s+Proceeding to internal job step number/gi;
  let match;
  while ((match = segmentRegex.exec(logText)) !== null) {
    segmentMarkers.push({ index: match.index, segmentId: 0 });
  }
  segmentMarkers.sort((a, b) => a.index - b.index);
  let segmentCounter = 0;
  for (let i = 0; i < segmentMarkers.length; i += 1) {
    if (i > 0 && segmentMarkers[i].index !== segmentMarkers[i - 1].index) {
      segmentCounter += 1;
    }
    segmentMarkers[i].segmentId = segmentCounter;
  }

  const stepMarkers = [];
  const routeMarkers = [];
  const titleMarkers = [];
  const stepRegex = /^\s*(?:Molecule\s*-\s*)?Step\s*([0-9]+)\s*(?:-\s*)?[^\n\r]*$/gim;
  while ((match = stepRegex.exec(logText)) !== null) {
    const text = shortenContextLabel(normalizeSpaces(match[0]));
    if (!text || /\\/.test(text)) continue;
    const segmentId = findSegmentIdForIndex(segmentMarkers, match.index);
    stepMarkers.push({
      index: match.index,
      stepNumber: Number.parseInt(match[1], 10),
      stepLabel: text,
      segmentId,
    });
    titleMarkers.push({
      index: match.index,
      titleLabel: text,
      segmentId,
    });
  }

  const routeRegex = /^\s*#p\s+([^\n\r]+)/gim;
  while ((match = routeRegex.exec(logText)) !== null) {
    const segmentId = findSegmentIdForIndex(segmentMarkers, match.index);
    routeMarkers.push({
      index: match.index,
      routeLabel: shortenContextLabel(normalizeSpaces(match[1])),
      segmentId,
    });
  }

  const framedTitleRegex = /-{8,}\s*\r?\n\s*([^\r\n\\][^\r\n]*)\s*\r?\n\s*-{8,}/gm;
  while ((match = framedTitleRegex.exec(logText)) !== null) {
    const titleRaw = normalizeSpaces(match[1]);
    if (!titleRaw || titleRaw.length < 5) continue;
    if (isGenericGaussianSectionTitle(titleRaw)) continue;
    if (/^(?:charge|structure from|recover|redundant internal coordinates)/i.test(titleRaw)) {
      continue;
    }
    if (
      /(?:quartic centrifugal distortion constants|sextic centrifugal distortion constants|average coordinates and mean square amplitudes|optimized parameters|input orientation|asymmetric top reduction|constants in the)/i.test(
        titleRaw
      )
    ) {
      continue;
    }
    const titleLabel = shortenContextLabel(titleRaw);
    const titleIndex = match.index + match[0].indexOf(match[1]);
    titleMarkers.push({
      index: titleIndex,
      titleLabel,
      segmentId: findSegmentIdForIndex(segmentMarkers, titleIndex),
    });
  }
  titleMarkers.sort((a, b) => a.index - b.index);

  return {
    segmentMarkers,
    stepMarkers,
    routeMarkers,
    titleMarkers,
  };
}

function extractContextInfo(logText, index, resolver = null) {
  if (resolver) {
    const segmentId = findSegmentIdForIndex(resolver.segmentMarkers || [], index);
    const title = findLastBeforeWithGroup(resolver.titleMarkers, index, null, segmentId);
    const step = findLastBeforeWithGroup(resolver.stepMarkers, index, null, segmentId);
    const route = findLastBeforeWithGroup(resolver.routeMarkers, index, null, segmentId);
    const stepLabel = step ? step.stepLabel : "";
    const stepNumber = step && Number.isFinite(step.stepNumber) ? step.stepNumber : null;
    const routeLabel = route ? route.routeLabel : "";
    const titleLabel = title ? title.titleLabel : "";
    const isGenericTitle = isGenericGaussianSectionTitle(titleLabel);
    const preferredTitleLabel = isGenericTitle ? "" : titleLabel;
    let contextLabel = "Gaussian block";
    if (stepLabel && preferredTitleLabel && !/^step\s*\d+/i.test(preferredTitleLabel)) {
      contextLabel = `${stepLabel} | ${preferredTitleLabel}`;
    } else if (stepLabel && routeLabel) {
      contextLabel = `${stepLabel} | ${routeLabel}`;
    } else if (stepLabel) {
      contextLabel = stepLabel;
    } else if (preferredTitleLabel) {
      contextLabel = preferredTitleLabel;
    } else if (routeLabel) {
      contextLabel = routeLabel;
    }
    const contextCore =
      stepNumber !== null
        ? `step-${stepNumber}`
        : preferredTitleLabel
          ? `title-${normalizeContextKey(preferredTitleLabel)}`
          : routeLabel
            ? `route-${normalizeContextKey(routeLabel)}`
            : "global";
    const groupKey = `seg-${segmentId}|${contextCore}`;
    return {
      segmentId,
      titleLabel,
      stepLabel,
      stepNumber,
      routeLabel,
      contextLabel,
      groupKey,
    };
  }

  const start = Math.max(0, index - 22000);
  const windowText = logText.slice(start, index + 1);
  let stepLabel = "";
  let stepNumber = null;
  let match;
  const stepRegex = /(?:Molecule\s*-\s*)?Step\s*([0-9]+)\s*(?:-\s*)?[^\n\r]*/gi;
  while ((match = stepRegex.exec(windowText)) !== null) {
    stepNumber = Number.parseInt(match[1], 10);
    stepLabel = shortenContextLabel(normalizeSpaces(match[0]));
  }

  let routeLabel = "";
  const routeRegex = /^\s*#p\s+([^\n\r]+)/gim;
  while ((match = routeRegex.exec(windowText)) !== null) {
    routeLabel = shortenContextLabel(normalizeSpaces(match[1]));
  }

  let contextLabel = "Gaussian block";
  if (stepLabel && routeLabel) contextLabel = `${stepLabel} | ${routeLabel}`;
  else if (stepLabel) contextLabel = stepLabel;
  else if (routeLabel) contextLabel = routeLabel;

  const groupKey = stepNumber !== null
    ? `step-${stepNumber}`
    : routeLabel
      ? `route-${routeLabel.toLowerCase().slice(0, 70)}`
      : "global";

  return {
    stepLabel,
    stepNumber,
    routeLabel,
    contextLabel,
    groupKey,
  };
}

function isGenericGaussianSectionTitle(value) {
  const text = normalizeSpaces(value).toLowerCase();
  if (!text) return true;
  if (text.length < 4) return true;
  if (text.includes("!")) return true;
  return /(?:thermochemistry|rotational constants|dipole moment|data for electric dipole|input orientation|normal coordinates|force constants|harmonic frequencies|vibrational analysis|frequencies --|anharmonic frequencies|zero-point correction|sum of electronic and|molecular mass|point group|moments of inertia|cartesian coordinates|job cpu time|archive entry)/i.test(
    text
  );
}

function findSegmentIdForIndex(segmentMarkers, index) {
  if (!Array.isArray(segmentMarkers) || !segmentMarkers.length) return 0;
  let out = segmentMarkers[0].segmentId || 0;
  for (let i = 0; i < segmentMarkers.length; i += 1) {
    const marker = segmentMarkers[i];
    if (!marker || !Number.isFinite(marker.index)) continue;
    if (marker.index > index) break;
    out = Number.isFinite(marker.segmentId) ? marker.segmentId : out;
  }
  return out;
}

function normalizeContextKey(value) {
  return String(value || "")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "")
    .slice(0, 80);
}

function normalizeSpaces(value) {
  return String(value || "")
    .replace(/[\u0000-\u001f]+/g, " ")
    .replace(/\\0,1\\.*$/g, "")
    .replace(/\\/g, " ")
    .replace(/\s+/g, " ")
    .trim();
}

function shortenContextLabel(value, maxLen = 92) {
  const text = String(value || "").trim();
  if (text.length <= maxLen) return text;
  return `${text.slice(0, maxLen - 3)}...`;
}

function extractAllRotationalConstantsOccurrences(logText) {
  const regex =
    /Rotational constants\s*\((MHZ|GHZ)\)\s*:\s*([-\d+.DEde]+)\s+([-\d+.DEde]+)\s+([-\d+.DEde]+)/gi;
  const events = [];
  let match;
  while ((match = regex.exec(logText)) !== null) {
    const unit = String(match[1] || "").toUpperCase();
    const factor = unit === "GHZ" ? 1000 : 1;
    const A = parseSciNumber(match[2]);
    const B = parseSciNumber(match[3]);
    const C = parseSciNumber(match[4]);
    if (!Number.isFinite(A) || !Number.isFinite(B) || !Number.isFinite(C)) continue;
    events.push({
      index: match.index,
      rotationalConstantsMHz: {
        A: A * factor,
        B: B * factor,
        C: C * factor,
      },
    });
  }
  return events;
}

function extractAllInputOrientationOccurrences(logText, resolver = null) {
  const marker = "Input orientation:";
  const events = [];
  let from = 0;
  while (from < logText.length) {
    const idx = logText.indexOf(marker, from);
    if (idx < 0) break;
    const atoms = tryExtractInputOrientationAtomsAt(logText, idx);
    if (atoms && atoms.length) {
      const context = extractContextInfo(logText, idx, resolver);
      events.push({
        index: idx,
        groupKey: context.groupKey,
        segmentId: context.segmentId,
        atoms,
      });
    }
    from = idx + marker.length;
  }
  return events;
}

function tryExtractInputOrientationAtomsAt(logText, startIndex) {
  try {
    return extractInputOrientationAtomsAt(logText, startIndex);
  } catch (error) {
    return null;
  }
}

function extractInputOrientationAtomsAt(logText, startIndex) {
  const sectionText = logText.slice(startIndex, Math.min(logText.length, startIndex + 8000));
  const lines = sectionText.split(/\r?\n/);
  const separatorLines = [];
  for (let i = 0; i < lines.length; i += 1) {
    if (/^\s*-{5,}\s*$/.test(lines[i])) {
      separatorLines.push(i);
      if (separatorLines.length >= 3) break;
    }
  }
  if (separatorLines.length < 3) {
    throw new Error('Unable to parse separators in the "Input orientation" section.');
  }
  const dataStart = separatorLines[1] + 1;
  const dataEnd = separatorLines[2];
  const atoms = [];
  for (let i = dataStart; i < dataEnd; i += 1) {
    const line = lines[i];
    if (!line || !line.trim()) continue;
    const match = line.match(
      new RegExp(
        `^\\s*\\d+\\s+(\\d+)\\s+\\d+\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`
      )
    );
    if (!match) continue;
    const atomicNumber = Number.parseInt(match[1], 10);
    const x = parseSciNumber(match[2]);
    const y = parseSciNumber(match[3]);
    const z = parseSciNumber(match[4]);
    if (
      !Number.isFinite(atomicNumber) ||
      atomicNumber <= 0 ||
      !Number.isFinite(x) ||
      !Number.isFinite(y) ||
      !Number.isFinite(z)
    ) {
      continue;
    }
    atoms.push({
      z: atomicNumber,
      coords: [x, y, z],
    });
  }
  if (!atoms.length) {
    throw new Error('No atomic rows found in "Input orientation" section.');
  }
  return atoms;
}

function extractAllPointGroupOccurrences(logText, resolver = null) {
  const events = [];
  const framework = /Framework Group\s*:?\s*([^\r\n]+)/gi;
  const full = /Full point group\s+([A-Za-z0-9]+)/gi;
  let match;
  while ((match = framework.exec(logText)) !== null) {
    const value = normalizeSpaces(match[1]).split(/\s+/)[0];
    const context = extractContextInfo(logText, match.index, resolver);
    events.push({
      index: match.index,
      pointGroup: value,
      groupKey: context.groupKey,
      segmentId: context.segmentId,
    });
  }
  while ((match = full.exec(logText)) !== null) {
    const value = normalizeSpaces(match[1]).split(/\s+/)[0];
    const context = extractContextInfo(logText, match.index, resolver);
    events.push({
      index: match.index,
      pointGroup: value,
      groupKey: context.groupKey,
      segmentId: context.segmentId,
    });
  }
  events.sort((a, b) => a.index - b.index);
  return events;
}

function extractAllDipoleOccurrences(logText, resolver = null) {
  const events = [];
  let from = 0;
  while (from < logText.length) {
    const idx = logText.indexOf(DIPOLE_HEADER, from);
    if (idx < 0) break;
    const dipole = tryExtractDipoleFromSectionAt(logText, idx);
    if (dipole) {
      const context = extractContextInfo(logText, idx, resolver);
      events.push({
        index: idx,
        dipole,
        groupKey: context.groupKey,
        segmentId: context.segmentId,
      });
    }
    from = idx + DIPOLE_HEADER.length;
  }
  return events;
}

function tryExtractDipoleFromSectionAt(logText, idx) {
  try {
    const sectionText = logText.slice(idx, Math.min(logText.length, idx + 2200));
    const x = extractDipoleDebyeComponent(sectionText, "x");
    const y = extractDipoleDebyeComponent(sectionText, "y");
    const z = extractDipoleDebyeComponent(sectionText, "z");
    return { a: x, b: y, c: z };
  } catch (error) {
    return null;
  }
}

function findLastBefore(events, index) {
  if (!Array.isArray(events) || !events.length) return null;
  let best = null;
  for (let i = 0; i < events.length; i += 1) {
    const event = events[i];
    if (event.index > index) break;
    best = event;
  }
  return best;
}

function findLastBeforeWithGroup(events, index, groupKey, segmentId = null) {
  if (!Array.isArray(events) || !events.length) return null;
  let best = null;
  for (let i = 0; i < events.length; i += 1) {
    const event = events[i];
    if (event.index > index) break;
    if (groupKey && event.groupKey !== groupKey) continue;
    if (segmentId !== null && Number.isFinite(segmentId)) {
      if (!Number.isFinite(event.segmentId) || event.segmentId !== segmentId) continue;
    }
    best = event;
  }
  return best;
}

function findNearestOccurrence(events, index) {
  if (!Array.isArray(events) || !events.length) return null;
  let best = null;
  let bestDistance = Number.POSITIVE_INFINITY;
  for (let i = 0; i < events.length; i += 1) {
    const event = events[i];
    const distance = Math.abs(event.index - index);
    if (distance < bestDistance) {
      best = event;
      bestDistance = distance;
    }
  }
  return best;
}

function findNearestByGroup(events, index, groupKey, segmentId = null) {
  if (!Array.isArray(events) || !events.length) return null;
  let best = null;
  let bestDistance = Number.POSITIVE_INFINITY;
  for (let i = 0; i < events.length; i += 1) {
    const event = events[i];
    if (groupKey && event.groupKey !== groupKey) continue;
    if (segmentId !== null && Number.isFinite(segmentId)) {
      if (!Number.isFinite(event.segmentId) || event.segmentId !== segmentId) continue;
    }
    const distance = Math.abs(event.index - index);
    if (distance < bestDistance) {
      best = event;
      bestDistance = distance;
    }
  }
  return best;
}

function extractMHzFromSectionOptional(sectionText, labelPattern) {
  if (!sectionText) return "";
  const regex = new RegExp(
    `^\\s*${labelPattern}\\s*:\\s*(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "m"
  );
  const match = sectionText.match(regex);
  if (!match) return "";
  return normalizeFortranNotation(match[2]);
}

function extractHzAsMHzFromSection(sectionText, labelPattern) {
  if (!sectionText) return "";
  const regex = new RegExp(
    `^\\s*${labelPattern}\\s*:\\s*(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "m"
  );
  const match = sectionText.match(regex);
  if (!match) return "";

  const hzValue = parseSciNumber(match[2]);
  if (!Number.isFinite(hzValue)) return "";
  return toFortranScientific(hzValue / 1e6);
}

function tryExtractDVibFromLog(logText) {
  try {
    return extractDVibFromLog(logText);
  } catch (error) {
    return null;
  }
}

function extractDVibFromLog(logText) {
  const dvibRegex =
    /Rotational Constants\s*\(in MHz\)[\s\S]*?Ae=\s*([-\d+.DEde]+)\s*A00=\s*([-\d+.DEde]+)[\s\S]*?Be=\s*([-\d+.DEde]+)\s*B00=\s*([-\d+.DEde]+)[\s\S]*?Ce=\s*([-\d+.DEde]+)\s*C00=\s*([-\d+.DEde]+)/gi;
  let match;
  let lastMatch = null;

  while ((match = dvibRegex.exec(logText)) !== null) {
    lastMatch = match;
  }

  if (!lastMatch) {
    throw new Error("Unable to find the A00/Ae rotational constants block in MHz.");
  }

  const Ae = parseSciNumber(lastMatch[1]);
  const A00 = parseSciNumber(lastMatch[2]);
  const Be = parseSciNumber(lastMatch[3]);
  const B00 = parseSciNumber(lastMatch[4]);
  const Ce = parseSciNumber(lastMatch[5]);
  const C00 = parseSciNumber(lastMatch[6]);

  if (
    !Number.isFinite(Ae) ||
    !Number.isFinite(A00) ||
    !Number.isFinite(Be) ||
    !Number.isFinite(B00) ||
    !Number.isFinite(Ce) ||
    !Number.isFinite(C00)
  ) {
    throw new Error("Unable to compute DVib values from the A00/Ae section.");
  }

  return {
    A: A00 - Ae,
    B: B00 - Be,
    C: C00 - Ce,
  };
}

function tryExtractPenultimateRotationalConstantsMHz(logText) {
  try {
    return extractPenultimateRotationalConstantsMHz(logText);
  } catch (error) {
    return null;
  }
}

function extractPenultimateRotationalConstantsMHz(logText) {
  const regex =
    /Rotational constants\s*\((?:MHZ|GHZ)\)\s*:\s*([-\d+.DEde]+)\s+([-\d+.DEde]+)\s+([-\d+.DEde]+)/gi;
  let match;
  const matches = [];

  while ((match = regex.exec(logText)) !== null) {
    matches.push(match);
  }

  if (!matches.length) {
    throw new Error('Unable to find "Rotational constants" occurrences in the Gaussian log.');
  }

  const best = matches.length >= 2 ? matches[matches.length - 2] : matches[matches.length - 1];
  const fullLineStart = best[0].match(/\(GHZ\)/i) ? 1000 : 1;
  const a = parseSciNumber(best[1]);
  const b = parseSciNumber(best[2]);
  const c = parseSciNumber(best[3]);
  if (!Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(c)) {
    throw new Error('Unable to parse rotational constants from "Rotational constants" line.');
  }

  return {
    A: a * fullLineStart,
    B: b * fullLineStart,
    C: c * fullLineStart,
  };
}

function extractPointGroup(logText) {
  // Prefer the last "Framework Group" entry (closest to the end of the log).
  const lines = logText.split(/\r?\n/);
  for (let i = lines.length - 1; i >= 0; i -= 1) {
    const line = lines[i];
    if (!/Framework Group/i.test(line)) continue;
    const withColon = line.match(/Framework Group\s*:\s*([^\r\n]+)/i);
    if (withColon && withColon[1].trim()) {
      return withColon[1].trim();
    }
    const withoutColon = line.match(/Framework Group\s+([^\r\n]+)/i);
    if (withoutColon && withoutColon[1].trim()) {
      return withoutColon[1].trim();
    }
  }

  // Fallback if "Framework Group" is not present.
  const regex = /Full point group\s+([A-Za-z0-9]+)/gi;
  let match;
  let pointGroup = null;
  while ((match = regex.exec(logText)) !== null) {
    pointGroup = match[1];
  }
  return pointGroup;
}

function tryExtractDipoleFromLog(logText) {
  try {
    return extractDipoleFromLog(logText);
  } catch (error) {
    return null;
  }
}

function extractDipoleFromLog(logText) {
  const idx = logText.lastIndexOf(DIPOLE_HEADER);
  if (idx < 0) {
    throw new Error("Unable to find the electric dipole moment section in the Gaussian log.");
  }

  const sectionText = logText.slice(idx, Math.min(logText.length, idx + 2000));
  const x = extractDipoleDebyeComponent(sectionText, "x");
  const y = extractDipoleDebyeComponent(sectionText, "y");
  const z = extractDipoleDebyeComponent(sectionText, "z");

  return { a: x, b: y, c: z };
}

function extractDipoleDebyeComponent(sectionText, axisLabel) {
  const regex = new RegExp(
    `^\\s*${axisLabel}\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`,
    "mi"
  );
  const match = sectionText.match(regex);
  if (!match) {
    throw new Error(
      `Unable to parse dipole component along ${axisLabel} from the Gaussian log.`
    );
  }
  return normalizeFortranNotation(match[2]);
}

function tryExtractLastInputOrientationAtoms(logText) {
  try {
    return extractLastInputOrientationAtoms(logText);
  } catch (error) {
    return null;
  }
}

function extractLastInputOrientationAtoms(logText) {
  const marker = "Input orientation:";
  const idx = logText.lastIndexOf(marker);
  if (idx < 0) {
    throw new Error('Unable to find "Input orientation" in the Gaussian log.');
  }

  const sectionText = logText.slice(idx, Math.min(logText.length, idx + 8000));
  const lines = sectionText.split(/\r?\n/);
  const separatorLines = [];

  for (let i = 0; i < lines.length; i += 1) {
    if (/^\s*-{5,}\s*$/.test(lines[i])) {
      separatorLines.push(i);
      if (separatorLines.length >= 3) break;
    }
  }

  if (separatorLines.length < 3) {
    throw new Error('Unable to parse separators in the last "Input orientation" section.');
  }

  const dataStart = separatorLines[1] + 1;
  const dataEnd = separatorLines[2];
  const atoms = [];

  for (let i = dataStart; i < dataEnd; i += 1) {
    const line = lines[i];
    if (!line || !line.trim()) continue;
    const match = line.match(
      new RegExp(
        `^\\s*\\d+\\s+(\\d+)\\s+\\d+\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s+(${NUMBER_TOKEN_PATTERN})\\s*$`
      )
    );
    if (!match) continue;

    const atomicNumber = Number.parseInt(match[1], 10);
    const x = parseSciNumber(match[2]);
    const y = parseSciNumber(match[3]);
    const z = parseSciNumber(match[4]);
    if (
      !Number.isFinite(atomicNumber) ||
      atomicNumber <= 0 ||
      !Number.isFinite(x) ||
      !Number.isFinite(y) ||
      !Number.isFinite(z)
    ) {
      continue;
    }

    atoms.push({
      z: atomicNumber,
      coords: [x, y, z],
    });
  }

  if (!atoms.length) {
    throw new Error('No atomic rows found in the last "Input orientation" section.');
  }
  return atoms;
}

function buildAtomicNumberXYZText(atoms) {
  const lines = [];
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    lines.push(
      `${atom.z} ${formatNumber(atom.coords[0])} ${formatNumber(atom.coords[1])} ${formatNumber(atom.coords[2])}`
    );
  }
  return `${lines.join("\n")}\n`;
}

function parseSciNumber(value) {
  const parsed = Number.parseFloat(String(value || "").replace(/D/gi, "E"));
  return Number.isFinite(parsed) ? parsed : null;
}

function normalizeFortranNotation(value) {
  return String(value || "").trim().replace(/e/gi, "D");
}

function toFortranScientific(value, fractionDigits = 10) {
  if (!Number.isFinite(value)) return "";
  const normalized = Math.abs(value) < 5e-30 ? 0 : value;
  return normalized.toExponential(fractionDigits).replace(/e/i, "D");
}

function parseFchkGeometry(text) {
  const atomicField = parseFchkNumberField(text, "Atomic numbers");
  const coordField = parseFchkNumberField(text, "Current cartesian coordinates");
  const atomCount = atomicField.values.length;
  if (!atomCount) {
    throw new Error("No atomic numbers found in FCHK geometry.");
  }
  if (coordField.values.length !== atomCount * 3) {
    throw new Error(
      `FCHK coordinate count mismatch: expected ${atomCount * 3} values, found ${coordField.values.length}.`
    );
  }

  const atoms = [];
  for (let i = 0; i < atomCount; i += 1) {
    const zValue = Math.round(atomicField.values[i]);
    if (!Number.isFinite(zValue) || zValue <= 0) {
      throw new Error(`Invalid atomic number in FCHK at index ${i + 1}.`);
    }
    const xBohr = coordField.values[3 * i + 0];
    const yBohr = coordField.values[3 * i + 1];
    const zBohr = coordField.values[3 * i + 2];
    atoms.push({
      label: String(zValue),
      z: zValue,
      coords: [
        xBohr * BOHR_TO_ANGSTROM,
        yBohr * BOHR_TO_ANGSTROM,
        zBohr * BOHR_TO_ANGSTROM,
      ],
    });
  }
  return { atoms };
}

function parseFchkDipoleMomentOptional(text) {
  try {
    const field = parseFchkNumberField(text, "Dipole Moment");
    if (!Array.isArray(field.values) || field.values.length < 3) return null;
    return [field.values[0], field.values[1], field.values[2]];
  } catch (error) {
    return null;
  }
}

function parseFchkNumberField(text, headerLabel) {
  const source = String(text || "");
  const escapedHeader = escapeRegExp(headerLabel);
  const headerRegex = new RegExp(`${escapedHeader}\\s+[IR]\\s+N=\\s*(\\d+)`, "i");
  const headerMatch = headerRegex.exec(source);
  if (!headerMatch) {
    throw new Error(`Unable to find "${headerLabel} ... N=..." in the selected FCHK file.`);
  }
  const expectedCount = Number.parseInt(headerMatch[1], 10);
  if (!Number.isFinite(expectedCount) || expectedCount <= 0) {
    throw new Error(`Invalid entry count for "${headerLabel}" in FCHK file.`);
  }
  const tail = source.slice(headerMatch.index + headerMatch[0].length);
  const numberRegex = new RegExp(NUMBER_TOKEN_PATTERN, "g");
  const tokens = tail.match(numberRegex) || [];
  if (tokens.length < expectedCount) {
    throw new Error(
      `FCHK section "${headerLabel}" is truncated: expected ${expectedCount} values, found ${tokens.length}.`
    );
  }
  const values = [];
  for (let i = 0; i < expectedCount; i += 1) {
    const parsed = parseSciNumber(tokens[i]);
    if (!Number.isFinite(parsed)) {
      throw new Error(`Invalid numeric value in FCHK section "${headerLabel}" at position ${i + 1}.`);
    }
    values.push(parsed);
  }
  return {
    count: expectedCount,
    values,
  };
}

function escapeRegExp(text) {
  return String(text || "").replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function parseFchkCartesianForceConstants(text) {
  const source = String(text || "");
  const headerRegex = /Cartesian Force Constants\s+R\s+N=\s*(\d+)/i;
  const headerMatch = headerRegex.exec(source);
  if (!headerMatch) {
    throw new Error('Unable to find "Cartesian Force Constants ... N=..." in the selected FCHK file.');
  }
  const lowerTriangularCount = Number.parseInt(headerMatch[1], 10);
  if (!Number.isFinite(lowerTriangularCount) || lowerTriangularCount <= 0) {
    throw new Error("Invalid lower-triangular count in FCHK Cartesian Force Constants header.");
  }

  const dimension = inferLowerTriangularDimension(lowerTriangularCount);
  const tail = source.slice(headerMatch.index + headerMatch[0].length);
  const valueRegex = new RegExp(NUMBER_TOKEN_PATTERN, "g");
  const tokens = tail.match(valueRegex) || [];
  if (tokens.length < lowerTriangularCount) {
    throw new Error(
      `FCHK Cartesian Force Constants section is truncated: expected ${lowerTriangularCount} values, found ${tokens.length}.`
    );
  }

  const values = [];
  for (let i = 0; i < lowerTriangularCount; i += 1) {
    const parsed = parseSciNumber(tokens[i]);
    if (!Number.isFinite(parsed)) {
      throw new Error(`Invalid numeric value in FCHK Cartesian Force Constants at position ${i + 1}.`);
    }
    values.push(parsed);
  }

  const matrix = Array.from({ length: dimension }, () => Array(dimension).fill(0));
  let cursor = 0;
  for (let i = 0; i < dimension; i += 1) {
    for (let j = 0; j <= i; j += 1) {
      const value = values[cursor];
      matrix[i][j] = value;
      matrix[j][i] = value;
      cursor += 1;
    }
  }
  return {
    lowerTriangularCount,
    dimension,
    matrix,
  };
}

function inferLowerTriangularDimension(lowerTriangularCount) {
  const delta = 1 + 8 * lowerTriangularCount;
  const root = Math.sqrt(delta);
  const raw = (root - 1) / 2;
  const dimension = Math.round(raw);
  if (!Number.isFinite(root) || Math.abs(raw - dimension) > 1e-8 || dimension <= 0) {
    throw new Error(
      `Invalid FCHK lower-triangular count (${lowerTriangularCount}): unable to infer square matrix dimension.`
    );
  }
  return dimension;
}

function cloneAtoms(atoms) {
  if (!Array.isArray(atoms)) return [];
  const out = [];
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    if (!atom || !Array.isArray(atom.coords) || atom.coords.length < 3) continue;
    out.push({
      label: atom.label,
      z: atom.z,
      coords: [Number(atom.coords[0]), Number(atom.coords[1]), Number(atom.coords[2])],
    });
  }
  return out;
}

function cloneMatrix(matrix) {
  if (!Array.isArray(matrix)) return [];
  return matrix.map((row) => (Array.isArray(row) ? row.slice() : []));
}

function isApproximatelySymmetricMatrix(matrix, tolerance = 1e-8) {
  if (!Array.isArray(matrix) || !matrix.length) return false;
  const n = matrix.length;
  for (let i = 0; i < n; i += 1) {
    if (!Array.isArray(matrix[i]) || matrix[i].length !== n) return false;
    for (let j = i + 1; j < n; j += 1) {
      const a = Number(matrix[i][j]);
      const b = Number(matrix[j][i]);
      if (!Number.isFinite(a) || !Number.isFinite(b)) return false;
      if (Math.abs(a - b) > tolerance) return false;
    }
  }
  return true;
}

function matrixTrace(matrix) {
  if (!Array.isArray(matrix) || !matrix.length) return Number.NaN;
  const n = matrix.length;
  let sum = 0;
  for (let i = 0; i < n; i += 1) {
    if (!Array.isArray(matrix[i]) || matrix[i].length !== n) return Number.NaN;
    const value = Number(matrix[i][i]);
    if (!Number.isFinite(value)) return Number.NaN;
    sum += value;
  }
  return sum;
}

function matrixFrobeniusNorm(matrix) {
  if (!Array.isArray(matrix) || !matrix.length) return Number.NaN;
  let sumSquares = 0;
  for (let i = 0; i < matrix.length; i += 1) {
    const row = matrix[i];
    if (!Array.isArray(row)) return Number.NaN;
    for (let j = 0; j < row.length; j += 1) {
      const value = Number(row[j]);
      if (!Number.isFinite(value)) return Number.NaN;
      sumSquares += value * value;
    }
  }
  return Math.sqrt(sumSquares);
}

function buildFchkDerivedPayloads(bundle, representation) {
  const rawAtoms = normalizeAtomsForQuartic(bundle.atoms);
  if (!rawAtoms.length) {
    throw new Error(`Unable to build FCHK geometry payload for "${bundle.fileName}".`);
  }
  const oriented = orientGeometryForRepresentation(rawAtoms, representation);
  const rotatedAtoms = oriented.atoms;
  const rotationalConstantsMHz = estimateRotationalConstantsMHz(rotatedAtoms);
  const dipole = formatFchkDipoleForRepresentation(bundle.dipoleMomentAu, oriented.rotation3);
  const sourceIndex = Date.now();
  const geometryVariant = {
    contextLabel: `FCHK geometry (rotated ${representation})`,
    sourceIndex,
    inputRepresentation: normalizeWatsonRepresentationToken(representation) || "Ir",
    pointGroup: "C1",
    rotationalConstantsMHz,
    dipole,
    dpcs3Atoms: rotatedAtoms.map((atom) => ({
      z: atom.z,
      coords: [atom.coords[0], atom.coords[1], atom.coords[2]],
    })),
    dpcs3XYZText: buildAtomicNumberXYZText(rotatedAtoms),
  };

  const harmonic = computeHarmonicQuarticFromHessian({
    hessianCartesian: bundle.hessianMatrix,
    geometryAtoms: rawAtoms,
    representation,
  });
  const coriolisLabel =
    harmonic && harmonic.coriolisModel && harmonic.coriolisModel.enabled
      ? " + Coriolis(harmonic)"
      : "";
  const quarticVariant = {
    contextLabel: `FCHK quartic harmonic${coriolisLabel} (rotated ${representation})`,
    sourceIndex: sourceIndex + 1,
    asym: { ...harmonic.asymFormatted },
    sym: { ...harmonic.symFormatted },
    inputRepresentation: normalizeWatsonRepresentationToken(representation) || "Ir",
  };

  return { geometry: geometryVariant, quartic: quarticVariant };
}

function orientGeometryForRepresentation(atoms, representation) {
  const center = computeCenterOfMass(atoms);
  const inertia = computeInertiaTensor(atoms, center);
  const principal = diagonalizeSymmetricTensor(inertia);
  const principalAtoms = atoms.map((atom) => ({
    label: atom.label,
    z: atom.z,
    coords: rotateToBasis(subtractVec(atom.coords, center), principal.basis),
  }));
  const rotatedAtoms = applyRepresentation(principalAtoms, representation);
  const rotation3 = buildRepresentationRotationMatrix(principal.basis, representation);
  return {
    atoms: rotatedAtoms,
    rotation3,
  };
}

function formatFchkDipoleForRepresentation(dipoleMomentAu, rotation3) {
  if (!Array.isArray(dipoleMomentAu) || dipoleMomentAu.length < 3) {
    return { a: "", b: "", c: "" };
  }
  if (!Array.isArray(rotation3) || rotation3.length !== 3) {
    return { a: "", b: "", c: "" };
  }
  const rotated = applyRotation3ToVector(
    [Number(dipoleMomentAu[0]), Number(dipoleMomentAu[1]), Number(dipoleMomentAu[2])],
    rotation3
  );
  const toDebyeText = (value) =>
    Number.isFinite(value) ? toFortranScientific(value * AU_TO_DEBYE, 10) : "";
  return {
    a: toDebyeText(rotated[0]),
    b: toDebyeText(rotated[1]),
    c: toDebyeText(rotated[2]),
  };
}

function formatInputDipoleForRepresentation(dipoleVectorDebye, rotation3) {
  if (!Array.isArray(dipoleVectorDebye) || dipoleVectorDebye.length < 3) {
    return null;
  }
  if (!Array.isArray(rotation3) || rotation3.length !== 3) {
    return null;
  }
  const vector = [
    Number(dipoleVectorDebye[0]),
    Number(dipoleVectorDebye[1]),
    Number(dipoleVectorDebye[2]),
  ];
  if (!Number.isFinite(vector[0]) || !Number.isFinite(vector[1]) || !Number.isFinite(vector[2])) {
    return null;
  }
  const rotated = applyRotation3ToVector(vector, rotation3);
  const formatted = {
    a: Number.isFinite(rotated[0]) ? toFortranScientific(rotated[0], 10) : "",
    b: Number.isFinite(rotated[1]) ? toFortranScientific(rotated[1], 10) : "",
    c: Number.isFinite(rotated[2]) ? toFortranScientific(rotated[2], 10) : "",
  };
  if (!formatted.a || !formatted.b || !formatted.c) {
    return null;
  }
  return formatted;
}

function applyRotation3ToVector(vector, rotation3) {
  return [
    rotation3[0][0] * vector[0] + rotation3[0][1] * vector[1] + rotation3[0][2] * vector[2],
    rotation3[1][0] * vector[0] + rotation3[1][1] * vector[1] + rotation3[1][2] * vector[2],
    rotation3[2][0] * vector[0] + rotation3[2][1] * vector[1] + rotation3[2][2] * vector[2],
  ];
}

function estimateRotationalConstantsMHz(atoms) {
  const center = computeCenterOfMass(atoms);
  const inertia = computeInertiaTensor(atoms, center);
  const principal = diagonalizeSymmetricTensor(inertia);
  const moments = principal.moments;
  const factor = 505379.0094;
  const toMHz = (moment) =>
    Number.isFinite(moment) && moment > 1e-14 ? factor / moment : Number.NaN;
  return {
    A: toMHz(moments[0]),
    B: toMHz(moments[1]),
    C: toMHz(moments[2]),
  };
}

function parseNumpyMatrixText(text) {
  const raw = String(text || "").replace(/\r/g, "\n");
  const withoutComments = raw.replace(/#.*$/gm, "");
  const normalized = withoutComments
    .replace(/\b(?:np\.)?array\s*\(/gi, "(")
    .replace(/,\s*dtype\s*=\s*[^)\]]+/gi, "")
    .replace(/\]\s*,\s*\[/g, "]\n[")
    .replace(/\]\s*\[/g, "]\n[")
    .replace(/;/g, "\n");

  const numberRegex = new RegExp(NUMBER_TOKEN_PATTERN, "g");
  const lines = normalized
    .split(/\n+/)
    .map((line) => line.trim())
    .filter(Boolean);
  const rows = [];
  for (let i = 0; i < lines.length; i += 1) {
    const tokens = lines[i].match(numberRegex);
    if (!tokens || !tokens.length) continue;
    const row = [];
    for (let j = 0; j < tokens.length; j += 1) {
      const value = parseSciNumber(tokens[j]);
      if (!Number.isFinite(value)) {
        throw new Error(`Invalid numeric token in Hessian row ${rows.length + 1}.`);
      }
      row.push(value);
    }
    rows.push(row);
  }

  if (!rows.length) {
    throw new Error("Unable to parse any numeric rows from the Hessian input.");
  }

  let matrix = rows;
  if (rows.length === 1 && rows[0].length > 1) {
    const side = Math.round(Math.sqrt(rows[0].length));
    if (side > 1 && side * side === rows[0].length) {
      matrix = [];
      for (let i = 0; i < side; i += 1) {
        matrix.push(rows[0].slice(i * side, (i + 1) * side));
      }
    }
  }

  const width = matrix[0].length;
  for (let i = 0; i < matrix.length; i += 1) {
    if (matrix[i].length !== width) {
      throw new Error("Hessian rows have different lengths.");
    }
  }
  if (matrix.length !== width) {
    throw new Error(`Hessian must be square. Parsed shape: ${matrix.length} x ${width}.`);
  }
  return matrix;
}

function parseHessianXYZBundleText(text) {
  const source = String(text || "").replace(/\r/g, "\n");
  const sections = parseBracketedSections(source);
  if (!sections.xyz || !sections.hessian) {
    throw new Error(
      'Invalid Hessian+XYZ file: expected both [XYZ] and [HESSIAN] sections.'
    );
  }
  const parsedXYZ = parseXYZ(sections.xyz);
  const hessianMatrix = parseNumpyMatrixText(sections.hessian);
  const expectedSize = parsedXYZ.atoms.length * 3;
  if (hessianMatrix.length !== expectedSize) {
    throw new Error(
      `Hessian/XYZ size mismatch: ${parsedXYZ.atoms.length} atoms require a ${expectedSize}x${expectedSize} Hessian, found ${hessianMatrix.length}x${hessianMatrix.length}.`
    );
  }
  const pointGroup = parsePointGroupSection(sections.point_group || sections.pointgroup || "");
  const dipoleVectorDebye = parseDipoleSection(
    sections.dipole || sections.dipole_debye || ""
  );
  return {
    atoms: parsedXYZ.atoms,
    hessianMatrix,
    pointGroup: pointGroup || null,
    dipoleVectorDebye: dipoleVectorDebye || null,
  };
}

function parseBracketedSections(text) {
  const source = String(text || "");
  const matches = [];
  const headerRegex = /^\s*\[(XYZ|HESSIAN|POINT_GROUP|POINTGROUP|DIPOLE|DIPOLE_DEBYE)\]\s*$/gim;
  let match;
  while ((match = headerRegex.exec(source)) !== null) {
    const rawKey = String(match[1] || "").toUpperCase();
    const key = rawKey
      .replace(/[^A-Z0-9]+/g, "_")
      .replace(/^_+|_+$/g, "")
      .toLowerCase();
    matches.push({
      key,
      start: match.index,
      contentStart: headerRegex.lastIndex,
    });
  }
  if (!matches.length) {
    return {};
  }
  const out = {};
  for (let i = 0; i < matches.length; i += 1) {
    const current = matches[i];
    const next = matches[i + 1];
    const contentEnd = next ? next.start : source.length;
    out[current.key] = source.slice(current.contentStart, contentEnd).trim();
  }
  return out;
}

function parsePointGroupSection(text) {
  const cleaned = String(text || "")
    .split(/\r?\n/)
    .map((line) => line.replace(/#.*$/, "").trim())
    .filter(Boolean);
  if (!cleaned.length) return null;
  return normalizePointGroupToken(cleaned[0]);
}

function parseDipoleSection(text) {
  const body = String(text || "");
  if (!body.trim()) return null;

  const labeled = {};
  const labeledRegex = new RegExp(
    `\\b([abcxyz])\\b\\s*(?:=|:)?\\s*(${NUMBER_TOKEN_PATTERN})`,
    "gi"
  );
  let labeledMatch;
  while ((labeledMatch = labeledRegex.exec(body)) !== null) {
    const axis = String(labeledMatch[1] || "").toLowerCase();
    const value = parseSciNumber(labeledMatch[2]);
    if (Number.isFinite(value)) {
      labeled[axis] = value;
    }
  }
  if (
    Number.isFinite(labeled.a) &&
    Number.isFinite(labeled.b) &&
    Number.isFinite(labeled.c)
  ) {
    return [labeled.a, labeled.b, labeled.c];
  }
  if (
    Number.isFinite(labeled.x) &&
    Number.isFinite(labeled.y) &&
    Number.isFinite(labeled.z)
  ) {
    return [labeled.x, labeled.y, labeled.z];
  }

  const numericTokens = body.match(new RegExp(NUMBER_TOKEN_PATTERN, "g")) || [];
  if (numericTokens.length >= 3) {
    const x = parseSciNumber(numericTokens[0]);
    const y = parseSciNumber(numericTokens[1]);
    const z = parseSciNumber(numericTokens[2]);
    if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
      return [x, y, z];
    }
  }
  return null;
}

function normalizePointGroupToken(value) {
  const token = String(value || "")
    .trim()
    .replace(/\s+/g, "")
    .replace(/[^A-Za-z0-9]/g, "")
    .toUpperCase();
  return token || null;
}

function buildTauCanonicalIndexLookup() {
  const lookup = {};
  for (let i = 0; i < TAU_COMPONENT_KEYS.length; i += 1) {
    const canonical = TAU_COMPONENT_KEYS[i].slice(1);
    const equivalents = expandTauEquivalentIndices(canonical);
    for (let j = 0; j < equivalents.length; j += 1) {
      lookup[equivalents[j]] = canonical;
    }
  }
  return lookup;
}

function expandTauEquivalentIndices(indices) {
  const seed = normalizeTauIndexToken(indices);
  if (!seed) return [];
  const out = new Set([seed]);
  const queue = [seed];
  while (queue.length) {
    const value = queue.shift();
    const variants = [
      `${value[1]}${value[0]}${value[2]}${value[3]}`,
      `${value[0]}${value[1]}${value[3]}${value[2]}`,
      `${value[2]}${value[3]}${value[0]}${value[1]}`,
    ];
    for (let i = 0; i < variants.length; i += 1) {
      const candidate = variants[i];
      if (!out.has(candidate)) {
        out.add(candidate);
        queue.push(candidate);
      }
    }
  }
  return Array.from(out);
}

function normalizeTauIndexToken(value) {
  const compact = String(value || "")
    .toLowerCase()
    .replace(/[^abc]/g, "");
  if (compact.length !== 4) return null;
  return compact;
}

function canonicalTauComponentKeyFromToken(value) {
  const normalized = normalizeTauIndexToken(value);
  if (!normalized) return null;
  const canonical = TAU_CANONICAL_INDEX_LOOKUP[normalized] || null;
  return canonical ? `t${canonical}` : null;
}

function parseTauPayloadText(text) {
  const source = String(text || "").replace(/\r/g, "\n");
  const tau = parseTauComponentsFromText(source);
  const sigma = parseTauSigmaFromText(source);
  const missing = TAU_COMPONENT_KEYS.filter((key) => !Number.isFinite(tau[key]));
  if (missing.length) {
    throw new Error(
      `Tau file is missing required components: ${missing.join(", ")}.`
    );
  }
  return {
    inputRepresentation: parseTauRepresentationFromText(source) || "Ir",
    tau,
    rotationalConstantsMHz: parseTauRotationalConstantsFromText(source),
    moments: parseTauMomentsFromText(source),
    sigma,
    sigmaSource: Number.isFinite(sigma) ? "sigma" : "none",
  };
}

function parseTauComponentsFromText(text) {
  const source = String(text || "");
  const out = {};
  for (let i = 0; i < TAU_COMPONENT_KEYS.length; i += 1) {
    out[TAU_COMPONENT_KEYS[i]] = null;
  }
  const patterns = [
    new RegExp(
      `\\b(?:taup?|tau|t)\\s*[_\\-\\s]*(?:\\(|\\[)?([abc]{4})(?:\\)|\\])?\\s*(?:=|:)?\\s*(${NUMBER_TOKEN_PATTERN})`,
      "gi"
    ),
    new RegExp(`\\b([abc]{4})\\b\\s*(?:=|:)?\\s*(${NUMBER_TOKEN_PATTERN})`, "gi"),
  ];
  for (let p = 0; p < patterns.length; p += 1) {
    const regex = patterns[p];
    regex.lastIndex = 0;
    let match;
    while ((match = regex.exec(source)) !== null) {
      const key = canonicalTauComponentKeyFromToken(match[1]);
      if (!key) continue;
      const parsed = parseSciNumber(match[2]);
      if (!Number.isFinite(parsed)) continue;
      out[key] = parsed;
    }
  }
  return out;
}

function parseTauRepresentationFromText(text) {
  const source = String(text || "");
  const taggedMatch = source.match(
    /\b(?:input_representation|output_representation)\b\s*[:=]\s*([Ii]{1,3}\s*[RrLl])\b/i
  );
  if (taggedMatch) {
    return normalizeWatsonRepresentationToken(taggedMatch[1]);
  }
  const explicitMatch = source.match(
    /\b(?:representation|rep|rappresentazione)\b\s*[:=]\s*([Ii]{1,3}\s*[RrLl])\b/i
  );
  if (explicitMatch) {
    return normalizeWatsonRepresentationToken(explicitMatch[1]);
  }
  return parseLastReductionRepresentationFromText(source);
}

function normalizeWatsonRepresentationToken(token) {
  const compact = String(token || "").replace(/\s+/g, "");
  const match = compact.match(/^(I{1,3})([RrLl])$/);
  if (!match) return null;
  const order = match[1].toUpperCase();
  const hand = String(match[2] || "").toLowerCase();
  if (order !== "I" && order !== "II" && order !== "III") {
    return null;
  }
  return `${order}${hand}`;
}

function parseLastReductionRepresentationFromText(text) {
  const source = String(text || "");
  const regex = /\b(?:Asymmetric|Symmetric)\s+Reduction\s+([Ii]{1,3}\s*[RrLl])\b/gi;
  let match;
  let last = null;
  while ((match = regex.exec(source)) !== null) {
    const normalized = normalizeWatsonRepresentationToken(match[1]);
    if (normalized) {
      last = normalized;
    }
  }
  return last;
}

function resolveLogSectionInputRepresentation(logText, sectionIndex, sectionText) {
  const source = String(logText || "");
  const safeIndex = Number.isFinite(sectionIndex)
    ? Math.max(0, Math.min(source.length, Math.trunc(sectionIndex)))
    : source.length;
  const fromRotorSymmetry = findLastRotorSymmetryRepresentationBefore(source, safeIndex);
  if (fromRotorSymmetry) return fromRotorSymmetry;

  const fromSection = parseTauRepresentationFromText(sectionText);
  if (fromSection) return fromSection;

  const lookbackStart = Math.max(0, safeIndex - 12000);
  const boundaryIndex = findLastLogBoundaryMarker(source, safeIndex);
  const contextStart = boundaryIndex >= 0 ? Math.max(lookbackStart, boundaryIndex) : lookbackStart;
  const preceding = source.slice(contextStart, safeIndex);
  return parseLastReductionRepresentationFromText(preceding) || "Ir";
}

function findLastRotorSymmetryRepresentationBefore(logText, beforeIndex) {
  const source = String(logText || "");
  const safeBefore = Number.isFinite(beforeIndex)
    ? Math.max(0, Math.min(source.length, Math.trunc(beforeIndex)))
    : source.length;
  const searchArea = source.slice(0, safeBefore);
  const regex =
    /Analysis of the Rotor Symmetry[\s\S]{0,3000}?Representation\s*:\s*([Ii]{1,3}\s*[RrLl])\s*Representation\b/gi;
  let match;
  let last = null;
  while ((match = regex.exec(searchArea)) !== null) {
    const normalized = normalizeWatsonRepresentationToken(match[1]);
    if (normalized) {
      last = normalized;
    }
  }
  return last;
}

function findLastLogBoundaryMarker(text, beforeIndex) {
  const source = String(text || "");
  const safeBefore = Number.isFinite(beforeIndex)
    ? Math.max(0, Math.min(source.length, Math.trunc(beforeIndex)))
    : source.length;
  const markers = [
    "Normal termination of Gaussian",
    "Link1:  Proceeding to internal job step number",
    "Link1: Proceeding to internal job step number",
    "Initial command:",
  ];
  let last = -1;
  for (let i = 0; i < markers.length; i += 1) {
    const idx = source.lastIndexOf(markers[i], safeBefore);
    if (idx > last) {
      last = idx;
    }
  }
  return last;
}

function parseTauRotationalConstantsFromText(text) {
  const A = parseLabeledSciNumber(text, "A(?:_?MHz)?");
  const B = parseLabeledSciNumber(text, "B(?:_?MHz)?");
  const C = parseLabeledSciNumber(text, "C(?:_?MHz)?");
  if (!Number.isFinite(A) || !Number.isFinite(B) || !Number.isFinite(C)) {
    return null;
  }
  return { a: A, b: B, c: C };
}

function parseTauMomentsFromText(text) {
  const Ia = parseLabeledSciNumber(text, "I[_\\s]*a");
  const Ib = parseLabeledSciNumber(text, "I[_\\s]*b");
  const Ic = parseLabeledSciNumber(text, "I[_\\s]*c");
  if (!Number.isFinite(Ia) || !Number.isFinite(Ib) || !Number.isFinite(Ic)) {
    return null;
  }
  return { a: Ia, b: Ib, c: Ic };
}

function parseTauSigmaFromText(text) {
  const sigmaLabels = [
    "sigma",
    "s[_\\s]*asym",
    "asymmetry(?:[_\\s]*parameter)?",
  ];
  for (let i = 0; i < sigmaLabels.length; i += 1) {
    const value = parseLabeledSciNumber(text, sigmaLabels[i]);
    if (Number.isFinite(value)) return value;
  }
  const inlineMatch = String(text || "").match(
    new RegExp(`(?:^|\\n)\\s*s\\s*[:=]\\s*(${NUMBER_TOKEN_PATTERN})(?=\\s*(?:$|[#!;]))`, "im")
  );
  if (!inlineMatch) return null;
  const parsed = parseSciNumber(inlineMatch[1]);
  return Number.isFinite(parsed) ? parsed : null;
}

function parseLabeledSciNumber(text, labelPattern) {
  const regex = new RegExp(
    `(?:^|\\n)\\s*(?:${labelPattern})\\b\\s*[:=]\\s*(${NUMBER_TOKEN_PATTERN})`,
    "im"
  );
  const match = String(text || "").match(regex);
  if (!match) return null;
  const parsed = parseSciNumber(match[1]);
  return Number.isFinite(parsed) ? parsed : null;
}

function normalizeTauRepresentation(value) {
  const raw = String(value || "").trim();
  if (!raw) return "Ir";
  if (TAU_REP_TO_RIGHT_HANDED[raw]) {
    return TAU_REP_TO_RIGHT_HANDED[raw];
  }
  const keys = Object.keys(TAU_REP_TO_RIGHT_HANDED);
  for (let i = 0; i < keys.length; i += 1) {
    if (keys[i].toLowerCase() === raw.toLowerCase()) {
      return TAU_REP_TO_RIGHT_HANDED[keys[i]];
    }
  }
  return "Ir";
}

function buildTauAxisRemap(inputRepresentation, outputRepresentation) {
  const from = normalizeTauRepresentation(inputRepresentation);
  const to = normalizeTauRepresentation(outputRepresentation);
  const fromMap = TAU_RIGHT_HANDED_AXIS_MAP[from] || TAU_RIGHT_HANDED_AXIS_MAP.Ir;
  const toMap = TAU_RIGHT_HANDED_AXIS_MAP[to] || TAU_RIGHT_HANDED_AXIS_MAP.Ir;
  const principalToInputAxis = {};
  principalToInputAxis[fromMap.a] = "a";
  principalToInputAxis[fromMap.b] = "b";
  principalToInputAxis[fromMap.c] = "c";
  return {
    inputRightHanded: from,
    outputRightHanded: to,
    axisOutToIn: {
      a: principalToInputAxis[toMap.a],
      b: principalToInputAxis[toMap.b],
      c: principalToInputAxis[toMap.c],
    },
  };
}

function getTauComponentValueByIndices(tau, indexToken) {
  const key = canonicalTauComponentKeyFromToken(indexToken);
  if (!key) return null;
  const value = tau ? tau[key] : null;
  return Number.isFinite(value) ? value : null;
}

function rotateTauComponentsForRepresentation(tau, inputRepresentation, outputRepresentation) {
  const remap = buildTauAxisRemap(inputRepresentation, outputRepresentation);
  const rotated = {};
  for (let i = 0; i < TAU_COMPONENT_KEYS.length; i += 1) {
    const key = TAU_COMPONENT_KEYS[i];
    const outputIndex = key.slice(1);
    const mappedInputIndex =
      remap.axisOutToIn[outputIndex[0]] +
      remap.axisOutToIn[outputIndex[1]] +
      remap.axisOutToIn[outputIndex[2]] +
      remap.axisOutToIn[outputIndex[3]];
    const value = getTauComponentValueByIndices(tau, mappedInputIndex);
    if (!Number.isFinite(value)) {
      throw new Error(`Unable to rotate tau component "${key}" to representation ${outputRepresentation}.`);
    }
    rotated[key] = value;
  }
  return {
    tau: rotated,
    remap,
  };
}

function rotateAxisTripletByRepresentation(values, inputRepresentation, outputRepresentation) {
  if (!values || !Number.isFinite(values.a) || !Number.isFinite(values.b) || !Number.isFinite(values.c)) {
    return null;
  }
  const mapping = REP_MAP[normalizeFullWatsonRepresentationToken(outputRepresentation)] || REP_MAP.Ir;
  return {
    a: values[mapping[0]],
    b: values[mapping[1]],
    c: values[mapping[2]],
  };
}

const QUARTIC_A_PAYLOAD_KEYS = ["deltaJ", "deltaJK", "deltaK", "smallDeltaJ", "smallDeltaK"];
const QUARTIC_S_PAYLOAD_KEYS = ["DJ", "DJK", "DK", "d1", "d2"];
const SEXTIC_A_PAYLOAD_KEYS = ["PHIN", "PHINK", "PHIKN", "PHIK", "phiN", "phiNK", "phiK"];
const SEXTIC_S_PAYLOAD_KEYS = ["HN", "HNK", "HKN", "HK", "h1", "h2", "h3"];

function normalizeFullWatsonRepresentationToken(token, fallback = "Ir") {
  return normalizeWatsonRepresentationToken(token) || fallback;
}

function inverseFullWatsonRepresentationToken(token, fallback = "Ir") {
  const rep = normalizeFullWatsonRepresentationToken(token, fallback);
  const axes = REP_MAP[rep] || REP_MAP[fallback] || REP_MAP.Ir;
  const inverse = [];
  for (let i = 0; i < axes.length; i += 1) {
    inverse[AXIS_INDEX[axes[i]]] = PRINCIPAL_AXIS_LABELS[i];
  }
  const inverseKey = inverse.join("");
  const reps = Object.keys(REP_MAP);
  for (let i = 0; i < reps.length; i += 1) {
    if (REP_MAP[reps[i]].join("") === inverseKey) return reps[i];
  }
  return fallback;
}

function rotateRepresentationABC(values, inputRepresentation, outputRepresentation) {
  if (!values) return null;
  const from = REP_MAP[normalizeFullWatsonRepresentationToken(inputRepresentation)] || REP_MAP.Ir;
  const to = REP_MAP[normalizeFullWatsonRepresentationToken(outputRepresentation)] || REP_MAP.Ir;
  const A = Number(values.A);
  const B = Number(values.B);
  const C = Number(values.C);
  if (!Number.isFinite(A) || !Number.isFinite(B) || !Number.isFinite(C)) return null;
  const principal = {};
  principal[from[0]] = A;
  principal[from[1]] = B;
  principal[from[2]] = C;
  return {
    A: principal[to[0]],
    B: principal[to[1]],
    C: principal[to[2]],
  };
}

function rotateRepresentationTripletObject(values, inputRepresentation, outputRepresentation) {
  if (!values || typeof values !== "object") return values;
  const from = REP_MAP[normalizeFullWatsonRepresentationToken(inputRepresentation)] || REP_MAP.Ir;
  const to = REP_MAP[normalizeFullWatsonRepresentationToken(outputRepresentation)] || REP_MAP.Ir;
  const principal = {};
  principal[from[0]] = values.a ?? "";
  principal[from[1]] = values.b ?? "";
  principal[from[2]] = values.c ?? "";
  return {
    a: principal[to[0]] ?? "",
    b: principal[to[1]] ?? "",
    c: principal[to[2]] ?? "",
  };
}

function transposeMatrix(matrix) {
  if (!Array.isArray(matrix) || !matrix.length) return [];
  return matrix[0].map((_, column) => matrix.map((row) => Number(row[column])));
}

function multiplyMatrices(left, right) {
  const rows = Array.isArray(left) ? left.length : 0;
  const inner = rows ? left[0].length : 0;
  const cols = Array.isArray(right) && right.length ? right[0].length : 0;
  const out = Array.from({ length: rows }, () => Array(cols).fill(0));
  for (let i = 0; i < rows; i += 1) {
    for (let k = 0; k < inner; k += 1) {
      const leftValue = Number(left[i][k]);
      for (let j = 0; j < cols; j += 1) {
        out[i][j] += leftValue * Number(right[k][j]);
      }
    }
  }
  return out;
}

function multiplyMatrixVector(matrix, vector) {
  return matrix.map((row) =>
    row.reduce((sum, value, index) => sum + Number(value) * Number(vector[index]), 0)
  );
}

function invertSquareMatrix(matrix) {
  const n = Array.isArray(matrix) ? matrix.length : 0;
  if (!n || !Array.isArray(matrix[0]) || matrix[0].length !== n) {
    throw new Error("Matrix must be square.");
  }
  const augmented = matrix.map((row, i) => {
    const base = row.map((value) => Number(value));
    const identity = Array.from({ length: n }, (_, j) => (i === j ? 1 : 0));
    return base.concat(identity);
  });
  for (let col = 0; col < n; col += 1) {
    let pivotRow = col;
    let pivotAbs = Math.abs(augmented[pivotRow][col]);
    for (let row = col + 1; row < n; row += 1) {
      const candidateAbs = Math.abs(augmented[row][col]);
      if (candidateAbs > pivotAbs) {
        pivotRow = row;
        pivotAbs = candidateAbs;
      }
    }
    if (!(pivotAbs > 1e-14)) {
      throw new Error("Matrix is singular.");
    }
    if (pivotRow !== col) {
      const temp = augmented[col];
      augmented[col] = augmented[pivotRow];
      augmented[pivotRow] = temp;
    }
    const pivot = augmented[col][col];
    for (let j = 0; j < 2 * n; j += 1) {
      augmented[col][j] /= pivot;
    }
    for (let row = 0; row < n; row += 1) {
      if (row === col) continue;
      const factor = augmented[row][col];
      if (Math.abs(factor) < 1e-18) continue;
      for (let j = 0; j < 2 * n; j += 1) {
        augmented[row][j] -= factor * augmented[col][j];
      }
    }
  }
  return augmented.map((row) => row.slice(n));
}

function pseudoInverseFullRowRank(matrix) {
  const transpose = transposeMatrix(matrix);
  return multiplyMatrices(transpose, invertSquareMatrix(multiplyMatrices(matrix, transpose)));
}

function pseudoInverseFullColumnRank(matrix) {
  const transpose = transposeMatrix(matrix);
  return multiplyMatrices(invertSquareMatrix(multiplyMatrices(transpose, matrix)), transpose);
}

function parseFortranVectorFromPayload(payload, keys) {
  if (!payload || typeof payload !== "object") return null;
  const out = [];
  for (let i = 0; i < keys.length; i += 1) {
    const parsed = parseSciNumber(payload[keys[i]]);
    if (!Number.isFinite(parsed)) return null;
    out.push(parsed);
  }
  return out;
}

function formatFortranVectorToPayload(values, keys) {
  const out = {};
  for (let i = 0; i < keys.length; i += 1) {
    out[keys[i]] = toFortranScientific(Number(values[i]), 10);
  }
  return out;
}

function computeTauSigmaFromRotationalConstantsMHz(values) {
  if (!values || !Number.isFinite(values.a) || !Number.isFinite(values.b) || !Number.isFinite(values.c)) {
    return null;
  }
  const denom = values.b - values.c;
  if (Math.abs(denom) < 1e-16) return 0;
  return (2 * values.a - values.b - values.c) / denom;
}

function computeTauSigmaFromMoments(values) {
  if (!values || !Number.isFinite(values.a) || !Number.isFinite(values.b) || !Number.isFinite(values.c)) {
    return null;
  }
  if (values.a <= 0 || values.b <= 0 || values.c <= 0) return null;
  const invA = 1 / values.a;
  const invB = 1 / values.b;
  const invC = 1 / values.c;
  const denom = invB - invC;
  if (Math.abs(denom) < 1e-16) return 0;
  return (2 * invA - invB - invC) / denom;
}

function buildQuarticFromTauComponents({ tau, sigma = 0 }) {
  const taaaa = Number(tau.taaaa);
  const tbbbb = Number(tau.tbbbb);
  const tcccc = Number(tau.tcccc);
  const taabb = Number(tau.taabb);
  const tbbcc = Number(tau.tbbcc);
  const tccaa = Number(tau.tccaa);
  const tabab = Number(tau.tabab);
  const tbcbc = Number(tau.tbcbc);
  const tcaca = Number(tau.tcaca);

  if (
    !Number.isFinite(taaaa) ||
    !Number.isFinite(tbbbb) ||
    !Number.isFinite(tcccc) ||
    !Number.isFinite(taabb) ||
    !Number.isFinite(tbbcc) ||
    !Number.isFinite(tccaa) ||
    !Number.isFinite(tabab) ||
    !Number.isFinite(tbcbc) ||
    !Number.isFinite(tcaca)
  ) {
    throw new Error("Invalid tau component set: all 9 independent components are required.");
  }

  const DJ = -(1 / 32) * (3 * tbbbb + 3 * tcccc + 2 * (tbbcc + tbcbc));
  const DK = DJ - 0.25 * (taaaa - (taabb + tabab) - (tccaa + tcaca));
  const DJK = -DJ - DK + 0.25 * taaaa;
  const R5 = -(1 / 32) * (tbbbb - (tcccc - 2 * taabb) + 2 * tabab + 2 * (tccaa + tcaca));
  const R6 = (1 / 64) * (tbbbb + tcccc - 2 * (tbbcc + tbcbc));
  const dJ = -(1 / 16) * (tbbbb - tcccc);
  const s = Number.isFinite(sigma) ? sigma : 0;

  const asymMHz = {
    deltaJ: DJ - 2 * R6,
    deltaJK: DJK + 12 * R6,
    deltaK: DK - 10 * R6,
    smallDeltaJ: dJ,
    smallDeltaK: 2 * R5 - 4 * s * R6,
  };
  const symMHz = {
    DJ,
    DJK,
    DK,
    d1: -dJ,
    d2: R6,
  };
  return {
    sigmaUsed: s,
    asymMHz,
    symMHz,
    asymFormatted: {
      deltaJ: toFortranScientific(asymMHz.deltaJ, 10),
      deltaJK: toFortranScientific(asymMHz.deltaJK, 10),
      deltaK: toFortranScientific(asymMHz.deltaK, 10),
      smallDeltaJ: toFortranScientific(asymMHz.smallDeltaJ, 10),
      smallDeltaK: toFortranScientific(asymMHz.smallDeltaK, 10),
    },
    symFormatted: {
      DJ: toFortranScientific(symMHz.DJ, 10),
      DJK: toFortranScientific(symMHz.DJK, 10),
      DK: toFortranScientific(symMHz.DK, 10),
      d1: toFortranScientific(symMHz.d1, 10),
      d2: toFortranScientific(symMHz.d2, 10),
    },
  };
}

function buildTauDerivedPayloads(bundle, representation) {
  if (!bundle || !bundle.tau) {
    throw new Error("Invalid tau bundle payload.");
  }
  const inputRep = bundle.inputRepresentation || "Ir";
  const rotation = rotateTauComponentsForRepresentation(bundle.tau, inputRep, representation);
  const rotatedRotationalConstantsMHz = rotateAxisTripletByRepresentation(
    bundle.rotationalConstantsMHz,
    inputRep,
    representation
  );
  const rotatedMoments = rotateAxisTripletByRepresentation(
    bundle.moments,
    inputRep,
    representation
  );
  let sigma =
    inputRep === representation && Number.isFinite(bundle.sigma)
      ? bundle.sigma
      : computeTauSigmaFromRotationalConstantsMHz(rotatedRotationalConstantsMHz);
  let sigmaSource =
    inputRep === representation && Number.isFinite(bundle.sigma)
      ? "sigma"
      : Number.isFinite(sigma)
        ? "ABC"
        : "none";
  if (!Number.isFinite(sigma)) {
    sigma = computeTauSigmaFromMoments(rotatedMoments);
    sigmaSource = Number.isFinite(sigma) ? "moments" : "none";
  }
  if (!Number.isFinite(sigma) && Number.isFinite(bundle.sigma)) {
    sigma = bundle.sigma;
    sigmaSource = "sigma";
  }
  if (!Number.isFinite(sigma)) {
    sigma = 0;
    sigmaSource = "default-0";
  }

  const quartic = buildQuarticFromTauComponents({
    tau: rotation.tau,
    sigma,
  });

  const inputRightHanded = normalizeTauRepresentation(inputRep);
  const outputRightHanded = normalizeTauRepresentation(representation);
  const sigmaLabel =
    sigmaSource === "ABC"
      ? "A/B/C"
      : sigmaSource === "moments"
        ? "Ia/Ib/Ic"
        : sigmaSource === "sigma"
          ? "sigma input"
          : "default s=0";
  const contextLabel = `Tau quartic rotated (${inputRep}→${representation}; qcent ${inputRightHanded}→${outputRightHanded}; s from ${sigmaLabel})`;
  const finalContextLabel =
    inputRep === representation
      ? `Tau quartic input (${inputRep}; s from ${sigmaLabel})`
      : contextLabel;

  return {
    quartic: {
      contextLabel: finalContextLabel,
      sourceIndex: Date.now(),
      asym: { ...quartic.asymFormatted },
      sym: { ...quartic.symFormatted },
      tau: { ...rotation.tau },
      sigma: quartic.sigmaUsed,
      inputRepresentation: inputRep,
    },
  };
}

function buildSexticCanonicalIndexLookup() {
  const lookup = {};
  for (let i = 0; i < SEXTIC_INERTIAL_KEYS.length; i += 1) {
    const canonical = SEXTIC_INERTIAL_KEYS[i].slice(1);
    const permutations = expandSexticEquivalentIndices(canonical);
    for (let j = 0; j < permutations.length; j += 1) {
      lookup[permutations[j]] = canonical;
    }
  }
  return lookup;
}

function expandSexticEquivalentIndices(indices) {
  const seed = normalizeSexticInertialIndexToken(indices);
  if (!seed) return [];
  const out = new Set();
  const arr = seed.split("");
  const permute = (start) => {
    if (start >= arr.length) {
      out.add(arr.join(""));
      return;
    }
    const seen = new Set();
    for (let i = start; i < arr.length; i += 1) {
      const value = arr[i];
      if (seen.has(value)) continue;
      seen.add(value);
      const tmp = arr[start];
      arr[start] = arr[i];
      arr[i] = tmp;
      permute(start + 1);
      arr[i] = arr[start];
      arr[start] = tmp;
    }
  };
  permute(0);
  return Array.from(out);
}

function normalizeSexticInertialIndexToken(value) {
  const compact = String(value || "")
    .toLowerCase()
    .replace(/[^abc]/g, "");
  if (compact.length !== 3) return null;
  return compact;
}

function canonicalSexticInertialKeyFromToken(value) {
  const normalized = normalizeSexticInertialIndexToken(value);
  if (!normalized) return null;
  const canonical = SEXTIC_CANONICAL_INDEX_LOOKUP[normalized] || null;
  return canonical ? `e${canonical}` : null;
}

function parseSexticPayloadText(text) {
  const source = String(text || "").replace(/\r/g, "\n");
  const inertial = parseSexticInertialComponentsFromText(source);
  const sigma = parseTauSigmaFromText(source);
  const missing = SEXTIC_INERTIAL_KEYS.filter((key) => !Number.isFinite(inertial[key]));
  if (missing.length) {
    throw new Error(
      `Sextic file is missing required inertial components: ${missing.join(", ")}.`
    );
  }
  return {
    inputRepresentation: parseTauRepresentationFromText(source) || "Ir",
    inertial,
    rotationalConstantsMHz: parseTauRotationalConstantsFromText(source),
    moments: parseTauMomentsFromText(source),
    sigma,
    sigmaSource: Number.isFinite(sigma) ? "sigma" : "none",
  };
}

function parseSexticInertialComponentsFromText(text) {
  const source = String(text || "");
  const out = {};
  for (let i = 0; i < SEXTIC_INERTIAL_KEYS.length; i += 1) {
    out[SEXTIC_INERTIAL_KEYS[i]] = null;
  }
  const patterns = [
    new RegExp(
      `\\b(?:eta|e)\\s*[_\\-\\s]*(?:\\(|\\[)?([abc]{3})(?:\\)|\\])?\\s*(?:=|:)?\\s*(${NUMBER_TOKEN_PATTERN})`,
      "gi"
    ),
    new RegExp(`\\b([abc]{3})\\b\\s*(?:=|:)?\\s*(${NUMBER_TOKEN_PATTERN})`, "gi"),
  ];
  for (let p = 0; p < patterns.length; p += 1) {
    const regex = patterns[p];
    regex.lastIndex = 0;
    let match;
    while ((match = regex.exec(source)) !== null) {
      const key = canonicalSexticInertialKeyFromToken(match[1]);
      if (!key) continue;
      const parsed = parseSciNumber(match[2]);
      if (!Number.isFinite(parsed)) continue;
      out[key] = parsed;
    }
  }
  return out;
}

function getSexticInertialValueByIndices(inertial, indexToken) {
  const key = canonicalSexticInertialKeyFromToken(indexToken);
  if (!key) return null;
  const value = inertial ? inertial[key] : null;
  return Number.isFinite(value) ? value : null;
}

function rotateSexticInertialComponents(inertial, inputRepresentation, outputRepresentation) {
  const remap = buildTauAxisRemap("Ir", inverseFullWatsonRepresentationToken(outputRepresentation));
  const rotated = {};
  for (let i = 0; i < SEXTIC_INERTIAL_KEYS.length; i += 1) {
    const key = SEXTIC_INERTIAL_KEYS[i];
    const outputIndex = key.slice(1);
    const mappedInputIndex =
      remap.axisOutToIn[outputIndex[0]] +
      remap.axisOutToIn[outputIndex[1]] +
      remap.axisOutToIn[outputIndex[2]];
    const value = getSexticInertialValueByIndices(inertial, mappedInputIndex);
    if (!Number.isFinite(value)) {
      throw new Error(
        `Unable to rotate sextic inertial component "${key}" to representation ${outputRepresentation}.`
      );
    }
    rotated[key] = value;
  }
  return {
    inertial: rotated,
    remap,
  };
}

function getRotationalConstantForSextic(values, key) {
  if (!values || typeof values !== "object") return null;
  const direct = Number(values[key]);
  if (Number.isFinite(direct)) return direct;
  const lower = Number(values[String(key).toLowerCase()]);
  return Number.isFinite(lower) ? lower : null;
}

function hasSexticProjectionInputs(quarticT, rotationalConstantsMHz) {
  const tKeys = ["T220", "T040", "T202", "T022", "T004"];
  const hasT =
    quarticT &&
    tKeys.every((key) => Number.isFinite(Number(quarticT[key])));
  const hasABC =
    Number.isFinite(getRotationalConstantForSextic(rotationalConstantsMHz, "A")) &&
    Number.isFinite(getRotationalConstantForSextic(rotationalConstantsMHz, "B")) &&
    Number.isFinite(getRotationalConstantForSextic(rotationalConstantsMHz, "C"));
  return Boolean(hasT && hasABC);
}

function buildSexticFromInertialComponents({
  inertial,
  sigma = 0,
  quarticT = null,
  rotationalConstantsMHz = null,
}) {
  const eaaa = Number(inertial.eaaa);
  const eaab = Number(inertial.eaab);
  const eaac = Number(inertial.eaac);
  const eabb = Number(inertial.eabb);
  const eabc = Number(inertial.eabc);
  const eacc = Number(inertial.eacc);
  const ebbb = Number(inertial.ebbb);
  const ebbc = Number(inertial.ebbc);
  const ebcc = Number(inertial.ebcc);
  const eccc = Number(inertial.eccc);
  if (
    !Number.isFinite(eaaa) ||
    !Number.isFinite(eaab) ||
    !Number.isFinite(eaac) ||
    !Number.isFinite(eabb) ||
    !Number.isFinite(eabc) ||
    !Number.isFinite(eacc) ||
    !Number.isFinite(ebbb) ||
    !Number.isFinite(ebbc) ||
    !Number.isFinite(ebcc) ||
    !Number.isFinite(eccc)
  ) {
    throw new Error("Invalid sextic inertial set: all 10 independent components are required.");
  }

  const s = Number.isFinite(sigma) ? sigma : 0;
  const sigma1 = Math.abs(s) > 1.0e-16 ? 1 / s : 0;

  if (hasSexticProjectionInputs(quarticT, rotationalConstantsMHz)) {
    const A = getRotationalConstantForSextic(rotationalConstantsMHz, "A");
    const B = getRotationalConstantForSextic(rotationalConstantsMHz, "B");
    const C = getRotationalConstantForSextic(rotationalConstantsMHz, "C");
    const T220 = Number(quarticT.T220);
    const T040 = Number(quarticT.T040);
    const T202 = Number(quarticT.T202);
    const T022 = Number(quarticT.T022);
    const T004 = Number(quarticT.T004);

    const Phi600 = (5 / 16) * (eaaa + ebbb) + (1 / 8) * (eaab + eabb);
    const Phi420 = (3 / 4) * (eaac + ebbc) + (1 / 4) * eabc - 3 * Phi600;
    const Phi240 = eacc + ebcc - 2 * Phi420 - 3 * Phi600;
    const Phi060 = eccc - Phi240 - Phi420 - Phi600;
    const Phi402 = (15 / 64) * (eaaa - ebbb) + (1 / 32) * (eaab - eabb);
    const Phi222 = 0.5 * (eaac - ebbc) - 2 * Phi402;
    const Phi042 = 0.5 * (eacc - ebcc) - Phi222 - Phi402;
    const Phi204 = (3 / 32) * (eaaa + ebbb) - (1 / 16) * (eaab + eabb);
    const Phi024 = (1 / 8) * (eaac + ebbc - eabc) - Phi204;
    const Phi006 = (1 / 64) * (eaaa - ebbb) - (1 / 32) * (eaab - eabb);

    const B200 = 0.5 * (A + B) - 4 * T004;
    const B020 = C - B200 + 6 * T004;
    const B002 = 0.25 * (A - B);
    if (Math.abs(B002) < 1.0e-16 || Math.abs(B020) < 1.0e-16) {
      throw new Error("Sextic projection is singular for this asymmetric-top representation.");
    }

    const PHIN = Phi600 + 2 * Phi204;
    const PHINK =
      Phi420 -
      12 * Phi204 +
      2 * Phi024 +
      16 * s * Phi006 +
      (8 * T022 * T004) / B002;
    const PHIKN = Phi240 + (10 / 3) * Phi420 - 30 * Phi204 - (10 / 3) * PHINK;
    const PHIK = Phi060 - (7 / 3) * Phi420 + 28 * Phi204 + (7 / 3) * PHINK;
    const phiN = Phi402 + Phi006;
    const phiK =
      Phi042 +
      (4 / 3) * s * Phi024 +
      ((32 / 3) * s * s + 9) * Phi006 +
      (4 * (T040 + (s * T022) / 3 - 2 * (s * s - 2) * T004) * T004) / B002;
    const phiNK =
      Phi222 +
      4 * s * Phi204 -
      10 * Phi006 +
      (2 * (T220 - 2 * s * T202 - 4 * T004) * T004) / B002;

    const rho = 0.25 * T022 * sigma1;
    const muDen = 2 * s * s + 27 / 16;
    const mu =
      (s * Phi042 - (9 / 8) * Phi024) / muDen +
      ((-2 * s * T040 + (s * s + 3) * T022 - 5 * s * T004) * T022) /
        (B020 * muDen);
    const nu = (3 / 16) * mu * sigma1 + (1 / 8) * Phi024 * sigma1 + (T004 * T022) / B020;
    const lambda =
      5 * nu * sigma1 +
      0.5 * Phi222 * sigma1 +
      ((-0.5 * T220 * sigma1 + T202 - T022 * sigma1 * sigma1 - 2 * T004 * sigma1) *
        T022) /
        B020;

    const asymMHz = {
      PHIN,
      PHINK,
      PHIKN,
      PHIK,
      phiN,
      phiNK,
      phiK,
    };
    const symMHz = {
      HN: Phi600 - lambda,
      HNK: Phi420 + 6 * lambda - 3 * mu,
      HKN: Phi240 - 5 * lambda + 10 * mu,
      HK: Phi060 - 7 * mu,
      h1: Phi402 - nu,
      h2: Phi204 + 0.5 * lambda,
      h3: Phi006 + nu,
    };

    return {
      sigmaUsed: s,
      sigma1,
      intermediates: {
        Phi600,
        Phi420,
        Phi240,
        Phi060,
        Phi402,
        Phi222,
        Phi042,
        Phi204,
        Phi024,
        Phi006,
        B200,
        B020,
        B002,
        rho,
        mu,
        nu,
        lambda,
      },
      asymMHz,
      symMHz,
      asymFormatted: {
        PHIN: toFortranScientific(asymMHz.PHIN, 10),
        PHINK: toFortranScientific(asymMHz.PHINK, 10),
        PHIKN: toFortranScientific(asymMHz.PHIKN, 10),
        PHIK: toFortranScientific(asymMHz.PHIK, 10),
        phiN: toFortranScientific(asymMHz.phiN, 10),
        phiNK: toFortranScientific(asymMHz.phiNK, 10),
        phiK: toFortranScientific(asymMHz.phiK, 10),
      },
      symFormatted: {
        HN: toFortranScientific(symMHz.HN, 10),
        HNK: toFortranScientific(symMHz.HNK, 10),
        HKN: toFortranScientific(symMHz.HKN, 10),
        HK: toFortranScientific(symMHz.HK, 10),
        h1: toFortranScientific(symMHz.h1, 10),
        h2: toFortranScientific(symMHz.h2, 10),
        h3: toFortranScientific(symMHz.h3, 10),
      },
    };
  }

  const omegaJ = (1 / 64) * (5 * (ebbb + eccc) - 2 * (ebbc + ebcc));
  const omegaK = (1 / 16) * (eaaa - eabb - eacc);
  const omegaJK = -(1 / 32) * (eaaa - 3 * (eabb + eacc) + 2 * eabc);
  const omegaKJ = (1 / 32) * (ebbb + eccc - 2 * eabc);
  const etaJ = (1 / 32) * (ebbb - eccc);
  const etaJK = -(1 / 16) * (ebbc - ebcc);
  const etaK = (1 / 16) * (eabb - eacc);

  const asymMHz = {
    PHIN: omegaJ,
    PHINK: omegaJK,
    PHIKN: omegaKJ,
    PHIK: omegaK,
    phiN: etaJ,
    phiNK: etaJK,
    phiK: etaK,
  };

  const omegaCoupling = omegaJK + omegaKJ;
  const symMHz = {
    HN: omegaJ - omegaCoupling / 3,
    HNK: omegaJK - omegaCoupling / 3,
    HKN: omegaKJ - omegaCoupling / 3,
    HK: omegaK + (2 * omegaCoupling) / 3,
    h1: etaJ,
    h2: etaK + s * etaJK,
    h3: etaK - s * etaJK,
  };

  return {
    sigmaUsed: s,
    intermediates: {
      omegaJ,
      omegaJK,
      omegaKJ,
      omegaK,
      etaJ,
      etaJK,
      etaK,
    },
    asymMHz,
    symMHz,
    asymFormatted: {
      PHIN: toFortranScientific(asymMHz.PHIN, 10),
      PHINK: toFortranScientific(asymMHz.PHINK, 10),
      PHIKN: toFortranScientific(asymMHz.PHIKN, 10),
      PHIK: toFortranScientific(asymMHz.PHIK, 10),
      phiN: toFortranScientific(asymMHz.phiN, 10),
      phiNK: toFortranScientific(asymMHz.phiNK, 10),
      phiK: toFortranScientific(asymMHz.phiK, 10),
    },
    symFormatted: {
      HN: toFortranScientific(symMHz.HN, 10),
      HNK: toFortranScientific(symMHz.HNK, 10),
      HKN: toFortranScientific(symMHz.HKN, 10),
      HK: toFortranScientific(symMHz.HK, 10),
      h1: toFortranScientific(symMHz.h1, 10),
      h2: toFortranScientific(symMHz.h2, 10),
      h3: toFortranScientific(symMHz.h3, 10),
    },
  };
}

function buildSexticDerivedPayloads(bundle, representation) {
  if (!bundle || !bundle.inertial) {
    throw new Error("Invalid sextic bundle payload.");
  }
  const inputRep = bundle.inputRepresentation || "Ir";
  const rotation = rotateSexticInertialComponents(bundle.inertial, inputRep, representation);
  const rotatedRotationalConstantsMHz = rotateAxisTripletByRepresentation(
    bundle.rotationalConstantsMHz,
    inputRep,
    representation
  );
  const rotatedMoments = rotateAxisTripletByRepresentation(
    bundle.moments,
    inputRep,
    representation
  );

  let sigma =
    inputRep === representation && Number.isFinite(bundle.sigma)
      ? bundle.sigma
      : computeTauSigmaFromRotationalConstantsMHz(rotatedRotationalConstantsMHz);
  let sigmaSource =
    inputRep === representation && Number.isFinite(bundle.sigma)
      ? "sigma"
      : Number.isFinite(sigma)
        ? "ABC"
        : "none";
  if (!Number.isFinite(sigma)) {
    sigma = computeTauSigmaFromMoments(rotatedMoments);
    sigmaSource = Number.isFinite(sigma) ? "moments" : "none";
  }
  if (!Number.isFinite(sigma) && Number.isFinite(bundle.sigma)) {
    sigma = bundle.sigma;
    sigmaSource = "sigma";
  }
  if (!Number.isFinite(sigma)) {
    sigma = 0;
    sigmaSource = "default-0";
  }

  const sextic = buildSexticFromInertialComponents({
    inertial: rotation.inertial,
    sigma,
  });

  const inputRightHanded = normalizeTauRepresentation(inputRep);
  const outputRightHanded = normalizeTauRepresentation(representation);
  const sigmaLabel =
    sigmaSource === "ABC"
      ? "A/B/C"
      : sigmaSource === "moments"
        ? "Ia/Ib/Ic"
        : sigmaSource === "sigma"
          ? "sigma input"
          : "default s=0";
  const contextLabel = `Sextic inertial rotated (${inputRep}→${representation}; sextic ${inputRightHanded}→${outputRightHanded}; s from ${sigmaLabel})`;
  const finalContextLabel =
    inputRep === representation
      ? `Sextic inertial input (${inputRep}; s from ${sigmaLabel})`
      : contextLabel;

  return {
    sextic: {
      contextLabel: finalContextLabel,
      sourceIndex: Date.now(),
      asymSextic: { ...sextic.asymFormatted },
      symSextic: { ...sextic.symFormatted },
      inertial: { ...rotation.inertial },
      sigma: sextic.sigmaUsed,
      intermediates: { ...sextic.intermediates },
      inputRepresentation: inputRep,
    },
  };
}

function computeHarmonicQuarticFromHessian({
  hessianCartesian,
  geometryAtoms,
  representation = "Ir",
}) {
  const mathLib = getMathLibrary();
  if (!mathLib) {
    throw new Error("math.js is required to diagonalize the Hessian.");
  }
  const atoms = normalizeAtomsForQuartic(geometryAtoms);
  if (!atoms.length) {
    throw new Error("No valid atoms available for harmonic quartic computation.");
  }
  const expectedSize = atoms.length * 3;
  if (!isSquareMatrix(hessianCartesian, expectedSize)) {
    throw new Error(
      `Hessian size mismatch: expected ${expectedSize} x ${expectedSize} for ${atoms.length} atoms.`
    );
  }

  const center = computeCenterOfMass(atoms);
  const inertia = computeInertiaTensor(atoms, center);
  const principalFrame = diagonalizeSymmetricTensor(inertia);
  const principalAtoms = atoms.map((atom) => ({
    label: atom.label,
    z: atom.z,
    coords: rotateToBasis(subtractVec(atom.coords, center), principalFrame.basis),
  }));
  const orientedAtoms = applyRepresentation(principalAtoms, representation);
  const orientedCenter = computeCenterOfMass(orientedAtoms);
  const centeredOrientedAtoms = orientedAtoms.map((atom) => ({
    label: atom.label,
    z: atom.z,
    coords: subtractVec(atom.coords, orientedCenter),
  }));

  const repRotation = buildRepresentationRotationMatrix(principalFrame.basis, representation);
  const rotatedHessian = transformCartesianHessian(hessianCartesian, repRotation);

  const massesAmu = centeredOrientedAtoms.map((atom) => {
    const mass = atomicMass[atom.z];
    if (!Number.isFinite(mass) || mass <= 0) {
      throw new Error(`Atomic mass not available for Z=${atom.z}.`);
    }
    return mass;
  });
  const massesAu = massesAmu.map((mass) => mass * AMU_TO_AU);
  const coordsBohr = centeredOrientedAtoms.map((atom) => [
    atom.coords[0] * ANGSTROM_TO_BOHR,
    atom.coords[1] * ANGSTROM_TO_BOHR,
    atom.coords[2] * ANGSTROM_TO_BOHR,
  ]);
  const momentsAmuAng2 = computeInertiaMomentsFromCoordinates(centeredOrientedAtoms, massesAmu);
  const momentsAu = computeInertiaMomentsBohr(coordsBohr, massesAu);

  const massWeighted = buildMassWeightedHessian(rotatedHessian, massesAu);
  const eig = mathLib.eigs(mathLib.matrix(massWeighted));
  const eigenValues = toArray(eig.values);
  const eigenVectors = toArray(eig.vectors);
  if (
    !Array.isArray(eigenValues) ||
    !Array.isArray(eigenVectors) ||
    eigenVectors.length !== expectedSize
  ) {
    throw new Error("Mass-weighted Hessian diagonalization failed.");
  }

  const modePairs = [];
  for (let col = 0; col < expectedSize; col += 1) {
    const vector = new Array(expectedSize);
    for (let row = 0; row < expectedSize; row += 1) {
      vector[row] = asReal(eigenVectors[row][col]);
    }
    modePairs.push({
      lambda: asReal(eigenValues[col]),
      vector,
    });
  }

  const zeroModes = isLikelyLinearFromMoments(momentsAu) ? 5 : 6;
  const expectedModeCount = Math.max(0, expectedSize - zeroModes);
  const orderedByAbs = modePairs.slice().sort((a, b) => Math.abs(a.lambda) - Math.abs(b.lambda));
  let vibModes = orderedByAbs.slice(zeroModes).filter((mode) => mode.lambda > 1e-12);
  if (vibModes.length > expectedModeCount && expectedModeCount > 0) {
    vibModes = vibModes
      .slice()
      .sort((a, b) => b.lambda - a.lambda)
      .slice(0, expectedModeCount);
  }
  if (vibModes.length < expectedModeCount) {
    throw new Error(
      `Not enough vibrational modes after removing translational/rotational modes (${vibModes.length}/${expectedModeCount}).`
    );
  }
  const coriolisPairs = HARMONIC_CORIOLIS_ENABLED
    ? buildHarmonicCoriolisPairs(vibModes, atoms.length)
    : [];

  const modeDerivatives = [];
  for (let modeIndex = 0; modeIndex < vibModes.length; modeIndex += 1) {
    const mode = vibModes[modeIndex];
    const deriv = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    for (let atomIndex = 0; atomIndex < massesAu.length; atomIndex += 1) {
      const massAu = massesAu[atomIndex];
      const rootMass = Math.sqrt(massAu);
      const base = atomIndex * 3;
      const r = coordsBohr[atomIndex];
      const u = [
        mode.vector[base + 0] / rootMass,
        mode.vector[base + 1] / rootMass,
        mode.vector[base + 2] / rootMass,
      ];
      const dotRu = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
      for (let alpha = 0; alpha < 3; alpha += 1) {
        for (let beta = 0; beta < 3; beta += 1) {
          const kronecker = alpha === beta ? 1 : 0;
          deriv[alpha][beta] +=
            massAu *
            (2 * dotRu * kronecker - (r[alpha] * u[beta] + u[alpha] * r[beta]));
        }
      }
    }
    modeDerivatives.push({
      lambda: mode.lambda,
      deriv,
    });
  }

  for (let axis = 0; axis < 3; axis += 1) {
    if (!Number.isFinite(momentsAu[axis]) || momentsAu[axis] <= 0) {
      throw new Error("Invalid inertia moment encountered while computing quartic constants.");
    }
  }

  const tPrefactorMHz =
    -(219474.6313705 * CM1_TO_MHZ) / (2 * TWO_PI * C_AU);
  const coriolisContribution = (a, b, c, d) => {
    if (!HARMONIC_CORIOLIS_ENABLED || a !== b || c !== d) {
      return 0;
    }
    let sum = 0;
    for (let i = 0; i < coriolisPairs.length; i += 1) {
      const pair = coriolisPairs[i];
      if (pair.lambdaSum <= HARMONIC_CORIOLIS_MIN_DENOM) continue;
      sum += (pair.zeta[a] * pair.zeta[c]) / pair.lambdaSum;
    }
    return HARMONIC_CORIOLIS_SCALE * sum;
  };
  const tComponent = (a, b, c, d) => {
    let sumWatson = 0;
    for (let i = 0; i < modeDerivatives.length; i += 1) {
      const item = modeDerivatives[i];
      sumWatson += (item.deriv[a][b] * item.deriv[c][d]) / item.lambda;
    }
    const sumCoriolis = coriolisContribution(a, b, c, d);
    const totalSum = sumWatson + sumCoriolis;
    const inertiaProduct = momentsAu[a] * momentsAu[b] * momentsAu[c] * momentsAu[d];
    return (tPrefactorMHz * totalSum) / inertiaProduct;
  };

  const taaaa = tComponent(0, 0, 0, 0);
  const tbbbb = tComponent(1, 1, 1, 1);
  const tcccc = tComponent(2, 2, 2, 2);
  const taabb = tComponent(0, 0, 1, 1);
  const tabab = tComponent(0, 1, 0, 1);
  const tccaa = tComponent(2, 2, 0, 0);
  const tcaca = tComponent(2, 0, 2, 0);
  const tbbcc = tComponent(1, 1, 2, 2);
  const tbcbc = tComponent(1, 2, 1, 2);

  const DJ = -(1 / 32) * (3 * tbbbb + 3 * tcccc + 2 * (tbbcc + tbcbc));
  const DK = DJ - 0.25 * (taaaa - (taabb + tabab) - (tccaa + tcaca));
  const DJK = -DJ - DK + 0.25 * taaaa;
  const R5 = -(1 / 32) * (tbbbb - (tcccc - 2 * taabb) + 2 * tabab + 2 * (tccaa + tcaca));
  const R6 = (1 / 64) * (tbbbb + tcccc - 2 * (tbbcc + tbcbc));
  const dJ = -(1 / 16) * (tbbbb - tcccc);

  const invA = 1 / momentsAu[0];
  const invB = 1 / momentsAu[1];
  const invC = 1 / momentsAu[2];
  const denom = invB - invC;
  const s = Math.abs(denom) < 1e-16 ? 0 : (2 * invA - invB - invC) / denom;

  const asymMHz = {
    deltaJ: DJ - 2 * R6,
    deltaJK: DJK + 12 * R6,
    deltaK: DK - 10 * R6,
    smallDeltaJ: dJ,
    smallDeltaK: 2 * R5 - 4 * s * R6,
  };
  const symMHz = {
    DJ,
    DJK,
    DK,
    d1: -dJ,
    d2: R6,
  };
  const asymFormatted = {
    deltaJ: toFortranScientific(asymMHz.deltaJ, 10),
    deltaJK: toFortranScientific(asymMHz.deltaJK, 10),
    deltaK: toFortranScientific(asymMHz.deltaK, 10),
    smallDeltaJ: toFortranScientific(asymMHz.smallDeltaJ, 10),
    smallDeltaK: toFortranScientific(asymMHz.smallDeltaK, 10),
  };
  const symFormatted = {
    DJ: toFortranScientific(symMHz.DJ, 10),
    DJK: toFortranScientific(symMHz.DJK, 10),
    DK: toFortranScientific(symMHz.DK, 10),
    d1: toFortranScientific(symMHz.d1, 10),
    d2: toFortranScientific(symMHz.d2, 10),
  };

  return {
    representation,
    modeCount: modeDerivatives.length,
    expectedModeCount,
    momentsAmuAng2,
    momentsAu,
    coriolisModel: {
      enabled: HARMONIC_CORIOLIS_ENABLED,
      scale: HARMONIC_CORIOLIS_SCALE,
      pairCount: coriolisPairs.length,
    },
    asymMHz,
    symMHz,
    asymFormatted,
    symFormatted,
  };
}

function buildHarmonicCoriolisPairs(vibModes, atomCount) {
  if (!Array.isArray(vibModes) || !vibModes.length || !Number.isFinite(atomCount) || atomCount <= 0) {
    return [];
  }
  const pairs = [];
  for (let i = 0; i < vibModes.length; i += 1) {
    const modeI = vibModes[i];
    for (let j = i + 1; j < vibModes.length; j += 1) {
      const modeJ = vibModes[j];
      const zeta = [0, 0, 0];
      for (let atom = 0; atom < atomCount; atom += 1) {
        const base = atom * 3;
        const liX = modeI.vector[base + 0];
        const liY = modeI.vector[base + 1];
        const liZ = modeI.vector[base + 2];
        const ljX = modeJ.vector[base + 0];
        const ljY = modeJ.vector[base + 1];
        const ljZ = modeJ.vector[base + 2];
        zeta[0] += liY * ljZ - liZ * ljY;
        zeta[1] += liZ * ljX - liX * ljZ;
        zeta[2] += liX * ljY - liY * ljX;
      }
      const lambdaSum = modeI.lambda + modeJ.lambda;
      pairs.push({
        modeI: i,
        modeJ: j,
        zeta,
        lambdaSum,
      });
    }
  }
  return pairs;
}

function normalizeAtomsForQuartic(atoms) {
  if (!Array.isArray(atoms)) return [];
  const out = [];
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    if (!atom || !Array.isArray(atom.coords) || atom.coords.length < 3) continue;
    let z = Number.isFinite(atom.z) ? Number(atom.z) : null;
    if (!Number.isFinite(z) || z <= 0) {
      z = lookupAtomicNumberByLabel(atom.label);
    }
    if (!Number.isFinite(z) || z <= 0) continue;
    const x = Number(atom.coords[0]);
    const y = Number(atom.coords[1]);
    const zCoord = Number(atom.coords[2]);
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(zCoord)) continue;
    out.push({
      label: atom.label || String(z),
      z: Number(z),
      coords: [x, y, zCoord],
    });
  }
  return out;
}

function isSquareMatrix(matrix, expectedSize = null) {
  if (!Array.isArray(matrix) || !matrix.length) return false;
  const width = matrix[0] && Array.isArray(matrix[0]) ? matrix[0].length : 0;
  if (!width || matrix.length !== width) return false;
  if (expectedSize !== null && matrix.length !== expectedSize) return false;
  for (let i = 0; i < matrix.length; i += 1) {
    if (!Array.isArray(matrix[i]) || matrix[i].length !== width) return false;
    for (let j = 0; j < width; j += 1) {
      if (!Number.isFinite(matrix[i][j])) return false;
    }
  }
  return true;
}

function getMathLibrary() {
  if (typeof math !== "undefined" && math && typeof math.eigs === "function") {
    return math;
  }
  return null;
}

function buildRepresentationRotationMatrix(principalBasis, representation) {
  const mapping = REP_MAP[representation] || REP_MAP.Ir;
  const xAxis = AXIS_INDEX[mapping[1]];
  const yAxis = AXIS_INDEX[mapping[2]];
  const zAxis = AXIS_INDEX[mapping[0]];
  return [
    principalBasis[xAxis].slice(),
    principalBasis[yAxis].slice(),
    principalBasis[zAxis].slice(),
  ];
}

function transpose3x3(matrix) {
  return [
    [matrix[0][0], matrix[1][0], matrix[2][0]],
    [matrix[0][1], matrix[1][1], matrix[2][1]],
    [matrix[0][2], matrix[1][2], matrix[2][2]],
  ];
}

function multiply3x3(a, b) {
  return [
    [
      a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0],
      a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1],
      a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2],
    ],
    [
      a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0],
      a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1],
      a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2],
    ],
    [
      a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0],
      a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1],
      a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2],
    ],
  ];
}

function transformCartesianHessian(hessian, rotation3) {
  const size = hessian.length;
  const out = Array.from({ length: size }, () => Array(size).fill(0));
  const atomCount = Math.floor(size / 3);
  const rotationT = transpose3x3(rotation3);

  for (let i = 0; i < atomCount; i += 1) {
    for (let j = 0; j < atomCount; j += 1) {
      const bi = i * 3;
      const bj = j * 3;
      const block = [
        [hessian[bi + 0][bj + 0], hessian[bi + 0][bj + 1], hessian[bi + 0][bj + 2]],
        [hessian[bi + 1][bj + 0], hessian[bi + 1][bj + 1], hessian[bi + 1][bj + 2]],
        [hessian[bi + 2][bj + 0], hessian[bi + 2][bj + 1], hessian[bi + 2][bj + 2]],
      ];
      const rotated = multiply3x3(multiply3x3(rotation3, block), rotationT);
      for (let r = 0; r < 3; r += 1) {
        for (let c = 0; c < 3; c += 1) {
          out[bi + r][bj + c] = rotated[r][c];
        }
      }
    }
  }
  return out;
}

function buildMassWeightedHessian(hessian, massesAu) {
  const size = hessian.length;
  const out = Array.from({ length: size }, () => Array(size).fill(0));
  for (let i = 0; i < size; i += 1) {
    const mi = massesAu[Math.floor(i / 3)];
    const sqrtMi = Math.sqrt(mi);
    for (let j = 0; j < size; j += 1) {
      const mj = massesAu[Math.floor(j / 3)];
      const sqrtMj = Math.sqrt(mj);
      out[i][j] = hessian[i][j] / (sqrtMi * sqrtMj);
    }
  }
  return out;
}

function computeInertiaMomentsFromCoordinates(atoms, massesAmu) {
  let ixx = 0;
  let iyy = 0;
  let izz = 0;
  for (let i = 0; i < atoms.length; i += 1) {
    const coord = atoms[i].coords;
    const m = massesAmu[i];
    ixx += m * (coord[1] * coord[1] + coord[2] * coord[2]);
    iyy += m * (coord[0] * coord[0] + coord[2] * coord[2]);
    izz += m * (coord[0] * coord[0] + coord[1] * coord[1]);
  }
  return [ixx, iyy, izz];
}

function computeInertiaMomentsBohr(coordsBohr, massesAu) {
  let ixx = 0;
  let iyy = 0;
  let izz = 0;
  for (let i = 0; i < coordsBohr.length; i += 1) {
    const coord = coordsBohr[i];
    const m = massesAu[i];
    ixx += m * (coord[1] * coord[1] + coord[2] * coord[2]);
    iyy += m * (coord[0] * coord[0] + coord[2] * coord[2]);
    izz += m * (coord[0] * coord[0] + coord[1] * coord[1]);
  }
  return [ixx, iyy, izz];
}

function isLikelyLinearFromMoments(moments) {
  if (!Array.isArray(moments) || moments.length < 3) return false;
  const maxMoment = Math.max(moments[0], moments[1], moments[2]);
  const minMoment = Math.min(moments[0], moments[1], moments[2]);
  if (!Number.isFinite(maxMoment) || maxMoment <= 0) return false;
  return minMoment / maxMoment < 1e-4;
}

function escapeHtml(text) {
  return String(text || "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

function renderWMSOutputRepresentationWarning(lines) {
  if (!ui.wmsOutputRepresentationWarning) return;
  const entries = Array.isArray(lines) ? lines.filter(Boolean) : [];
  if (!entries.length) {
    ui.wmsOutputRepresentationWarning.classList.remove("is-visible");
    ui.wmsOutputRepresentationWarning.innerHTML = "";
    return;
  }
  ui.wmsOutputRepresentationWarning.innerHTML = entries
    .map((line) => `<div class="representation-warning__line">${escapeHtml(line)}</div>`)
    .join("");
  ui.wmsOutputRepresentationWarning.classList.add("is-visible");
}

function loadExternalScriptOnce(src) {
  return new Promise((resolve, reject) => {
    if (typeof document === "undefined") {
      reject(new Error("Cannot load Pyodide dynamically outside a browser."));
      return;
    }
    const existing = document.querySelector(`script[data-wms-dynamic-src="${src}"]`);
    if (existing) {
      existing.addEventListener("load", () => resolve(), { once: true });
      existing.addEventListener("error", () => reject(new Error(`Unable to load ${src}`)), {
        once: true,
      });
      if (existing.dataset.loaded === "true") resolve();
      return;
    }
    const script = document.createElement("script");
    script.src = src;
    script.async = true;
    script.dataset.wmsDynamicSrc = src;
    script.addEventListener(
      "load",
      () => {
        script.dataset.loaded = "true";
        resolve();
      },
      { once: true }
    );
    script.addEventListener("error", () => reject(new Error(`Unable to load ${src}`)), {
      once: true,
    });
    document.head.appendChild(script);
  });
}

async function ensureLoadPyodideFunction() {
  if (typeof loadPyodide === "function") return;
  loadingOverlay.setMessage("Pyodide locale non trovato: scarico il motore Python dal CDN…");
  await loadExternalScriptOnce(PYODIDE_CDN_SCRIPT_URL);
  if (typeof loadPyodide !== "function") {
    throw new Error("Pyodide loader is not available after loading the CDN script.");
  }
}

async function ensureDistortionPyodide() {
  if (distortionPyodide) return distortionPyodide;

  const bootstrap = async (indexURL) => {
    await ensureLoadPyodideFunction();
    loadingOverlay.setMessage(
      indexURL === PYODIDE_CDN_INDEX_URL
        ? "Caricamento Pyodide dal CDN…"
        : "Caricamento Pyodide locale…"
    );
    if (!distortionPyodideReady) {
      distortionPyodideReady = loadPyodide({ indexURL });
    }
    const py = await distortionPyodideReady;
    if (!py.__wmsPrepDistortionLoaded) {
      loadingOverlay.setMessage("Preparazione dei moduli WMS-Prep…");
      const modules = await Promise.all(
        PYODIDE_DISTORTION_MODULES.map(async (entry) => {
          const response = await fetch(entry.url, { cache: "no-cache" });
          if (!response.ok) {
            throw new Error(`Unable to load Pyodide helper: ${entry.url}`);
          }
          return {
            fileName: entry.fileName,
            text: await response.text(),
          };
        })
      );
      for (let i = 0; i < modules.length; i += 1) {
        py.FS.writeFile(modules[i].fileName, modules[i].text);
      }
      py.runPython("import importlib; importlib.invalidate_caches(); import wms_prep_bridge");
      py.__wmsPrepDistortionLoaded = true;
    }
    distortionPyodide = py;
    return py;
  };

  try {
    return await bootstrap(PYODIDE_INDEX_URL);
  } catch (error) {
    distortionPyodide = null;
    distortionPyodideReady = null;
    if (PYODIDE_INDEX_URL === PYODIDE_CDN_INDEX_URL) {
      throw error;
    }
    try {
      loadingOverlay.setMessage("Pyodide locale non disponibile: ritento dal CDN…");
      return await bootstrap(PYODIDE_CDN_INDEX_URL);
    } catch (fallbackError) {
      distortionPyodide = null;
      distortionPyodideReady = null;
      throw fallbackError;
    }
  }
}

async function runPyodideDistortionTransform(payload) {
  loadingOverlay.show("Preparazione del calcolo delle costanti di distorsione…");
  try {
    await waitForNextPaint();
    const py = await ensureDistortionPyodide();
    loadingOverlay.setMessage("Conversione delle costanti di distorsione…");
    const inputJson = JSON.stringify(payload || {});
    py.globals.set("__wms_prep_payload_json", inputJson);
    const resultJson = await py.runPythonAsync(`
import json
from wms_prep_bridge import process_gaussian_log_payload
json.dumps(process_gaussian_log_payload(json.loads(__wms_prep_payload_json)))
`);
    return JSON.parse(String(resultJson || "{}"));
  } finally {
    try {
      if (distortionPyodide) {
        distortionPyodide.runPython("globals().pop('__wms_prep_payload_json', None)");
      }
    } catch (error) {
      // no-op
    }
    loadingOverlay.hide();
  }
}

function rotateAtomsBetweenRepresentations(atoms, inputRepresentation, outputRepresentation) {
  if (!Array.isArray(atoms) || !atoms.length) return [];
  const from = REP_MAP[normalizeFullWatsonRepresentationToken(inputRepresentation)] || REP_MAP.Ir;
  const to = REP_MAP[normalizeFullWatsonRepresentationToken(outputRepresentation)] || REP_MAP.Ir;
  return atoms.map((atom) => {
    const coords = Array.isArray(atom && atom.coords) ? atom.coords : [];
    const principal = {
      [from[1]]: Number(coords[0]),
      [from[2]]: Number(coords[1]),
      [from[0]]: Number(coords[2]),
    };
    return {
      ...atom,
      coords: [
        principal[to[1]],
        principal[to[2]],
        principal[to[0]],
      ],
    };
  });
}

function resolveRequestedOutputRepresentation(parsed, fallback = "Ir") {
  syncWMSOutputRepresentationControl();
  return normalizeFullWatsonRepresentationToken(
    ui.wmsOutputRepresentation ? ui.wmsOutputRepresentation.value : fallback,
    fallback
  );
}

async function assembleParsedForRepresentation(parsed, targetRepresentation) {
  const targetRep = normalizeFullWatsonRepresentationToken(targetRepresentation, "Ir");
  const geometryRep = normalizeFullWatsonRepresentationToken(
    parsed.geometryRepresentation || targetRep,
    targetRep
  );
  const dvibRep = normalizeFullWatsonRepresentationToken(
    parsed.dvibRepresentation || parsed.geometryRepresentation || targetRep,
    targetRep
  );
  const rotatedGeometryAtoms = rotateAtomsBetweenRepresentations(
    parsed.dpcs3Atoms,
    geometryRep,
    targetRep
  );
    const fixedDvib = {
      A: parsed.dvib ? parsed.dvib.A : 0,
      B: parsed.dvib ? parsed.dvib.B : 0,
      C: parsed.dvib ? parsed.dvib.C : 0,
    };

    const output = {
      ...parsed,

      // ABC sempre stampati nel loro ordine fisso originale
      rotationalConstantsMHz: parsed.rotationalConstantsMHz,

      // la geometria può continuare a ruotare se ti serve per XYZ / quadrupolo / ecc.
      geometryRepresentation: targetRep,

      // Unrotated DVib
      dvibRepresentation: dvibRep,
      dvib: {
        A: fixedDvib.A,
        B: fixedDvib.B,
        C: fixedDvib.C,
      },

      // Unrotated dipole
      dipole: parsed.dipole,

      dpcs3Atoms: rotatedGeometryAtoms,
      dpcs3XYZText: rotatedGeometryAtoms.length ? buildAtomicNumberXYZText(rotatedGeometryAtoms) : "",

      asym: parsed.asym ? { ...parsed.asym } : { ...EMPTY_QUARTIC_ASYM },
      sym: parsed.sym ? { ...parsed.sym } : { ...EMPTY_QUARTIC_SYM },
      asymSextic: parsed.asymSextic ? { ...parsed.asymSextic } : { ...EMPTY_SEXTIC_ASYM },
      symSextic: parsed.symSextic ? { ...parsed.symSextic } : { ...EMPTY_SEXTIC_SYM },

      // queste invece restano legate alla rappresentazione finale
      quarticRepresentation: targetRep,
      sexticRepresentation: targetRep,
      quadrupoleRepresentation: targetRep,

      quadrupole: {
        nuclei: [],
        wmsFields: { ...EMPTY_QUADRUPOLE_WMS_FIELDS },
      },

      assemblyWarnings: buildInitialWMSOutputWarnings(parsed, targetRep),
    };

  const transformed = await runPyodideDistortionTransform({
    outputRepresentation: targetRep,
    quartic: parsed.quarticRaw,
    sextic: parsed.sexticRaw,
    quadrupole: parsed.quadrupoleRaw,
  });

  if (transformed && transformed.quartic) {
    output.asym = transformed.quartic.asym || { ...EMPTY_QUARTIC_ASYM };
    output.sym = transformed.quartic.sym || { ...EMPTY_QUARTIC_SYM };
  }
  if (transformed && transformed.sextic) {
    output.asymSextic = transformed.sextic.asymSextic || { ...EMPTY_SEXTIC_ASYM };
    output.symSextic = transformed.sextic.symSextic || { ...EMPTY_SEXTIC_SYM };
  }
  if (transformed && transformed.quadrupole) {
    output.quadrupole = {
      nuclei: Array.isArray(transformed.quadrupole.nuclei)
        ? transformed.quadrupole.nuclei
        : [],
      wmsFields: {
        ...EMPTY_QUADRUPOLE_WMS_FIELDS,
        ...((transformed.quadrupole && transformed.quadrupole.wmsFields) || {}),
      },
    };
  }
  if (Array.isArray(transformed && transformed.warnings)) {
    output.assemblyWarnings.push(...transformed.warnings.filter(Boolean));
  }

  return output;
}

async function refreshWMSInputOutput() {
  const refreshToken = ++state.outputRefreshToken;
  if (!state.logSources.length) {
    state.combinedParsedData = null;
    state.wmsInputText = "";
    state.catrefInputText = "";
    state.catrefBundles = [];
    state.catrefVarText = "";
    state.catrefIntText = "";
    state.dpcs3XYZText = "";
    ui.wmsOutput.value = "";
    if (ui.catrefOutput) ui.catrefOutput.value = "";
    if (ui.dpcs3XYZOutput) ui.dpcs3XYZOutput.value = "";
    renderWMSOutputRepresentationWarning([]);
    ui.downloadWmsInputBtn.disabled = true;
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
    if (ui.openBdpcs3Btn) ui.openBdpcs3Btn.disabled = true;
    return;
  }

  const combined = combineConnectedSourcesForWMS();
  state.combinedParsedData = combined;
  if (!combined) {
    state.wmsInputText = "";
    state.catrefInputText = "";
    state.catrefVarText = "";
    state.catrefIntText = "";
    state.dpcs3XYZText = "";
    ui.wmsOutput.value = "# Connect at least one Geometry source to generate WMS-Rot input.\n";
    if (ui.catrefOutput) {
      ui.catrefOutput.value = "# Connect at least one Geometry source to generate CATREF input.\n";
    }
    if (ui.dpcs3XYZOutput) ui.dpcs3XYZOutput.value = "";
    renderWMSOutputRepresentationWarning([]);
    ui.downloadWmsInputBtn.disabled = true;
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
    if (ui.openBdpcs3Btn) ui.openBdpcs3Btn.disabled = true;
    return;
  }

  state.dpcs3XYZText = combined.dpcs3XYZText || "";
  if (ui.dpcs3XYZOutput) ui.dpcs3XYZOutput.value = state.dpcs3XYZText;
  const representationWarnings = buildInitialWMSOutputWarnings(
    combined,
    resolveRequestedOutputRepresentation(combined)
  );

  let rep = "Ir";
  let assembled = combined;
  try {
    rep = resolveRequestedOutputRepresentation(combined);
    assembled = await assembleParsedForRepresentation(combined, rep);
    if (refreshToken !== state.outputRefreshToken) return;
    state.combinedParsedData = assembled;
    renderWMSOutputRepresentationWarning(assembled.assemblyWarnings || representationWarnings);
    const text = buildWMSInputText(assembled, rep);
    state.wmsInputText = text;
    ui.wmsOutput.value = text;
    ui.downloadWmsInputBtn.disabled = false;
    state.dpcs3XYZText = assembled.dpcs3XYZText || "";
    if (ui.dpcs3XYZOutput) ui.dpcs3XYZOutput.value = state.dpcs3XYZText;
  } catch (error) {
    if (refreshToken !== state.outputRefreshToken) return;
    state.wmsInputText = "";
    state.catrefInputText = "";
    state.catrefBundles = [];
    state.catrefVarText = "";
    state.catrefIntText = "";
    renderWMSOutputRepresentationWarning(representationWarnings);
    ui.wmsOutput.value = `# Unable to build WMS-Rot input: ${error instanceof Error ? error.message : String(error)}\n`;
    if (ui.catrefOutput) {
      ui.catrefOutput.value = "# Unable to build CATREF input until the WMS-Rot input is available.\n";
    }
    ui.downloadWmsInputBtn.disabled = true;
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
    if (ui.openBdpcs3Btn) {
      ui.openBdpcs3Btn.disabled = !Array.isArray(combined.dpcs3Atoms) || !combined.dpcs3Atoms.length;
    }
    return;
  }

  try {
    const catrefSet = buildCATREFInputSet(assembled, rep, {
      baseName: buildLogOutputBaseName("catref"),
    });
    if (refreshToken !== state.outputRefreshToken) return;
    const defaultBundle =
      catrefSet.byReduction[DEFAULT_CATREF_REDUCTION] || catrefSet.bundles[0] || null;
    state.catrefInputText = catrefSet.bundleText;
    state.catrefBundles = catrefSet.bundles.slice();
    state.catrefVarText = defaultBundle ? defaultBundle.varText : "";
    state.catrefIntText = defaultBundle ? defaultBundle.intText : "";
    if (ui.catrefOutput) {
      ui.catrefOutput.value = catrefSet.bundleText;
    }
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = false;
  } catch (error) {
    if (refreshToken !== state.outputRefreshToken) return;
    state.catrefInputText = "";
    state.catrefBundles = [];
    state.catrefVarText = "";
    state.catrefIntText = "";
    if (ui.catrefOutput) {
      ui.catrefOutput.value = `# Unable to build CATREF input: ${error instanceof Error ? error.message : String(error)}\n`;
    }
    if (ui.downloadCatrefInputBtn) ui.downloadCatrefInputBtn.disabled = true;
  }
  if (ui.openBdpcs3Btn) {
    ui.openBdpcs3Btn.disabled = !Array.isArray(combined.dpcs3Atoms) || !combined.dpcs3Atoms.length;
  }
}

function combineConnectedSourcesForWMS() {
  const geometrySource = getConnectedSource("geometry");
  if (!geometrySource) return null;
  const quarticSource = getConnectedSource("quartic");
  const sexticSource = getConnectedSource("sextic");
  const dvibSource = getConnectedSource("dvib");
  const quadrupoleSource = getConnectedSource("quadrupole");
  return combineParsedComponents({
    geometry: geometrySource ? geometrySource.payload : null,
    quartic: quarticSource ? quarticSource.payload : null,
    sextic: sexticSource ? sexticSource.payload : null,
    dvib: dvibSource ? dvibSource.payload : null,
    quadrupole: quadrupoleSource ? quadrupoleSource.payload : null,
  });
}

function getConnectedSource(type) {
  const sourceId = state.connections[type];
  if (!sourceId) return null;
  return state.sourceById[sourceId] || null;
}

function combineParsedComponents({ geometry, quartic, sextic, dvib, quadrupole }) {
  if (!geometry) return null;
  return {
    pointGroup: geometry.pointGroup || "CS",
    rotationalConstantsMHz: geometry.rotationalConstantsMHz || { A: NaN, B: NaN, C: NaN },
    geometryRepresentation: geometry.inputRepresentation || null,
    dvibRepresentation: dvib ? dvib.inputRepresentation || null : null,
    dvib: dvib || { A: 0, B: 0, C: 0 },
    dipole: geometry.dipole || { a: "", b: "", c: "" },
    dpcs3Atoms: Array.isArray(geometry.dpcs3Atoms) ? geometry.dpcs3Atoms : [],
    dpcs3XYZText: geometry.dpcs3XYZText || "",
    quarticRepresentation: quartic ? quartic.inputRepresentation || null : null,
    quarticRaw: quartic
      ? {
          tauPrime: quartic.tauPrime ? { ...quartic.tauPrime } : null,
          sigma: Number.isFinite(quartic.sigma) ? quartic.sigma : null,
          rotationalConstantsMHz: quartic.rotationalConstantsMHz
            ? { ...quartic.rotationalConstantsMHz }
            : null,
          moments: quartic.moments ? { ...quartic.moments } : null,
          contextLabel: quartic.contextLabel || "",
          inputRepresentation: quartic.inputRepresentation || null,
        }
      : null,
    asym: quartic ? quartic.asym || { ...EMPTY_QUARTIC_ASYM } : { ...EMPTY_QUARTIC_ASYM },
    sym: quartic ? quartic.sym || { ...EMPTY_QUARTIC_SYM } : { ...EMPTY_QUARTIC_SYM },
    sexticRepresentation: sextic ? sextic.inputRepresentation || null : null,
    sexticRaw: sextic
      ? {
          inertial: sextic.inertial ? { ...sextic.inertial } : null,
          sigma: Number.isFinite(sextic.sigma) ? sextic.sigma : null,
          rotationalConstantsMHz: sextic.rotationalConstantsMHz
            ? { ...sextic.rotationalConstantsMHz }
            : null,
          moments: sextic.moments ? { ...sextic.moments } : null,
          contextLabel: sextic.contextLabel || "",
          inputRepresentation: sextic.inputRepresentation || null,
        }
      : null,
    asymSextic: sextic
      ? sextic.asymSextic || { ...EMPTY_SEXTIC_ASYM }
      : { ...EMPTY_SEXTIC_ASYM },
    symSextic: sextic
      ? sextic.symSextic || { ...EMPTY_SEXTIC_SYM }
      : { ...EMPTY_SEXTIC_SYM },
    quadrupoleRepresentation: quadrupole ? quadrupole.inputRepresentation || null : null,
    quadrupoleRaw: quadrupole
      ? {
          nuclei: Array.isArray(quadrupole.nuclei)
            ? quadrupole.nuclei.map((nucleus) => ({
                ...nucleus,
                chi: nucleus && nucleus.chi ? { ...nucleus.chi } : {},
              }))
            : [],
          contextLabel: quadrupole.contextLabel || "",
          inputRepresentation: quadrupole.inputRepresentation || null,
        }
      : null,
    quadrupole: {
      nuclei: [],
      wmsFields: { ...EMPTY_QUADRUPOLE_WMS_FIELDS },
    },
  };
}

function collectConnectedRepresentations(parsed) {
  if (!parsed || typeof parsed !== "object") return [];
  const candidates = [
    parsed.quarticRepresentation,
    parsed.sexticRepresentation,
    parsed.dvibRepresentation,
    parsed.quadrupoleRepresentation,
    parsed.geometryRepresentation,
  ];
  const representations = [];
  for (let i = 0; i < candidates.length; i += 1) {
    const normalized = normalizeWatsonRepresentationToken(candidates[i]);
    if (!normalized || representations.includes(normalized)) continue;
    representations.push(normalized);
  }
  return representations;
}

function resolveConnectedRepresentation(parsed, fallback = "Ir") {
  const normalizedFallback = normalizeWatsonRepresentationToken(fallback) || "Ir";
  const representations = collectConnectedRepresentations(parsed);
  if (!representations.length) {
    return normalizedFallback;
  }
  return representations[0];
}

function buildFlowGeneratorRepresentationNote() {
  const combined = combineConnectedSourcesForWMS();
  const targetRep = resolveRequestedOutputRepresentation(combined, "Ir");
  if (!combined) {
    return `final representation = ${targetRep}`;
  }
  const all = collectConnectedRepresentations(combined);
  if (!all.length) {
    return `final representation = ${targetRep} | source representation unavailable`;
  }
  if (all.length === 1) {
    return `final representation = ${targetRep} | source representation = ${all[0]}`;
  }
  return `final representation = ${targetRep} | source representations = ${all.join(", ")}`;
}

function buildWMSInputText(parsed, representation) {
  const outputParsed = parsed;
  const quadrupoleFields =
    outputParsed &&
    outputParsed.quadrupole &&
    outputParsed.quadrupole.wmsFields &&
    typeof outputParsed.quadrupole.wmsFields === "object"
      ? outputParsed.quadrupole.wmsFields
      : EMPTY_QUADRUPOLE_WMS_FIELDS;
  const lines = [];
  lines.push("#ROTATIONAL");
  lines.push("rotor_type = asymmetric");
  lines.push(`representation = ${representation}`);
  lines.push("Watson Reduction = A");
  lines.push(`Point Group = ${outputParsed.pointGroup || "CS"}`);
  lines.push("T_K = 4");
  lines.push(`A_MHz = ${formatCompactNumber(outputParsed.rotationalConstantsMHz.A)}`);
  lines.push(`B_MHz = ${formatCompactNumber(outputParsed.rotationalConstantsMHz.B)}`);
  lines.push(`C_MHz = ${formatCompactNumber(outputParsed.rotationalConstantsMHz.C)}`);
  lines.push(`DVibA_MHz= ${formatCompactNumber(outputParsed.dvib.A)}`);
  lines.push(`DVibB_MHz= ${formatCompactNumber(outputParsed.dvib.B)}`);
  lines.push(`DVibC_MHz= ${formatCompactNumber(outputParsed.dvib.C)}`);
  lines.push(`Dipole_a_D = ${outputParsed.dipole.a}`);
  lines.push(`Dipole_b_D = ${outputParsed.dipole.b}`);
  lines.push(`Dipole_c_D = ${outputParsed.dipole.c}`);
  lines.push(`I_NUC = ${quadrupoleFields.I_NUC}`);
  lines.push(`I_NUC_1 = ${quadrupoleFields.I_NUC_1}`);
  lines.push(`I_NUC_2 = ${quadrupoleFields.I_NUC_2}`);
  for (let i = 0; i < QUAD_CART_KEYS.length; i += 1) {
    const key = QUAD_CART_KEYS[i];
    lines.push(`${key} = ${quadrupoleFields[key]}`);
  }
  for (let i = 0; i < QUAD_CART_KEYS.length; i += 1) {
    const key = QUAD_CART_KEYS[i];
    lines.push(`${key}_1 = ${quadrupoleFields[`${key}_1`]}`);
  }
  for (let i = 0; i < QUAD_CART_KEYS.length; i += 1) {
    const key = QUAD_CART_KEYS[i];
    lines.push(`${key}_2 = ${quadrupoleFields[`${key}_2`]}`);
  }
  lines.push(`DELTA J_MHz = ${outputParsed.asym.deltaJ}`);
  lines.push(`DELTA JK_MHz = ${outputParsed.asym.deltaJK}`);
  lines.push(`DELTA K_MHz = ${outputParsed.asym.deltaK}`);
  lines.push(`delta J_MHz = ${outputParsed.asym.smallDeltaJ}`);
  lines.push(`delta K_MHz = ${outputParsed.asym.smallDeltaK}`);
  lines.push(`DJ_MHz = ${outputParsed.sym.DJ}`);
  lines.push(`DJK_MHz = ${outputParsed.sym.DJK}`);
  lines.push(`DK_MHz = ${outputParsed.sym.DK}`);
  lines.push(`d1_MHz = ${outputParsed.sym.d1}`);
  lines.push(`d2_MHz = ${outputParsed.sym.d2}`);
  lines.push(`PHI N_MHz = ${outputParsed.asymSextic.PHIN}`);
  lines.push(`PHI NK_MHz = ${outputParsed.asymSextic.PHINK}`);
  lines.push(`PHI KN_MHz = ${outputParsed.asymSextic.PHIKN}`);
  lines.push(`PHI K_MHz = ${outputParsed.asymSextic.PHIK}`);
  lines.push(`phi N_MHz = ${outputParsed.asymSextic.phiN}`);
  lines.push(`phi NK_MHz = ${outputParsed.asymSextic.phiNK}`);
  lines.push(`phi K_MHz = ${outputParsed.asymSextic.phiK}`);
  lines.push(`H N_MHz = ${outputParsed.symSextic.HN}`);
  lines.push(`H NK_MHz = ${outputParsed.symSextic.HNK}`);
  lines.push(`H KN_MHz = ${outputParsed.symSextic.HKN}`);
  lines.push(`H K_MHz = ${outputParsed.symSextic.HK}`);
  lines.push(`h1_MHz = ${outputParsed.symSextic.h1}`);
  lines.push(`h2_MHz = ${outputParsed.symSextic.h2}`);
  lines.push(`h3_MHz = ${outputParsed.symSextic.h3}`);
  return `${lines.join("\n")}\n`;
}

function buildCATREFInputFiles(parsed, representation, options = {}) {
  const outputParsed = parsed;
  const baseName = sanitizeCATREFBaseName(options.baseName || "gaussian_catref");
  const reduction = normalizeCATREFReduction(options.reduction || DEFAULT_CATREF_REDUCTION);
  const pointGroup = normalizePointGroupToken(outputParsed && outputParsed.pointGroup) || "C1";
  const temperatureK = normalizeCATREFNumericOption(
    options.temperatureK,
    DEFAULT_CATREF_TEMPERATURE_K
  );
  const str0 = normalizeCATREFNumericOption(options.str0, DEFAULT_CATREF_STR0);
  const str1 = normalizeCATREFNumericOption(options.str1, DEFAULT_CATREF_STR1);
  const freqMaxGHz = normalizeCATREFNumericOption(
    options.freqMaxGHz,
    DEFAULT_CATREF_FREQ_MAX_GHZ
  );
  const constants = buildCATREFConstantsFromParsed(outputParsed, reduction);
  const rotorInfo = inferCATREFRotorInfo(constants.A, constants.B, constants.C);
  const qrot = computeCATREFPartitionFunction({
    rotorType: rotorInfo.rotorType,
    rotationalConstantsMHz: constants,
    pointGroup,
    temperatureK,
  });
  const jMax = normalizeCATREFIntegerOption(
    options.jMax,
    estimateCATREFJMax({
      rotorType: rotorInfo.rotorType,
      rotationalConstantsMHz: constants,
      freqMaxGHz,
    })
  );
  const bundle = {
    baseName,
    representation: normalizeWatsonRepresentationToken(representation) || "Ir",
    reduction,
    pointGroup,
    rotorType: rotorInfo.rotorType,
    asymmetryClass: rotorInfo.asymmetryClass,
    kappa: rotorInfo.kappa,
    nvibSign: rotorInfo.nvibSign,
    qrot,
    jMax,
    str0,
    str1,
    freqMaxGHz,
    temperatureK,
    constants,
    varFileName: `${baseName}.var`,
    intFileName: `${baseName}.int`,
  };
  bundle.varText = buildCATREFVarText(bundle);
  bundle.intText = buildCATREFIntText(bundle);
  bundle.bundleText = buildCATREFBundleText(bundle);
  return bundle;
}

function buildCATREFInputSet(parsed, representation, options = {}) {
  const baseName = sanitizeCATREFBaseName(options.baseName || "gaussian_catref");
  const bundles = CATREF_REDUCTIONS.map((reduction) =>
    buildCATREFInputFiles(parsed, representation, {
      ...options,
      baseName: `${baseName}_${reduction}`,
      reduction,
    })
  );
  const byReduction = {};
  for (let i = 0; i < bundles.length; i += 1) {
    byReduction[bundles[i].reduction] = bundles[i];
  }
  return {
    baseName,
    bundles,
    byReduction,
    bundleText: buildCATREFInputSetText(bundles),
  };
}

function buildCATREFBundleText(bundle) {
  const lines = [];
  lines.push(`WMS-Prep CATREF export`);
  lines.push(`Representation: ${bundle.representation}`);
  lines.push(`Reduction: ${bundle.reduction}`);
  lines.push(`Rotor type: ${bundle.rotorType}${bundle.asymmetryClass ? ` (${bundle.asymmetryClass})` : ""}`);
  if (Number.isFinite(bundle.kappa)) {
    lines.push(`Ray kappa: ${bundle.kappa.toFixed(6)}`);
  }
  lines.push(`Point group: ${bundle.pointGroup}`);
  lines.push(`Qrot(${formatCompactNumber(bundle.temperatureK)} K): ${formatCompactNumber(bundle.qrot)}`);
  lines.push("");
  lines.push(`===== ${bundle.varFileName} =====`);
  lines.push(bundle.varText.trimEnd());
  lines.push("");
  lines.push(`===== ${bundle.intFileName} =====`);
  lines.push(bundle.intText.trimEnd());
  lines.push("");
  return lines.join("\n");
}

function buildCATREFInputSetText(bundles) {
  return `${bundles
    .map((bundle) => buildCATREFBundleText(bundle).trimEnd())
    .join("\n\n")}\n`;
}

function buildCATREFVarText(bundle) {
  const title = String(bundle.baseName || "gaussian_catref");
  const constants = bundle.constants || {};
  const rotorType = String(bundle.rotorType || "asymmetric");
  const reduction = normalizeCATREFReduction(bundle.reduction);
  const reductionCode = reduction.toLowerCase();
  const parameterLines = [];

  if (rotorType === "linear") {
    parameterLines.push(buildCATREFVarParameterLine(100, constants.B, "B"));
    parameterLines.push(buildCATREFVarParameterLine(200, -constants.DJ, "-D"));
    parameterLines.push(buildCATREFVarParameterLine(300, constants.HJ, "H"));
  } else if (reduction === "S") {
    const d2Id = Number.isFinite(constants.d2Id) ? constants.d2Id : 50000;
    parameterLines.push(buildCATREFVarParameterLine(10000, constants.A, "A"));
    parameterLines.push(buildCATREFVarParameterLine(20000, constants.B, "B"));
    parameterLines.push(buildCATREFVarParameterLine(30000, constants.C, "C"));
    parameterLines.push(buildCATREFVarParameterLine(200, constants.DJ, "DJ"));
    parameterLines.push(buildCATREFVarParameterLine(1100, constants.DJK, "DJK"));
    parameterLines.push(buildCATREFVarParameterLine(2000, constants.DK, "DK"));
    parameterLines.push(buildCATREFVarParameterLine(40100, constants.dJ, "d1"));
    parameterLines.push(buildCATREFVarParameterLine(d2Id, constants.dK, d2Id === 41000 ? "d2" : "d2"));
    parameterLines.push(buildCATREFVarParameterLine(300, constants.HJ, "H N"));
    parameterLines.push(buildCATREFVarParameterLine(1200, constants.HJK, "H NK"));
    parameterLines.push(buildCATREFVarParameterLine(2100, constants.HKJ, "H KN"));
    parameterLines.push(buildCATREFVarParameterLine(3000, constants.HK, "H K"));
    parameterLines.push(buildCATREFVarParameterLine(50100, constants.h1, "h 1"));
    parameterLines.push(buildCATREFVarParameterLine(51000, constants.h2, "h 2"));
    parameterLines.push(buildCATREFVarParameterLine(60000, constants.h3, "h 3"));
  } else {
    parameterLines.push(buildCATREFVarParameterLine(10000, constants.A, "A"));
    parameterLines.push(buildCATREFVarParameterLine(20000, constants.B, "B"));
    parameterLines.push(buildCATREFVarParameterLine(30000, constants.C, "C"));
    parameterLines.push(buildCATREFVarParameterLine(200, constants.DJ, "DJ"));
    parameterLines.push(buildCATREFVarParameterLine(1100, constants.DJK, "DJK"));
    parameterLines.push(buildCATREFVarParameterLine(2000, constants.DK, "DK"));
    parameterLines.push(buildCATREFVarParameterLine(40100, constants.dJ, "dJ"));
    parameterLines.push(buildCATREFVarParameterLine(41000, constants.dK, "dK"));
    parameterLines.push(buildCATREFVarParameterLine(300, constants.HJ, "HJ"));
    parameterLines.push(buildCATREFVarParameterLine(1200, constants.HJK, "HJK"));
    parameterLines.push(buildCATREFVarParameterLine(2100, constants.HKJ, "HKJ"));
    parameterLines.push(buildCATREFVarParameterLine(3000, constants.HK, "HK"));
    parameterLines.push(buildCATREFVarParameterLine(40200, constants.h1, "h1"));
    parameterLines.push(buildCATREFVarParameterLine(41100, constants.h2, "h2"));
    parameterLines.push(buildCATREFVarParameterLine(42000, constants.h3, "h3"));
  }

  const npar = parameterLines.length;
  const lines = [title];
  if (rotorType === "linear") {
    lines.push(
      `${String(npar).padStart(4)}   31    4    0    1.0000E-010    1.0000E+006    1.0000E+000 1.0000000000`
    );
    lines.push("l   -1    1    0    0    0    1    1    1         0   -1    0");
  } else {
    lines.push(
      `${String(npar).padStart(4)}  ${String(Math.max(50, bundle.jMax + 25)).padStart(3)}  0  0  0.0000E+000  1.0000E+006  -1.0000E+000  1.0000000000`
    );
    lines.push(
      `${reductionCode}   1  ${bundle.nvibSign > 0 ? "1" : "-1"}  0  ${String(bundle.jMax).padStart(2)}  0  1  1  1  0  0  0`
    );
  }
  lines.push(...parameterLines);
  return `${lines.join("\n")}\n`;
}

function buildCATREFIntText(bundle) {
  const lines = [];
  lines.push(String(bundle.baseName || "gaussian_catref"));
  lines.push(
    `${String(DEFAULT_CATREF_FLAGS)}  ${String(DEFAULT_CATREF_TAG)}  ${formatCATREFIntNumber(bundle.qrot, 4)}  0  ${bundle.jMax}  ${formatCATREFIntNumber(bundle.str0, 1)}  ${formatCATREFIntNumber(bundle.str1, 1)}  ${formatCATREFIntNumber(bundle.freqMaxGHz, 3)}  ${formatCATREFIntNumber(bundle.temperatureK, 3)}`
  );

  const dipoles = bundle.constants && bundle.constants.dipole ? bundle.constants.dipole : { a: 0, b: 0, c: 0 };
  const rotorType = String(bundle.rotorType || "asymmetric");
  const primaryLinearDipole = pickPrimaryDipoleEntry(dipoles);
  const entries = rotorType === "linear"
    ? [{ id: "001", value: primaryLinearDipole.value, axis: primaryLinearDipole.axis }]
    : [
        { id: "001", value: dipoles.a, axis: "a" },
        { id: "002", value: dipoles.b, axis: "b" },
        { id: "003", value: dipoles.c, axis: "c" },
      ];

  for (let i = 0; i < entries.length; i += 1) {
    const entry = entries[i];
    lines.push(
      `${String(entry.id).padStart(6)}  ${formatCATREFIntNumber(entry.value, 6)}  / mu(${entry.axis}) /`
    );
  }
  return `${lines.join("\n")}\n`;
}

function buildCATREFVarParameterLine(id, value, label, uncertainty = 0) {
  const idText = String(id).padStart(13);
  const valueText = toFortranScientific(Number.isFinite(value) ? value : 0, 12).padStart(22);
  const uncertaintyText = toFortranScientific(
    Number.isFinite(uncertainty) ? uncertainty : 0,
    8
  ).padStart(18);
  return `${idText} ${valueText} ${uncertaintyText} / ${label} /`;
}

function buildCATREFConstantsFromParsed(parsed, reduction) {
  const asym = parsed && parsed.asym ? parsed.asym : {};
  const sym = parsed && parsed.sym ? parsed.sym : {};
  const asymSextic = parsed && parsed.asymSextic ? parsed.asymSextic : {};
  const symSextic = parsed && parsed.symSextic ? parsed.symSextic : {};
  const dipole = parsed && parsed.dipole ? parsed.dipole : {};
  const A = Number(parsed && parsed.rotationalConstantsMHz && parsed.rotationalConstantsMHz.A);
  const B = Number(parsed && parsed.rotationalConstantsMHz && parsed.rotationalConstantsMHz.B);
  const C = Number(parsed && parsed.rotationalConstantsMHz && parsed.rotationalConstantsMHz.C);
  const reductionKey = normalizeCATREFReduction(reduction);
  const constants = {
    A,
    B,
    C,
    dipole: {
      a: parseCATREFNumeric(dipole.a, 0),
      b: parseCATREFNumeric(dipole.b, 0),
      c: parseCATREFNumeric(dipole.c, 0),
    },
    d2Id: 50000,
  };
  if (reductionKey === "S") {
    constants.DJ = parseCATREFNumeric(sym.DJ, parseCATREFNumeric(asym.deltaJ, 0));
    constants.DJK = parseCATREFNumeric(sym.DJK, parseCATREFNumeric(asym.deltaJK, 0));
    constants.DK = parseCATREFNumeric(sym.DK, parseCATREFNumeric(asym.deltaK, 0));
    constants.dJ = parseCATREFNumeric(sym.d1, parseCATREFNumeric(asym.smallDeltaJ, 0));
    constants.dK = parseCATREFNumeric(sym.d2, parseCATREFNumeric(asym.smallDeltaK, 0));
    constants.HJ = parseCATREFNumeric(symSextic.HN, parseCATREFNumeric(asymSextic.PHIN, 0));
    constants.HJK = parseCATREFNumeric(symSextic.HNK, parseCATREFNumeric(asymSextic.PHINK, 0));
    constants.HKJ = parseCATREFNumeric(symSextic.HKN, parseCATREFNumeric(asymSextic.PHIKN, 0));
    constants.HK = parseCATREFNumeric(symSextic.HK, parseCATREFNumeric(asymSextic.PHIK, 0));
    constants.h1 = parseCATREFNumeric(symSextic.h1, parseCATREFNumeric(asymSextic.phiN, 0));
    constants.h2 = parseCATREFNumeric(symSextic.h2, parseCATREFNumeric(asymSextic.phiNK, 0));
    constants.h3 = parseCATREFNumeric(symSextic.h3, parseCATREFNumeric(asymSextic.phiK, 0));
  } else {
    constants.DJ = parseCATREFNumeric(asym.deltaJ, parseCATREFNumeric(sym.DJ, 0));
    constants.DJK = parseCATREFNumeric(asym.deltaJK, parseCATREFNumeric(sym.DJK, 0));
    constants.DK = parseCATREFNumeric(asym.deltaK, parseCATREFNumeric(sym.DK, 0));
    constants.dJ = parseCATREFNumeric(asym.smallDeltaJ, parseCATREFNumeric(sym.d1, 0));
    constants.dK = parseCATREFNumeric(asym.smallDeltaK, parseCATREFNumeric(sym.d2, 0));
    constants.HJ = parseCATREFNumeric(asymSextic.PHIN, parseCATREFNumeric(symSextic.HN, 0));
    constants.HJK = parseCATREFNumeric(asymSextic.PHINK, parseCATREFNumeric(symSextic.HNK, 0));
    constants.HKJ = parseCATREFNumeric(asymSextic.PHIKN, parseCATREFNumeric(symSextic.HKN, 0));
    constants.HK = parseCATREFNumeric(asymSextic.PHIK, parseCATREFNumeric(symSextic.HK, 0));
    constants.h1 = parseCATREFNumeric(asymSextic.phiN, parseCATREFNumeric(symSextic.h1, 0));
    constants.h2 = parseCATREFNumeric(asymSextic.phiNK, parseCATREFNumeric(symSextic.h2, 0));
    constants.h3 = parseCATREFNumeric(asymSextic.phiK, parseCATREFNumeric(symSextic.h3, 0));
  }
  return constants;
}

function inferCATREFRotorInfo(A, B, C) {
  const close = (x, y) =>
    Number.isFinite(x) &&
    Number.isFinite(y) &&
    Math.abs(x - y) <= Math.max(1e-9, 1e-9 * Math.max(Math.abs(x), Math.abs(y)));
  const finiteA = Number.isFinite(A);
  const finiteB = Number.isFinite(B);
  const finiteC = Number.isFinite(C);
  let rotorType = "asymmetric";

  if (finiteB && finiteC && close(B, C) && (!finiteA || A > Math.max(B, C) * 1e4)) {
    rotorType = "linear";
  } else if (finiteA && finiteB && finiteC && close(A, B) && close(B, C)) {
    rotorType = "spherical";
  } else if (finiteA && finiteB && finiteC && close(B, C) && A > B) {
    rotorType = "prolate";
  } else if (finiteA && finiteB && finiteC && close(A, B) && B > C) {
    rotorType = "oblate";
  }

  const kappa = computeRayAsymmetryKappa(A, B, C);
  let asymmetryClass = "";
  if (rotorType === "asymmetric" && Number.isFinite(kappa)) {
    if (kappa > 0.05) asymmetryClass = "quasi-oblate";
    else if (kappa < -0.05) asymmetryClass = "quasi-prolate";
    else asymmetryClass = "intermediate";
  } else if (rotorType === "prolate") {
    asymmetryClass = "prolate";
  } else if (rotorType === "oblate") {
    asymmetryClass = "oblate";
  } else if (rotorType === "linear") {
    asymmetryClass = "linear";
  } else if (rotorType === "spherical") {
    asymmetryClass = "spherical";
  }

  return {
    rotorType,
    asymmetryClass,
    kappa,
    nvibSign: rotorType === "oblate" || (rotorType === "asymmetric" && Number.isFinite(kappa) && kappa > 0)
      ? -1
      : 1,
  };
}

function computeRayAsymmetryKappa(A, B, C) {
  if (!Number.isFinite(A) || !Number.isFinite(B) || !Number.isFinite(C)) return Number.NaN;
  const denom = A - C;
  if (Math.abs(denom) < 1e-12) return Number.NaN;
  return (2 * B - A - C) / denom;
}

function computeCATREFPartitionFunction({
  rotorType,
  rotationalConstantsMHz,
  pointGroup,
  temperatureK,
}) {
  const sigma = getRotationalSymmetryNumber(pointGroup);
  const T = Number.isFinite(temperatureK) && temperatureK > 0 ? temperatureK : DEFAULT_CATREF_TEMPERATURE_K;
  const AHz = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.A, 0)) * 1e6;
  const BHz = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.B, 0)) * 1e6;
  const CHz = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.C, 0)) * 1e6;
  if (String(rotorType || "").toLowerCase() === "linear") {
    const bHz = firstFinitePositive(BHz, CHz, AHz, 0);
    if (bHz <= 0) return 1.0;
    return Math.max(1.0, (BOLTZMANN_CONSTANT_SI * T) / (PLANCK_CONSTANT_SI * bHz * Math.max(1, sigma)));
  }
  if (AHz <= 0 || BHz <= 0 || CHz <= 0) return 1.0;
  const qrot =
    Math.sqrt(Math.PI) *
    Math.pow((BOLTZMANN_CONSTANT_SI * T) / PLANCK_CONSTANT_SI, 1.5) /
    (Math.sqrt(AHz * BHz * CHz) * Math.max(1, sigma));
  return Math.max(1.0, qrot);
}

function getRotationalSymmetryNumber(pointGroup) {
  const pg = normalizePointGroupToken(pointGroup);
  if (!pg) return 1;
  if (pg === "CINFV") return 1;
  if (pg === "DINFH") return 2;
  if (pg === "C1" || pg === "CS" || pg === "CI") return 1;
  if (pg === "TD") return 12;
  if (pg === "OH") return 24;
  if (pg === "IH") return 60;
  let match = pg.match(/^D(\d+)$/i);
  if (match) return Math.max(1, 2 * Number.parseInt(match[1], 10));
  match = pg.match(/^C(\d+)$/i);
  if (match) return Math.max(1, Number.parseInt(match[1], 10));
  match = pg.match(/^S(\d+)$/i);
  if (match) {
    const n = Number.parseInt(match[1], 10);
    return Number.isFinite(n) && n > 0 ? Math.max(1, Math.floor(n / 2)) : 1;
  }
  return 1;
}

function estimateCATREFJMax({ rotorType, rotationalConstantsMHz, freqMaxGHz }) {
  const freqMaxMHz = Math.max(1, normalizeCATREFNumericOption(freqMaxGHz, DEFAULT_CATREF_FREQ_MAX_GHZ) * 1000);
  const A = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.A, 0));
  const B = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.B, 0));
  const C = Math.abs(parseCATREFNumeric(rotationalConstantsMHz.C, 0));
  let scaleMHz = firstFinitePositive(C, B, A, 0);
  if (String(rotorType || "").toLowerCase() === "linear") {
    scaleMHz = firstFinitePositive(B, C, A, 0);
  }
  if (scaleMHz <= 0) return DEFAULT_CATREF_MIN_JMAX;
  const estimate = Math.ceil(freqMaxMHz / (2 * scaleMHz)) + 5;
  return Math.max(
    DEFAULT_CATREF_MIN_JMAX,
    Math.min(DEFAULT_CATREF_MAX_JMAX, estimate)
  );
}

function normalizeCATREFReduction(value) {
  const raw = String(value || DEFAULT_CATREF_REDUCTION).trim().toUpperCase();
  return raw === "S" ? "S" : "A";
}

function sanitizeCATREFBaseName(value) {
  const raw = String(value || "gaussian_catref").trim();
  const withoutExt = raw.replace(/\.[^.]+$/, "");
  const safe = withoutExt
    .replace(/\s+/g, "_")
    .replace(/[^A-Za-z0-9_.-]/g, "_")
    .replace(/_+/g, "_")
    .replace(/^_+|_+$/g, "");
  return safe || "gaussian_catref";
}

function parseCATREFNumeric(value, fallback = 0) {
  if (typeof value === "number") {
    return Number.isFinite(value) ? value : fallback;
  }
  const parsed = parseSciNumber(value);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function normalizeCATREFNumericOption(value, fallback) {
  return Number.isFinite(value) ? Number(value) : fallback;
}

function normalizeCATREFIntegerOption(value, fallback) {
  const source = Number.isFinite(value) ? value : fallback;
  const rounded = Math.round(source);
  return Number.isFinite(rounded) ? rounded : Math.round(fallback);
}

function firstFinitePositive(...values) {
  for (let i = 0; i < values.length; i += 1) {
    const value = Number(values[i]);
    if (Number.isFinite(value) && value > 0) return value;
  }
  return 0;
}

function pickPrimaryDipoleEntry(dipoles) {
  const candidates = [
    { axis: "a", value: Number(dipoles && dipoles.a) },
    { axis: "b", value: Number(dipoles && dipoles.b) },
    { axis: "c", value: Number(dipoles && dipoles.c) },
  ];
  for (let i = 0; i < candidates.length; i += 1) {
    if (Number.isFinite(candidates[i].value) && candidates[i].value > 0) {
      return candidates[i];
    }
  }
  return { axis: "a", value: 0 };
}

function formatCATREFIntNumber(value, fractionDigits = 3) {
  if (!Number.isFinite(value)) return "";
  const sanitized = Math.abs(value) < 5e-13 ? 0 : value;
  return sanitized.toFixed(fractionDigits).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

function formatCompactNumber(value) {
  if (!Number.isFinite(value)) return "";
  const sanitized = Math.abs(value) < 5e-13 ? 0 : value;
  return sanitized.toFixed(6).replace(/(?:\.0+|(\.\d+?)0+)$/, "$1");
}

function downloadWMSInput() {
  if (!state.wmsInputText) return;
  const fileName = buildLogOutputFileName("input.txt");
  downloadTextFile(state.wmsInputText, fileName);
}

function downloadCATREFInput() {
  if (Array.isArray(state.catrefBundles) && state.catrefBundles.length) {
    for (let i = 0; i < state.catrefBundles.length; i += 1) {
      const bundle = state.catrefBundles[i];
      if (bundle && bundle.varText && bundle.varFileName) {
        downloadTextFile(bundle.varText, bundle.varFileName);
      }
      if (bundle && bundle.intText && bundle.intFileName) {
        downloadTextFile(bundle.intText, bundle.intFileName);
      }
    }
    return;
  }
  if (!state.catrefVarText || !state.catrefIntText) return;
  downloadTextFile(state.catrefVarText, buildLogOutputFileName("catref.var"));
  downloadTextFile(state.catrefIntText, buildLogOutputFileName("catref.int"));
}

function onOpenInBDPCS3Clicked() {
  try {
    openInBDPCS3();
  } catch (error) {
    setStatus(error instanceof Error ? error.message : String(error), true);
  }
}

function openInBDPCS3() {
  if (!state.combinedParsedData || !Array.isArray(state.combinedParsedData.dpcs3Atoms)) return;
  const xyzText = buildBDPCS3XYZFileText(state.combinedParsedData.dpcs3Atoms);
  const fileName = buildLogOutputFileName("dpcs3.xyz");

  try {
    localStorage.setItem(
      BDPCS3_PRELOAD_STORAGE_KEY,
      JSON.stringify({
        xyz: xyzText,
        fileName,
        createdAt: Date.now(),
      })
    );
  } catch (error) {
    throw new Error("Unable to store XYZ data for BDPCS3 preload.");
  }

  const targetUrl = "/BDPCS3/index.html";
  const opened = window.open(targetUrl, "_blank", "noopener");
  if (!opened) {
    throw new Error("Popup blocked while opening BDPCS3.");
  }
}

function buildBDPCS3XYZFileText(atoms) {
  const lines = [];
  lines.push(String(atoms.length));
  lines.push("DPCS3 geometry from WMS-Prep");
  for (let i = 0; i < atoms.length; i += 1) {
    const atom = atoms[i];
    lines.push(
      `${atom.z} ${formatNumber(atom.coords[0])} ${formatNumber(atom.coords[1])} ${formatNumber(atom.coords[2])}`
    );
  }
  return `${lines.join("\n")}\n`;
}

function downloadTextFile(content, fileName) {
  const blob = new Blob([content], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = fileName;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function buildOutputFileName(rep, extension) {
  const base = state.fileName ? state.fileName.replace(/\.[^.]+$/, "") : "representation_fixed";
  return `${base}_${rep}.${extension}`;
}

function buildLogOutputBaseName(suffixBase) {
  const base = state.logFileName ? state.logFileName.replace(/\.[^.]+$/, "") : "gaussian";
  return `${base}_${suffixBase}`;
}

function buildLogOutputFileName(suffix) {
  return `${buildLogOutputBaseName("").replace(/_$/, "")}_${suffix}`;
}

function setStatus(message, isError = false) {
  if (!ui.status) return;
  ui.status.textContent = message;
  ui.status.style.color = isError ? "#b10000" : "#0b4e9b";
}

const REPRESENTATION_FIXER_EXPORTS = {
  parseXYZ,
  parseNumpyMatrixText,
  parseHessianXYZBundleText,
  parseFchkGeometry,
  parseFchkCartesianForceConstants,
  parseFchkDipoleMomentOptional,
  parseTauPayloadText,
  parseTauComponentsFromText,
  rotateTauComponentsForRepresentation,
  rotateTauPrimeComponentsForRepresentation,
  detectInputRepresentationFromPrincipalBasis,
  buildQuarticFromTauComponents,
  normalizeTauRepresentation,
  parseSexticPayloadText,
  rotateSexticInertialComponents,
  buildSexticFromInertialComponents,
  parseGaussianLogComponents,
  parseGaussianLogForWMSInput,
  createSourcesFromParsedLog,
  buildBDPCS3DerivedGeometryPayload,
  autoConnectFlowSources,
  combineParsedComponents,
  resolveConnectedRepresentation,
  assembleParsedForRepresentation,
  buildWMSInputText,
  buildCATREFInputFiles,
  buildCATREFInputSet,
  buildCATREFVarText,
  buildCATREFIntText,
  computeHarmonicQuarticFromHessian,
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = REPRESENTATION_FIXER_EXPORTS;
}
