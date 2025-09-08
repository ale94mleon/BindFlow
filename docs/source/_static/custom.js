// === highlight.js theme switcher ===
function setHighlightTheme() {
  const old = document.getElementById("hljs-theme");
  if (old) old.remove();

  const themeLink = document.createElement("link");
  themeLink.rel = "stylesheet";
  themeLink.id = "hljs-theme";

  if (document.documentElement.getAttribute("data-theme") === "dark") {
    themeLink.href =
      "https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/atom-one-dark.min.css";
  } else {
    themeLink.href =
      "https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css";
  }

  document.head.appendChild(themeLink);
}

// === main logic ===
document.addEventListener("DOMContentLoaded", function () {
  // apply theme once
  setHighlightTheme();
  // watch for theme toggles (pydata-sphinx-theme / furo etc.)
  const observer = new MutationObserver(() => setHighlightTheme());
  observer.observe(document.documentElement, {
    attributes: true,
    attributeFilter: ["data-theme"],
  });

  // === DOM elements ===
  const preview = document.getElementById("env-preview");
  const downloadBtn = document.getElementById("download-env");
  const styleButtons = document.querySelectorAll(
    "#installer-config button[data-style]"
  );
  const featureCheckboxes = document.querySelectorAll(
    "#installer-config input[data-feature]"
  );

  // === Base environments ===
  const basePinnedEnv = {
    name: "BindFlow",
    channels: ["conda-forge"],
    conda_dependencies: [
      "python=3.10.14",
      "pip",
      "bioconda::snakemake=7.32.2",
      "openmm=8.0.0",
      "openmmforcefields>=0.11.2",
      "ambertools>=22.5",
      "pdbfixer=1.9",
      "openff-toolkit=0.14.2",
      "pyyaml=6.0",
      "parmed=4.1.0",
      "pandas=2.0.3",
      "numpy=1.26.4",
    ],
    pip_dependencies: [
      "alchemlyb==2.0.0",
      "mdanalysis==2.5.0",
      // "mdrestraintsgenerator>=0.2.1",
      "git+https://github.com/ale94mleon/MDRestraintsGenerator@51a543afc8ee2b0b14188995412b7865103203a6",
      "pymbar==4.0.1",
      "rdkit==2023.3.2",
      "toff==0.2.0"
    ],
  };

  const baseRelaxedEnv = {
    name: "BindFlow",
    channels: ["conda-forge"],
    conda_dependencies: [
      // Need this python version for the current MDRestraintsGenerator
      "python >=3.8,<3.12",
      "pip",
      "bioconda::snakemake",
      "openmm",
      "openmmforcefields",
      "ambertools",
      "pdbfixer",
      "openff-toolkit",
      "pyyaml",
      "parmed",
      "pandas",
      "numpy"
    ],
    pip_dependencies: [
      "alchemlyb",
      "mdanalysis",
      // "mdrestraintsgenerator>=0.2.1",
      "git+https://github.com/ale94mleon/MDRestraintsGenerator@51a543afc8ee2b0b14188995412b7865103203a6",
      "pymbar",
      "rdkit",
      "toff>=0.2.0"
    ],
  };

  // === Feature dependencies ===
  // Pinned features (exact versions)
  const featureDepsPinned = {
    mmpbsa: { conda: ["mpi4py=3.1.5"], pip: [] },
    espaloma: { conda: ["espaloma=0.3.2"], pip: [] },
  };

  // Relaxed features (looser versions)
  const featureDepsRelaxed = {
    mmpbsa: { conda: ["mpi4py=3.1.5"], pip: [] },
    espaloma: { conda: ["espaloma>=0.3.2"], pip: [] }, // same version can be relaxed if desired
  };

  // === State ===
  let state = {
    style: null,
    features: new Set(),
  };

  // Deep copy
  function clone(obj) {
    return JSON.parse(JSON.stringify(obj));
  }

  // Build YAML string
  function envToYAML(env) {
    let yaml = `name: ${env.name}\nchannels:\n`;
    env.channels.forEach((c) => {
      yaml += `  - ${c}\n`;
    });

    yaml += "dependencies:\n";
    env.conda_dependencies.forEach((dep) => {
      yaml += `  - ${dep}\n`;
    });

    if (env.pip_dependencies && env.pip_dependencies.length > 0) {
      yaml += "  - pip:\n";
      env.pip_dependencies.forEach((dep) => {
        yaml += `    - ${dep}\n`;
      });
    }

    return yaml;
  }

  // Update preview and download link
  function updateEnv() {
    if (!state.style) {
      preview.textContent = "Select options above...";
      downloadBtn.style.display = "none";
      return;
    }

    let env =
      state.style === "pinned" ? clone(basePinnedEnv) : clone(baseRelaxedEnv);

    // Select the right featureDeps object
    const currentFeatureDeps =
      state.style === "pinned" ? featureDepsPinned : featureDepsRelaxed;

    // Add feature dependencies
    state.features.forEach((f) => {
      if (currentFeatureDeps[f]) {
        env.conda_dependencies.push(...currentFeatureDeps[f].conda);
        env.pip_dependencies.push(...currentFeatureDeps[f].pip);
      }
    });

    const yaml = envToYAML(env);

    // reset preview content
    preview.innerHTML = "";
    const codeEl = document.createElement("code");
    codeEl.className = "language-yaml";
    codeEl.textContent = yaml;
    preview.appendChild(codeEl);

    // apply highlight
    hljs.highlightElement(codeEl);

    // prepare download
    downloadBtn.style.display = "inline-block";
    const blob = new Blob([yaml], { type: "text/yaml" });
    const url = URL.createObjectURL(blob);
    downloadBtn.setAttribute("href", url);
    downloadBtn.setAttribute("download", "environment.yml");
  }

  // === Event listeners ===
  // Style buttons
  styleButtons.forEach((btn) => {
    btn.addEventListener("click", () => {
      styleButtons.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      state.style = btn.getAttribute("data-style");
      updateEnv();
    });
  });

  // Feature checkboxes
  featureCheckboxes.forEach((cb) => {
    cb.addEventListener("change", () => {
      const f = cb.getAttribute("data-feature");
      if (cb.checked) {
        state.features.add(f);
      } else {
        state.features.delete(f);
      }
      updateEnv();
    });
  });
});