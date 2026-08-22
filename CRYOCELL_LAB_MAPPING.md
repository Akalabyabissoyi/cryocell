# CryoCell Lab — Integration Mapping & Phased Plan

This is the required agent deliverable (Agent-Instructions §17, Visualisation §30).
It is written after reading `CRYOCELL_AGENT_INSTRUCTIONS.md`,
`CRYOCELL_VISUALISATION_SPECIFICATION.md`, and `REFERENCES.md` in full, and
inspecting the existing engine (`cryocell/model.py`), analysis layer
(`cryocell/analysis.py`), and GUI (`cryocell/app.py`, `cellview.py`).

**Budget note:** this integration was intentionally scoped to a bounded, additive
foundation to stay within the requested ~20% weekly-usage ceiling. It does **not**
implement the general-biology layers or the new "CryoCell Lab" GUI. It lays the
persistent-cell + provenance architecture and honestly exposes what the calibrated
physics already supports. The rest is planned below.

## 1. What was done this pass (additive; nothing existing removed)
- **Folder merge:** `~/cryocell` is the single canonical copy; `~/cryocel` and
  `~/cryoce` were aligned to it and can now be deleted (they are byte-identical).
- **`cryocell/virtualcell.py` (new):** a persistent `CryoCell` object with identity,
  state, provenance history, `save_state / restore / clone / branch / compare /
  replay / export`, `state_hash`, the **5 perturbation methods**, and the **15
  signature methods**. Cryobiology signatures wrap the real engine; general-biology
  signatures return explicit `{"supported": False}`.
- **`tests/test_virtualcell.py` (new):** proves the mandated invariants —
  signature immutability, state-transition on perturbation, clone/branch parent
  preservation, honest-unsupported, provenance completeness, calibration anchor
  held (hMSC 1 °C/min within 30–42%). All pass.
- **Untouched:** `model.py`, `analysis.py`, `app.py`, `cellview.py`, the calibration
  (hMSC/DMSO optimum ~1 °C/min ≈ 35%), and every existing GUI workflow.

## 2. REFERENCES.md mechanisms → implementation modules
The mandatory mechanisms (Agent-Instructions §10) are **already implemented** in the
existing engine — this pass exposes them through the new API, it did not re-derive them.

| # | Mechanism | Module / status |
|---|---|---|
| 1 | Kedem–Katchalsky water/CPA transport | `model.py` transport block · get_membrane_transport_state |
| 2 | CPA-specific σ, Ea | `CPAS` library · exposed |
| 3 | Volume-dependent membrane area + reserve | `model.py` A_inst |
| 4–6 | Non-ideal virial thermodynamics; composition Tg′/eutectic | `model.py` dT_from_osm, B_eff, tg_run |
| 7–9 | Extra/intra-cellular ice; propagation; IIF magnitude | nucleation + `k_prop` + supercool scaling |
| 10 | Mitochondrial osmometer/CPA/ice/ΔΨm/ROS/MPT | mito block · get_organelle_state |
| 11 | Delayed-onset death (CIDOC) | ROS→S_72 · get_cell_fate(72h) |
| 12 | Vitrification / glass fraction / mobility | WLF `mobility` · glass in get_thermal_ice_state |
| 13 | Cold-dependent cytoskeleton | MT/actin/IF cold + tensegrity |
| 14 | Protein hydration/water-replacement/denaturation | protein-stability block · prot |
| 15–16 | CPA + additive combinations; membrane/IRI protection | `ADDITIVES` · get_membrane_transport_state |
| 17–18 | Cytoskeletal tensegrity; mechanical→mito coupling | tensegrity + mech_mito · get_osmotic_mechanical_state |
| 19 | Cell-line robustness + stated limits | phenotype params · CELL_LINES + warnings |
| 20–21 | Numerical stability; validation limits | `analysis.stability_check`; documented as literature_prior |

## 3. Canonical 20-method API (support status)
**15 signatures** — cryobiology = real; biology = honest `unsupported`:

| Signature | Status |
|---|---|
| get_thermal_ice_state | ✅ simulated |
| get_osmotic_mechanical_state | ✅ simulated (incl. cytoskeletal element loads) |
| get_membrane_transport_state | ✅ simulated |
| get_organelle_state | ✅ simulated |
| get_cell_morphology | ✅ features (image = GUI render, flagged) |
| get_stress_response | ✅ partial (ROS/mech/tox/protein; ER/DNA/ISR unsupported) |
| get_metabolic_state | ✅ partial (ATP/ΔΨm/ROS; full flux network unsupported) |
| get_cell_fate | ✅ simulated (24 h + 72 h, endpoints flagged) |
| get_protein_state | ⚠️ bulk native-fraction only; per-protein unsupported |
| get_genomic / epigenomic / transcriptomic | ⛔ unsupported (honest) |
| get_protein_interactions / _complexes | ⛔ unsupported (honest) |
| get_pathway_state | ⛔ unsupported (honest) |

**5 perturbations** — all implemented, mapping to engine `Params`:
apply_chemical_perturbation · apply_physical_protocol · apply_environmental_perturbation ·
apply_genetic_perturbation (mechanism-flag knockouts) · apply_biological_perturbation
(adhesion + phenotype). Compatibility aliases: apply_cpa, apply_thermal_protocol,
get_intracellular_ice, get_mitochondrial_state.

## 4. Honest limits (must not be overstated)
- **Not a stateful step integrator.** The engine runs a complete 9-phase protocol per
  call. "Stateful perturbation" = protocol composition + re-simulation (same as the GUI).
  True resumable mid-protocol handoff is a future refactor (Phase 2).
- **Single calibration anchor** (hMSC + DMSO). HeLa/A549/Jurkat contexts are
  `literature_prior`, flagged, never reported as cell-type-validated.
- **No biology layers.** Genome/epigenome/transcriptome/proteome/interactions/pathways
  are not modelled; their methods return `unsupported`. No pathway/interaction value is
  ever fabricated.
- **No new "CryoCell Lab" GUI** in this pass. The existing GUI (illustrative/
  fluorescence/phase-contrast render, scientific labels, cryo-signatures, cell-line and
  additive controls, video export) is preserved and untouched.

## 5. Phased plan for the remainder (each a separate budgeted unit)
- **Phase 2 — engine & state:** resumable mid-protocol integration; versioned
  save-file migration; conservation-check unit tests (water/CPA/solute/energy).
- **Phase 3 — pathway knowledge layer:** versioned, evidence-tagged pathway graph
  (Reactome/GO where licensing permits) wired to the physical drivers the engine
  already computes (cooling→membrane→ion→ΔΨm→ROS→apoptosis). Cross-talk edges typed
  per §6.4. Read-only, uncertainty-tagged.
- **Phase 4 — proteome / interactions:** bulk→per-protein where data supports it;
  interaction/complex retrieval with source-typed confidence (never a "probability"
  unless calibrated).
- **Phase 5 — CryoCell Lab GUI:** the AIDO-inspired workspace (experiment thread,
  command bar, branching, evidence workspace) built additively around the existing
  cell view, per the visualisation spec's MVP (§28) first.

## 6. Suggested next experiments (ranked by uncertainty reduction)
1. Measure the target line's **water permeability Lp** (dominant uncertainty; the
   `Measure next?` tool already ranks this).
2. **24 h AND 72 h recovery** (the model predicts delayed loss — validates CIDOC).
3. **Intracellular-ice fraction vs cooling rate** (anchors the nucleation parameters
   for a new cell type before trusting the U-curve).
