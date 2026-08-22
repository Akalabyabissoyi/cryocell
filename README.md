# CryoCell

A virtual human mesenchymal stem cell for testing cryoprotectants in silico.

Desktop application. Multi-compartment biophysics, a deformable cell that
crenates and blebs because the mechanics say it should, an extracellular ice
field that closes in from every direction, the full interphase organelle
complement annotated with Human Protein Atlas data, and a post-thaw
mechanotransduction module that produces apoptosis and necrosis rather than
assuming them.

## Install and run

```bash
pip install -r requirements.txt
python3 run_cryocell.py
```

Python 3.10+. On macOS everything is pip-installable; no system Qt needed.

## What is on screen

**Left** — every protocol and cell parameter as a live control. Changing one
re-runs the model (debounced ~0.2 s).

**Centre** — the cell.

* **scroll** to zoom, 0.35× to 14×. Past ~3.5× the plasma membrane resolves
  into two leaflets with CPA partitioned into the headgroup region and pores
  opening where the MD regime says they should.
* **shift-drag** (or middle-drag) to pan.
* **hover** any compartment for its live simulated state and its Human Protein
  Atlas annotation — gene count, share of the proteome, how much of it is
  multi-localised.
* **drag** the cell to deform it, or drag an organelle to displace it. The soft
  body responds and recovers; tension, Piezo1 and ROCK respond with it.
* **compartment checkboxes** strip the view down, the way a confocal channel
  picker does.

**Right** — outcome, live plots, and the analysis tabs:

* **Population** — 120 individual cells with sampled volume, Lp, Ps, Vb, sterol
  and nucleation-site density, run across all cores. Returns a distribution and
  a variance attribution, not a single number.
* **Perturbation screen** — disables each model component in turn and ranks
  what it costs. Most rows carry a bench analogue.
* **Measure next?** — sweeps every uncertain parameter across how badly it is
  known and tells you which measurement would most reduce the uncertainty.
* **Stability** — re-integrates at half and double step size and reports drift.
* **HPA reference** — the Cell Atlas numbers the inspector draws on.

## The ice field

The cell does not sit in a slab between two ice sheets. It sits in an irregular
unfrozen pocket bounded by the growth fronts of many grains, and that pocket
closes in from every direction as the ice fraction rises. Grain boundaries are
the Voronoi tessellation of the grain centres, so the ice reads as a
polycrystal. The pocket is elongated rather than isotropic, because a cell is
not compressible: when the intergranular channel narrows below the cell
diameter the cell deforms into a longer pocket rather than being crushed.

## The deformable cell

The membrane is a 160-vertex contour with tension springs whose **rest length
comes from the membrane's own contour length**. Membrane area is conserved when
a cell dehydrates; enclosed area is not — so the contour has nowhere to go but
buckle. Crenation is a consequence of the osmotic module, not a drawing
decision, and DMSO makes it worse because area per lipid grows +0.85%/mol%.

Bending stiffness sets the fold wavelength and has a narrow useful window
(swept numerically: above ~0.030 the membrane stops expressing slack as folds
at all). The actin cortex uses a dead-band, because restoring pointwise irons
the folds flat and restoring on a smoothed profile lets spikes grow into a
star. Blebs nucleate stochastically at a rate set by the ROCK/pMLC index.

## Biophysics (literature-grounded)

The transport and thermodynamics core is built from published cryobiology, not
convenience approximations — see `REFERENCES.md` for the citation behind every
term. In brief:

* **Kedem–Katchalsky 3-parameter transport** — water and CPA fluxes with a real
  reflection coefficient σ and solvent-drag coupling (Kedem & Katchalsky 1958;
  Kleinhans 1998), on a **volume-dependent membrane area**, with CPA-specific
  permeation activation energies.
* **Non-ideal solution thermodynamics** — freezing-point depression via the
  osmotic virial equation (Elliott et al. 2007), with **composition-specific
  Tg′ and eutectic** per CPA, so e.g. glycerol arrests far warmer than DMSO.
* **Intercellular ice propagation** — junction-coupled, autocatalytic, so
  adherent monolayers and spheroids suffer far more intracellular ice than a
  dissociated suspension (Irimia & Karlsson 2002; Acker & McGann 2000). IIF
  lethality scales with the supercooling at nucleation.
* **Delayed-onset death (CIDOC)** — a ROS variable accumulates on rewarming and
  drives a day-1-to-3 apoptotic wave, so the outcome now reports a settled 72 h
  viability (S_72), not just the immediate window (Baust 2000, 2009).

## Provenance

Every parameter is tiered **measured / estimated / fitted** in the project
documents. Three coefficients are outright fits with nothing behind them: the
mobility exponent for damage chemistry, the grain-growth constant, and the
mechanical-squeeze coefficient. The **Measure next?** tab will tell you where
those sit relative to Lp and ATP turnover in terms of what to calibrate first.

Subcellular data: Human Protein Atlas Cell Atlas
(https://www.proteinatlas.org/humanproteome/subcellular), retrieved 16 August
2026 — 13,603 genes localised experimentally across 49 subcellular locations.
Volume fractions used for the drawing are conventional mammalian values, not
HPA data; the HPA counts proteins, not volume, and this is flagged in the app.

Cryobiology parameters and the calibration table are documented in
`cryocell-virtual-cell-model-spec.md`, `cryocell-whole-cell-architecture.md`
and `cryocell-mechanotransduction-module.md` in the project.

## Validation

The model reproduces, without being fitted to any of them: the inverted-U
survival curve against cooling rate with an optimum at ~1 °C/min; hMSC in
suspension at 38% against Heng's measured 39.8%; the biphasic Y-27632
dose-response including the fall above 30 µM; z-VAD-fmk as the strongest single
intervention; and the protective effect of a 4 °C post-thaw hold.

## Not a validated predictor

A mechanistic model with real parameters. Use it to rank protocols and reason
about mechanism; calibrate before quoting a number.

## Tests

```bash
pip install -r requirements.txt pytest
pytest
```

The suite checks the persistent-cell invariants (signature immutability,
provenance completeness, clone/branch parenting) and pins the calibration
anchor: hMSC + DMSO at 1 °C/min stays within 30–42% 24 h survival (35.4% at
head).

## License

MIT — see [LICENSE](LICENSE). If you use CryoCell in your research, citation
metadata is in [CITATION.cff](CITATION.cff).
