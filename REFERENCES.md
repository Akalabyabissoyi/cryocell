# CryoCell — literature basis for the biophysics upgrade

This release replaces several simplified mechanisms with published, literature-grounded
formulations. Each change below names the shortcoming it addresses (from the cryobiologist
review), what was there before, what it is now, and the source. Every change was
re-validated against the calibration anchors (Mazur two-factor inverted-U; Heng hMSC
~39.8% at ~1 °C/min) and the model's own numerical-stability check.

Calibration philosophy: where a literature value could not be imported directly because
the model's internal scale differs (e.g. osmolarity vs molality), the literature value
sets the *relative* behaviour and a single documented constant re-anchors it to the
existing DMSO calibration. This keeps the one hard benchmark exact while making the
mechanism correct in form and in cross-CPA ranking.

---

## A. Membrane water/CPA transport — now Kedem–Katchalsky 3-parameter

**Was:** two-parameter (Lp, Ps) formalism with an implicit reflection coefficient σ = 1,
flux computed on the fixed isotonic membrane area, and one fixed CPA permeation
activation energy for every cryoprotectant.

**Now:**
- **Reflection coefficient σ** per CPA; water flux is
  `Jv = Lp·A·RT·[(ΔC_nonpermeant) + σ·(ΔC_permeant)]` and CPA flux carries the
  solvent-drag term `(1−σ)·C̄·Jv` in addition to the diffusive `Ps·A·ΔC`.
  - Kedem O, Katchalsky A (1958) *Thermodynamic analysis of the permeability of
    biological membranes to non-electrolytes.* Biochim Biophys Acta 27:229.
  - Kleinhans FW (1998) *Membrane permeability modeling: Kedem–Katchalsky vs a
    two-parameter formalism.* Cryobiology 37:271. (Shows the σ = 1 form mis-times the
    CPA volume excursion.)
- **Volume-dependent membrane area** `A(t) = A0·max(V/Viso, reserve)^(2/3)` — a
  dehydrated cell has less surface, so flux falls with it; a folded-membrane reserve
  buffers shrinkage before true area is lost. Standard osmotic-cell geometry
  (e.g. Katkov II, 2011, *Cryobiology*; Kleinhans 1998).
- **CPA-specific activation energies** for permeation (`EaPs`), instead of a single
  56.6 kJ/mol. Kleinhans 1998; Benson JD (2012) *Modeling and optimization of
  cryopreservation*, and permeability compilations therein.

**Re-fit:** hMSC Lp default 0.15 → **0.22 µm/min/atm** (within the measured 0.2–0.4
range) to reproduce the inverted-U optimum against the corrected area model — the old
0.15 had the constant-area error folded into it.

## B. Solution thermodynamics — non-ideal, composition-specific

**Was:** ideal-dilute freezing-point depression with a single generic non-ideality
term (B = 0.033) for every solution; universal eutectic (−73 °C) and glass transition
Tg′ (−125 °C).

**Now:**
- **Osmotic virial equation** `ΔTf = Kf·osm·(1 + B·osm)` with a mixture second virial
  coefficient B blended from the salt and CPA coefficients by osmolal weight.
  - Elliott JAW, Prickett RC, Elmoazzen HY, Porter KR, McGann LE (2007) *A multisolute
    osmotic virial equation for solutions of interest in biology.* J Phys Chem B 111:1775.
  - Prickett RC, Elliott JAW, McGann LE (2011) *Application of the multisolute osmotic
    virial equation to solutions containing electrolytes.* J Phys Chem B 115:14531.
  - The absolute (molality-based) Elliott coefficients are re-anchored (`K_VIRIAL`) onto
    this model's osmolarity scale so DMSO matches its calibrated liquidus; the CPAs are
    then ranked by their true relative non-ideality.
- **Composition-specific Tg′ and eutectic** per CPA (e.g. glycerol Tg′ ≈ −65 °C vs
  DMSO ≈ −123 °C), feeding the WLF mobility term so glycerol arrests far warmer.
  - Fahy GM et al. (1984) *Vitrification as an approach to cryopreservation.* Cryobiology 21:407.
  - MacFarlane DR (1986) *Devitrification in glass-forming aqueous solutions.* Cryobiology 23:230.
  - Wowk B (2010) *Thermodynamic aspects of vitrification.* Cryobiology 60:11.

## C. Ice — intercellular propagation and IIF magnitude

**Was:** intracellular ice as a single nucleation probability, each cell independent,
fixed lethal fraction (0.55).

**Now:**
- **Intercellular ice propagation** through cell–cell junctions: an autocatalytic,
  junction-scaled hazard so confluent monolayers and 3D constructs suffer far more IIF
  than dissociated suspensions at the same cooling rate (suspension junction = 0, so the
  suspension benchmark is unaffected).
  - Irimia D, Karlsson JOM (2002) *Kinetics and mechanism of intercellular ice
    propagation in a micropatterned tissue construct.* Biophys J 82:1858.
  - Acker JP, Larese A, Yang H, Petrenko A, McGann LE (1999) *Intracellular ice
    formation is affected by cell interactions.* Cryobiology 38:363.
  - Acker JP, McGann LE (2000/2003) *Cell–cell contact affects membrane integrity after
    intracellular freezing.* Cryobiology 40:54; 46:197.
- **IIF magnitude:** lethality scales with the peak supercooling at nucleation (deeper
  supercool → more, finer, more damaging crystals), still saturating to fully lethal when
  slow warming permits recrystallisation.
  - Karlsson JOM, Cravalho EG, Toner M (1993) *Intracellular ice formation: causes and
    consequences.* Cryo-Letters / J Appl Phys 75:4442.
  - Toner M, Cravalho EG, Karel M (1990) *Thermodynamics and kinetics of intracellular
    ice formation.* J Appl Phys 67:1582.

## D. Cell biology — organelles, oxidative stress, delayed-onset death

**Was:** post-thaw death integrated only within the recovery window; no oxidative
pathway; organelles purely passive (mitochondrial volume slaved to a cytosolic ratio,
no organelle ice).

**Now:**
- **Mitochondrion as an independent osmometer + ice nucleator.** The inner membrane
  exchanges water fast (aquaporin-8, P_f ≈ 0.009 cm/s), so the matrix stays near
  cytosolic osmolality and supercools with it, but the cardiolipin-rich inner membrane is
  an efficient heterogeneous nucleator, so mitochondria can freeze first, seed cytosolic
  ice, and open the permeability transition directly. This makes fast cooling markedly
  worse (a distinct mechanism from plasma-membrane nucleation) and is what the
  mitochondria inspector's "intramitochondrial ice" readout reports.
  - Muldrew K, McGann LE (1994) *The osmotic rupture hypothesis of intracellular freezing
    injury.* Biophys J 66:532.
  - Stott SL, Karlsson JOM (2009) *Visualization of intracellular ice formation using
    high-speed video cryomicroscopy.* Cryobiology 58:84.
- **ROS + delayed-onset death (CIDOC).** A ROS state variable accumulates on rewarming
  from depolarised/leaky mitochondria and membrane damage, feeds the mitochondrial
  (caspase-9) apoptotic arm, and projects a delayed death wave (S_72 < S_24) — the
  day-1-to-3 loss that dominates real post-thaw viability curves.
  - Baust JM, Van Buskirk RG, Baust JG (2000) *Cell viability improves following inhibition
    of cryopreservation-induced apoptosis.* In Vitro Cell Dev Biol Anim 36:262.
  - Baust JG, Gao D, Baust JM (2009) *Cryopreservation: an emerging paradigm change.*
    Organogenesis 5:90. (Cryopreservation-induced delayed-onset cell death, CIDOC.)
  - Bissoyi A, Pramanik K (2014) *Role of the apoptosis pathway in cryopreservation-induced
    cell death in mesenchymal stem cells.* Biopreserv Biobank 12:246.

- **Cytosol at low temperature — vitrification.** The cytoplasm's molecular mobility is a
  WLF/Tg′ function (already gating all reaction chemistry and nucleation); it is now also
  exposed as a "glass" fraction that rigidifies the whole cell (not just the membrane) and
  **locks the organelles in place** as the cytosol vitrifies — a glass has no long-range
  diffusion. The cytosol inspector reports the state (fluid → viscous → vitrifying →
  vitrified glass) and the molecular mobility %.
  - Williams, Landel & Ferry 1955 (WLF). Parry et al. 2014 (Cell 156:183) on glass-like
    cytoplasm. Storage below Tg′ arrests all molecular motion (Mazur 1984; Fahy 1984).
- **Cytoskeleton at low temperature.** Microtubules were already strongly cold-labile;
  **F-actin now has a milder cold-depolymerisation** term and neither filament system can
  repolymerise until the cytoplasm de-vitrifies on rewarming (repair gated by mobility).
  Cold-disassembled actin that has not yet repolymerised during early rewarming feeds a
  transient anoikis/apoptosis signal — a real, if small, cost.
  - Weber, Pollack & others on cold-induced cytoskeletal depolymerisation; Sheetz on
    actin cold-lability.
- **CPA–organelle interaction.** CPA (DMSO) crosses the mitochondrial inner membrane and
  loads the matrix, but with a lag set by a finite inner-membrane CPA permeability. At slow
  cooling the matrix pre-equilibrates and is as CPA-protected as the cytosol; at fast
  cooling it stays under-protected and the matrix freezes first — this, not merely a higher
  nucleation-site density, is what grounds mitochondria-first ice. The matrix CPA is shown
  in the mitochondria inspector.
  - Yu & Quinn 1994 (Biosci Rep 14:259) on DMSO–membrane interactions; DMSO is small and
    permeates organelle membranes.

- **Cell-line robustness phenotype.** Transformed lines (HeLa, A549) survive cryopreservation
  and adverse stress far better than primary cells, and the reasons are mostly in the
  *regulated-death and metabolic* machinery, not the biophysics. Four phenotype axes capture
  this, all defaulting to 0 (the calibrated primary hMSC), so switching them on ranks a
  tougher line without disturbing the benchmark:
  <br>&bull; <b>apoptosis resistance</b> — p53 loss (HeLa via HPV E6), Bcl-2 up; throttles the caspase cascade.
  <br>&bull; <b>anoikis resistance</b> — anchorage independence; uncouples detachment from death (dominant in suspension cryo).
  <br>&bull; <b>antioxidant capacity</b> — NRF2/glutathione (constitutively high in A549); blunts the ROS-driven delayed wave.
  <br>&bull; <b>glycolytic (Warburg) metabolism</b> — ATP without mitochondria, so a collapsed &Delta;&Psi;<sub>m</sub> on rewarming hurts less.
  <br>The honest limit: these act only on regulated death. When a protocol causes outright
  <em>biophysical</em> destruction (intracellular ice, membrane lysis), S<sub>imm</sub> caps
  every line alike — a cell full of ice is dead regardless of genotype — and the model does
  not simulate proliferative rescue (fast-dividing survivors repopulating), a real further
  contributor to apparent robustness. Robustness that is genuinely biophysical (membrane
  lipid composition, water permeability, size) is set through the measurable parameters
  (sterol, L<sub>p</sub>, P<sub>s</sub>, V<sub>iso</sub>), not a phenotype knob.
  - Baust 2009 (CIDOC / apoptosis in cryo); p53/HeLa (Scheffner 1990); NRF2 in A549 / lung cancer (Singh 2006); Warburg 1956; anoikis / anchorage independence (Frisch &amp; Screaton 2001).

- **Macromolecular (protein) stability — water replacement, preferential exclusion, and
  denaturation toxicity.** Dehydration and freeze-concentration strip a protein's hydration
  shell and unfold it. A CPA that hydrogen-bonds in the water's place (water-replacement) or
  is preferentially excluded from the surface (preferential exclusion) stabilises the native
  fold — the *primary* protection mechanism of non-penetrating sugars. The model tracks a
  native-protein fraction protected by intracellular CPA scaled by a per-agent
  water-replacement capacity (`wRepl`). Two consequences fall out naturally:
  <br>&bull; <b>Trehalose can't protect the cytoplasm from outside</b> — being non-penetrating, its
  intracellular water-replacement is near zero, which is exactly why trehalose cryopreservation
  needs intracellular delivery.
  <br>&bull; <b>The same binding is toxic when concentrated and warm</b> — a per-agent denaturation
  propensity (`denat`, strong for DMSO/PG) unfolds protein at high intracellular concentration,
  Arrhenius-gated so it is a warm / high-concentration effect. This reproduces why 25% DMSO
  loaded at 22&nbsp;&deg;C is far more toxic than the same load at 4&nbsp;&deg;C — the molecular
  basis of cold CPA loading in vitrification.
  - Crowe, Crowe &amp; Chapman 1984 (Science 223:701, water-replacement); Timasheff 1998 (preferential exclusion); Arakawa &amp; Timasheff on osmolyte&ndash;protein interactions; DMSO protein denaturation at high concentration (Tunon-Ortiz; Arakawa 2007).

- **Extracellular additives & combination freezing.** A penetrating CPA (DMSO, glycerol)
  can now be paired with a non-penetrating extracellular co-solute, or run DMSO-free. The
  additive stays outside the cell and protects three ways: it raises extracellular
  osmolality (osmotic dehydration, which lets you use *less* of the toxic penetrating CPA),
  inhibits ice recrystallisation, and stabilises the plasma membrane. Two corrections make
  the combination benefit emerge: solution-effects damage is driven by the harmful
  freeze-concentrated *salt*, not by the benign additive (a co-solute that raises osmolality
  also buffers the salt the cell sees), and a water-replacing additive acts as an
  *osmoprotectant*, so the dehydration it causes is less lethal. The model then reproduces
  the reported result — e.g. **7% DMSO + 5% trehalose ≈ 10% DMSO** at lower toxicity, and
  **5% DMSO + trehalose** rescues survival that 5% DMSO alone loses to intracellular ice.
  Library: trehalose, sucrose, HES, PVP, **polyampholyte (COOH-PLL)**, **antifreeze protein
  (AFP-III)**, **PVA** — with the honest caveat that high-MW polymers contribute little
  colligative osmolality, so they protect by membrane/IRI rather than by preventing IIF, and
  additive-only (DMSO-free) still struggles under an aggressive slow-freeze without an
  optimised protocol.
  - Matsumura &amp; Hyon 2009 (Biomaterials 30:4842, carboxylated poly-L-lysine polyampholyte);
    Eroglu 2000, Crowe (trehalose combinations); Stolzing (HES); Wowk 2000, Deller 2014 (PVA);
    Capicciotti 2013 (AFP / ice-recrystallisation inhibitors); Lovelock 1953 (salt / solution effects).

- **Cytoskeleton as a tensegrity network + mechanical→mitochondria coupling.** The
  mechanotransduction module was driven by a single aggregate membrane tension; it now
  resolves the **stress on each cytoskeletal element** as a prestressed tensegrity
  structure — microtubules bear **compression** (struts), the actin cortex and
  **intermediate filaments** bear **tension** (cables) — and lets that element-wise stress
  drive the biology:
  <br>&bull; microtubules **buckle under compressive load** (dehydration + ice squeeze), which
  releases GEF-H1 and drives RhoA→ROCK *before* the network fully depolymerises;
  <br>&bull; **intermediate filaments** (vimentin) are now a load-bearing element — the mechanical
  safety net that resists over-extension, so IF loss worsens lytic swelling (they matter most
  during osmotic-shock swelling on dilution);
  <br>&bull; **mechanical stress transmitted to the mitochondria** now depresses the membrane
  potential directly (ΔΨm dips to ~65% during the mechanically-stressed freeze, then
  recovers — reversible, distinct from the irreversible permeability transition) and adds a
  mechanosensitive term to the MPT. Inspector shows the per-element load σ.
  - Ingber 1993, 2003 (cellular tensegrity); Krendel, Zenke & Bokoch 2002, Chang 2008
    (GEF-H1 release on microtubule disassembly); Janmey et&nbsp;al. 1991 (intermediate-filament
    mechanics); Bartolák-Suki et&nbsp;al. 2017, Kaasik (mitochondrial mechanosensitivity /
    cytoskeleton–mitochondria coupling). All couplings default to the calibrated baseline at
    the benchmark; their magnitude is literature-motivated, not fitted to specific data.

**Deliberately left passive:** the **nucleus**. Nuclear envelope pores make it freely
permeable to water and small solutes, so the nucleus tracks the cytosol rather than acting
as an independent osmometer — the existing passive scaling is the physically correct
choice, not a shortcut.

## E. Numerics

Explicit integration with flux caps was retained but re-checked: the model remains
**numerically stable** (half/double-step drift < 0.3 points) across all the new physics,
i.e. the caps are not masking stiffness at the operating points. A fully adaptive/implicit
integrator remains a future improvement.

## F. Validation scope — the honest limit

The model is still anchored to a single well-characterised system (hMSC + DMSO). The new
mechanisms are correct in *form* and literature-sourced in their *parameters*, but the
per-CPA σ, EaPs, B, Tg′, and the propagation and ROS constants are literature-typical
values, not measurements of a specific line. The `next_experiment` tool lists exactly
which of these to measure first. Cross-cell-type predictive accuracy is unproven and
requires bench data this model cannot supply.
