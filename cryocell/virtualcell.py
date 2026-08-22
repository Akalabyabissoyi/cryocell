"""
virtualcell.py — persistent virtual-cell wrapper (CryoCell Lab foundation).

An ADDITIVE layer over the existing, calibrated model.py engine. It wraps the
cryobiology simulation in a persistent, observable, perturbable cell object with
provenance and branching, per CRYOCELL_AGENT_INSTRUCTIONS.md. It changes nothing
in model.py, app.py, or the GUI — the existing workflows are untouched.

Scope honesty (mandatory per the spec):
  * The underlying engine runs a COMPLETE 9-phase protocol per call; it is not
    yet a resumable mid-protocol integrator. So "stateful perturbation" here means
    protocol/context COMPOSITION followed by re-simulation — the same thing the
    GUI already does. True resumable mid-protocol handoff is a future refactor.
  * The general-biology layers (genome, epigenome, transcriptome, proteome,
    protein interactions/complexes, pathways) are NOT modelled. Their signature
    methods return an explicit `{"supported": False, ...}` result. They never
    fabricate pathway activities, expression, or interaction scores.
  * Every readout carries a source_type and calibration_status. Literature-
    motivated parameters are never reported as cell-type-validated.
"""
from __future__ import annotations
import hashlib
import json
import time
import uuid
from dataclasses import asdict

from .model import Params, simulate, CPAS, ADDITIVES

MODEL_VERSION = "cryocell-2026.08"
PARAM_SET_VERSION = "hMSC-DMSO-calibration-1"

# Calibrated / literature-parameterised cell contexts (map to model params).
# hMSC is the calibration anchor; others carry the documented robustness
# phenotype and are literature_prior, NOT cell-type-validated.
CELL_LINES = {
    "hMSC":   dict(Viso=1800, cyto=1.35, sterol=25),
    "HeLa":   dict(Viso=1400, cyto=1.00, sterol=30, apop_resist=0.75,
                   anoikis_resist=0.80, glycolytic=0.60, antioxidant=0.35),
    "A549":   dict(Viso=1100, cyto=0.85, sterol=28, apop_resist=0.40,
                   anoikis_resist=0.55, glycolytic=0.55, antioxidant=0.85),
    "Jurkat": dict(Viso=700,  cyto=0.90, sterol=27, apop_resist=0.30,
                   anoikis_resist=0.50, glycolytic=0.55, antioxidant=0.40),
}

_CPA_NAMES = {"dmso": "dmso", "glycerol": "glycerol", "ethylene glycol": "eg",
              "eg": "eg", "propylene glycol": "pg", "pg": "pg", "trehalose": "tre"}
_ADD_NAMES = {"trehalose": "tre", "sucrose": "suc", "hes": "hes",
              "hydroxyethyl starch": "hes", "pvp": "pvp", "polyampholyte": "pll",
              "cooh-pll": "pll", "afp": "afp", "antifreeze protein": "afp", "pva": "pva"}


def _unsupported(readout, reason):
    return {"supported": False, "source_type": "unsupported",
            "readout": readout, "reason": reason,
            "note": "This layer is not modelled; no value is fabricated."}


class CryoCell:
    """A persistent virtual cell: observe (read-only), perturb (state-changing),
    save / restore / clone / branch, with full provenance."""

    # ------------------------------------------------------------------ init
    def __init__(self, cell_type="hMSC", initial_temperature_c=22.0, **context):
        self.id = "cell_" + uuid.uuid4().hex[:12]
        self.cell_type = cell_type
        self.context = dict(context)
        self.warnings = []
        base = CELL_LINES.get(cell_type)
        if base is None:
            self.warnings.append(
                f"'{cell_type}' is not a calibrated context; using hMSC baseline "
                f"parameters. Treat all outputs as extrapolation.")
            base = CELL_LINES["hMSC"]
        self.P = Params(**{k: v for k, v in base.items() if k in Params.__annotations__})
        self.P.T_add = float(initial_temperature_c)
        self.history = []
        self.parent = None            # (parent_id, from_state_name)
        self._states = {}             # name -> params dict snapshot
        self._result = None           # cached (out, res, events)
        self._sim_hash = None
        self._record("create", dict(cell_type=cell_type,
                                    initial_temperature_c=initial_temperature_c, **context),
                     None, self.state_hash())

    # -------------------------------------------------------------- internals
    def _params_dict(self):
        d = asdict(self.P)
        d["ko"] = dict(d.get("ko") or {})
        return d

    def state_hash(self):
        blob = json.dumps(self._params_dict(), sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()[:16]

    def _ensure(self):
        h = self.state_hash()
        if self._result is None or self._sim_hash != h:
            self._result = simulate(self.P)
            self._sim_hash = h
        return self._result

    def _cal(self):
        return "calibrated" if self.cell_type in CELL_LINES else "literature_prior"

    def _record(self, op, params, before, after, extra=None):
        rec = dict(operation=op, parameters=params, state_before=before,
                   state_after=after,
                   timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   model_version=MODEL_VERSION, parameter_set_version=PARAM_SET_VERSION,
                   pathway_database_versions={}, reference_keys=["REFERENCES.md"],
                   uncertainty={"basis": "literature-motivated where noted; "
                                         "anchored to hMSC/DMSO calibration"},
                   validity_flags=list(self.warnings),
                   calibration_status=self._cal())
        if extra:
            rec.update(extra)
        self.history.append(rec)
        return rec

    def _transition_summary(self):
        _, R, _ = self._ensure()
        return {"S_imm_pct": round(R["S_imm"] * 100, 1),
                "S_24h_pct": round(R["S_24"] * 100, 1),
                "S_72h_pct": round(R["S_72"] * 100, 1),
                "P_iif_pct": round(R["P_iif"] * 100, 1),
                "min_volume_frac": round(R["minV"], 3)}

    def _perturb(self, op, updates, params_log):
        before = self.state_hash()
        for k, v in updates.items():
            if k == "ko":
                cur = dict(self.P.ko or {}); cur.update(v); self.P.ko = cur
            elif k in Params.__annotations__:
                setattr(self.P, k, v)
        after = self.state_hash()
        self._ensure()
        rec = self._record(op, params_log, before, after)
        return {"operation": op, "state_before": before, "state_after": after,
                "changed": before != after, "summary": self._transition_summary(),
                "provenance": rec}

    # =================================================================== #
    #  FIVE PERTURBATION METHODS (state-changing)                         #
    # =================================================================== #
    def apply_chemical_perturbation(self, agents=None, solution="DPBS",
                                    temperature_c=None, exposure_duration_min=None,
                                    addition_schedule="single_step", **kw):
        upd = {}
        for a in (agents or []):
            name = (a.get("name") or "").strip().lower()
            comp = a.get("compartment", "extracellular")
            conc = a.get("concentration_percent_vv",
                         a.get("concentration_percent_wv", a.get("concentration")))
            if comp == "extracellular_additive" or name in _ADD_NAMES and name not in _CPA_NAMES:
                key = _ADD_NAMES.get(name)
                if key:
                    upd["additive"] = key
                    if conc is not None:
                        upd["add_conc"] = float(conc)
            elif name in _CPA_NAMES:
                upd["cpa_key"] = _CPA_NAMES[name]
                if conc is not None:
                    upd["conc_pct"] = float(conc)
            else:
                self.warnings.append(f"Unknown agent '{name}' ignored (not in CPA/additive library).")
        if temperature_c is not None:
            upd["T_add"] = float(temperature_c)
        if exposure_duration_min is not None:
            upd["hold_min"] = float(exposure_duration_min)
        upd["add_steps"] = 1 if addition_schedule == "single_step" else 3
        return self._perturb("apply_chemical_perturbation", upd,
                             dict(agents=agents, solution=solution,
                                  temperature_c=temperature_c,
                                  exposure_duration_min=exposure_duration_min,
                                  addition_schedule=addition_schedule))

    def apply_physical_protocol(self, protocol):
        p = dict(protocol or {})
        m = {"cooling_rate_c_per_min": "CR", "warming_rate_c_per_min": "WR",
             "end_temperature_c": "T_store", "storage_temperature_c": "T_store",
             "storage_time_days": "days", "dilution_temperature_c": "T_dil",
             "recovery_temperature_c": "T_recover", "recovery_duration_h": "recover_h",
             "dilution": "dilution"}
        upd = {m[k]: v for k, v in p.items() if k in m}
        nuc = p.get("nucleation")
        if isinstance(nuc, dict) and nuc.get("temperature_c") is not None:
            upd["T_seed"] = float(nuc["temperature_c"])
        return self._perturb("apply_physical_protocol", upd, p)

    def apply_environmental_perturbation(self, temperature_c=None, recovery_temperature_c=None,
                                         recovery_duration_h=None, **kw):
        upd = {}
        if temperature_c is not None:
            upd["T_recover"] = float(temperature_c)
        if recovery_temperature_c is not None:
            upd["T_recover"] = float(recovery_temperature_c)
        if recovery_duration_h is not None:
            upd["recover_h"] = float(recovery_duration_h)
        return self._perturb("apply_environmental_perturbation", upd,
                             dict(temperature_c=temperature_c,
                                  recovery_temperature_c=recovery_temperature_c,
                                  recovery_duration_h=recovery_duration_h))

    def apply_genetic_perturbation(self, target=None, mode="knockout", **kw):
        # supported knockouts map to the engine's mechanism-isolation flags
        supported = {"aqp", "cpaperm", "serca", "mpt", "oxphos", "sterol", "msc",
                     "mt", "actin", "intf", "gel", "nucl", "iri", "rock", "piezo",
                     "casp", "calpain", "fa", "yap"}
        upd = {}
        if target in supported and mode in ("knockout", "knockdown"):
            upd["ko"] = {target: True}
        else:
            self.warnings.append(
                f"Genetic target '{target}' mode '{mode}' not supported as a "
                f"mechanism flag; no change applied.")
        return self._perturb("apply_genetic_perturbation", upd,
                             dict(target=target, mode=mode))

    def apply_biological_perturbation(self, adhesion=None, cell_line_phenotype=None, **kw):
        upd = {}
        if adhesion in ("suspension", "adherent", "sheared", "spheroid"):
            upd["adhesion"] = adhesion
        for k in ("apop_resist", "anoikis_resist", "antioxidant", "glycolytic"):
            if cell_line_phenotype and k in cell_line_phenotype:
                upd[k] = float(cell_line_phenotype[k])
        return self._perturb("apply_biological_perturbation", upd,
                             dict(adhesion=adhesion, cell_line_phenotype=cell_line_phenotype))

    # =================================================================== #
    #  FIFTEEN SIGNATURE METHODS (read-only)                              #
    # =================================================================== #
    def _sim(self, source_type="mechanistically_simulated"):
        _, R, _ = self._ensure()
        return R, {"source_type": source_type, "calibration_status": self._cal(),
                   "model_version": MODEL_VERSION,
                   "applicability": ("calibration anchor" if self.cell_type == "hMSC"
                                     else "extrapolation from hMSC/DMSO calibration"),
                   "warnings": list(self.warnings)}

    # --- cryobiology signatures (real, from the calibrated engine) ---
    def get_thermal_ice_state(self, **kw):
        R, meta = self._sim()
        return {**meta, "temperature_c": self.P.T_store,
                "extracellular_ice_fraction": R["fIce"],
                "P_intracellular_ice": R["P_iif"],
                "intracellular_ice_amount": R["iifAmount"],
                "P_mitochondrial_ice": R["P_mito"],
                "supercooling_peak_C": R["supercoolPeak"],
                "grain_size_um": R["grain"], "grain_max_um": R["grainMax"],
                "glass_fraction": R["glass"], "recrystallisation_D": R["D_recry_ice"],
                "intercellular_propagation_coupling": R["junction"],
                "units": {"temperature": "C", "grain": "um"}}

    def get_osmotic_mechanical_state(self, **kw):
        R, meta = self._sim()
        return {**meta, "min_volume_frac": R["minV"], "max_volume_frac": R["maxV"],
                "membrane_tension": R["tension"],
                "cytoskeletal_load": {"microtubule_compression_sigma": R["sigMT"],
                                      "intermediate_filament_tension_sigma": R["sigIF"],
                                      "intermediate_filament_integrity": R["intf"],
                                      "microtubule_integrity": R["mt"],
                                      "actin_integrity": R["actin"]},
                "ice_squeeze_min_channel_per_cell": R["squeezeMin"],
                "osmotic_injury_D": R["D_osm"], "mechanical_injury_D": R["D_mech"],
                "units": {"volume": "fraction of isotonic", "sigma": "dimensionless load"}}

    def get_membrane_transport_state(self, **kw):
        R, meta = self._sim()
        cpa = CPAS[self.P.cpa_key]
        return {**meta, "Lp_um_per_min_per_atm_25C": self.P.lp,
                "Ps_um_per_s_25C": self.P.ps,
                "reflection_coefficient_sigma": cpa["sigma"],
                "Ps_activation_energy_kJ_per_mol": cpa["EaPs"] / 1e3,
                "membrane_area_scaling": "V^(2/3) with folded reserve",
                "membrane_phase": {"gel_fraction": R["memb"]["gel"],
                                   "fluid_fraction": R["memb"]["fluid"],
                                   "Tm_C": R["memb"]["Tm"]},
                "defect_pore_index": R["memb"].get("domainLeak"),
                "membrane_integrity_D": R["D_mem"],
                "formulation": "Kedem-Katchalsky 3-parameter (see REFERENCES.md)"}

    def get_organelle_state(self, **kw):
        R, meta = self._sim()
        return {**meta,
                "mitochondria": {"matrix_volume_rel": R["Vmito_final"],
                                 "membrane_potential_dPsi": R["dPsi"],
                                 "permeability_transition_frac": R["mptFrac"],
                                 "matrix_CPA_M": R["matrixCPA"],
                                 "ice_probability": R["P_mito"],
                                 "treatment": "active osmometer + nucleator (REFERENCES.md)"},
                "nucleus": {"chromatin_condensation": R["chromCond"],
                            "treatment": "passively coupled to cytosol (documented choice)"},
                "endoplasmic_reticulum": {"lumenal_Ca_remaining": R["caER"],
                                          "cytosolic_Ca": R["caCyt"]},
                "cytoskeleton": self.get_osmotic_mechanical_state()["cytoskeletal_load"]}

    def get_cell_morphology(self, modality="schematic", **kw):
        R, meta = self._sim()
        img_ok = modality in ("schematic",)
        return {**meta, "modality": modality,
                "image_supported": img_ok,
                "features": {"relative_volume_min": R["minV"],
                             "relative_volume_max": R["maxV"],
                             "membrane_slack_proxy": max(0.0, R["memb"]["APL"] - 1),
                             "intracellular_ice_area_fraction_proxy": R["P_iif"] * R["iifAmount"],
                             "gel_membrane_fraction": R["memb"]["gel"]},
                "note": ("Quantitative features only; brightfield/confocal/Raman image "
                         "channels are the GUI's render, not a validated microscopy model."
                         if not img_ok else "schematic feature set")}

    def get_stress_response(self, **kw):
        R, meta = self._sim()
        return {**meta, "oxidative_ROS_peak": R["rosPeak"],
                "delayed_death_projection_D": R["D_delayed"],
                "cytosolic_calcium": R["caCyt"], "mechanical_tension": R["tension"],
                "CPA_toxicity_D": R["D_tox"], "protein_native_fraction": R["prot"],
                "membrane_stress_D": R["D_mem"],
                "cold_shock": "captured via WLF mobility + cold-labile cytoskeleton",
                "ER_DNA_inflammatory_ISR": _unsupported(
                    "ER/DNA-damage/inflammatory/integrated-stress-response",
                    "not modelled as discrete pathways in the current engine")}

    def get_metabolic_state(self, **kw):
        R, meta = self._sim()
        return {**meta, "source_type": "mechanistically_simulated (partial)",
                "ATP_rel_resting": R["atp"], "ATP_min_rel": R["atpMin"],
                "mitochondrial_membrane_potential": R["dPsi"],
                "oxphos_glycolysis_split": {"glycolytic_bias": self.P.glycolytic},
                "redox_ROS": R["rosPeak"],
                "full_metabolic_network": _unsupported(
                    "glycolysis/TCA/PPP/lipid/amino-acid/nucleotide fluxes",
                    "only ATP balance, ΔΨm and ROS are modelled; full flux network not implemented")}

    def get_cell_fate(self, recovery_duration_h=24, **kw):
        R, meta = self._sim()
        if recovery_duration_h <= 30:
            surv, tp = R["S_24"], "24h (modelled)"
        else:
            surv, tp = R["S_72"], "72h (delayed-death projection)"
        return {**meta, "recovery_time_point": tp,
                "membrane_integrity_survival": R["S_imm"],
                "survival": surv, "S_24h": R["S_24"], "S_72h": R["S_72"],
                "functional_recovery": R["F_rec"],
                "live_fraction": R["liveFrac"], "apoptotic_fraction": R["apopFrac"],
                "necrotic_fraction": R["necrFrac"],
                "endpoints_modelled": ["viability", "apoptosis", "necrosis",
                                       "delayed-onset death", "functional recovery"],
                "endpoints_not_modelled": ["ferroptosis", "senescence",
                                           "proliferation competence"]}

    # --- general-biology signatures (explicitly unsupported, never fabricated) ---
    def get_genomic_state(self, **kw):
        return _unsupported("genomic state", "genome layer not implemented")

    def get_epigenomic_state(self, **kw):
        return _unsupported("epigenomic state", "epigenome layer not implemented")

    def get_transcriptomic_state(self, **kw):
        return _unsupported("transcriptomic state", "transcriptome layer not implemented")

    def get_protein_state(self, **kw):
        # honest partial: the engine has a bulk native-protein fraction only
        R, meta = self._sim()
        return {**meta, "source_type": "mechanistically_simulated (bulk only)",
                "native_protein_fraction": R["prot"],
                "per_protein_abundance_structure_localisation_PTM": _unsupported(
                    "per-protein proteome", "only a bulk native-fraction is modelled")}

    def get_protein_interactions(self, **kw):
        return _unsupported("protein interactions", "interaction layer not implemented")

    def get_protein_interaction_complexes(self, **kw):
        return _unsupported("protein complexes", "structure/complex layer not implemented")

    def get_pathway_state(self, **kw):
        return _unsupported("pathway state",
                            "versioned pathway knowledge layer not implemented; "
                            "no pathway activity is fabricated")

    # =================================================================== #
    #  PERSISTENCE: save / restore / clone / branch / compare / export   #
    # =================================================================== #
    def save_state(self, name):
        self._states[name] = self._params_dict()
        self._record("save_state", dict(name=name), self.state_hash(), self.state_hash())
        return name

    def restore(self, name):
        snap = self._states.get(name)
        if snap is None:
            raise KeyError(f"no saved state '{name}'")
        before = self.state_hash()
        self.P = Params(**{k: v for k, v in snap.items() if k in Params.__annotations__})
        self.P.ko = dict(snap.get("ko") or {})
        self._result = None
        self._record("restore", dict(name=name), before, self.state_hash())
        return self.state_hash()

    def clone(self, from_state=None):
        child = CryoCell.__new__(CryoCell)
        child.id = "cell_" + uuid.uuid4().hex[:12]
        child.cell_type = self.cell_type
        child.context = dict(self.context)
        child.warnings = list(self.warnings)
        src = self._states.get(from_state, self._params_dict()) if from_state else self._params_dict()
        child.P = Params(**{k: v for k, v in src.items() if k in Params.__annotations__})
        child.P.ko = dict(src.get("ko") or {})
        child.history = []
        child.parent = (self.id, from_state)
        child._states = dict(self._states)
        child._result = None
        child._sim_hash = None
        child._record("clone", dict(parent=self.id, from_state=from_state),
                      None, child.state_hash())
        return child

    branch = clone   # alias: a branch is a clone with a recorded parent

    def compare(self, other):
        a, _ = self._sim()
        b, _ = other._sim()
        keys = ["S_imm", "S_24", "S_72", "F_rec", "P_iif", "D_osm", "D_tox",
                "D_mem", "D_iif", "minV", "dPsi", "rosPeak"]
        return {"self": self.id, "other": other.id,
                "delta": {k: round((a[k] - b[k]) * (100 if k.startswith(("S_", "D_", "P_")) else 1), 3)
                          for k in keys}}

    def replay(self):
        return list(self.history)

    def export(self):
        _, R, _ = self._ensure()
        return {"cell_id": self.id, "cell_type": self.cell_type,
                "model_version": MODEL_VERSION, "parameter_set_version": PARAM_SET_VERSION,
                "parameters": self._params_dict(), "state_hash": self.state_hash(),
                "outcome": self._transition_summary(), "history": self.history,
                "saved_states": list(self._states), "warnings": self.warnings,
                "reference_file": "REFERENCES.md",
                "calibration_status": self._cal()}

    # =================================================================== #
    #  Backwards-compatibility aliases (§14)                             #
    # =================================================================== #
    def apply_cpa(self, name="dmso", concentration=10.0, temperature_c=22.0, **kw):
        return self.apply_chemical_perturbation(
            agents=[{"name": name, "concentration": concentration}],
            temperature_c=temperature_c, **kw)

    def apply_thermal_protocol(self, cooling_rate_c_per_min=1.0, **kw):
        return self.apply_physical_protocol({"cooling_rate_c_per_min": cooling_rate_c_per_min, **kw})

    def get_intracellular_ice(self):
        return self.get_thermal_ice_state()["P_intracellular_ice"]

    def get_mitochondrial_state(self):
        return self.get_organelle_state()["mitochondria"]
