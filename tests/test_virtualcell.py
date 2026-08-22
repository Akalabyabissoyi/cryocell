"""Invariant tests for the persistent virtual-cell layer (spec §15, §27).
Run: python tests/test_virtualcell.py    (from the repo root, with deps on path)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cryocell.virtualcell import CryoCell


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    assert cond, name


def main():
    print("Signature immutability (read-only methods must not change state):")
    cell = CryoCell(cell_type="hMSC", initial_temperature_c=22.0)
    for meth in ["get_thermal_ice_state", "get_osmotic_mechanical_state",
                 "get_membrane_transport_state", "get_organelle_state",
                 "get_cell_morphology", "get_stress_response", "get_metabolic_state",
                 "get_cell_fate", "get_pathway_state", "get_protein_interactions"]:
        before = cell.state_hash()
        getattr(cell, meth)()
        check(f"{meth} leaves state_hash unchanged", cell.state_hash() == before)

    print("\nPerturbations create a new state + history entry:")
    before_h, before_hist = cell.state_hash(), len(cell.history)
    r = cell.apply_chemical_perturbation(
        agents=[{"name": "DMSO", "concentration": 10.0, "compartment": "extracellular"}],
        temperature_c=22.0, exposure_duration_min=10.0)
    check("chemical perturbation changed state_hash", cell.state_hash() != before_h)
    check("history grew", len(cell.history) > before_hist)
    check("transition summary present", "S_24h_pct" in r["summary"])

    print("\nPhysical protocol runs and yields survival:")
    cell.apply_physical_protocol({"cooling_rate_c_per_min": 1.0, "end_temperature_c": -196,
                                  "warming_rate_c_per_min": 400})
    fate = cell.get_cell_fate(recovery_duration_h=24)
    check("S_24 in a sane range 0..1", 0.0 <= fate["survival"] <= 1.0)
    check("hMSC 1C/min optimum near calibration (30-42%)", 0.30 <= fate["S_24h"] <= 0.42)

    print("\nClone / branch preserves the parent:")
    parent_h = cell.state_hash()
    child = cell.clone()
    child.apply_physical_protocol({"cooling_rate_c_per_min": 100.0})
    check("parent state unchanged after child perturbation", cell.state_hash() == parent_h)
    check("child records parent id", child.parent[0] == cell.id)
    check("child differs from parent", child.state_hash() != cell.state_hash())

    print("\nUnsupported biology layers are honest (never fabricated):")
    for meth in ["get_genomic_state", "get_transcriptomic_state",
                 "get_pathway_state", "get_protein_interactions"]:
        out = getattr(cell, meth)()
        check(f"{meth} returns supported=False", out.get("supported") is False)

    print("\nProvenance + export completeness:")
    exp = cell.export()
    for k in ("cell_id", "model_version", "parameter_set_version", "state_hash",
              "history", "calibration_status", "reference_file"):
        check(f"export contains '{k}'", k in exp)
    rec = cell.history[-1]
    for k in ("operation", "state_before", "state_after", "timestamp",
              "calibration_status", "reference_keys"):
        check(f"history entry has '{k}'", k in rec)

    print("\nBackwards-compatibility aliases work:")
    c2 = CryoCell("HeLa")
    c2.apply_cpa("dmso", 10.0)
    check("apply_cpa alias runs", c2.get_intracellular_ice() is not None)
    check("get_mitochondrial_state alias runs", "membrane_potential_dPsi" in c2.get_mitochondrial_state())
    check("HeLa flagged/calibrated appropriately", c2._cal() == "calibrated")

    print("\nALL INVARIANT TESTS PASSED.")


def test_invariants():
    """pytest entry point: runs every invariant check via assert."""
    main()


if __name__ == "__main__":
    main()
