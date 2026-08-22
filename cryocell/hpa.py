"""
Subcellular compartment reference data from the Human Protein Atlas Cell Atlas.

Numbers are quoted as published at https://www.proteinatlas.org/humanproteome/subcellular
(retrieved 16 August 2026). The Cell Atlas annotates 49 subcellular locations and
provides experimental localisation for 13,603 genes (67% of protein-coding genes),
rising to 17,062 genes (85%) when predictions are included.

`multi` is the fraction of that compartment's proteome that is also found elsewhere —
the single most useful number here, because it says how much of a compartment is
private to it. Mitochondria are the most self-contained major compartment (53%
shared); the cytosol is the least (83% shared).
"""

HPA_TOTALS = {
    "genes_experimental": 13603,
    "genes_with_prediction": 17062,
    "pct_experimental": 67,
    "pct_with_prediction": 85,
    "locations_annotated": 49,
}

# key -> (display name, genes, % of protein-coding genes, fraction multi-localised, note)
HPA = {
    "nucleoplasm":  ("Nucleoplasm",          6880, 34, 0.70,
                     "Largest annotated compartment. Transcription, RNA processing, DNA repair."),
    "cytosol":      ("Cytosol",              5341, 26, 0.83,
                     "Most promiscuous compartment; shares heavily with nucleus and plasma membrane."),
    "vesicles":     ("Vesicles",             2488, 13, 0.72,
                     "Secretory, transport and endocytic vesicles."),
    "plasma_mem":   ("Plasma membrane",      2467, 12, 0.84,
                     "Signalling, transport, adhesion. The compartment CPAs act on first."),
    "nucleoli":     ("Nucleoli",             1478,  7, 0.66,
                     "rRNA synthesis and processing. Fibrillar centre 325, rim 153."),
    "golgi":        ("Golgi apparatus",      1279,  6, 0.79,
                     "Glycosylation and sorting."),
    "mitochondria": ("Mitochondria",         1132,  6, 0.53,
                     "The most self-contained major compartment: only 53% shared."),
    "centrosome":   ("Centrosome",            716,  4, 0.90,
                     "Microtubule organising centre. Centriolar satellites 228."),
    "microtubules": ("Microtubules",          356,  3, 0.89,
                     "Cold-labile. Depolymerise below ~12 C, releasing GEF-H1 onto RhoA."),
    "er":           ("Endoplasmic reticulum", 576,  3, 0.60,
                     "Protein folding and the cell's main Ca2+ store."),
    "actin":        ("Actin filaments",       377,  2, 0.85,
                     "Includes focal adhesion sites (148 genes)."),
    "focal_adh":    ("Focal adhesion sites",  148,  1, 0.85,
                     "Integrin-based mechanosensors; vinculin-positive."),
    "nuclear_mem":  ("Nuclear membrane",      295,  1, 0.72,
                     "Lined by nuclear lamins; pores make it leaky to small solutes."),
    "lipid_drop":   ("Lipid droplets",         40,  0, 0.60, "Neutral lipid stores."),
    "peroxisomes":  ("Peroxisomes",            24,  0, 0.60, "Beta-oxidation, ROS handling."),
    "lysosomes":    ("Lysosomes",              20,  0, 0.60, "Acidic degradative compartment."),
    "endosomes":    ("Endosomes",              16,  0, 0.60, "Endocytic sorting."),
    "interm_fil":   ("Intermediate filaments", None, None, None,
                     "Not separately quantified in the Cell Atlas; lamins are the nuclear class."),
}

# How much of the cell's volume each compartment occupies, and how many discrete
# objects of it to draw. Volume fractions are conventional mammalian-cell values,
# NOT from the HPA (the HPA counts proteins, not volume) -- flagged as such.
GEOMETRY = {
    #  key            vol frac   count   source of the volume fraction
    "nucleus":       (0.100,  1,  "conventional mammalian value"),
    "mitochondria":  (0.055, 22,  "conventional; 4-8% in MSC"),
    "er":            (0.080,  1,  "conventional (rough + smooth)"),
    "golgi":         (0.020,  1,  "conventional"),
    "lysosomes":     (0.010, 12,  "conventional"),
    "peroxisomes":   (0.005,  9,  "conventional"),
    "endosomes":     (0.008, 10,  "conventional"),
    "lipid_drop":    (0.004,  6,  "conventional"),
    "vesicles":      (0.010, 34,  "conventional"),
    "centrosome":    (0.001,  1,  "one MTOC per interphase cell"),
}

def summary_line(key):
    e = HPA.get(key)
    if not e:
        return None
    name, genes, pct, multi, note = e
    if genes is None:
        return f"{name} — not separately quantified in the Cell Atlas"
    return (f"{name}: {genes:,} genes ({pct}% of protein-coding), "
            f"{multi*100:.0f}% also localise elsewhere")
