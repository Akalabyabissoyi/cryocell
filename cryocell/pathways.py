"""
Reactome pathway membership for the cell-stress pathways (data-driven, replacing
hand-mapped gene counts). Gene lists are the human participants of each Reactome
pathway, fetched from the Reactome ContentService.

Source: Reactome (https://reactome.org), retrieved 2026-08-26.
Licence: Reactome data is released under CC0 1.0 (public domain) and is bundled
here in full. Reactome: Milacic et al. 2024, Nucleic Acids Res 52:D672.
"""
REACTOME_RETRIEVED = "2026-08-26"
REACTOME_LICENSE = "CC0 1.0 (public domain)"

REACTOME = {
    "oxid": dict(id="R-HSA-3299685", name='Detoxification of reactive oxygen species', count=37,
               genes=['AQP8', 'ATOX1', 'ATP7A', 'CAT', 'CCS', 'CYBA', 'CYBB', 'CYCS', 'ERO1A', 'GPX1', 'GPX2', 'GPX3', 'GPX5', 'GPX6', 'GPX7', 'GPX8', 'GSR', 'GSTP1', 'NCF1', 'NCF2', 'NCF4', 'NOX4', 'NOX5', 'NUDT2', 'P4HB', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX5', 'PRDX6', 'SOD1', 'SOD2', 'SOD3', 'TXN', 'TXN2', 'TXNRD1', 'TXNRD2']),
    "apop": dict(id="R-HSA-109606", name='Intrinsic pathway for apoptosis', count=55,
               genes=['AKT1', 'AKT2', 'AKT3', 'APAF1', 'APIP', 'AVEN', 'BAD', 'BAK1', 'BAX', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BID', 'BMF', 'C1QBP', 'CARD8', 'CASP3', 'CASP7', 'CASP8', 'CASP9', 'CDKN2A', 'CYCS', 'DIABLO', 'DYNLL1', 'DYNLL2', 'E2F1', 'GSDMD', 'GSDME', 'GZMB', 'MAPK1', 'MAPK3', 'MAPK8', 'NMT1', 'PMAIP1', 'PPP1R13B', 'PPP3CC', 'PPP3R1', 'SEPTIN4', 'SFN', 'STAT3', 'TFDP1', 'TFDP2', 'TP53', 'TP53BP2', 'TP63', 'TP73', 'UACA', 'XIAP', 'YWHAB', 'YWHAE', 'YWHAG', 'YWHAH', 'YWHAQ', 'YWHAZ']),
    "prot": dict(id="R-HSA-381119", name='Unfolded protein response (UPR)', count=93,
               genes=['ACADVL', 'ADD1', 'ARFGAP1', 'ASNS', 'ATF3', 'ATF4', 'ATF6', 'ATF6B', 'ATP6V0D1', 'CALR', 'CCL2', 'CEBPB', 'CEBPG', 'CREB3', 'CREB3L1', 'CREB3L2', 'CREB3L3', 'CREB3L4', 'CREBRF', 'CTDSP2', 'CUL7', 'CXCL8', 'CXXC1', 'DCP2', 'DCSTAMP', 'DCTN1', 'DDIT3', 'DDX11', 'DIS3', 'DNAJB11', 'DNAJB9', 'DNAJC3', 'EDEM1', 'EIF2AK3', 'EIF2S1', 'EIF2S2', 'EIF2S3', 'ERN1', 'EXOSC1', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC5', 'EXOSC6', 'EXOSC7', 'EXOSC8', 'EXOSC9', 'EXTL1', 'EXTL2', 'EXTL3', 'FKBP14', 'GFPT1', 'GOSR2', 'GSK3A', 'HDGF', 'HERPUD1', 'HSP90B1', 'HSPA5', 'HYOU1', 'IGFBP1', 'KDELR3', 'KHSRP', 'KLHDC3', 'LMNA', 'MBTPS1', 'MBTPS2', 'MYDGF', 'NFYA', 'NFYB', 'NFYC', 'PARN', 'PDIA5', 'PDIA6', 'PLA2G4B', 'PPP2R5B', 'PREB', 'SEC31A', 'SERP1', 'SHC1', 'SRPRA', 'SRPRB', 'SSR1', 'SULT1A3', 'SYVN1', 'TATDN2', 'TLN1', 'TPP1', 'TSPYL2', 'WFS1', 'WIPI1', 'XBP1', 'YIF1A', 'ZBTB17']),
    "mech": dict(id="R-HSA-5627117", name='RHO GTPases activate ROCKs', count=19,
               genes=['CFL1', 'LIMK1', 'LIMK2', 'MYH10', 'MYH11', 'MYH14', 'MYH9', 'MYL12B', 'MYL6', 'MYL9', 'PAK1', 'PPP1CB', 'PPP1R12A', 'PPP1R12B', 'RHOA', 'RHOB', 'RHOC', 'ROCK1', 'ROCK2']),
}
