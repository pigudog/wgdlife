"""Cohort definitions.
"""

# The 12 carcinoma types of the article: epithelial carcinomas with >= 10 WGD+ tumours in each of the TCGA
# primaries, the MSK-IMPACT primaries and MSK-MET.
CANCER_TYPES_12 = (
    "Non-Small Cell Lung Cancer",
    "Breast Cancer",
    "Bladder Cancer",
    "Endometrial Cancer",
    "Prostate Cancer",
    "Head and Neck Cancer",
    "Colorectal Cancer",
    "Esophagogastric Cancer",
    "Ovarian Cancer",
    "Renal Cell Carcinoma",
    "Hepatobiliary Cancer",
    "Pancreatic Cancer",
)

# Current HGNC symbol (MSK-IMPACT target list) -> the older spelling used in the MSK-IMPACT release MAF and in the
# TCGA MC3 / PanCanAtlas tables.
TARGET_OLD_SYMBOLS = {
    "ABRAXAS1": "FAM175A", "CCNQ": "FAM58A", "COP1": "RFWD2", "ELOC": "TCEB1", "EZHIP": "CXorf67",
    "H1-2": "HIST1H1C", "H2BC5": "HIST1H2BD", "H3-3A": "H3F3A", "H3-3B": "H3F3B", "H3-4": "HIST3H3",
    "H3-5": "H3F3C", "H3C1": "HIST1H3A", "H3C2": "HIST1H3B", "H3C3": "HIST1H3C", "H3C4": "HIST1H3D",
    "H3C6": "HIST1H3E", "H3C7": "HIST1H3F", "H3C8": "HIST1H3G", "H3C10": "HIST1H3H", "H3C11": "HIST1H3I",
    "H3C12": "HIST1H3J", "H3C13": "HIST2H3D", "KMT5A": "SETD8", "MRE11": "MRE11A", "NSD2": "WHSC1",
    "NSD3": "WHSC1L1", "PAK5": "PAK7", "PRKN": "PARK2", "TENT5C": "FAM46C",
}

OLD2NEW_SYMBOLS = {old: new for new, old in TARGET_OLD_SYMBOLS.items()}

# MAF Variant_Classification values that are not non-synonymous (the mutation filter of every SNV analysis)
SILENT_CLASSES = frozenset({"Silent", "3'UTR", "5'UTR", "Intron", "IGR", "5'Flank", "3'Flank", "RNA", "Splice_Region"})

def current_symbol(s):
    """Map a Series of gene symbols from the release spelling to the current HGNC symbol (others unchanged)."""
    return s.map(OLD2NEW_SYMBOLS).fillna(s)
