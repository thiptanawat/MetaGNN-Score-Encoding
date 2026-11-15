"""Shared loaders and comparison helpers for the identifier audit of the second reconstruction.

Paths are resolved from the repository root (three levels above this file) or from environment
variables that name an input outside the repository:

  MAIN_MODEL_JSON    the Recon3D JSON of the main analysis (default: data/raw/Recon3D.json, the raw input
                     that data/raw/README.md describes)
  SECOND_ROOT_JSON   the converted Human-GEM model written by second/build_root.py, that is
                     <second root>/data/raw/Recon3D.json (required)
  HUMANGEM_MODEL_DIR the model/ directory of a Human-GEM checkout at tag v2.1.1 (required; the audit reads
                     metabolites.tsv and reactions.tsv, and check_yml.py reads Human-GEM.yml)

The build records are read from second/ and the main medium module and metabolite map from src/ and
manifests/ of the repository. RDKit is used only to turn Human-GEM InChI strings into InChIKeys.
"""
import ast, csv, json, os, re
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]


def _required(name):
    if name not in os.environ:
        raise SystemExit("set the environment variable %s (see the docstring of audit_lib.py)" % name)
    return Path(os.environ[name])


MAIN_JSON = Path(os.environ.get("MAIN_MODEL_JSON", REPO / "data/raw/Recon3D.json"))
SECOND_JSON = _required("SECOND_ROOT_JSON")
HG_DIR = _required("HUMANGEM_MODEL_DIR")
SEC = REPO / "second"
MEDIUM_PY = REPO / "src/medium_v2.py"
MAPFILE = REPO / "manifests/metabolite_map.tsv"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def as_list(v):
    if v is None or v == "":
        return []
    if isinstance(v, list):
        return [str(x) for x in v]
    return [str(v)]


def split_ids(s):
    return [t.strip() for t in re.split(r"[;,]", s or "") if t.strip()]


def norm_hmdb(x):
    m = re.fullmatch(r"HMDB0*(\d+)", x.strip())
    return "HMDB%07d" % int(m.group(1)) if m else x.strip()


def parse_formula(f):
    if not f:
        return None
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def heavy(f):
    d = parse_formula(f)
    if d is None:
        return None
    return {k: v for k, v in d.items() if k != "H"}


def heavy_str(d):
    return "" if d is None else "".join("%s%d" % (k, d[k]) for k in sorted(d))


def hcount(f):
    d = parse_formula(f)
    return None if d is None else d.get("H", 0)


def inchikey(inchi):
    from rdkit.Chem import inchi as ci
    if not inchi or not inchi.startswith("InChI="):
        return ""
    try:
        return ci.InchiToInchiKey(inchi) or ""
    except Exception:
        return ""


def key_parts(k):
    """(skeleton block, stereo/isotope block, protonation letter) of a standard InChIKey."""
    m = re.fullmatch(r"([A-Z]{14})-([A-Z]{8})([A-Z])([A-Z])-([A-Z])", k or "")
    return (m.group(1), m.group(2), m.group(5)) if m else None


def norm_name(n):
    n = (n or "").lower()
    n = re.sub(r"\b(d|l|dl|r|s|sn|alpha|beta)-", "", n)
    return re.sub(r"[^a-z0-9]+", " ", n).strip()


def name_tokens(n):
    return re.findall(r"[a-z0-9]+", norm_name(n))


def stereo_tags(n):
    return sorted(set(re.findall(r"(?:^|[\s(\[,])((?:D|L|DL|R|S|sn)(?:-|\)))", n or "")))


def medium_lists():
    tree = ast.parse(MEDIUM_PY.read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") in ("BASAL", "CULTURE", "SERUM"):
            out[node.targets[0].id] = ast.literal_eval(node.value)
    return out


class Models:
    def __init__(self):
        self.main = json.load(open(MAIN_JSON))
        self.second = json.load(open(SECOND_JSON))
        self.mm = {m["id"]: m for m in self.main["metabolites"]}
        self.mr = {r["id"]: r for r in self.main["reactions"]}
        self.sm = {m["id"]: m for m in self.second["metabolites"]}
        self.sr = {r["id"]: r for r in self.second["reactions"]}
        self.hgt = {r["mets"]: r for r in read_tsv(HG_DIR / "metabolites.tsv")}
        self.hgr = {r["rxns"]: r for r in read_tsv(HG_DIR / "reactions.tsv")}
        self.exmap = read_tsv(SEC / "exchange_map.tsv")
        self.ex_by_id = {r["pipeline_id"]: r for r in self.exmap}
        self.ex_by_met = {}
        for r in self.exmap:
            self.ex_by_met.setdefault(r["humangem_metabolite"], []).append(r)
        self._key = {}

    def hg_key(self, mid):
        if mid not in self._key:
            self._key[mid] = inchikey(self.hgt.get(mid, {}).get("metInChI", ""))
        return self._key[mid]

    def main_info(self, mid):
        m = self.mm.get(mid)
        if not m:
            return None
        a = m.get("annotation", {})
        return dict(id=mid, name=m["name"], formula=m.get("formula") or "", charge=m.get("charge"),
                    kegg=as_list(a.get("kegg.compound")), chebi=as_list(a.get("chebi")),
                    hmdb=[norm_hmdb(x) for x in as_list(a.get("hmdb"))], keys=as_list(a.get("inchi_key")),
                    mnx=as_list(a.get("metanetx.chemical")), bigg=as_list(a.get("bigg.metabolite")))

    def hg_info(self, mid):
        m = self.sm.get(mid)
        t = self.hgt.get(mid)
        if not m or not t:
            return None
        return dict(id=mid, name=m["name"], formula=m.get("formula") or "", charge=m.get("charge"),
                    kegg=split_ids(t["metKEGGID"]), chebi=split_ids(t["metChEBIID"]),
                    hmdb=[norm_hmdb(x) for x in split_ids(t["metHMDBID"])], keys=[self.hg_key(mid)] if self.hg_key(mid) else [],
                    inchi=t["metInChI"], smiles=t["metSmiles"], mnx=split_ids(t["metMetaNetXID"]),
                    bigg=split_ids(t["metBiGGID"]))


def classify(a, b):
    """Compare main-model info a and Human-GEM info b. Returns (class, detail dict)."""
    d = {}
    ha, hb = heavy(a["formula"]), heavy(b["formula"])
    d["heavy_main"], d["heavy_hg"] = heavy_str(ha), heavy_str(hb)
    d["heavy_equal"] = "" if (ha is None or hb is None) else str(ha == hb)
    d["H_main"], d["H_hg"] = hcount(a["formula"]), hcount(b["formula"])
    d["charge_main"], d["charge_hg"] = a["charge"], b["charge"]
    d["chebi_overlap"] = ";".join(sorted(set(a["chebi"]) & set(b["chebi"])))
    d["hmdb_overlap"] = ";".join(sorted(set(a["hmdb"]) & set(b["hmdb"])))
    d["kegg_overlap_main_hg"] = ";".join(sorted(set(a["kegg"]) & set(b["kegg"])))
    d["name_norm_equal"] = str(norm_name(a["name"]) == norm_name(b["name"]))
    if ha is not None and hb is not None and not ha and not hb:
        # no heavy atoms (the proton): InChI layers carry only charge and isotope information
        ok = (d["kegg_overlap_main_hg"] or d["hmdb_overlap"]) and a["formula"] == b["formula"]
        d["protonation_differs"] = str(a["charge"] != b["charge"])
        d["note"] = "no heavy atoms; InChI layers differ only in charge/isotope layers (Human-GEM InChI: %s)" % b.get("inchi", "")
        return ("SAME_PROTONATION_DIFFERS" if ok else "CANNOT_TELL"), d
    pa = [key_parts(k) for k in a["keys"] if key_parts(k)]
    pb = [key_parts(k) for k in b["keys"] if key_parts(k)]
    d["key_main"] = ";".join(a["keys"]); d["key_hg"] = ";".join(b["keys"])
    d["stereo_tags_main"] = "".join(stereo_tags(a["name"])); d["stereo_tags_hg"] = "".join(stereo_tags(b["name"]))
    d["protonation_differs"] = ""
    if ha is not None and hb is not None and ha != hb:
        return "DIFFERENT_COMPOUND", d
    if pa and pb:
        best = None
        for x in pa:
            for y in pb:
                if x[0] == y[0]:
                    best = (x, y); break
            if best: break
        if not best:
            # Connectivity hash differs although the heavy atoms agree. Treat as a tautomer or ionization drawing
            # difference only when the stereo hash and at least one database identifier agree.
            ids = d["chebi_overlap"] or d["hmdb_overlap"] or d["kegg_overlap_main_hg"]
            if ids and any(x[1] == y[1] for x in pa for y in pb):
                d["protonation_differs"] = "True"
                d["note"] = "InChIKey connectivity block differs (tautomer or ionization drawing); stereo hash and identifiers agree"
                return "SAME_PROTONATION_DIFFERS", d
            kdis = a["kegg"] and b["kegg"] and not (set(a["kegg"]) & set(b["kegg"]))
            if kdis:
                return "DIFFERENT_COMPOUND", d
            d["missing"] = "main-model InChI layers (only a hashed key is stored) and a shared identifier"
            d["note"] = "InChIKey connectivity block differs although heavy atoms and formula agree and no KEGG identifiers conflict; tautomer drawing versus constitutional isomer cannot be decided"
            return "CANNOT_TELL", d
        x, y = best
        prot = (x[2] != y[2]) or (a["charge"] != b["charge"]) or (d["H_main"] != d["H_hg"])
        d["protonation_differs"] = str(bool(prot))
        UNS = "UHFFFAOY"
        if x[1] != y[1]:
            if x[1] == UNS or y[1] == UNS:
                return "STEREO_UNSPECIFIC_ON_ONE_SIDE", d
            return "STEREO_MISMATCH", d
        return ("SAME_PROTONATION_DIFFERS" if prot else "SAME"), d
    miss = []
    if not pa: miss.append("main-model inchi_key")
    if not pb: miss.append("Human-GEM InChI")
    d["missing"] = "+".join(miss)
    return "CANNOT_TELL", d

def build_hg_ex(M):
    """Extracellular Human-GEM metabolites that have a one-metabolite boundary reaction in the second root."""
    hg_ex = {}
    for r in M.second["reactions"]:
        if len(r["metabolites"]) == 1:
            mid = next(iter(r["metabolites"]))
            if M.sm[mid]["compartment"] == "e":
                hg_ex.setdefault(mid, []).append(r["id"])
    return hg_ex, {mid: M.hg_info(mid) for mid in hg_ex}


def alternatives(M, hg_ex, hg_e_info, mapped_id, main_inf, src_kegg):
    h = hg_e_info[mapped_id]
    sk = set(h["kegg"]) | set(main_inf["kegg"] if main_inf else []) | set(split_ids(src_kegg))
    mk = [key_parts(k) for k in (h["keys"] + (main_inf["keys"] if main_inf else [])) if key_parts(k)]
    skel = {p[0] for p in mk}
    hn = norm_name(h["name"]); ht = name_tokens(h["name"])
    base = re.sub(r"__[A-Za-z]+$", "", (h["bigg"] or [""])[0])
    res = []
    for mid, inf in hg_e_info.items():
        if mid == mapped_id: continue
        tags_s, tags_w = [], []
        if set(inf["kegg"]) & sk: tags_s.append("kegg:" + ",".join(sorted(set(inf["kegg"]) & sk)))
        if set(inf["chebi"]) & set(h["chebi"]): tags_s.append("chebi")
        if main_inf and set(inf["chebi"]) & set(main_inf["chebi"]): tags_s.append("chebi_main")
        ik = {key_parts(k)[0] for k in inf["keys"] if key_parts(k)}
        if ik & skel: tags_s.append("inchikey_skeleton")
        exact = {(p[0], p[1]) for p in (key_parts(k) for k in inf["keys"]) if p}
        if exact & {(p[0], p[1]) for p in mk}: tags_s.append("inchikey_skeleton_and_stereo_equal")
        if set(inf["mnx"]) & set(h["mnx"]): tags_s.append("metanetx")
        if norm_name(inf["name"]) == hn and hn: tags_s.append("name_normalized")
        ib = re.sub(r"__[A-Za-z]+$", "", (inf["bigg"] or [""])[0])
        if base and ib == base and "__" in ((inf["bigg"] or [""])[0] + (h["bigg"] or [""])[0]): tags_s.append("bigg_stem")
        it = name_tokens(inf["name"])
        def contig(a, b):
            return len(a) > 0 and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
        if not tags_s and len("".join(ht)) >= 5 and (contig(ht, it) or contig(it, ht)) and len("".join(it)) >= 5:
            tags_w.append("name_token_overlap")
        if tags_s or tags_w:
            res.append((mid, inf, ";".join(tags_s), ";".join(tags_w)))
    return res

