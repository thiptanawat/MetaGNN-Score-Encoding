"""Cross-check name, formula and charge of the metabolites used in the audit between Human-GEM.yml (read line by line)
and the converted COBRA JSON of the second root. Usage: python -I check_yml.py OUTFILE"""
import sys, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_lib import *
M = Models()
need = set()
for r in read_tsv(SEC / "target_map.tsv"): need.add(r["humangem_metabolite"])
hg_ex, _ = build_hg_ex(M)
need |= set(hg_ex)
cur, rec, n, bad = None, {}, 0, []
with open(HG_DIR / "Human-GEM.yml") as fh:
    for line in fh:
        m = re.match(r"\s+- id: (MAM\d+[a-z])\s*$", line)
        if m:
            cur = m.group(1); rec = {"id": cur}; continue
        if cur:
            m = re.match(r"\s+- (name|formula|charge|compartment): (.*)$", line)
            if m:
                rec[m.group(1)] = m.group(2).strip().strip('"')
                if m.group(1) == "charge" or "formula" in rec and "charge" in rec:
                    pass
            if re.match(r"\s+- metFrom", line) or line.startswith("- reactions") or line.startswith("- genes"):
                if cur in need and cur in M.sm:
                    j = M.sm[cur]; n += 1
                    ok = (rec.get("name") == j["name"] and rec.get("formula", "") == (j.get("formula") or "")
                          and float(rec.get("charge", 0)) == float(j.get("charge") or 0))
                    if not ok: bad.append((cur, rec, j["name"], j.get("formula"), j.get("charge")))
                cur = None
open(sys.argv[1], "w").write("compared\t%d\nmismatching\t%d\n" % (n, len(bad)) + "".join(repr(b) + "\n" for b in bad[:50]))
print(n, len(bad), bad[:3])
