"""Download the MSK-IMPACT 50K release files used by the examples (Zenodo record 18445440; CC BY-NC-ND 4.0).
"""
import argparse
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

RECORD = "18445440"
BASE = f"https://zenodo.org/records/{RECORD}/files/"
FILES = {   # name: (bytes, md5) as listed by the Zenodo record
    "clinical":  ("data_clinical_sample.txt",     9_864_808,     "627040924f7f88ee571b88749b830cbd"),
    "segments":  ("data_facets_gene_level.txt",   3_534_793_826, "775f174f3442a7cba7e151548bed1d2d"),
    "mutations": ("data_mutations_ccf.maf.txt",   620_354_633,   "1249ca252de6938148679edf230177cc"),
}

def md5sum(path: Path, block: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()

def download(url: str, dest: Path, size: int) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "wgdlife/1.0"})
    with urllib.request.urlopen(req) as r, open(tmp, "wb") as f:
        done = 0
        for chunk in iter(lambda: r.read(1 << 22), b""):
            f.write(chunk)
            done += len(chunk)
            print(f"\r  {dest.name}: {done / 1e6:,.0f} / {size / 1e6:,.0f} MB", end="", flush=True)
    print()
    tmp.replace(dest)

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default=None, help="comma list of keys to fetch (clinical, segments, mutations); default all")
    ap.add_argument("--data", default=os.environ.get("WGDLIFE_DATA", "data"), help="data directory (default $WGDLIFE_DATA or data/)")
    a = ap.parse_args()
    dest_dir = Path(a.data) / "msk_impact_50k"
    dest_dir.mkdir(parents=True, exist_ok=True)
    keys = [k for k in a.only.split(",") if k] if a.only else list(FILES)
    for k in keys:
        name, size, md5 = FILES[k]
        dest = dest_dir / name
        if dest.exists() and dest.stat().st_size == size and md5sum(dest) == md5:
            print(f"  {name}: present, MD5 verified")
            continue
        print(f"  {name}: downloading {size / 1e6:,.0f} MB from Zenodo record {RECORD}")
        download(BASE + name + "?download=1", dest, size)
        got = md5sum(dest)
        if got != md5:
            dest.unlink()
            sys.exit(f"MD5 mismatch for {name}: {got} != {md5} (file removed; retry)")
        print(f"  {name}: MD5 verified")
    print(f"done: {dest_dir}")

if __name__ == "__main__":
    main()
