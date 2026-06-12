#!/bin/bash
# Download + arrange the two integration datasets from NCBI GEO so the v2.1
# preprocess builders work first-try. Run where there is internet + storage —
# RECOMMENDED on Shirokane (a login node is fine; this is just I/O, no compute):
#
#     bash scripts/download_new_datasets.sh
#
# It populates:
#     data/gse298212/<sample>/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}
#     data/gse218855/<the processed .h5ad and/or 10x files>
#
# Then run jobs/shirokane/run_v21_new_datasets_preprocess.sh.
#
# If GEO's file names differ from what the organizer expects, it still downloads
# everything into data/<acc>/_raw/ and prints an `ls` so we can finish the layout.
set -uo pipefail
ROOT="${TRAJ_PROJECT_ROOT:-$(pwd)}"
cd "$ROOT"
mkdir -p data/gse298212 data/gse218855

FTP () {  # $1 = accession -> echoes the GEO suppl URL
  local acc="$1"; local pre="${acc:0:$((${#acc}-3))}nnn"
  echo "https://ftp.ncbi.nlm.nih.gov/geo/series/${pre}/${acc}/suppl/"
}

dl_suppl () {  # $1 = accession, $2 = dest dir : mirror the whole suppl/ folder
  local acc="$1" dest="$2"; local url; url="$(FTP "$acc")"
  echo ">>> downloading ${acc} supplementary from ${url}"
  mkdir -p "${dest}/_raw"
  # -r recursive, -np no-parent, -nd flatten, -e robots=off; retry + continue.
  wget -r -np -nd -e robots=off --tries=3 --continue --timeout=60 \
       -P "${dest}/_raw" "${url}" 2>&1 | tail -5 || {
    echo "!! wget failed for ${acc}. On Shirokane try: module load wget or use curl."
    return 1; }
  echo ">>> ${acc} raw files:"; ls -lh "${dest}/_raw" | tail -20
}

# ── GSE298212 (human blood; 10x per sample) ──────────────────────────────────
dl_suppl GSE298212 data/gse298212
echo ">>> extracting any RAW tar for GSE298212"
for t in data/gse298212/_raw/*RAW.tar; do [ -f "$t" ] && tar -xf "$t" -C data/gse298212/_raw/; done
# Organize GSM*_<sample>_{barcodes,features/genes,matrix} into per-sample 10x dirs.
python3 - <<'PY'
import os, re, glob, shutil
raw = "data/gse298212/_raw"
files = glob.glob(os.path.join(raw, "*"))
# group by the sample token between the GSM id and the 10x role
role = {"barcodes":"barcodes.tsv.gz","features":"features.tsv.gz","genes":"features.tsv.gz","matrix":"matrix.mtx.gz"}
samples = {}
for f in files:
    b = os.path.basename(f).lower()
    m = re.search(r"(gsm\d+)[_\.]?(.*?)[_\.](barcodes|features|genes|matrix)", b)
    if not m: continue
    sample = (m.group(2) or m.group(1)).strip("_.") or m.group(1)
    samples.setdefault(sample, {})[role[m.group(3)]] = f
made = 0
for sample, parts in samples.items():
    if len(parts) >= 2:
        d = os.path.join("data/gse298212", sample); os.makedirs(d, exist_ok=True)
        for dst, src in parts.items():
            shutil.copy(src, os.path.join(d, dst)); made += 1
print(f"[organize] arranged {len(samples)} sample dirs ({made} files) under data/gse298212/")
print("[organize] sample dirs:", sorted([s for s in samples if len(samples[s])>=2]))
PY
echo ">>> data/gse298212 layout:"; ls -R data/gse298212 | grep -v _raw | head -40

# ── GSE218855 (mouse; processed .h5ad expected) ──────────────────────────────
dl_suppl GSE218855 data/gse218855
echo ">>> extracting any RAW tar for GSE218855"
for t in data/gse218855/_raw/*RAW.tar; do [ -f "$t" ] && tar -xf "$t" -C data/gse218855/_raw/; done
# Surface any .h5ad to data/gse218855/ where the builder looks.
for h in data/gse218855/_raw/*.h5ad data/gse218855/_raw/*.h5ad.gz; do
  [ -f "$h" ] || continue
  case "$h" in *.gz) gunzip -kf "$h"; h="${h%.gz}";; esac
  cp -f "$h" "data/gse218855/$(basename "$h")"
done
echo ">>> data/gse218855 contents:"; ls -lh data/gse218855 | grep -v _raw | head

echo
echo "============================================================"
echo "DONE. Verify the layout above, then:"
echo "  qsub jobs/shirokane/run_v21_new_datasets_preprocess.sh"
echo "If GSE298212 sample dirs are not named PBMC_EPC/S1D1/S1D3/S1D6/S1D8,"
echo "edit SAMPLE_TO_ABS_DAY in scripts/build_gse298212_chemical_input.py."
echo "If GSE218855 shipped only 10x (no .h5ad), paste 'ls -R data/gse218855/_raw'"
echo "and I'll point the builder at it."
echo "============================================================"
