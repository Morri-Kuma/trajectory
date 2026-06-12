# Deletion Proposal — awaiting your approval

**Date:** 2026-06-09 · **Verification pass:** 2026-06-09
**Status:** ⛔ **PROPOSAL ONLY — nothing has been deleted.** No file will be
removed or renamed without your explicit approval.

This itemises cleanup candidates from the audit (`docs/reports/project_audit_report.md` §6),
**now corrected with a verification pass** (grep for references, dir diffs,
submodule check). All listed paths are **gitignored** (`data/`, `logs/`,
`models/`, and the named `benchmark` intermediates), so removal does not rewrite
history — it frees local disk and reduces "what is current" confusion.

> **Correction vs the first draft:** the GSE178325 duplicate is **`rna_seq/`**,
> *not* `rna_seq_10x/`. The build script
> `scripts/build_gse178325_raw_full_gene_input.py` reads `.../rna_seq_10x`, and
> that folder also holds the QC outputs (`qc_gse178325/`). `rna_seq/` is the
> older, more-deeply-nested layout that **no script references**. Deleting
> `rna_seq_10x/` would have broken the GSE178325 pipeline — keep it.

**Legend — Category:** duplicate / obsolete / broken / temporary / unsafe-naming.

---

## Tier A — high-value cleanup (large, verified non-essential)

| # | Path | Category | Verified finding | Recommend | Risk |
|---|---|---|---|---|---|
| A1 | `benchmark/archive/` (**~32 GB**) | obsolete | `grep` of `result_manifest.yaml`, configs, reports, registry → **no current asset references it**. Contains only `README.md` + `rerun_backups/` (already gitignored). | **Delete** (or cold-archive externally first) | Low |
| A2 | `data/gse178325_human/rna_seq/` (**~2.6 GB**) | duplicate | The **older nested** copy (`GSM*/<sample>/...`). `rna_seq_10x/` is the flattened canonical one the builder uses + has QC. **No script references plain `rna_seq/`.** | **Delete `rna_seq/`; KEEP `rna_seq_10x/`** | Low (data re-downloadable from GEO if ever needed) |

## Tier B — small, safe, or clearly stray

| # | Path | Category | Verified finding | Recommend | Risk |
|---|---|---|---|---|---|
| B1 | `logs/.git/` | broken | **No `.gitmodules`** → it is a stray nested repo, not a submodule | Delete the inner `.git` | Low |
| B2 | `benchmark/ground_truth/_expanded_official_backup/` (~38 MB) | obsolete | Backup of superseded providers | Delete (keep one cold copy) | Low |
| B3 | `benchmark/results/_invalid_lineage_label_space/` (~372 KB) | broken | Name marks it invalid | Delete after a final `grep` for references | Low |
| B4 | `benchmark/results/smoke/` (0 B) | temporary | Empty | Delete | None |
| B5 | `**/__pycache__/`, `.pytest_cache/` | temporary | Regenerated build caches | Delete | None |

## Tier C — verify/transform, not blind delete

| # | Path | Category | Verified finding | Recommend | Risk |
|---|---|---|---|---|---|
| C1 | `data/trajectory_raw_data_audit_report_v2.md` (265 lines) vs `docs/reports/…_v2.md` (290 lines) | duplicate | They differ; the **`docs/reports/` copy is longer and canonical** (referenced by the repo structure). The `data/` copy is the stale one. | Move/delete the `data/` copy; keep `docs/reports/` | Low–medium |
| C2 | `data/gse230659(human/` (literal `(`) | unsafe-naming | Parenthesis breaks naive shell/glob | **Rename** → `data/gse230659_human/` and update references (not a deletion) | Medium — update any hard-coded path |

---

## Recommended approval set (my advice)

- **Approve A1, A2, and all of Tier B now** — verified safe, ~**34.6 GB** reclaimed (A1 32 GB + A2 2.6 GB + small B items), zero loss of tracked code, current results, or benchmark inputs.
- **A1:** I suggest a quick external `tar` copy first *only if* you ever revisit old reruns; otherwise straight delete.
- **C1/C2:** do them as a separate, careful step (diff-then-keep for C1; rename-and-update-paths for C2) — not in the bulk delete.

## What I will do on approval
For each approved item: small Tier B items I will remove and report freed space;
the large Tier A items I will (your choice) either delete or print a `tar`/move
command for cold storage; Tier C I will diff/rename without deleting anything
unconfirmed. **Until you approve, nothing is removed or renamed.**
