#!/usr/bin/env bash
# Prepare a clean, citable release commit + tag for the scTimeBench conditions benchmark.
#
# WHY THIS SCRIPT: the repo's git index was wedged by a stale .git/index.lock left by a
# crashed git operation on 2026-06-01, which is the root cause of the ~800 mixed
# changes/deletions `git status` showed. The Cowork sandbox mount cannot delete files
# inside .git, so this final step must run on your own machine, where you have full
# permissions. It is safe and reversible (see UNDO at the bottom).
#
# Usage:   bash scripts/prepare_release.sh [TAG]
# Example: bash scripts/prepare_release.sh v1.0.0-rc1
set -euo pipefail

TAG="${1:-v1.0.0-rc1}"
MAXBYTES=$((5*1024*1024))   # abort if any staged file exceeds 5 MB

cd "$(git rev-parse --show-toplevel 2>/dev/null || echo .)"
echo "==> repo: $(pwd)"

# 1) clear the stale lock (root cause of the messy index)
if [ -f .git/index.lock ]; then
  echo "==> removing stale .git/index.lock ($(stat -c '%y' .git/index.lock 2>/dev/null || true))"
  rm -f .git/index.lock
fi

# 2) sanitize the index: unstage everything, working tree untouched
echo "==> git reset (mixed; working tree untouched)"
git reset -q

# 3) re-apply .gitignore to already-tracked files, then stage the true working state
echo "==> re-evaluating .gitignore across all paths"
git rm -r --cached -q . >/dev/null
git add -A

# 4) SAFETY: refuse to commit if any large file slipped into the index
echo "==> checking staged blob sizes (limit 5 MB)"
BIG=$(git diff --cached --name-only | while read -r f; do
  [ -f "$f" ] || continue
  sz=$(git cat-file -s ":$f" 2>/dev/null || echo 0)
  [ "${sz:-0}" -gt "$MAXBYTES" ] && printf '%s\t%s\n' "$sz" "$f"
done || true)
if [ -n "$BIG" ]; then
  echo "!! Large files staged (add them to .gitignore or 'git rm --cached', then re-run):"
  echo "$BIG" | sort -rn | awk -F'\t' '{printf "   %.1f MB  %s\n", $1/1048576, $2}'
  exit 1
fi

# 5) show what will be committed (summary only)
echo "==> staged summary:"
git diff --cached --shortstat
TESTS=$(grep -rh "def test_" tests/*.py 2>/dev/null | wc -l | tr -d ' ')

# 6) commit + annotated tag
git commit -q -m "Release ${TAG}: reviewer-revision r2 (OOD gate integrated, planned-extension framing, multi-seed provenance, docs sync)

- Geometric OOD gate integrated into src/pipeline.py (per-state + accept/reject decision); k/quantile sensitivity sweep added.
- scIMF/PI-SDE/Squidiff scoped strictly as planned extensions (VendoringRequired); excluded from all results.
- Conclusions framed as a silver-reference benchmark (not biological ground truth).
- Multi-seed CI labeled a derived artifact; manifest + config SHA-256 checksums + run logs; HPC rerun recipe documented.
- Embedding-evaluator name collision resolved; threshold-free lineage-robustness module added.
- Docs synced (README/REPRODUCIBILITY/TODO): ${TESTS} unit tests; large regenerable artifacts git-ignored."
git tag -a "${TAG}" -m "scTimeBench conditions benchmark — ${TAG} (reviewer revision r2)"

echo
echo "==> DONE. Created commit + tag ${TAG}"
git --no-pager log --oneline -1
echo "Tag: $(git describe --tags --exact-match 2>/dev/null || echo "${TAG}")"
echo
echo "Next (optional): push with tags ->  git push origin HEAD --tags"
echo "Archive a citable snapshot -> create a GitHub release from ${TAG}, or upload to Zenodo/OSF."
echo
echo "UNDO if needed:  git tag -d ${TAG} && git reset --soft HEAD~1"
