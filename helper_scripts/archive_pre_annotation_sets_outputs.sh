#!/bin/bash
# Title: Archive post-08 outputs written before versioned annotation sets
# Date: 2026-11-22
# Summary: One-off move of the unversioned script 08/09/10/helper/poster
# outputs (old `cell_type` column, no annotation stamp) into
#   data/xenium/_pre_annotation_sets_2026-11-22/<same relative path>
# so that only versioned outputs (.../<cohort>/vN/...) remain in place.
# Versioned folders (v1, v2, ...) are never moved. Nothing is deleted.
#
# Run from the repository root (locally or on Phoenix). Dry run by default:
#   bash helper_scripts/archive_pre_annotation_sets_outputs.sh          # list only
#   bash helper_scripts/archive_pre_annotation_sets_outputs.sh --apply  # move
# Undo: move each archived path back to the same path without the archive prefix.

set -euo pipefail

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1
COHORT=min15_area95_selected_resolutions
ARCHIVE=data/xenium/_pre_annotation_sets_2026-11-22
[[ -d data/xenium ]] || { echo "Run from the repository root (no data/xenium here)"; exit 1; }

targets=()
# Old per-sample script 08 objects: data/xenium/processed/vbct/<sample>/cell_type_annotation/
for dir in data/xenium/processed/vbct/*/cell_type_annotation; do
    [[ "$dir" == data/xenium/processed/vbct/cell_type_annotation ]] && continue
    [[ -d "$dir" ]] && targets+=("$dir")
done
# Unversioned entries inside each step's cohort folder (everything except vN/).
for kind in output processed; do
    for step in cell_type_annotation cell_type_proportions spatial_analysis poster; do
        base="data/xenium/$kind/vbct/$step/$COHORT"
        [[ -d "$base" ]] || continue
        for entry in "$base"/*; do
            [[ -e "$entry" ]] || continue
            [[ "$(basename "$entry")" =~ ^v[0-9]+$ ]] && continue
            targets+=("$entry")
        done
    done
done

if [[ ${#targets[@]} -eq 0 ]]; then
    echo "Nothing to archive."
    exit 0
fi

echo "Archive: $ARCHIVE"
for path in "${targets[@]}"; do
    dest="$ARCHIVE/${path#data/xenium/}"
    printf '%8s  %s -> %s\n' "$(du -sh "$path" | cut -f1)" "$path" "$dest"
    if [[ -e "$dest" ]]; then
        echo "  destination exists; refusing to overwrite" >&2
        exit 1
    fi
done

if [[ $APPLY -eq 0 ]]; then
    echo "Dry run: nothing moved. Re-run with --apply to move ${#targets[@]} paths."
    exit 0
fi

for path in "${targets[@]}"; do
    dest="$ARCHIVE/${path#data/xenium/}"
    mkdir -p "$(dirname "$dest")"
    mv "$path" "$dest"
done

cat > "$ARCHIVE/README.md" <<'EOF'
# Outputs from before versioned annotation sets (archived 2026-11-22)

Script 08, 09, 10, spatial-helper and poster outputs written before annotation
sets were introduced (branch `feature/annotation-sets`). They use the labels of
annotation set v1 (`config/annotations/vbct/min15_area95_selected_resolutions/v1/`)
but the old single `cell_type` column, carry no annotation stamp, and sit at
unversioned paths, so current code does not read them. Paths below this folder
mirror their original location under `data/xenium/`.

- `processed/vbct/<sample>/cell_type_annotation/`: old script 08 annotated H5ADs
  (read by the conference-figure config `config/conference_figures/`, repointed here).
- `{processed,output}/vbct/<step>/min15_area95_selected_resolutions/`: old tables
  and figures (script 08 counts, script 09 proportions, script 10 results and
  enrichment tables, cross-sample spatial helpers, poster helpers).

v1 was regenerated on Phoenix by script 08 (job 16939096) with identical counts.
Context: `context/tasks/cell_type_annotation.md`.
EOF
echo "Moved ${#targets[@]} paths; wrote $ARCHIVE/README.md"
