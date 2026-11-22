"""
Title: Shared code for the Xenium BANKSY workflow (script 08 onwards)
Date: 2026-11-22
Summary: Modules imported by the numbered pipeline scripts, the helpers in
helper_scripts/, and tools/. Not a pipeline step itself.

    cohorts          cohort configs (clustering facts) and --input-root/--output-root
    annotation_sets  versioned cell-type annotation sets: loading, validation,
                     derived output paths; command line for new/check/diff/release
                     (python -m xenium_workflow.annotation_sets ...)
    spatial_tables   script 10 neighbourhood enrichment table format

Moved here from the repository root on 2026-11-22 (workflow_config.py was
renamed cohorts.py).
"""
