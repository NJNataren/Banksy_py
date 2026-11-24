"""
Title: Dead-citation checker for project documentation
Created: 2026-11-20
Last modified: 2026-11-24
Summary: Find file and directory paths cited in agent/handoff documentation
(backtick spans and fenced code blocks) and report those that no longer exist,
so agents do not follow instructions that point at renamed or moved scripts,
configs or folders. Dead citations in history sections, on lines that say
something was renamed/superseded/archived, or anywhere in a doc kept under an
`archive/` folder, are warnings; dead citations in current sections are errors
(exit status 1).
Inputs: Markdown files (default: AGENTS.md, README_*.md, context/**/*.md,
hpc/slurm/README.md in the repository).
Outputs: Report on stdout (or JSON with --json). Writes nothing.
Usage: python tools/check_doc_citations.py [files ...] [--include-data] [--all]
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GLOBS = ["AGENTS.md", "README_*.md", "context/**/*.md", "hpc/slurm/README.md"]

# Extensions that mark a token as a file reference even without a slash.
CODE_EXTS = (".py", ".json", ".sl", ".R", ".r", ".md", ".sh", ".ipynb", ".yml", ".yaml", ".ini", ".toml")
DATA_EXTS = (".h5ad", ".csv", ".tsv", ".png", ".pdf", ".svg", ".pkl", ".gz", ".parquet", ".zarr", ".log")

# Generated outputs that legitimately differ between machines. They are only
# checked with --include-data, and then reported as "not local", not dead.
DATA_PREFIXES = ("data/", "figures/", "logs/", "hpc/raw_data", "hpc/output")

# HPC paths cannot be checked from here without connecting (which is forbidden).
REMOTE_PREFIXES = ("/scratchdata1/", "/hpcfs/", "/uofaresstor/")

# Directories skipped when looking up a bare filename anywhere in the repo.
INDEX_SKIP = {".git", "data", "figures", "logs", "renv", "__pycache__", ".ipynb_checkpoints", "node_modules"}

# Headings that mark dated or historical material.
HISTORY_HEADING_RE = re.compile(
    r"history|session|update|archive|supersed|changelog|previous|legacy|log\b|handoff|snapshot|renamed", re.I)

# How many lines back a cited folder still counts as the base for a bare
# filename (e.g. "Directory: `x/`" followed by "Files:" and a list of names).
RECENT_DIR_WINDOW = 15
# Line wording that marks an intentional mention of an old name.
HISTORY_LINE_RE = re.compile(
    r"supersed|archiv|renamed|moved (from|to)|formerly|previously|legacy|deprecated|removed|"
    r"no longer|old(er)? |was \S+ (until|before)|replaced", re.I)

# Wording that marks a path as proposed rather than existing (design notes).
PLANNED_RE = re.compile(r"propos|planned|possible|would |could |future|suggest|candidate|draft|not yet|to be created", re.I)

INLINE_CODE_RE = re.compile(r"`([^`\n]+)`")

# Local mount points for HPC storage; see the HPC rule in AGENTS.md.
MOUNT_RE = re.compile(r"^(~|/home/[^/]+)/(mnt|gadi_mount)(/|$)")

# First path segments that exist in the repository; a slash token starting
# with anything else is treated as prose rather than a repo path.
TOP_LEVEL = {p.name for p in REPO_ROOT.iterdir()} | {"Banksy_py"}
TOKEN_STRIP = "\"'(),;:[]{}<>!?"


def build_index():
    """Map each filename in the repository (minus data/output dirs) to its paths."""
    index = {}
    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dirnames[:] = [d for d in dirnames if d not in INDEX_SKIP]
        for name in filenames + dirnames:
            index.setdefault(name, []).append(Path(dirpath, name).relative_to(REPO_ROOT))
    return index


def iter_code_text(lines):
    """Yield (line_number, text, heading, context) for code spans and code-block lines.

    `context` is the prose line introducing a fenced block (e.g. "Possible
    package layout:"), or the line itself for inline spans.
    """
    heading, in_fence, intro, last_prose = "", False, "", ""
    for number, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("```"):
            in_fence = not in_fence
            intro = last_prose if in_fence else ""
            continue
        if in_fence:
            yield number, line, heading, intro
            continue
        if stripped.startswith("#"):
            heading = stripped.lstrip("#").strip()
        if stripped:
            last_prose = stripped
        for span in INLINE_CODE_RE.findall(line):
            yield number, span, heading, line


def candidate_tokens(text):
    """Split a code span or command into tokens that look like file paths."""
    for raw in re.split(r"\s+", text):
        # `--export=ALL,CONFIG_DIR=path,X=1` and `--config=path` forms: check
        # each comma-separated value after its last `=`.
        parts = raw.split(",") if "=" in raw else [raw]
        for part in parts:
            token = part.rsplit("=", 1)[-1].strip(TOKEN_STRIP).rstrip(".")
            if not token or token.startswith("-") or "://" in token or "@" in token:
                continue
            if "/" in token or token.endswith(CODE_EXTS + DATA_EXTS):
                yield token


def classify(token, doc_path, index, include_data):
    """Return (status, detail) for one cited path.

    Status is one of ok, DEAD, ARCHIVED, REMOTE, SKIP, NOTLOCAL.
    """
    # Placeholders, globs and shell variables cannot be resolved literally.
    if any(ch in token for ch in "<>{}*$") or "..." in token or token.startswith("."):
        return "SKIP", "template or variable"
    if token.startswith(REMOTE_PREFIXES):
        return "REMOTE", "HPC path, not checked"
    if token.startswith("/tmp/"):
        return "SKIP", "temporary file"
    if "path/to/" in token:
        return "SKIP", "template or variable"
    if MOUNT_RE.match(token):
        # A local view of the HPC: only checkable while the mount is active.
        mount_path = Path(token).expanduser()
        return ("ok", "") if mount_path.exists() else ("MOUNT", "mount path; mount not active or path absent")
    if token.startswith("~"):
        return ("ok", "") if Path(token).expanduser().exists() else ("DEAD", "home path not found")
    if token.startswith("/"):
        return ("ok", "") if Path(token).exists() else ("DEAD", "absolute path not found")
    rel = token[len("Banksy_py/"):] if token.startswith("Banksy_py/") else token
    is_data = rel.startswith(DATA_PREFIXES) or (rel.endswith(DATA_EXTS) and "/" not in rel)
    if is_data and not include_data:
        return "SKIP", "data/output path (use --include-data)"
    for base in (REPO_ROOT, doc_path.parent, REPO_ROOT.parent):
        if (base / rel).exists():
            return "ok", ""
    bare = rel.rstrip("/")
    if "/" not in bare:
        # A bare filename: accept it if it exists anywhere outside archive/.
        hits = index.get(bare, [])
        live = [h for h in hits if h.parts[0] != "archive"]
        if live:
            return "ok", str(live[0])
        if hits:
            return "ARCHIVED", f"only in {hits[0]}"
        if rel.endswith("/"):
            return "SKIP", "folder name without location (e.g. an output subfolder)"
        if not re.search(r"[_\-\d]", Path(bare).stem):
            return "SKIP", "generic placeholder name"
    elif not rel.endswith(CODE_EXTS + DATA_EXTS) and not rel.endswith("/") \
            and rel.split("/")[0] not in TOP_LEVEL:
        # Slash-separated words in prose (QC/inspection, PTMT/PC35) or branch names.
        return "SKIP", "not a repository path"
    elif rel.split("/")[0] not in TOP_LEVEL and not rel.startswith(DATA_PREFIXES):
        # Relative to some output folder the doc describes, not to the repo.
        return "SKIP", "relative to an unnamed folder"
    if is_data:
        return "NOTLOCAL", "generated output not present on this machine"
    return "DEAD", "not found"


def check_file(doc_path, index, include_data):
    """Return a list of citation results for one Markdown file."""
    results, seen = [], set()
    recent_dirs = []  # (line_number, resolved folder) for folders cited recently
    lines = doc_path.read_text(errors="replace").splitlines()
    # Archived docs (any `archive/` folder, e.g. context/tasks/dotplot/archive/) are
    # history throughout, whatever their original headings say.
    in_archive = "archive" in doc_path.relative_to(REPO_ROOT).parts[:-1]
    for number, text, heading, context in iter_code_text(lines):
        planned = bool(PLANNED_RE.search(heading)) or bool(PLANNED_RE.search(context))
        # An inline span that is one existing path containing spaces
        # (e.g. "poster v1 (1).pdf") is a single citation, not several tokens.
        whole = text.strip()
        if whole.startswith(("/", "~")) and " " in whole and Path(whole).expanduser().exists():
            continue
        for token in candidate_tokens(text):
            if (number, token) in seen:
                continue
            seen.add((number, token))
            # Judge history wording without the token itself, so a filename
            # such as `..._archive_labels.json` does not excuse itself.
            line_text = lines[number - 1].replace(token, "")
            history = (in_archive or bool(HISTORY_HEADING_RE.search(heading))
                       or bool(HISTORY_LINE_RE.search(line_text)))
            status, detail = classify(token, doc_path, index, include_data)
            folder = REPO_ROOT / token.removeprefix("Banksy_py/")
            if token.endswith("/") and folder.is_dir():
                recent_dirs.append((number, folder))
            elif status == "DEAD" and "/" not in token:
                # A bare filename listed under a recently cited folder.
                for _, folder in reversed([d for d in recent_dirs if number - d[0] <= RECENT_DIR_WINDOW]):
                    if (folder / token).exists():
                        status, detail = "ok", str(folder / token)
                        break
            # Old names in history sections are expected; keep them visible but
            # do not fail the check for them.
            if status in ("DEAD", "ARCHIVED") and history:
                status = "WARN-" + status
            elif status == "DEAD" and planned:
                status, detail = "PLANNED", "proposed in a design note; does not exist yet"
            results.append({"file": str(doc_path.relative_to(REPO_ROOT)), "line": number, "token": token,
                            "status": status, "detail": detail, "heading": heading})
    return results


def parse_args():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description="Report dead file/directory citations in project docs (read-only).")
    parser.add_argument("files", nargs="*", type=Path, help="Markdown files to check (default: project docs)")
    parser.add_argument("--include-data", action="store_true",
                        help="Also check data/, figures/ and logs/ paths (reported as NOTLOCAL if absent)")
    parser.add_argument("--all", action="store_true", help="Also list citations that resolved or were skipped")
    parser.add_argument("--json", action="store_true", help="Print results as JSON")
    return parser.parse_args()


def main():
    """Check documentation files and exit 1 if any current-section citation is dead."""
    args = parse_args()
    if args.files:
        files = [f.resolve() for f in args.files]
    else:
        files = sorted({p.resolve() for pattern in DEFAULT_GLOBS for p in REPO_ROOT.glob(pattern)})
    files = [f for f in files if f.is_file()]
    index = build_index()
    results = [r for f in files for r in check_file(f, index, args.include_data)]

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        shown = ("DEAD", "ARCHIVED", "WARN-DEAD", "WARN-ARCHIVED", "NOTLOCAL")
        for r in results:
            if args.all or r["status"] in shown:
                detail = f"  ({r['detail']})" if r["detail"] else ""
                print(f"{r['file']}:{r['line']}  {r['status']:<13} {r['token']}{detail}")
        counts = {}
        for r in results:
            counts[r["status"]] = counts.get(r["status"], 0) + 1
        print(f"\n{len(files)} files, {len(results)} citations: "
              + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
        print("DEAD/ARCHIVED = current section (fails); WARN-* = history section or line mentions a rename/supersession;\n"
              "PLANNED = proposed path in a design note; MOUNT = HPC mount path (mount inactive); use --all for the rest")

    if any(r["status"] in ("DEAD", "ARCHIVED") for r in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
