"""
Title: Task-file layout checker
Created: 2026-11-25
Last modified: 2026-11-25
Summary: Check task handoffs in context/tasks/ against the "Task files" rule in
AGENTS.md: required sections in the agreed order, one Next Step, no
"Session Update" sections, and the line limit. Files without both a
"## Current State" and a "## History" heading have not been converted to the
layout yet; they are
reported as NOT CONVERTED and only fail with --strict. Files under an
`archive/` folder are skipped.
Inputs: Markdown files (default: context/tasks/**/*.md, excluding archive/).
Outputs: Report on stdout. Exit status 1 if any ERROR is found. Writes nothing.
Usage: python tools/check_task_file.py [files ...] [--strict] [--quiet]

Keep SECTIONS and MAX_LINES in sync with the "Task files" rule in AGENTS.md.
"""

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TASKS_DIR = REPO_ROOT / "context" / "tasks"

# Canonical section order; True = required.
SECTIONS = [
    ("Goal", True),
    ("Current State", True),
    ("Decisions", True),
    ("Next Step", True),
    ("Later", False),
    ("Open Questions", False),
    ("Pitfalls", False),
    ("History", True),
]
MAX_LINES = 250
WARN_LINES = 225          # warn when close to the limit, so archiving happens first
NEXT_STEP_WARN_LINES = 25  # a long Next Step usually holds detail or several steps

H2_RE = re.compile(r"^##\s+(.+?)\s*$")
H3_RE = re.compile(r"^###\s+")
SESSION_RE = re.compile(r"session\s+(update|handoff|history|log)", re.I)
# Markers that a Later item is already finished.
DONE_RE = re.compile(r"~~|✅|\[x\]|\b(done|completed|finished)\b", re.I)


def canonical(heading):
    """Map a heading like 'Current State (updated 2026-11-24)' to its section name."""
    base = re.sub(r"\s*\(.*\)\s*$", "", heading).strip().rstrip(":")
    for name, _ in SECTIONS:
        if base.lower() == name.lower():
            return name
    return None


def check_file(path, strict):
    """Return (status, messages) where messages are (level, line, text)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    msgs = []

    in_fence = False
    headings = []  # (line_no, raw, canonical)
    for i, line in enumerate(lines, 1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = H2_RE.match(line)
        if m:
            headings.append((i, m.group(1), canonical(m.group(1))))

    # Converted files have both of these; some older files have "Current State" alone.
    found = {c for _, _, c in headings}
    if not {"Current State", "History"} <= found:
        level = "ERROR" if strict else "INFO"
        msgs.append((level, 0, "NOT CONVERTED: lacks '## Current State' and/or '## History'; convert when next worked on"))
        return msgs

    # Session logs.
    for ln, raw, _ in headings:
        if SESSION_RE.search(raw):
            msgs.append(("ERROR", ln, f"session-log heading '## {raw}': fold it into Current State/History"))

    # Duplicates, unknown sections, order.
    seen = {}
    order = [name for name, _ in SECTIONS]
    last_idx = -1
    for ln, raw, c in headings:
        if c is None:
            if not SESSION_RE.search(raw):
                msgs.append(("WARN", ln, f"extra section '## {raw}': not in the layout; consider moving it into Current State or a linked doc"))
            continue
        if c in seen:
            msgs.append(("ERROR", ln, f"duplicate '## {c}' (first at line {seen[c]})"))
            continue
        seen[c] = ln
        idx = order.index(c)
        if idx < last_idx:
            msgs.append(("ERROR", ln, f"'## {c}' is out of order; expected order: {', '.join(order)}"))
        last_idx = max(last_idx, idx)

    for name, required in SECTIONS:
        if required and name not in seen:
            msgs.append(("ERROR", 0, f"missing required section '## {name}'"))

    # Section bodies.
    bounds = {}
    for k, (ln, _, c) in enumerate(headings):
        end = headings[k + 1][0] - 1 if k + 1 < len(headings) else len(lines)
        if c and c not in bounds:
            bounds[c] = (ln, end)

    if "Next Step" in bounds:
        start, end = bounds["Next Step"]
        body = [l for l in lines[start:end] if l.strip()]
        if not body:
            msgs.append(("ERROR", start, "Next Step is empty"))
        subs = [l for l in lines[start:end] if H3_RE.match(l)]
        if len(subs) > 1:
            msgs.append(("WARN", start, f"Next Step has {len(subs)} subheadings: keep one current step and queue the rest under '## Later'"))
        if len(body) > NEXT_STEP_WARN_LINES:
            msgs.append(("WARN", start, f"Next Step is {len(body)} non-blank lines: keep the step itself short and move detail or queued steps elsewhere"))

    if "Later" in bounds:
        start, end = bounds["Later"]
        for j in range(start, end):
            l = lines[j]
            if l.lstrip().startswith(("-", "*", "1", "2", "3", "4", "5", "6", "7", "8", "9")) and DONE_RE.search(l):
                msgs.append(("WARN", j + 1, "Later item looks finished: remove it (record outcome in Current State/History)"))

    n = len(lines)
    if n > MAX_LINES:
        msgs.append(("ERROR", 0, f"{n} lines (limit {MAX_LINES}): archive superseded detail to archive/<task>_to_<date>.md"))
    elif n > WARN_LINES:
        msgs.append(("WARN", 0, f"{n} lines, close to the {MAX_LINES}-line limit: archive before adding much more"))

    return msgs


def default_files():
    return sorted(p for p in TASKS_DIR.rglob("*.md") if "archive" not in p.relative_to(TASKS_DIR).parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Usage:")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", type=Path, help="task files (default: all in context/tasks, excluding archive/)")
    ap.add_argument("--strict", action="store_true", help="treat unconverted files as errors")
    ap.add_argument("--quiet", action="store_true", help="print only ERROR lines")
    args = ap.parse_args()

    files = [p.resolve() for p in args.files] or default_files()
    n_err = n_warn = 0
    for path in files:
        if "archive" in path.parts:
            continue
        try:
            rel = path.relative_to(REPO_ROOT)
        except ValueError:
            rel = path
        msgs = check_file(path, args.strict)
        for level, ln, text in msgs:
            n_err += level == "ERROR"
            n_warn += level == "WARN"
            if args.quiet and level != "ERROR":
                continue
            print(f"{rel}:{ln}  {level}  {text}")
        if not msgs and not args.quiet:
            print(f"{rel}  OK")

    if not args.quiet:
        print(f"\n{len(files)} file(s) checked: {n_err} error(s), {n_warn} warning(s)")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
