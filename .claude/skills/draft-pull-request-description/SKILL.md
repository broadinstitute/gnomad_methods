---
name: draft-pull-request-description
description: Use this whenever the user wants to open a pull request (PR), draft a PR, or summarize their branch for merging.
---

## Objective

You are a Staff Bioinformatics Engineer preparing a comprehensive, professional, and concise Pull Request description for gnomad_methods, the shared Hail utility library that gnomAD pipelines and other repos import. Your goal is to produce a description thorough enough that the reviewer can understand:

- The intent and motivation for the change
- The specific modules or resources affected (sample QC, variant QC, frequency annotation, VDS/HT/MT operations, Dataproc cluster config, etc.)
- The design/architecture chosen and why (especially if Hail-specific tradeoffs were involved)
- The exact changes and any new output paths or resource names
- What tests were added or changed (check `tests/` and `.github/workflows/ci.yml`)
- Anything requiring specific reviewer attention (e.g. breaking changes to resource paths, schema changes to HTs, large compute cost implications)

## Workflow

### 1. Determine the base branch

Default to `main`. Ask the user only if they indicate this is a stacked branch. Use the chosen branch as `<base>` everywhere below.

### 2. Analyze the diff

Run `git diff <base>...HEAD` to understand the full set of changes. Then read the changed files in full — especially any new or modified Python scripts, resource definitions, or config changes.

### 3. Identify gnomAD-specific context

Check whether the diff touches any of the following and note it in the description:

- **Hail operations**: new/modified VDS, MT, HT read/write paths, aggregations, or annotations
- **Dataproc / cluster config**: cluster size, autoscaling policy, init scripts, `hailctl dataproc submit` invocations
- **Resource paths**: new or renamed GCS paths (especially under `gs://gnomad*/`)
- **Public API**: new, renamed, or removed public functions, and new parameters with their defaults. Say whether each change is additive or breaking for code that imports the library (see Project Overview in `CLAUDE.md`).
- **Privacy or release filters**: anything touching allele frequency output, sample suppression, or release HTs
- **Test coverage**: added/changed tests under `tests/`
- **Dependencies**: if `requirements.txt`, `requirements-dev.txt`, `setup.py`, or `pyproject.toml` changed, list each added or re-pinned library and check it against the existing pins — Hail, pandas, and numpy are the usual conflicts (see Dependencies in `CLAUDE.md`). A quick check is `pip install --dry-run -r requirements.txt` in a throwaway venv. If a new library needs an existing pin loosened or bumped, say so and suggest the pin change, which belongs in the same PR.
- **Downstream callers**: if the signature, return schema, field names, or defaults of a public function changed, grep `gnomad_qc` for callers, on its `main` branch and not only a local checkout. List the call sites that break and which pipeline version each belongs to (see Pull Requests in `CLAUDE.md`). Flag them in the description rather than fixing them here; the fix belongs in a coordinated `gnomad_qc` PR.
- **Repo layout**: if a subpackage was added, removed, renamed, or moved, check that this PR also updates the repo layout table in the README.

### 4. Check for AI code reviews

`CLAUDE.md` asks for `/review` on every PR, and for the description to say what
came of it. Gather what you can, then ask the user to confirm it:

1. **Find the latest `/review` run.** If `/review` ran earlier in this
   conversation, use its results. Otherwise look in `.code_review/`, where the
   skill writes its prompts: the newest `*CODE_REVIEW_COMPLETED` marker's name
   starts with the run's start time (`YYYYMMDD_HHMMSS`), the prompt files with
   the same prefix name its models (`…__opus.md`), and its `…_validate.md`
   prompts list the findings. The directory is shared by every branch in the
   checkout and doesn't record which branch was reviewed.
2. **Ask one question.** Show the run's start time, models, and findings, plus
   any commits made after it started
   (`git log --format='%h %cd %s' --date=format-local:%Y%m%d_%H%M%S <base>..HEAD`).
   Ask whether that run reviewed this branch, and whether those later commits
   only fix its findings.
3. **Fill in the disposition** of each finding (fixed here, deferred, or
   rejected as a false positive) from the branch's commits, and ask only about
   the findings they don't settle.

If you find no run, or the files aren't laid out as described, ask the user
which reviews they ran, with which models, and what came of the findings.

Then write the `/review` line of the AI code review section as exactly one of:

- `- [x] /review run — models: <models>`, followed by its findings and their
  disposition, if the user confirms it reviewed this branch with no changes
  after it except fixes for its findings.
- `- [ ] /review last run <date>, before the latest commits`, if it ran before
  later changes other than fixes for its findings. Offer to re-run it.
- `- [ ] /review not run`.

Then ask whether they ran any other AI reviews on this branch (for example
`/code-review`, or a review with another model or tool), and list each one they
confirm with what came of it. Don't list a review the user hasn't confirmed.

### 5. Draft the PR description

**Ask the user for the summary first.** Before drafting, ask them to describe in
a sentence or two what changed, why, and which tests cover it. The person opening
the PR should be able to explain it, and reviewers rely on that. Build the
Summary section from their answer: tidy the wording, but don't add claims they
didn't make. If the diff shows a significant change their summary leaves out,
point it out and ask rather than adding it yourself.

Use the following template — fill every section; do not leave placeholders:

```markdown
## Summary

<!-- 1-3 bullet points built from the user's own summary of what this PR does and why -->

## Changes

<!-- Bullet list of specific changes: functions added/modified, new CLI flags, resource paths changed, etc. -->

## Hail / Compute notes

<!-- Any cluster config, VDS/HT schema changes, estimated cost or runtime if significant. Write "None" if not applicable. -->

## Testing

<!-- How was this tested? Dataproc job IDs / log snippets if relevant, or unit test coverage. -->

## AI code review

<!-- Written in step 4: the /review line (checked only when the user confirms
     the run reviewed this branch and later commits only fix its findings) with
     its findings and their disposition, then any other AI reviews the user
     confirmed running. -->

## Reviewer notes

<!-- Anything requiring specific attention: breaking changes, schema migrations, large GCS writes, dependency changes (and whether they fit the existing pins), or design tradeoffs the reviewer should weigh in on. Write "None" if not applicable. -->
```

### 6. Merge strategy

Specify whether this PR should be **squash-merged** (preferred for single-topic changes) or **rebased** (preferred for a clean commit stack). Default to squash unless the commit history is intentionally structured.

### 7. Output

Print the complete filled-out markdown. Do not add UI screenshot instructions — this codebase has no frontend.
