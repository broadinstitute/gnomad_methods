---
name: draft-pull-request-description
description: Use this whenever the user wants to open a pull request (PR), draft a PR, or summarize their branch for merging.
---

## Objective

You are a Staff Bioinformatics Engineer preparing a comprehensive, professional, and concise Pull Request description for a gnomAD pipeline codebase. Your goal is to produce a description thorough enough that the reviewer can understand:

- The intent and motivation for the change
- The specific pipeline steps or resources affected (sample QC, variant QC, frequency annotation, VDS/HT/MT operations, Dataproc cluster config, etc.)
- The design/architecture chosen and why (especially if Hail-specific tradeoffs were involved)
- The exact changes and any new output paths or resource names
- What tests were added or changed (check `tests/` and `.github/workflows/ci.yml`)
- Anything requiring specific reviewer attention (e.g. breaking changes to resource paths, schema changes to HTs, large compute cost implications)

## Workflow

### 1. Determine the base branch

Default to `main`. Ask the user only if they indicate this is a stacked branch.

### 2. Analyze the diff

Run `git diff main...HEAD` to understand the full set of changes. Then read the changed files in full — especially any new or modified Python scripts, resource definitions, or config changes.

### 3. Identify gnomAD-specific context

Check whether the diff touches any of the following and note it in the description:

- **Hail operations**: new/modified VDS, MT, HT read/write paths, aggregations, or annotations
- **Dataproc / cluster config**: cluster size, autoscaling policy, init scripts, `hailctl dataproc submit` invocations
- **Resource paths**: new or renamed GCS paths (especially under `gs://gnomad*/`)
- **Pipeline step flags**: new CLI arguments, `--overwrite`, step-gating logic
- **Privacy or release filters**: anything touching allele frequency output, sample suppression, or release HTs
- **Test coverage**: added/changed tests under `tests/`
- **Dependencies**: if `requirements.txt`, `requirements-dev.txt`, `setup.py`, or `pyproject.toml` changed, list each added or re-pinned library and check it against the existing pins — Hail, pandas, and numpy are the usual conflicts (see Dependencies in `CLAUDE.md`). A quick check is `pip install --dry-run -r requirements.txt` in a throwaway venv. If a new library needs an existing pin loosened or bumped, say so and suggest the pin change, which belongs in the same PR.

### 4. Check for AI code reviews

Look for completion markers written by the `/review` skill, newest first:

```bash
ls .code_review/*CODE_REVIEW_COMPLETED 2>/dev/null | sort -r
```

`.code_review/` is shared by every branch in the checkout, and the review skill
records neither the branch nor the commit it reviewed, so the newest marker may
belong to another branch. Each marker's name starts with its run's timestamp
(`YYYYMMDD_HHMMSS`, local time). Use the first marker that passes both checks:

1. **It is newer than the branch's last commit**, so the review saw the current
   code. Compare its timestamp with
   `git log -1 --format=%cd --date=format-local:%Y%m%d_%H%M%S`.
2. **It reviewed this branch's files.** Compare the files its findings cite
   (`path:line` references in that run's reviewer outputs,
   `.code_review/<timestamp>_*`) with `git diff --name-only main...HEAD`. Most
   should be files this branch changes; a widely shared file such as
   `CLAUDE.md` or `README.md` appearing in both is not enough. If the run cites
   too little to tell, ask the user whether it was a review of this branch.

If no marker passes both, but one passes check 2, that run reviewed this branch
before its latest commits: don't check the box. Say when it ran and offer to
re-run `/review`.

For the marker that passes, identify which models were in the roster by listing prompt files that share its timestamp prefix:

```bash
ls .code_review/YYYYMMDD_HHMMSS_PROMPT__*.md 2>/dev/null | grep -v '_validate\.md$'  # detection
ls .code_review/YYYYMMDD_HHMMSS_PROMPT__*_validate.md 2>/dev/null                     # validation
```

The model name (`opus`, `sonnet`, `agy`) is the last `__`-separated token of a
detection file (`…__opus.md`) and the token just before `_validate.md` in a
validation file (`…_opus_validate.md`). Filter `_validate.md` out of the
detection list: when a review is split into subsets, the label itself contains
`__`, so a looser pattern also matches the validation prompts.

**Check both lists — the two rosters can differ.** A reviewer that completes
detection can still drop out of validation (an `agy` run that times out at the
validation call, a model whose CLI invocation fails partway), and the review
skill treats the two rosters as decided independently. Report the models that
actually *validated*, since that is what the confidence labels rest on. If the
two lists differ, say so rather than printing one merged roster — e.g.
"reviewed by opus, sonnet, agy; validated by opus, sonnet".

Use what you found to write the `/review` line of the AI code review section.
The template does not pre-fill it, so write exactly one of these:

- A marker passed both checks: `- [x] /review run — models: <validated models>`,
  followed by its findings and their disposition (fixed here / deferred /
  rejected as a false positive).
- Only a run from before the latest commits matched:
  `- [ ] /review last run <date>, before the latest commits`.
- Otherwise: `- [ ] /review not run`.

Then ask the user whether they ran any other AI reviews on this branch (for
example `/code-review`, or a review with another model or tool), and list each
one they confirm with what came of it. Don't list a review the user hasn't
confirmed.

### 5. Draft the PR description

Use the following template — fill every section; do not leave placeholders:

```markdown
## Summary

<!-- 1-3 bullet points on what this PR does and why -->

## Changes

<!-- Bullet list of specific changes: functions added/modified, new CLI flags, resource paths changed, etc. -->

## Hail / Compute notes

<!-- Any cluster config, VDS/HT schema changes, estimated cost or runtime if significant. Write "None" if not applicable. -->

## Testing

<!-- How was this tested? Dataproc job IDs / log snippets if relevant, or unit test coverage. -->

## AI code review

<!-- Written in step 4: the /review line (checked only for a review of this
     branch's latest commit) with its findings and their disposition, then any
     other AI reviews the user confirmed running. -->

## Reviewer notes

<!-- Anything requiring specific attention: breaking changes, schema migrations, large GCS writes, dependency changes (and whether they fit the existing pins), or design tradeoffs the reviewer should weigh in on. Write "None" if not applicable. -->
```

### 6. Merge strategy

Specify whether this PR should be **squash-merged** (preferred for single-topic changes) or **rebased** (preferred for a clean commit stack). Default to squash unless the commit history is intentionally structured.

### 7. Output

Print the complete filled-out markdown. Do not add UI screenshot instructions — this codebase has no frontend.
