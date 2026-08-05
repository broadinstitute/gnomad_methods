---
name: draft-pull-request-description
description: Use this whenever the user wants to open a pull request (PR), draft a PR, or summarize their branch for merging.
allowed-tools: bash, read_file
---

## Objective

You are a Staff Bioinformatics Engineer preparing a comprehensive, professional, and concise Pull Request description for a gnomAD pipeline codebase. Your goal is to produce a description thorough enough that the reviewer can understand:

- The intent and motivation for the change
- The specific pipeline steps or resources affected (sample QC, variant QC, frequency annotation, VDS/HT/MT operations, Dataproc cluster config, etc.)
- The design/architecture chosen and why (especially if Hail-specific tradeoffs were involved)
- The exact changes and any new output paths or resource names
- What tests were added or changed (check `gnomad_qc/tests/` and `.github/workflows/ci.yml`)
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
- **Test coverage**: added/changed tests under `gnomad_qc/tests/`

### 4. Check for a completed code review

Look for a completion marker written by the `/review` skill:

```bash
ls .code_review/*CODE_REVIEW_COMPLETED 2>/dev/null | sort | tail -1
```

If a marker exists, take its timestamp prefix (format `YYYYMMDD_HHMMSS`) and identify which models were in the roster by listing prompt files that share that prefix:

```bash
ls .code_review/YYYYMMDD_HHMMSS_PROMPT__*__*.md 2>/dev/null        # detection
ls .code_review/YYYYMMDD_HHMMSS_PROMPT__*_validate.md 2>/dev/null  # validation
```

The suffix before `.md` is the model name (`opus`, `sonnet`, `agy`).

**Check both lists — the two rosters can differ.** A reviewer that completes
detection can still drop out of validation (an `agy` run that times out at the
validation call, a model whose CLI invocation fails partway), and the review
skill treats the two rosters as decided independently. Report the models that
actually *validated*, since that is what the confidence labels rest on. If the
two lists differ, say so rather than printing one merged roster — e.g.
"reviewed by opus, sonnet, agy; validated by opus, sonnet".

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

## Code review

<!-- Emit exactly ONE of the two lines below, whichever is true, and delete the other.
     Never emit both: an unchecked "not run" sitting above a checked "run" reads as
     if no review happened, and a skimmer will believe the first line.
       run:     - [x] `/review` run — models: Opus, Sonnet
       not run: - [ ] `/review` not run
     If a review was run, follow the line with its findings and their disposition
     (fixed here / deferred / rejected as a false positive). -->
- [x] `/review` run — models: <!-- e.g. Opus, Sonnet -->

## Reviewer notes

<!-- Anything requiring specific attention: breaking changes, schema migrations, large GCS writes, or design tradeoffs the reviewer should weigh in on. Write "None" if not applicable. -->
```

### 6. Merge strategy

Specify whether this PR should be **squash-merged** (preferred for single-topic changes) or **rebased** (preferred for a clean commit stack). Default to squash unless the commit history is intentionally structured.

### 7. Output

Print the complete filled-out markdown. Do not add UI screenshot instructions — this codebase has no frontend.
