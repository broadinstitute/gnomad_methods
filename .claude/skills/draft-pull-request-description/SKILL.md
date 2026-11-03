---
name: draft-pull-request-description
description: Use this whenever the user wants to open a pull request (PR), draft a PR, or summarize their branch for merging.
---

## Objective

You are a Staff Bioinformatics Engineer preparing a concise, professional Pull Request description for gnomad_methods, the shared Hail utility library that gnomAD pipelines and other repos import. Your goal is a description that lets the reviewer understand:

- The intent and motivation for the change
- The specific modules or resources affected (sample QC, variant QC, frequency annotation, VDS/HT/MT operations, Dataproc cluster config, etc.)
- The design/architecture chosen and why (especially if Hail-specific tradeoffs were involved)
- The exact changes and any new output paths or resource names
- What tests were added or changed (check `tests/` and `.github/workflows/ci.yml`)
- Anything requiring specific reviewer attention (e.g. breaking changes to resource paths, schema changes to HTs, large compute cost implications)

## Workflow

### 1. Determine the base branch

Find the base with the `find-base-branch` skill (`.claude/skills/find-base-branch/SKILL.md`), and use the ref it reports as `<base>` everywhere below.

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
- **Dependencies**: if `requirements.txt`, `requirements-dev.txt`, `setup.py`, or `pyproject.toml` changed, follow Dependencies in `CLAUDE.md`, and list each added or re-pinned library and any pin change it needs.
- **Downstream callers**: if a public function's signature, return schema, field names, or defaults changed, do the `gnomad_qc` check under Pull Requests in `CLAUDE.md` (against `gnomad_qc`'s `main`, not only a local checkout), and list the call sites that break with the pipeline version each belongs to.
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
2. **Find what changed since the run.** The branch's reflog records which
   commit it pointed to when the run started:
   `git rev-parse <branch>@{YYYY-MM-DD HH:MM:SS}`, using the run's start time.
   `git diff --stat <that commit> HEAD` then shows what changed since. Compare
   content rather than commit dates: a rebase or `/clean-git-history` gives
   every commit a new date even when the content is unchanged. If the branch
   was rebased onto a newer base since, limit the diff to this branch's files
   (`-- $(git diff --name-only <base>...HEAD)`). If git warns that the log
   only goes back to a later date, the reflog can't tell you (the branch is
   newer than the run, or this is a different clone); list the commits dated
   after the run instead
   (`git log --format='%h %cd %s' --date=format-local:%Y%m%d_%H%M%S <base>..HEAD`)
   and say the dates may reflect a rewrite.
3. **Ask one question.** Show the run's start time, models, and findings, plus
   what changed since it started. Ask whether that run reviewed this branch,
   and, if anything changed, whether those changes only fix its findings.
4. **Fill in the disposition** of each finding (fixed here, deferred, or
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

**Keep it short.** Use at most two sentences per bullet, and leave out anything
the diff already makes obvious. Reviewers read the whole description, and a long
one hides the parts that need their attention.

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

### 6. Title and Changelog label

`.github/pull_request_template.md` asks for both, because the release notes are
generated from merged PRs and sorted by label (`.github/release.yml`):

- **Title:** a short, descriptive phrase that reads well as one line of the
  release notes. If the PR already exists and its title no longer covers the
  branch, suggest a new one.
- **Label:** one of the `Changelog:` labels defined in `.github/release.yml`,
  with a one-line reason. Changes that don't touch the published package, such
  as CI or contributor tooling, take `Changelog: ignore`.

### 7. Merge strategy

Specify whether this PR should be **squash-merged** (preferred for single-topic changes) or **rebased** (preferred for a clean commit stack). Default to squash unless the commit history is intentionally structured.

### 8. Output

Print the suggested title and label, then the complete filled-out markdown. Set
the title or label on the PR (`gh pr edit --title` / `--add-label`) only if the
user asks. Do not add UI screenshot instructions — this codebase has no frontend.
