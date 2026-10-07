---
name: commit
description: Use this when the user wants to commit staged changes. Drafts a commit message matching the repo's style, confirms with the user, then commits.
allowed-tools: bash, read_file
---

## Objective

Draft a concise commit message matching the style of recent commits in this repository, run the repo's pre-commit hooks on the staged changes, show the user what will be committed, confirm, then run the commit.

## Workflow

### 1. Check staged changes

Run `git diff --staged --stat` to see what is staged. If nothing is staged:
- Run `git status` to show what is available
- Tell the user what is unstaged and ask what to add before proceeding

### 2. Run pre-commit hooks

Always run the repo's pre-commit hooks on the staged files before reviewing, so the commit matches what CI enforces (black, isort, autopep8, pydocstyle, whitespace fixers, etc.). This runs even if the user has not installed the git hook.

Skip this step only if the repo has no `.pre-commit-config.yaml`.

Capture the staged file list first (only run hooks on staged files — never all files, to avoid reformatting unrelated code):
```bash
git diff --staged --name-only
```

Run the hooks scoped to those files:
```bash
pre-commit run --files <staged files>
```
- If `pre-commit` is not on PATH, try `python -m pre_commit run --files <staged files>`.
- If pre-commit is unavailable entirely, fall back to running the tools listed in `.pre-commit-config.yaml` directly on the staged files, using the versions pinned in `requirements-dev.txt` / the config's `rev:` (a throwaway venv is fine). For this repo that means: `black`, `isort --profile black --filter-files`, `autopep8 --exit-code --in-place`, `pydocstyle`. autopep8 needs `lib2to3`, so run it under Python ≤ 3.12, not 3.13.
- If neither pre-commit nor the individual tools can be run, tell the user the hooks could not run and why, and ask how to proceed — do NOT silently commit unformatted code.

Auto-fixing hooks (black, isort, autopep8, end-of-file-fixer, trailing-whitespace) modify the working tree, which unstages those edits. Re-stage the same paths and re-run until the hooks pass cleanly:
```bash
git add -- <staged files>
```
If a non-autofixable hook fails (e.g. pydocstyle, check-yaml), stop and surface the failure to the user rather than committing.

### 3. Review the diff

Run `git diff --staged` and read the full diff. Look for anything that shouldn't be there: debug prints, commented-out code, unintended files, or anything that looks like a secret. If something looks wrong, stop and tell the user before proceeding.

### 4. Read recent commit style

Run `git log --oneline -10` to understand the repo's commit message convention (capitalization, prefix style, length, etc.).

### 5. Draft a commit message

Match the observed style exactly. General rules unless the repo overrides them:
- Short imperative phrase, ≤ 72 characters
- No conventional commit prefixes (`feat:`, `fix:`, etc.) unless the repo uses them
- No trailing period
- Add a blank line and a body paragraph only when the change genuinely needs more context than the subject line can carry

### 6. Confirm with the user

Show:
- The list of staged files (`git diff --staged --name-only`)
- The draft commit message

Ask the user to confirm, request edits, or abort. Do not commit until they confirm.

### 7. Commit

Run:
```bash
git commit -m "$(cat <<'EOF'
<message here>

Assisted-by: ClaudeCode:<model-id>
EOF
)"
```

`<model-id>` is the model **you are actually running as** — e.g.
`claude-opus-5`, `claude-sonnet-5`. Do not copy a model name out of this
file or out of an earlier commit; the trailer records which model
assisted, so a stale name makes it useless. If the repo's existing
trailers use a different format, match the repo.

Do NOT push unless the user explicitly asks.
