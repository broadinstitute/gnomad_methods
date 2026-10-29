---
name: commit
description: Use this when the user wants to commit staged changes. Drafts a commit message matching the repo's style, confirms with the user, then commits.
---

## Objective

Draft a concise commit message matching the style of recent commits in this repository, run the repo's pre-commit hooks on the staged changes, show the user what will be committed, confirm, then run the commit.

## Workflow

### 1. Check staged changes

Run `git diff --staged --stat` to see what is staged. If nothing is staged:
- Run `git status` to show what is available
- Tell the user what is unstaged and ask what to add before proceeding

Whenever you stage files, name the paths. Never use `git add .` or `git add -A`:
they also pick up untracked files that don't belong in the commit, such as
scratch files or local notes.

### 2. Run pre-commit hooks

Always run the repo's pre-commit hooks on the staged files before reviewing, so the commit matches what CI enforces. This runs even if the user has not installed the git hook.

Capture the staged file list first (only run hooks on staged files — never all files, to avoid reformatting unrelated code). Leave out files staged for deletion: the hooks have nothing to check there, and `git add` on a deleted path fails and stages nothing.
```bash
git diff --staged --name-only --diff-filter=d
```
If that list is empty (everything staged is a deletion), skip the rest of this step: there is nothing for the hooks to check, and the commands below would run on every file instead of none.

Then check whether any of those files also has unstaged changes:
```bash
git diff --name-only -- <staged files>
```
A file listed here is only partly staged. The hooks reformat the whole file on disk, and re-staging it after they run (below) would also stage the changes the user deliberately left out. Tell the user which files are partly staged and ask how to proceed before running the hooks; do not re-stage those files without their OK.

Run the hooks scoped to the staged files:
```bash
pre-commit run --files <staged files>
```
- If `pre-commit` is not on PATH, try `python -m pre_commit run --files <staged files>`.
- If pre-commit is unavailable entirely, fall back to running the tools listed in `.pre-commit-config.yaml` directly on the staged files, using the versions pinned in `requirements-dev.txt` / the config's `rev:` (a throwaway venv is fine). Ensure Python version is compatible with the pre-commit tools. Pass each tool only the staged files that match its hook's `files:` and `exclude:` patterns, with the hook's `args:`; some hooks only cover part of the repo, and running them on other files reformats code CI leaves alone.
- If neither pre-commit nor the individual tools can be run, tell the user the hooks could not run and why, and ask how to proceed — do NOT silently commit unformatted code.

Auto-fixing hooks (formatters and whitespace fixers) modify the working tree, which unstages those edits. Re-stage the same paths and re-run until the hooks pass cleanly:
```bash
git add -- <staged files>
```
Then confirm nothing is left unstaged in those files (`git diff --name-only -- <staged files>` prints nothing); if it does, the re-add failed and the index still holds unformatted content.
If a hook that can't fix files fails (a docstring or YAML check, for example), stop and surface the failure to the user rather than committing.

### 3. Review the diff

Run `git diff --staged` and read the full diff. Look for anything that shouldn't be there: debug prints, commented-out code, unintended files, or anything that looks like a secret. If something looks wrong, stop and tell the user before proceeding.

**Check for tests.** Any new or modified function needs a test (see Testing in `CLAUDE.md`). If the staged diff adds or changes code under `gnomad/` but stages nothing under `tests/`, tell the user which functions look untested and ask whether to add tests before committing. Don't block on it: the tests may be in another commit on the branch, or the change may not alter behavior (docs, formatting, a rename).

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
`claude-opus-5-5`, `claude-sonnet-5-5`. Do not copy a model name out of this
file or out of an earlier commit; the trailer records which model
assisted, so a stale name makes it useless. If the repo's existing
trailers use a different format, match the repo.

Do NOT push unless the user explicitly asks.
