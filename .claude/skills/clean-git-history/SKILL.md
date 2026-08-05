---
name: clean-git-history
description: Use this when the user wants to clean up a messy branch, squash commits, or tidy history before opening a PR. Rewrites local history into clean atomic commits without losing any changes.
allowed-tools: bash, read_file
---

## Objective

Rewrite a noisy branch into a clean, atomic series of commits. Every change must be preserved exactly — verified by an empty diff against a backup branch.

Never use interactive rebase. Use `reset --soft` + deliberate restaging instead.

## Workflow

### 1. Identify the base branch

**Do not run any diff or log commands spanning the branch until the base is confirmed.**

1. Run `git log --oneline --graph --decorate -n 20` to see where the branch diverges.
2. Check for stacking: find the oldest commit unique to this branch and run `git branch -a --contains <hash>`. If another non-main branch contains it, this branch is stacked on it.
3. **Stop and ask the user:** "Is this branch stacked on another branch? I think the base is `[branch]` — confirm or correct me."
4. After the user confirms, verify: `git merge-base <confirmed-base> HEAD` should equal the tip of the base.

### 2. Understand the existing history (read-only, mutate nothing)

1. `git status` — must be clean. If not, stop and ask the user to commit or stash first.
2. `git log -p --reverse <base>..HEAD` — read every commit's diff and message.
3. `gh pr view` if a PR exists.
4. Write a **commit plan**: an ordered list of target commits, each with:
   - Title (matching this repo's commit style — check `git log --oneline -10`)
   - One-line rationale
   - Which original commits / files it absorbs

   Map every original commit into exactly one target commit.

5. Present the plan to the user and ask for explicit approval before touching anything.

### 3. Restructure

1. Create a backup: `git branch backup-$(git branch --show-current)-$(date +%Y-%m-%d--%H-%M)`. Tell the user the backup branch name.
2. `git reset --soft <base>` then `git reset` to unstage everything.
3. For each planned commit:
   - `git add <specific files>` — never `git add .` or `git add -A`
   - `git diff --staged` — read the full diff and confirm it matches the plan before committing
   - Commit with the planned title and an `Assisted-by: ClaudeCode:<model-id>`
     trailer, where `<model-id>` is the model **you are actually running as**
     (e.g. `claude-opus-5`). Do not copy the placeholder or a model name from
     an earlier commit — a stale name defeats the point of the trailer. Match
     the repo's existing trailer format if it differs.
     ```bash
     git commit -m "$(cat <<'EOF'
     <title>

     Assisted-by: ClaudeCode:<model-id>
     EOF
     )"
     ```

### 4. Verify

1. `git diff <backup-branch> HEAD` — this diff **must be empty**. If it is not, immediately run `git reset --hard <backup-branch>` to restore the original, tell the user what happened, and re-plan.
2. `git status` — must be clean.
3. Show the user the new `git log --oneline <base>..HEAD`.

### 5. Handoff

Do NOT push. Tell the user to run `git push --force-with-lease` themselves once they are satisfied. Offer to draft a PR description using the draft-pull-request-description skill.
