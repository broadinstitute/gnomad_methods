---
name: clean-git-history
description: Use this when the user wants to clean up a messy branch, squash commits, or tidy history before opening a PR. Rewrites local history into clean atomic commits without losing any changes.
---

## Objective

Rewrite a noisy branch into a clean, atomic series of commits. Every change must be preserved exactly — verified by an empty diff against a backup branch.

Never use interactive rebase. Use `reset --soft` + deliberate restaging instead.

## Workflow

### 1. Identify the base branch

**Do not run any diff or log commands spanning the branch until the base is confirmed.**

1. Run `git log --oneline --graph --decorate -n 20` to see where the branch diverges.
2. Find the base with the `find-base-branch` skill (`.claude/skills/find-base-branch/SKILL.md`), and have the user confirm it even if it is `main`: a rewrite on the wrong base is hard to spot afterwards.
3. Record the fork point: `git merge-base <base> HEAD`. If it is not the tip of the base, the base has moved on since the branch was cut. Tell the user, and ask whether to rebase onto the base first or to clean up the history in place on top of the fork point. If they want to rebase, stop here and let them do it (it may need conflicts resolved), then start this skill again.

### 2. Understand the existing history (read-only, mutate nothing)

1. `git status` — must be clean. If not, stop and ask the user to commit or stash first.
2. `git log -p --reverse <base>..HEAD` — read every commit's diff and message.
3. `gh pr view` if a PR exists.
4. Write a **commit plan**: an ordered list of target commits, each with:
   - Title (matching this repo's commit style — check `git log --oneline -10`)
   - One-line rationale
   - Which original commits / files it absorbs
   - Authorship (see Attribution in `CLAUDE.md`), from
     `git log --format='%h %an <%ae>%n%(trailers:key=Co-authored-by,valueonly)' <fork-point>..HEAD`:
     - **Author:** the dominant original author of the commits it absorbs. A
       rewrite otherwise makes whoever runs it the author of every commit.
     - **Co-authors:** the other original authors, plus the people in the
       absorbed commits' `Co-authored-by:` trailers. GitHub adds those when a
       reviewer's suggestion is applied, and a rewrite would otherwise drop
       them. Keep only people: leave out AI co-author lines (e.g.
       `Co-Authored-By: Claude ...`), since the rebuilt commit's `Assisted-by`
       trailer covers those.

   Map every original commit into exactly one target commit.

   Target commits must own disjoint sets of files: step 3 stages whole files, so one file's changes can't be split across commits that way. If two original commits edit the same file and belong in different target commits, either merge those target commits or plan to commit an intermediate version of the file in the earlier one. When the earlier target commit absorbs the file's earlier original commits, that version is exactly `git show <last-absorbed-commit>:<file>`; otherwise it has to be edited by hand. Say which in the plan, so the user approves it.

5. Present the plan to the user and ask for explicit approval before touching anything.

### 3. Restructure

1. Create a backup: `git branch backup-$(git branch --show-current)-$(date +%Y-%m-%d--%H-%M)`. Tell the user the backup branch name.
2. `git reset --soft <fork-point>` (the merge-base from step 1.3) then `git reset` to unstage everything. Never reset to the tip of the base: if the base has moved, the new commits would also undo every commit added to the base since the fork, and the check in step 4 would not catch it.
3. For each planned commit, stage and commit the way the `commit` skill
   (`.claude/skills/commit/SKILL.md`) says to: name the paths when staging
   (its "Check staged changes" section) and end the message with its
   `Assisted-by` trailer (its "Commit" section). Skip its
   other steps, such as running hooks and confirming each commit: the plan is
   already approved, and a hook that reformats a file would make the check in
   step 4 fail.
   - `git add <specific files>`
   - For a file the plan splits across commits: copy its final version to a
     temporary path outside the repo, write the intermediate version into the
     file (`git show <commit>:<file> > <file>`, or the hand-edited version),
     `git add` it, and copy the final version back after committing. This works
     for files the branch added, which are untracked after the reset.
   - `git diff --staged` — read the full diff and confirm it matches the plan before committing
   - `git commit --no-verify --author="<Name> <email>"` with the planned
     author and title, a `Co-authored-by:` line for each co-author the plan
     lists for this commit, and the trailer. The
     `--no-verify` matters: if the user installed the pre-commit git hook, it
     would reformat files during the commit and the result would no longer match
     the backup branch in step 4. This skill reorganizes the branch's content
     without changing it; formatting fixes belong in their own commit afterwards,
     via the `commit` skill.

### 4. Verify

1. `git diff <backup-branch> HEAD` — this diff **must be empty**. If it is not, immediately run `git reset --hard <backup-branch>` to restore the original, tell the user what happened, and re-plan.
2. `git status` — must be clean.
3. Show the user the new `git log --oneline <base>..HEAD`.

### 5. Handoff

Do NOT push. Tell the user to run `git push --force-with-lease` themselves once they are satisfied. Offer to draft a PR description using the draft-pull-request-description skill.
