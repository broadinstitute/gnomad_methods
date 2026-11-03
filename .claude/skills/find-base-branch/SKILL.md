---
name: find-base-branch
description: Use this to work out which branch the current branch is based on, for example before drafting a PR description or rewriting history. The draft-pull-request-description and clean-git-history skills use it.
---

## Objective

Find the ref this branch should be compared against: the branch it was cut from, at the version the PR is actually compared with. Getting this wrong puts other people's commits into a diff or a history rewrite.

## Workflow

1. **Fetch.** Run `git fetch` so the remote-tracking refs are current.
2. **Check for stacking.** Find the oldest commit not on `origin/main` (`git rev-list --reverse origin/main..HEAD | head -1`) and run `git branch -a --contains <hash> --no-contains HEAD`, which leaves out branches stacked on top of this one. Also ignore this branch's own remote-tracking ref (`origin/<this branch>`) and any `backup-*` branches left by `clean-git-history`. If another branch still contains the commit, this branch is stacked on it; otherwise the base is `main`.
3. **Compare local and remote.** Compare the local base with its remote-tracking ref (`origin/<base>`). If they point to the same commit, or only one of them exists, use that one. If they differ, show the user how (`git rev-list --left-right --count <base>...origin/<base>` counts the commits only on each side) and ask which to use: a stale local ref makes upstream commits look like part of this branch, and a local base with unpushed commits makes `origin/<base>` too old.
4. **Confirm.** Tell the user which base you found. Ask them to confirm it if it isn't `main`, or if the skill that called this one says to always confirm. Use the chosen ref (for example `origin/main`) as `<base>`.
