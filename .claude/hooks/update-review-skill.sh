#!/bin/sh
# SessionStart hook: keep a personal clone of Ben Weisburd's /review skill
# current, so contributors get upstream updates without doing anything.
#
# A personal skill in ~/.claude/skills/ takes precedence over the copy of the
# same skill in this repo, so the clone is what /review actually runs. Fails
# quietly when offline; prints only when there is something to tell the user.

dir="$HOME/.claude/skills/review"
repo="https://github.com/bw2/claude-code-review-skill.git"

if [ ! -e "$dir" ]; then
    if git clone -q --depth 1 "$repo" "$dir" >/dev/null 2>&1; then
        echo "Installed the /review skill from $repo into $dir."
    fi
elif git -C "$dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git -C "$dir" pull -q --ff-only >/dev/null 2>&1
else
    echo "$dir is not a git clone of $repo, so the /review skill was not updated."
fi

exit 0
