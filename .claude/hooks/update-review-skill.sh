#!/bin/sh
# SessionStart hook: keep a personal clone of Ben Weisburd's /review skill
# current, so contributors get upstream updates without doing anything.
#
# The clone is installed as a personal skill (~/.claude/skills/review), which
# makes /review available in every project; this repo ships no copy of it. The
# hook checks upstream at most once a day, is silent when offline or already
# current, and prints only when there is something to tell the user.

dir="$HOME/.claude/skills/review"
repo="https://github.com/bw2/claude-code-review-skill.git"

if [ ! -e "$dir" ]; then
    if git clone -q --depth 1 "$repo" "$dir" >/dev/null 2>&1; then
        echo "Installed the /review skill from $repo into $dir as a personal" \
            "skill, so it is available in every project. Restart Claude Code" \
            "(or start a new session) to use it."
    fi
    exit 0
fi

# Only update a directory that is itself a clone of the upstream repo, not a
# copied skill or a directory nested inside some other repository.
real_dir=$(cd "$dir" 2>/dev/null && pwd -P)
top=$(git -C "$dir" rev-parse --show-toplevel 2>/dev/null)
origin_url=$(git -C "$dir" remote get-url origin 2>/dev/null)
case "$origin_url" in
    *bw2/claude-code-review-skill*) ;;
    *) origin_url="" ;;
esac
if [ -z "$top" ] || [ "$top" != "$real_dir" ] || [ -z "$origin_url" ]; then
    echo "$dir is not a clone of $repo, so the /review skill was not updated."
    exit 0
fi

# Rate-limit the network round-trip to once a day.
if [ -n "$(find "$dir/.git/FETCH_HEAD" -mmin -1440 2>/dev/null)" ]; then
    exit 0
fi

# Offline: say nothing.
git -C "$dir" fetch -q origin >/dev/null 2>&1 || exit 0

if [ -n "$(git -C "$dir" status --porcelain --untracked-files=no)" ]; then
    echo "The /review skill clone at $dir has local changes, so it was not" \
        "updated to the latest upstream version."
    exit 0
fi
if ! git -C "$dir" merge -q --ff-only '@{u}' >/dev/null 2>&1; then
    echo "The /review skill clone at $dir could not be fast-forwarded to" \
        "upstream (local commits, or upstream rewrote history), so it was not" \
        "updated."
fi

exit 0
