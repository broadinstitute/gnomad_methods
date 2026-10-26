#!/bin/sh
# SessionStart hook: keep a personal clone of Ben Weisburd's /review skill
# current, so contributors get upstream updates without doing anything.
#
# The clone is installed as a personal skill (~/.claude/skills/review), which
# makes /review available in every project; this repo ships no copy of it. The
# hook checks upstream at most once a day, is silent when offline or already
# current, and prints only when there is something to tell the user. Set
# GNOMAD_SKIP_REVIEW_HOOK=1 to turn it off.

[ -n "$GNOMAD_SKIP_REVIEW_HOOK" ] && exit 0
[ -n "$HOME" ] || exit 0

dir="$HOME/.claude/skills/review"
repo="https://github.com/bw2/claude-code-review-skill.git"
stamp="$HOME/.claude/.review-skill-last-check"
partial="$HOME/.claude/.review-skill-clone"
reload_skills='"hookSpecificOutput": {"hookEventName": "SessionStart", "reloadSkills": true}'

# Plain stdout from a SessionStart hook only reaches Claude, so messages for
# the user go in the JSON systemMessage field.
json_escape() {
    printf '%s' "$*" | sed 's/\\/\\\\/g; s/"/\\"/g'
}
notify() {
    printf '{"systemMessage": "%s"}\n' "$(json_escape "$@")"
}

# Clone into a temporary directory and move it into place, so a clone cut off
# by the hook's timeout can't leave a partial directory behind.
if [ ! -e "$dir" ]; then
    rm -rf "$partial"
    if mkdir -p "$(dirname "$dir")" \
        && git clone -q --depth 1 "$repo" "$partial" >/dev/null 2>&1 \
        && mv "$partial" "$dir"; then
        printf '{"systemMessage": "%s", %s}\n' "$(json_escape \
            "Installed the /review skill from $repo into $dir as a personal" \
            "skill, so it is available in every project.")" "$reload_skills"
    else
        rm -rf "$partial"
    fi
    exit 0
fi

# Check upstream at most once a day.
today=$(date +%Y%m%d)
[ "$(cat "$stamp" 2>/dev/null)" = "$today" ] && exit 0

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
    echo "$today" > "$stamp"
    notify "$dir is not a clone of $repo, so this repo's hook does not update" \
        "the /review skill there. To get automatic updates, move that" \
        "directory aside; the next session clones a fresh copy."
    exit 0
fi

# Offline: say nothing, and try again next session.
git -C "$dir" fetch -q origin >/dev/null 2>&1 || exit 0
echo "$today" > "$stamp"

if [ -n "$(git -C "$dir" status --porcelain --untracked-files=no 2>/dev/null)" ]; then
    notify "The /review skill at $dir has local changes, so it was not" \
        "updated. Discard them, or delete the directory and the next session" \
        "clones a fresh copy."
    exit 0
fi
before=$(git -C "$dir" rev-parse HEAD 2>/dev/null)
if ! git -C "$dir" merge -q --ff-only '@{u}' >/dev/null 2>&1; then
    notify "The /review skill at $dir could not be fast-forwarded to upstream" \
        "(it has local commits, is on a detached HEAD, or upstream rewrote" \
        "history), so it was not updated. Delete the directory and the next" \
        "session clones a fresh copy."
elif [ "$(git -C "$dir" rev-parse HEAD 2>/dev/null)" != "$before" ]; then
    printf '{%s}\n' "$reload_skills"
fi

exit 0
