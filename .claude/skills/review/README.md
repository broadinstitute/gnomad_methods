# claude-code-review-skill

A [Claude Code](https://claude.com/claude-code) skill that runs code review using Claude under two different models — Opus and Sonnet — plus an optional third reviewer, [agy](https://antigravity.google/) (Google Antigravity). The reviewers in the roster independently review your code in parallel, merge the results into a deduplicated list, and then do a second validation pass to see if they agree with each other's findings. At the end, you are presented with a list of all findings along with the level of agreement between reviewers (`Agreed`, `Disputed`, `Unanimous`, etc — the exact tiers depend on whether agy is in the roster), and the proposed fixes.

agy is optional: the skill drops it automatically if the `agy` CLI isn't installed, and you can also ask to skip it explicitly (e.g. `/review --no-agy`) even when it is installed.

## Requirements

- [Claude Code](https://claude.com/claude-code) (`claude` CLI) — invoked twice per step, once pinned to Opus and once to Sonnet.
- [`agy`](https://antigravity.google/) (Google Antigravity CLI), installed and authenticated — **optional**.

## Install

Clone this repo directly into your personal skills directory:

```bash
git clone https://github.com/bw2/claude-code-review-skill.git ~/.claude/skills/review
```

Restart Claude Code (or start a new session) and the `/review` command is available in any project.

## Usage

```
/review
```

With no argument, reviews `git diff main...HEAD` if there is one, otherwise the whole current directory. You can also target something specific:

```
/review path/to/file.py path/to/other.py
/review --diff main...HEAD
/review focus on error handling in the auth module
/review --no-agy path/to/file.py
```

## How it works

1. **Setup** — creates a scratch `.code_review/` directory with one subdirectory per reviewer (`opus`, `sonnet`, `agy`). agy is included only if the CLI is found on PATH and you haven't asked to skip it.
2. **Detect** — each reviewer in the roster reads the target from disk independently and reports findings as structured JSON (file, lines, severity, claim, evidence, suggested fix).
3. **Consolidate** — duplicate findings across reviewers are merged; each unique finding gets a stable ID.
4. **Validate** — every reviewer in the roster re-reads every consolidated finding and votes `agree` / `disagree` / `unsure`, under a prompt that hides which model originally flagged it and shuffles item order to defeat position bias.
5. **Synthesize** — findings are grouped by severity, each with a confidence label derived purely from the vote split (tiers depend on whether agy is in the roster: `Unanimous`/`Agreed`/`Controversial`/`Rejected`/`Inconclusive` with 3 reviewers, `Agreed`/`Disputed`/`Rejected`/`Inconclusive` with 2), and a concrete fix.

Everything under `.code_review/` is scratch state for a single run — safe to delete between reviews.

## License

MIT
