# Working with an AI assistant on this codebase

Not gnomAD or Hail knowledge — this is about the tooling around it: how
to keep an assistant's output useful, what wastes tokens, and what
belongs in a project config file versus your own.

The rest of this tree tells an assistant *what is true about the data*.
This page is about *how to work with one*.

## Keeping output terse

The default failure mode is too much prose: restating the request,
narrating each step, docstrings that paraphrase the signature, comments
that restate the line above.

Stating a preference in a config file works, but only where the
assistant is already reading. **Put the rule next to the thing it
governs** — a verbosity rule inside the code-style section of a project's
`CLAUDE.md` lands better than the same sentence in a general preferences
list, because it's in context exactly when style is being decided. A
vague global "be concise" competes with everything else in the file and
gets diluted.

What's worth stating explicitly, because assistants get these wrong by
default:

- Comments explain **why**, never **what**.
- Don't add docstrings or comments to code you didn't touch.
- Don't restate the request before answering it.
- Delete a sentence; if nothing is lost, it shouldn't have been there.

If a rule is being ignored even when it's in context, that's the point
to reach for a **hook** — the only mechanism that enforces rather than
requests. Rules are advisory; hooks run.

## Using compute and tokens well

Two separate costs are easy to conflate: **your token spend** with the
assistant, and **the cluster time** it causes you to spend. The second
is usually the larger number.

**Cluster cost.** Everything in [Lazy evaluation and
materialization](lazy-eval-and-materialization.md#lazy-evaluation-and-materialization)
is a token-cost issue too — a stray `.count()` for logging costs a full
re-run of the upstream query. When reviewing generated pipeline code,
that entry's grep (`.count()`, `.show()`) is the highest-value check you
can do before submitting a job.

**Token cost.** What actually burns tokens, in rough order:

- **Re-reading whole files** to answer a question about one function.
  Ask for a targeted `grep`/search first; read the file only when the
  answer needs surrounding context.
- **Shell inspection loops** — `cat`-ing several files in sequence,
  re-running `ls` to re-orient. One combined command that prints exactly
  what's needed beats five round trips.
- **Permission round-trips.** Every prompt for a read-only command costs
  a turn. Allowlisting genuinely read-only commands in a project's
  `.claude/settings.json` removes that overhead — see the config
  section below.
- **Long outputs you don't read.** Piping to `head`, `tail`, or a
  filtered `grep` is cheaper than dumping a 5,000-line log into context.

Worth knowing: a suggested command you *reject* still costs the tokens
that proposed it. Rejecting a bad plan is not free, so it's cheaper to
be specific up front than to iterate through proposals.

> **To fill in.** This section is deliberately partial. Add what you
> find: model choice per task type, when a fresh session beats
> continuing a long one, whether long-context sessions are worth their
> cost on this codebase, and any measured before/after numbers.

## `CLAUDE.md` vs `CLAUDE.local.md`

Most repos here carry two config files, and the split is by
**portability**, not by secrecy.

| | `CLAUDE.md` | `CLAUDE.local.md` |
|---|---|---|
| tracked | **committed** — everyone gets it | **gitignored** — yours alone |
| holds | facts true for anyone in this repo | facts true only on your machine |
| examples | project structure, pipeline architecture, code conventions, cross-repo coupling, recurring gotchas, invocation patterns with placeholders | absolute paths, conda env names, your cluster names, GCP project IDs, scratch buckets, output-postfix conventions, personal workflow preferences |

**The test:** would this sentence be true for a new team member on a
different laptop? If yes it's `CLAUDE.md`; if it names your machine, it's
`CLAUDE.local.md`. When in doubt, prefer `CLAUDE.md` **with a
placeholder** (`<repo-root>`, `<cluster>`) over `CLAUDE.local.md` with a
concrete value — the former helps everyone.

**Confirm `CLAUDE.local.md` is actually ignored.** Being untracked is not
the same as being ignored: an untracked-but-unignored file is one
`git add -A` away from being committed, along with your paths and cluster
names. Check with `git check-ignore -v CLAUDE.local.md`; if it prints
nothing, add it to `.gitignore`.

### What makes these files work

- **Keep them tight.** These are reference cards, not changelogs. A
  bullet in the right existing section beats a new heading per finding.
  If a section sprawls, consolidate rather than append.
- **Record what you can't re-derive.** Non-obvious API behavior, a
  version gotcha, why a workaround exists, where a resource lives.
  Skip anything a 30-second `grep` would answer.
- **Put rules where they apply.** See the verbosity discussion above —
  placement affects whether a rule is followed.
- **Update them as you learn**, in the same session you learned it.
  The knowledge that never gets written down is the knowledge you
  acquired while busy.

### Bootstrapping your own

To create or refresh these files, hand your assistant this section and a
request along these lines:

> Read `knowledge/working-with-ai.md`, then look through this repo —
> structure, entry points, test setup, CI, and how jobs actually get
> run. Draft a `CLAUDE.md` covering what any contributor would need, and
> a `CLAUDE.local.md` for anything specific to my machine. Sort every
> fact by the portability test in that doc. Show me both before writing
> them, and don't include anything you haven't verified against the
> repo.

The last two clauses matter. Ask to see it first — these files are read
by every future session, so a wrong fact in one is expensive. And an
assistant asked to write a config file will otherwise cheerfully invent
plausible-looking conventions that this repo doesn't follow.

To refresh an existing pair, point at them instead: *"reconcile
`CLAUDE.md` against the current repo — flag anything stale or no longer
true, and propose cuts for anything I could re-derive with a grep."*
Staleness is the main failure mode of these files, and pruning is the
part nobody does.
