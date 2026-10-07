---
name: review
description: >
  Performs code review by launching Claude twice under two different models —
  Opus and Sonnet — plus (optionally) agy (Google Antigravity), all in
  parallel. Then, it merges their findings into a deduplicated list and
  launches a validation pass, asking each reviewer in the roster to agree or
  disagree with the others' findings. Then, it gives you a final list of
  findings so that you can choose which ones to fix.
---

# Review Skill

Two Claude reviewers — one pinned to **Opus**, one pinned to **Sonnet** — plus an optional third reviewer, **agy** (Google Antigravity), detect issues in parallel by reading local files directly. Every consolidated finding is then scored by every reviewer in this run's roster (each running locally with disk access) using a deliberately neutral prompt — no reviewer attribution, no meta-framing about peer review. The prompt removes explicit authorship signals. Findings are presented with confidence labels derived purely from validator votes.

agy is optional. The skill auto-detects whether the `agy` CLI is installed and drops it from the roster if not, and you can also ask to skip it explicitly (`/review --no-agy`, or "skip agy" / "without agy" anywhere in your instructions) even when it is installed. With agy in the roster there are 3 reviewers/validators (confidence tiers: Unanimous / Agreed / Controversial / Rejected / Inconclusive); without it there are 2 (Agreed / Disputed / Rejected / Inconclusive).

## Requirements

- `claude` CLI (this skill runs inside it, so it's already present) — invoked twice per step, once with `--model opus` and once with `--model sonnet`.
- `agy` CLI ([Google Antigravity](https://antigravity.google/)), installed and authenticated — **optional**. If missing (or explicitly skipped), the skill runs with the two Claude reviewers only.

## Workflow

> **Working-directory hygiene.** Every shell command in this skill assumes the project root as cwd. Because the harness's bash cwd can drift across tool calls (a prior `cd subdir` from earlier in the session can persist), each command below anchors itself with `cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)" && …` (or its parametrized form). Do not strip those anchors — they're what make the skill robust to a drifted cwd. The expression is `git rev-parse --show-toplevel` for git repos (finds the project root regardless of harness cwd), falling back to `pwd` for non-git projects. **Caveat for non-git projects:** the `pwd` fallback only works if the harness cwd happens to be the project root — if it has drifted, you must `cd` back to the project root yourself first, or replace the whole `$(…)` expression with a hardcoded absolute project-root path inside each command before running it.

### Step 1: Setup

```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)" && mkdir -p .code_review/opus .code_review/sonnet .code_review/agy
```

**Decide whether agy is in this run's roster.** Include it only if BOTH:
1. The user did not ask to skip it — scan the raw `/review` argument (before Step 2 parses it for a target) for `--no-agy`, or plain language such as "skip agy" / "without agy" / "don't use agy". If found, strip it from the argument string before Step 2 resolves the target, so it isn't mistaken for a path.
2. `command -v agy >/dev/null 2>&1` finds it on PATH.

If either check fails, proceed with the **opus + sonnet roster only** for the rest of this run. Tell the user which roster is in effect before launching anything (e.g. "agy not found on PATH — reviewing with opus + sonnet only" or "skipping agy per your request"). Unlike opus and sonnet — core to every run — agy is purely additive: if it's out, there's no stand-in to find, just a smaller, still-complete roster.

Generate a unique run id `XX` as the current timestamp in the format `YYYYMMDD_HHMMSS`:

```bash
XX=$(date +%Y%m%d_%H%M%S)
```

Use this same `XX` value as the prefix for every file the skill writes under `.code_review/` during this run (prompts, the completion marker, etc.) — capture it once and reuse it so all artifacts from a single run share a stable prefix.

### Step 1b: Decide whether to split into independent subsets

Large projects can exceed what the detection models can usefully attend to in a single review, even when each agent reads the files from disk itself — review quality degrades long before the context-window hard limit, and a single agent juggling 50 files reasons more shallowly than splitting the work across reviewers. When the target is large, the review is split into multiple subsets that are reviewed in parallel as independent pipelines, then merged at synthesis (Step 6).

**Resolve the file list and measure the target.** First resolve the user's argument the same way Step 2 will (explicit paths, `git diff main...HEAD`, or all files under a directory via `git ls-files`). Then compute:
- Total source lines: `wc -l <files>` (sum the `total` row).
- Total file count.
- Approximate total byte size: `wc -c <files>` summed, or `du -bc <files> | tail -1`.

**Threshold for splitting (moderate):**
- ≤ ~6,000 source lines **AND** ≤ ~40 files **AND** ≤ ~300 KB → **single subset called `all`**, proceed to Step 2 unchanged.
- Otherwise → **plan a split** before launching reviewers.

**Plan the split.** The goal is 2–5 subsets where each subset (a) fits comfortably under the threshold above and (b) is internally cohesive and externally as independent as possible from the other subsets. Apply these heuristics in order of preference:

1. **By top-level subdirectory** (`src/`, `scripts/`, `tests/`, language-specific dirs, etc.) — usually the cleanest natural boundary.
2. **By language / runtime** (e.g. Python build scripts vs. browser JS/HTML/CSS) — often correlates with module boundaries.
3. **By named feature module** when a feature has its own subdir and few cross-imports.
4. **By coarse import/require analysis** — quickly grep for cross-file references; group strongly-connected files together.

Sanity-check the split: each subset should be reviewable without constantly cross-referencing other subsets. If two candidate subsets are deeply interdependent (one imports many symbols from the other, shared types, etc.), merge them. When a cross-subset reference is unavoidable, you may include the referenced file in BOTH subsets if it is small (under a few hundred lines) and the dependency is non-trivial; otherwise note the dependency in a one-line `> Out-of-subset dependency: <file>` header inside the dependent subset's prompt so reviewers know what's intentionally out of scope.

**Oversize fallback.** If a single file or directory is itself larger than the threshold and cannot be cleanly split along module boundaries, review it as its own subset anyway. Append a one-line `> Note: this subset exceeds the recommended review size; reviewer attention may be degraded.` header to that subset's prompt so the user sees the caveat in the synthesis output. Do NOT mechanically chunk a file mid-way.

**Output of this step:** a list of subsets, each with:
- A short kebab-case name (e.g. `web-ui`, `build-scripts`, `audio-pipeline`).
- The file list.
- A one-sentence rationale for the boundary.

Print the planned splits to the user as a brief preview before launching reviewers (subset name → file count, short rationale) so they can interrupt if the boundaries look wrong.

**Threading subsets through Steps 2–5.** When more than one subset is planned, run Steps 2–5 **once per subset**, with these naming conventions:
- The Step 2 `<label>` becomes `<base-label>__<subset-name>` (e.g. `swedish_flashcards__web-ui`). All existing `XX_PROMPT__<label>__…` and `XX_PROMPT__<label>_…_validate.md` paths in Steps 2–5 keep their exact shape; only the `<label>` value carries the subset.
- The reviewer/validator subdirectories (`.code_review/opus`, `.code_review/sonnet`, `.code_review/agy`) are shared across subsets — each subset's CLI invocations use a distinct prompt-file path, so parallel runs from the same subdirectory don't collide as long as we rely on captured stdout (we do; no CLI writes files we depend on later).

**Concurrency cap.** Launching every (subset × reviewer) at once can trip rate limits and overload the machine. **Cap at ~6 concurrent CLI invocations when agy is in the roster (3 reviewers), ~4 when it isn't (2 reviewers).** Batch the work: e.g. with 3 subsets × 3 reviewers = 9 runs, launch 6 first and queue the remaining 3 to start as earlier ones finish. Use `run_in_background: true` with `Monitor` to wait, or chain sequentially after a first batch returns. The same cap applies to the Step 5 validation pass.

If only a single subset is planned (project is below threshold), the rest of the skill behaves exactly as it did before splitting was introduced.

### Step 2: Resolve the target and write the detection prompt

> **Per subset.** When Step 1b planned multiple subsets, perform this step once per subset using that subset's file list as the resolved `{TARGET}`. Express `{TARGET}` as an explicit file list (e.g. "the files `path/to/foo.py` and `path/to/bar.py` (paths are relative to the repo root)") so the agent reviews only the subset's files. Use `<label>` = `<base-label>__<subset-name>` in the prompt-file paths. Write all subsets' prompts before moving on to Step 3.

Parse the user's argument to `/review` (after the Step 1 skip-agy scan has already stripped `--no-agy`, if present):
- If it names files/directories or a git ref (e.g. `--diff main...HEAD`), use that.
- Otherwise default to `git diff main...HEAD` if there's a non-empty diff, else "the code in the current directory".

Resolve `{TARGET}` to a concrete description the agent can act on, e.g.:
- "the files `path/to/foo.py` and `path/to/bar.py` (paths are relative to the repo root)"
- "the changes in `git diff main...HEAD` (run the diff yourself to see them, then read the changed files in full)"
- "all source files under `./analysis/` (use `git ls-files analysis/` to enumerate, then read each)"

Write one detection prompt per reviewer in the active roster — `.code_review/XX_PROMPT__<label>__opus.md`, `.code_review/XX_PROMPT__<label>__sonnet.md`, and, if agy is in the roster, `.code_review/XX_PROMPT__<label>__agy.md` — using the template below with `<agent>` substituted to the agent's actual subdirectory name (`opus`, `sonnet`, or `agy`) in each file. The prompts are identical EXCEPT the agy prompt additionally gets the **agy narration-suppression block** (defined below the template) appended at its very end. (A single shared file would leave the literal `<agent>` placeholder in the prompt, which is confusing to the model even if path resolution via `../..` still works.)

````
Do a thorough code review of {TARGET}.

**Read the files from disk yourself.** You are running from a subdirectory two levels below the repo root; the repo root is `../..` (your cwd is `.code_review/<agent>/` under the repo root). Use `../..` to construct all file paths. Start by listing what's in scope (e.g. `git -C ../.. ls-files` or the explicit paths above, prefixed with `../../`), then read each file in full before forming opinions.

For each candidate issue:

1. Pinpoint the exact file and line(s).
2. Re-open the file and confirm the cited lines actually say what you claim.
3. Confirm the failure mode is reachable: trace at least one caller / test / config that exercises this code, and confirm the bad input can occur in practice. If you cannot, do not report the finding. Speculative "this would crash if X were None" concerns when every caller statically prevents X from being None are the dominant false-positive class — drop them. **Do not invent hypothetical inputs to argue reachability** — phrases like "if the test used X" or "if a caller passed Y" mean you have not found a real trigger. Cite a concrete input an existing caller / test / config actually produces, or drop the finding. A hypothetical that would-trigger-the-bug-if-it-existed is not evidence of reachability; it is a confession that you couldn't find one.

**Severity bar — report only:**
- `critical` — wrong output, data loss, security issue, or crash on a realistic input.
- `high` — likely bug under normal use, OR significant maintainability hazard.
- `medium` — verifiable quality issue (dead code with no caller; duplicated logic that has actually drifted; documentation that contradicts the implementation in a misleading way).

Do NOT report `low` / style / nits.

**Specifically DO NOT flag:**
- Missing type hints, formatting.
- Defensive checks for inputs that cannot occur given the call sites.
- "Consider extracting a helper" / suggested abstractions where existing code is clear and used in one place.
- Hypothetical security / untrusted-input concerns when the code is a one-shot script or runs in a trusted environment.
- Pre-existing issues outside the target.
- Anything a linter would catch.

**Output format — STRICT.** Before the JSON, you may include a brief `<scratch>…</scratch>` block listing the files you inspected and any candidates you ruled out — use this to free up exploration, not as part of the answer. The final answer must be a single JSON object inside one ```json fenced code block, emitted as the LAST ```json block in your response. **If multiple ```json blocks appear in your response (e.g. you echoed the example schema), only the LAST block is parsed.** Anything outside the parsed block (including the scratch) is discarded.

```json
{
  "summary": {"critical": 0, "high": 2, "medium": 1},
  "findings": [
    {
      "id": "F1",
      "file": "path/to/file.py",
      "lines": "42-45",
      "severity": "high",
      "title": "Short one-line title (≤ 100 chars)",
      "claim": "Self-contained 1-2 sentence statement of what is wrong, grounded in the code. Another reader who has not seen your other findings must be able to judge it from this field alone. Cite line content where useful.",
      "evidence": "Concrete reachability: which caller/test/config triggers it, and the input that produces the bad outcome. Cite file:line.",
      "suggested_fix": "Brief concrete fix, or empty string if not obvious"
    }
  ]
}
```

Field rules:
- `id`: `F1`, `F2`, … in listed order.
- `lines`: `42`, `42-45`, or `42,50,73`.
- `severity`: exactly one of `critical`, `high`, `medium`.
- `claim` and `evidence` will be shown verbatim to peer validators — keep them self-contained and grounded.
- If you find no issues, return `{"summary": {"critical": 0, "high": 0, "medium": 0}, "findings": []}`.
````

If the user passed extra instructions beyond a target, append them after the template (do not modify the output schema). **Do not** append a project file tree or file contents — each agent reads the disk itself.

**agy narration-suppression block (append to the agy prompt ONLY — both detection here and validation in Step 5).** agy's print mode streams step-by-step tool narration to stdout ("I will view file X", "I will search for Y", …) and can exhaust its `--print-timeout` budget on narration before it ever emits the final answer (observed: a detection run that printed dozens of "I will view …" lines then died on `Error: timed out waiting for response` with zero ```json blocks). agy has no CLI flag to suppress this, so instruct it in-prompt. Append this verbatim as the last lines of the agy prompt file:

```
**Suppress narration — IMPORTANT.** Do NOT print step-by-step progress narration, plans, or reasoning traces (e.g. "I will view file X", "I will search for Y", "Next I will check Z"). Read the files you need SILENTLY using your tools, then emit ONLY the final answer. Any such narration wastes your limited print-timeout budget and risks the run timing out before the answer is produced. Your entire stdout should be (optionally) the brief `<scratch>` block followed by the single final ```json block — nothing else.
```

(The reference to `<scratch>` matches the detection template; for the Step 5 validation prompt, which has no `<scratch>` block, drop that clause so the last sentence reads "Your entire stdout should be the single final ```json block — nothing else.")

### Step 3: Launch the active roster's detection agents in parallel

All runs use `run_in_background: true`. Each agent runs in its own subdirectory (`.code_review/opus`, `.code_review/sonnet`, and, if in the roster, `.code_review/agy`) two levels below the repo root, so it reads the repo from `../..`. **When multiple subsets exist, launch each subset's agents the same way, but obey the Step 1b concurrency cap — batch the remainder and start queued runs as earlier ones finish.**

**Opus — pinned model, fast non-interactive:**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/opus" && cat ../XX_PROMPT__<label>__opus.md | claude -p --model opus --output-format text
```

**Sonnet — pinned model, fast non-interactive:**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/sonnet" && cat ../XX_PROMPT__<label>__sonnet.md | claude -p --model sonnet --output-format text
```

**agy (only when in this run's roster) — permissions auto-approved, 10-minute print timeout (default is 5m, which can be short for review-sized prompts):**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/agy" && OUTPUT=$(cat ../XX_PROMPT__<label>__agy.md | agy -p --print-timeout 10m 2>&1); STATUS=$?; echo "$OUTPUT"; if [ -z "$(echo "$OUTPUT" | tr -d '[:space:]')" ] || [ "$STATUS" -ne 0 ]; then echo "AGY_FAILED"; fi
```

If agy prints `AGY_FAILED`, drop it from this run's roster and proceed with opus + sonnet only — no stand-in needed, since that pair is already a complete two-reviewer roster on its own. If opus or sonnet itself fails (non-zero exit, empty output) — unusual, since both run through the same `claude` CLI executing this skill — proceed with whichever survived (plus agy, if in the roster); do not fabricate output for the missing slot.

### Step 4: Consolidate findings

> **Per subset.** When Step 1b planned multiple subsets, perform consolidation independently per subset (each subset has its own set of agent JSON outputs). Use subset-local finding ids like `<subset-name>:X1`, `<subset-name>:X2`, … so the subset is obvious downstream. Cross-subset deduplication is unnecessary because subsets contain disjoint files; if a finding happens to span subsets (e.g. an out-of-subset dependency you exposed via the `> Out-of-subset dependency:` header), assign it to whichever subset's reviewers raised it and leave the other subset alone.

After every agent in the active roster returns, parse its JSON. **If the response contains multiple ```json fenced blocks (e.g. the model echoed the example schema), use the LAST block as the answer** — the templates instruct the model to emit its answer last, and any earlier block is an echo. If an agent returned prose, salvage what you can and flag the run as less reliable.

- `AGY_FAILED` → proceed with opus + sonnet findings only.
- If opus or sonnet failed with no counterpart to fall back on → proceed with the survivor(s).

Deduplicate across agents (same file, overlapping lines, same root cause → one merged finding). Be conservative — only merge when clearly the same issue. When merging, take the higher severity.

For each merged finding, assign a stable cross-val id `X1`, `X2`, … and record `file`, `lines`, `severity`, `title`, `claim`, `evidence`, `suggested_fix` (most self-contained variant when multiple agents flagged it). **`suggested_fix` must end up non-empty for every finding** — if every agent that flagged it left `suggested_fix` empty, write a concrete one-line fix yourself from the `claim`/`evidence` before moving on. Step 6 always prints a planned fix per finding, so there must be one to print.

The set of unique findings — regardless of who flagged each one — is the input to Step 5. **Every reviewer in the active roster validates every finding** (including findings it originally flagged); the prompt removes explicit authorship signals (no reviewer name, no ordering hint, no "you said" framing). Note that the verbatim `claim` text retains its original stylistic fingerprint, so this is a soft (not absolute) protection against self-recognition bias.

### Step 5: Peer validation pass

> **Per subset.** When Step 1b planned multiple subsets, run a separate validation pass per subset, scoring only that subset's consolidated findings. Validator prompts use `<label>` = `<base-label>__<subset-name>`. Across all subsets the same concurrency cap from Step 1b applies; batch validators if running many subsets at once. Findings from different subsets are never mixed in a single validator prompt — that would force a validator to chase paths outside the subset it was scoped to.

Every reviewer in the active roster (opus, sonnet, and agy if included) runs against **every** consolidated finding — including findings it originally flagged. The prompt is deliberately neutral and hides any framing about who proposed the issue. **All validators read from disk — no file contents are inlined**, since every CLI in this skill has disk access.

Write `.code_review/XX_PROMPT__<label>_<agent>_validate.md`:

````
For each item below: read the cited source code and answer — do you agree with the finding? Provide a reason grounded in the actual code.

Your cwd is `.code_review/<agent>/`. The source lives at `../..`. Read each cited file from disk before answering, prefixing paths with `../../`. Read freely within the target project directory and its subdirectories. Consult files outside the project only reluctantly, when a finding genuinely cannot be settled from inside the project.

Output: a single JSON array in one ```json fenced code block, in the same order as the items below. No preamble. If multiple ```json blocks appear in your response (e.g. you echoed the example), only the LAST block is parsed. Each entry: `{"issue_id": "...", "verdict": "agree" | "disagree" | "unsure", "reason": "<cite specific code>"}`. Use `unsure` only when the cited code does not let you reach a confident answer.

```json
[
  {"issue_id": "X1", "verdict": "agree",    "reason": "foo.py:42 calls bar() before bar is defined; bar is only set inside the `if cfg:` branch on line 47. Falsy cfg → NameError."},
  {"issue_id": "X2", "verdict": "disagree", "reason": "Line 98 has `if x is None: return`, so x is provably non-None at line 100."},
  {"issue_id": "X3", "verdict": "unsure",   "reason": "The cited code does not show enough of the surrounding control flow to decide."}
]
```

Items:

<for each consolidated finding, render in shuffled order, using ONLY the fields below — no reviewer names, no ordering hints, no framing about where the finding came from>
### X<N>
- file: <file>
- lines: <lines>
- finding: <claim>
- reachability: <evidence>
````

Per-agent prompt rules:
- **Shuffle the order of items independently for each validating agent.** This defeats position bias.
- Every reviewer in the active roster gets the same full list of consolidated findings (including ones it originally flagged at detection).
- **No reviewer attribution and no meta-framing.** Do not name any original finder; do not say "you wrote this" or "you flagged this"; do not reveal that this prompt follows an earlier *detection* pass or that another agent produced these claims. The agree/disagree question itself IS the task — `do you agree with the finding?` is fine and intentional; what's forbidden is exposing the multi-agent / peer-review workflow context.
- Include `file`, `lines`, `finding` (the `claim` text), and `reachability` (the `evidence` text) per item — the detection schema promises validators will see both. Severity and any other Step-4 metadata are still omitted to avoid anchoring the verdict.
- Do NOT inline file contents — every validator reads from disk.
- **agy validator prompt ONLY (when agy is in the roster):** append the **agy narration-suppression block** (defined in Step 2, using the no-`<scratch>` wording) as the last lines of `XX_PROMPT__<label>_agy_validate.md`, for the same timeout reason as detection. The opus and sonnet validator prompts do not get it.

**Launch the active roster's validators in parallel.** Each Claude validator uses the same pinned model as its own detection run — opus stays opus, sonnet stays sonnet. Pinning matters here: the point of splitting Claude into two slots is to compare Opus's and Sonnet's independent judgment on every finding, not just add a second inference pass at whatever the CLI's default happens to be.

**Opus validator:**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/opus" && cat ../XX_PROMPT__<label>_opus_validate.md | claude -p --model opus --output-format text
```

**Sonnet validator:**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/sonnet" && cat ../XX_PROMPT__<label>_sonnet_validate.md | claude -p --model sonnet --output-format text
```

**agy validator (only when in this run's roster):**
```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)/.code_review/agy" && OUTPUT=$(cat ../XX_PROMPT__<label>_agy_validate.md | agy -p --print-timeout 10m 2>&1); STATUS=$?; echo "$OUTPUT"; if [ -z "$(echo "$OUTPUT" | tr -d '[:space:]')" ] || [ "$STATUS" -ne 0 ]; then echo "AGY_FAILED"; fi
```

If the agy validator returns `AGY_FAILED`, drop it from validation for this run and proceed with the opus + sonnet verdicts using the two-slot tier scheme (Step 6) — this can happen even when agy succeeded at detection, since validation is a separate CLI call; treat detection and validation rosters as decided independently.

If opus or sonnet's validator run fails, proceed with the survivor(s); do NOT manufacture a verdict for the missing slot.

Parse each validator's JSON into per-issue verdicts. You now have, for each finding, one verdict per validator that actually ran (2 or 3).

### Step 6: Synthesis

> **Multiple subsets.** When Step 1b planned multiple subsets, the synthesis report combines all subsets into a single output. Re-number findings into a global `X1`, `X2`, … sequence and append `(subset: <subset-name>)` after each X-id label so the reader can see where each finding came from. Order severity tiers globally (all Critical findings first across all subsets, then High across all subsets, etc.). If any subset had its `> Note: this subset exceeds the recommended review size` header, surface that caveat once near the top of the report. If a subset returned zero findings at any tier, you do not need a per-subset placeholder — the tier-level `None.` rule still applies globally.

For each finding, compute a confidence label **purely from validator verdicts**. Who originally flagged the issue at detection is not used at synthesis — every active validator scored every finding.

Each validator returns one of `agree`, `disagree`, `unsure` per finding (those are the only valid verdicts — "unverified" is not a verdict).

**Pick the tier scheme once per run, based on how many validator slots actually returned a verdict** (not how many were originally planned — a mid-run `AGY_FAILED` at validation still counts as "didn't run"):

**Three validators ran (opus, sonnet, agy):** apply top-down, first match wins:
1. **Rejected** — 2 or more `disagree`, no `agree`.
2. **Unanimous** — all three returned `agree`.
3. **Controversial** — at least one `agree` AND at least one `disagree` (a genuine split among voters). Show every dissenting reason inline.
4. **Agreed** — 2 or more `agree`, no `disagree`.
5. **Inconclusive** — everything else (unsure-heavy, or too few opinions to decide).

**Two validators ran (opus, sonnet only — agy skipped, unavailable, or dropped mid-run):** apply top-down, first match wins:
1. **Rejected** — both slots ran and both returned `disagree`.
2. **Agreed** — both slots ran and both returned `agree`.
3. **Disputed** — both slots ran, one `agree` and one `disagree`. Show both reasons inline.
4. **Inconclusive** — everything else: any `unsure` involved, OR only one validator slot ran at all. A single surviving verdict is never enough on its own to reach Agreed/Rejected/Disputed — it needs corroboration from a second validator. Show the lone verdict and reason so the reader can weigh it themselves.

**Do NOT re-judge severity at synthesis.** The reviewers/validators saw the code; you're working from JSON. Your job is narrow: deduplicate, label confidence, present. Use the merged severity from Step 4. If you think a verdict is weak, surface that as a user-facing note next to the finding — do not silently downgrade or drop it.

**Present** grouped by severity tier. For each tier, list findings or write `None.` explicitly. Every bullet starts with the X-id in bold. The label lists each active validator's verdict in a fixed order (Opus, Sonnet, agy — dropping agy from the list entirely when it isn't in this run's roster) so the reader can see at a glance how each model voted — but do NOT mention who originally flagged the issue.

**Every finding gets a `Fix:` line — no exceptions.** Print it as a second line under the bullet, using the finding's (possibly synthesizer-authored, per Step 4) `suggested_fix`. This lets the reader decide apply/skip from the report alone, without re-opening each agent's raw JSON. For a `Rejected` finding, replace the fix text with a one-line note that no fix is planned since validators treat it as a false positive — do not invent a code fix for something the review concluded isn't real.

```
## Critical (must fix immediately)
- **X1** [Unanimous: Opus agree, Sonnet agree, agy agree] Description… (file:line)
  Fix: <concrete one-line fix>

## High (likely bugs / significant hazards)
- **X2** [Agreed: Opus agree, Sonnet agree] Description… (file:line) — agy not in this run's roster.
  Fix: <concrete one-line fix>
- **X3** [Controversial: Opus agree, Sonnet disagree — "<reason>", agy agree] Description… (file:line)
  Fix: <concrete one-line fix>
- **X4** [Disputed: Opus agree, Sonnet disagree — "<reason>"] Description… (file:line)
  Fix: <concrete one-line fix>

## Medium (quality issues)
- **X5** [Rejected: Opus disagree — "<reason>", Sonnet disagree — "<reason>", agy disagree — "<reason>"] Description… (file:line) — included for transparency; no validator agrees, treat as a likely false positive.
  Fix: None planned — treated as a false positive.
- **X6** [Inconclusive: Opus agree, agy unavailable] Description… (file:line) — agy's validator failed mid-run; only one validator vote available.
  Fix: <concrete one-line fix>
```

For Disputed / Controversial / Inconclusive items, always include each non-`agree` validator's short reason in quotes.

**Offer to fix.** Ask which findings to apply, by X-id (e.g. "apply X1 and X3").

### Step 7: Mark the review complete

Once the user has triaged the findings — either by responding to all or most of them (apply / dismiss / "skip the rest") or explicitly signalling they're done — write an empty marker file so future runs (and the user) can tell that this review was closed out:

```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)" && touch .code_review/XX__CODE_REVIEW_COMPLETED
```

`XX` is the same run id chosen in Step 1, so the marker sits alongside that run's prompts under `.code_review/`. Do NOT write this marker preemptively after just presenting the synthesis — it means the user has actually engaged with the findings, not merely seen them.

## Notes

- Large projects are split into independent subsets at Step 1b (~6k lines / ~40 files / ~300 KB threshold), with Steps 2–5 running per subset in parallel under a concurrency cap, and Step 6 merging all subsets into a single global report.
- The two Claude slots are **pinned** (`--model opus`, `--model sonnet`) in both detection (Step 3) and validation (Step 5) — the whole point of the split is to compare Opus's and Sonnet's independent judgment, not to add a second run at whatever the CLI's default model happens to be.
- agy is optional and purely additive. It's included only when the `agy` CLI is found on PATH AND the user hasn't asked to skip it (`--no-agy`, "skip agy", "without agy"). If it fails at detection or validation, it's simply dropped from the roster for that step — there's no stand-in, since opus + sonnet is already a complete two-reviewer roster on its own.
- **Every reviewer in the active roster validates every finding** (including ones it originally flagged at detection). The prompt removes explicit origin signals — no reviewer name, no ordering hint, no "you said" framing. Self-recognition via stylistic fingerprints in the verbatim `claim` text is a known residual risk (frontier evaluators recognize their own outputs above chance — Panickssery et al. 2024), so treat this as a soft, not absolute, protection.
- Claim order is shuffled independently per validator to defeat position bias.
- Synthesis labels are derived **purely from validator verdicts**, not from who originally flagged the finding. Detection only determines what gets into the candidate list; validation determines confidence.
- Validators are told to read freely inside the target project directory; consulting files outside the project is allowed only reluctantly when a claim genuinely cannot be settled from in-project code.
- The confidence-tier scheme (Step 6) is chosen per run based on how many validator slots actually returned a verdict — 3 slots use Unanimous/Agreed/Controversial/Rejected/Inconclusive; 2 slots use Agreed/Disputed/Rejected/Inconclusive. Detection and validation rosters are decided independently, so it's possible for agy to contribute detection findings but drop out of validation (or vice versa).
- If opus or sonnet fails with no counterpart to fall back on, proceed with the survivor(s) — never invent a verdict for the missing slot. A single surviving validator verdict on its own always lands in the Inconclusive tier (see Step 6) since it lacks corroboration.
- If a review is requested multiple times, write each prompt as if it's the first review — do not mention prior reviews.

