# gnomad_methods

Shared Hail utility library for gnomAD pipelines. **This is a library, not a pipeline** — it exposes APIs that other repos (`gnomad_qc`, `gnomad-constraint`, and others) import, so public-API changes ripple outward. Prefer additive changes over breaking renames.

## gnomAD & Hail institutional knowledge

`knowledge/` is a shared reference of gnomAD + Hail institutional knowledge maintained by the methods group: silent-fail bug classes, design tradeoffs, release conventions, and facts that are hidden in the source.

**Consult [`knowledge/README.md`](knowledge/README.md) before writing gnomAD- or Hail-touching code.** It is a routing index — read it, then load only the sub-document you need (missingness, Hail idioms, lazy evaluation, partitioning/shuffles, VDS, gnomAD-specific, biology, analysis patterns, Dataproc, release conventions, testing) rather than the whole tree.

It is also where new institutional knowledge belongs. The contribution rules live in that README; the two that matter most are **never add a fact you have not verified** (against the Hail source, the official docs, or a job you actually ran) and **propose entries to the user before writing them**.

## Verbosity: write less prose than you think you should

Applies to docstrings, comments, and the text of a PR description alike. The default failure mode of an AI assistant here is too much prose, not too little.

- **Comments explain *why*, never *what*.** If a line needs a comment to say what it does, rename something instead. A comment restating the code is worse than no comment: it goes stale and it dilutes the comments that matter.
- **Don't comment or docstring code you didn't touch.** A diff that adds explanation to surrounding lines is harder to review and buries the actual change.
- **Docstrings: one summary line, then only what a caller can't infer from the signature.** `:param ht: Input Table.` earns nothing. Spend the words on units, defaults that matter, filters already applied, and invariants the caller must satisfy.
- **No section headers, bullet summaries, or restatements of the request** in an explanation. Say the thing once.
- **When there's nothing to say, say nothing.**

The test: delete a sentence. If nothing is lost, it should not have been there.

## CLAUDE.md vs CLAUDE.local.md

`CLAUDE.md` is committed and holds facts true for anyone working in this repo. `CLAUDE.local.md` is gitignored and holds anything tied to one person's machine — absolute paths, conda env names, cluster names, GCP projects.

Decide by portability: if a new contributor on a different laptop would need it, it belongs here. When in doubt, prefer this file with a placeholder (`<repo-root>`, `<cluster>`) over the local file with a concrete value. See [`knowledge/working-with-ai.md`](knowledge/working-with-ai.md) for the longer version, including how to bootstrap both files.
