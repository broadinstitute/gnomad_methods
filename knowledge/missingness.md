# Missingness

Null and NaN are separate concepts in Hail and are handled
*inconsistently* between functions. Most of the silent-fail bugs in this
doc reduce to one of them propagating somewhere you didn't expect.

## `hl.case()` without `missing_false=True` silently drops rows

**`[Hail]`**

```python
case = hl.case()                            # BAD
case = hl.case(missing_false=True)          # GOOD
```

Without `missing_false=True`, a `.when(cond, val)` where `cond` evaluates
to null makes the entire case return null. When that null becomes a
`group_by` key (or feeds into another null-propagating op), Hail silently
drops the row.

Concrete failure: an expression like
```python
hl.case()
  .when(hl.is_missing(cv), "none")
  .when(cv.is_plp, "plp")
  ...
  .default("none")
```
returns null when `cv` is defined but `cv.is_plp` is null (edge case:
partial struct). Whole case → null → dropped from group_by. Symptom:
~184k trio-truth pairs missing from a stratification table with no error.

**Fix:** `hl.case(missing_false=True)` — treats null conditions as False,
falling through to the default. Same fix applies to
`.when(source.contains(tag), tag)` patterns where `source` might be null.

**Symptom checklist:**
- Aggregation total < input rowcount, no error
- One or more "large" cells in a stratification are missing entirely
- Off-diagonal cells populated at zero when they should have entries

---

## `hl.if_else(null_cond, a, b)` returns null, not `b`

**`[Hail]`**

Unlike SQL COALESCE or Python's ternary, Hail's `hl.if_else` propagates
null when the condition is null:
```python
is_v1_rarer = ht.v1_ann.ac <= ht.v2_ann.ac  # null if either .ac is null
rarer_ac = hl.if_else(is_v1_rarer, ht.v1_ann.ac, ht.v2_ann.ac)  # → null!
```

If a downstream analysis expects `rarer_ac` to hold *one of* v1 or v2
(like a min), you get all-nulls when either input is null. Mysteriously
empty buckets in stratifications.

**Fixes:**
- `hl.coalesce(v1.ac, v2.ac)` — if you just want first-defined.
- `hl.if_else(cond, ..., ..., missing_false=True)` — treat null
  condition as False.
- Coalesce inputs before the comparison: `hl.or_else(ht.v1_ann.ac, ...)`.

---

## `group_by` with null-string keys drops rows from labeled outputs

**`[Hail]`**

`hl.agg.group_by(k, ...)` where `k` is a null string: Hail preserves the
null-keyed items in the raw result dict (as `key=None`), but consumers
that build a matrix or filter by known labels (`if k not in labels:
continue`) silently drop them. Compounds with the `hl.case`/`hl.if_else`
issues above.

**Debug pattern:** print the raw `result.items()` keys before matrix
conversion; look for `None` and empty-string keys.

---

## `hl.explode` on a null array/set drops the row

**`[Hail]`**

`ht.explode("gene_id")` when `gene_id` is null (from a left-join miss)
silently drops that row. Not preserved with a null; not errored on; just
vanishes.

If you want to keep the row: `.filter(hl.is_defined(ht.gene_id))`
upstream to make the drop explicit, or coalesce to `hl.empty_array(...)`
— but note that exploding empty gives zero rows too. Use nullability
handling deliberately.

---

## `hl.set` / `hl.dict` operations on null return null, not empty

**`[Hail]` · Verified on Hail 0.2.134**

- `null_set.contains(x)` → null
- `null_set.difference(other)` → null
- `null_set.length()` → null
- `null_dict.get(k, default)` → null, **not** the default (see the
  `dict.get` entry below)

These null returns propagate. In a `filter()`, null → row dropped. In an
`agg_group_by` key, null → row dropped.

**Defensive pattern for potentially-null Set / Array fields:**
```python
source = hl.or_else(ht.v1_ann.source, hl.empty_set(hl.tstr))
gene_id = hl.or_else(ht.v1_ann.gene_id, hl.empty_array(hl.tstr))
```

---

## `dict.get(k, default)` covers an absent key, not a present-but-missing value

**`[Hail]` · Verified on Hail 0.2.134**

```python
hl.dict({"a": 1}).get("b", 5)                       # 5    -- key absent
hl.dict({"a": 1}).get(hl.missing(hl.tstr), 5)       # 5    -- key itself missing
hl.dict({"a": hl.missing(hl.tint32)}).get("a", 5)   # None -- value is missing
hl.missing(hl.tdict(hl.tstr, hl.tint32)).get("a", 5)  # None -- the dict is missing
```

The default fires on *lookup failure*, not on *null value*. A dict built
from a join whose right side had null fields hands those nulls straight
through your default, and a dict field that is itself null ignores the
default entirely. Use `hl.or_else(d.get(k), default)` when you want
"default on null" semantics too.

---

## `hl.parse_json` with a mismatched schema returns nulls, not an error

**`[Hail]` · Verified on Hail 0.2.134**

When you hand `hl.parse_json` a type that doesn't match the JSON, it
does not complain — it produces the closest thing it can:

```python
j = '{"gene":"BRCA1","impact":"HIGH","score":3}'

hl.parse_json(j, hl.tstruct(gene=hl.tstr, impact=hl.tstr, score=hl.tint32))
# Struct(gene='BRCA1', impact='HIGH', score=3)          <- correct

hl.parse_json(j, hl.tstruct(gene=hl.tstr, consequence=hl.tstr))
# Struct(gene='BRCA1', consequence=None)                <- field name wrong -> null

hl.parse_json(j, hl.tstruct(gene=hl.tstr, score=hl.tstr))
# Struct(gene='BRCA1', score='3')                       <- type wrong -> coerced

hl.parse_json(j, hl.tarray(hl.tstr))
# None                                                  <- wholly wrong -> all null
```

A misspelled or renamed field yields a **column of nulls** that looks
exactly like "this annotation is legitimately missing for these
variants," and a wrong scalar type is silently coerced. This is the
mechanism behind hand-written schemas for large nested JSON (VEP output
being the usual suspect) quietly losing fields after an upstream version
bump.

**Seen in production, and note there are two distinct causes.** gnomAD
v3 and v4 both shipped VEP structs with several fully-missing fields
(`ancestral`, `context`, `swissprot`, `trembl`, `uniparc`, the
`*.minimised` ones), and the two explanations are different problems that
present identically (*reported · project meeting*):

1. **The field was never produced.** `ancestral` and `context` come from
   the LOFTEE plugin; without that plugin supplying them there is nothing
   to parse.
2. **The declared structure didn't match the JSON**, so parsing returned
   nulls for the fields that didn't line up — the mechanism above.

Both look like "this annotation is missing for every variant." The fix
in that case was to drop the dead fields from the release, but you can
only make that call once you know which of the two you're looking at.

**Rule of thumb:** after parsing JSON with an explicit type, assert
definedness on a field you *know* is populated —
`ht.aggregate(hl.agg.count_where(hl.is_defined(ht.parsed.some_field)))`
— rather than trusting that the parse "worked." Canary #3 in
[Testing](testing.md#canary-tests) covers this shape.

---

## `hl.max` drops missing but propagates NaN — and the aggregator disagrees

**`[Hail]` · Verified on Hail 0.2.134**

| expression | result |
|---|---|
| `hl.max([1, NA, 3])` | `3` — missing dropped (`filter_missing=True` is the default) |
| `hl.max([1, NA, 3], filter_missing=False)` | missing |
| `hl.max([NA, NA])`, `hl.max(hl.empty_array(...))` | missing |
| `hl.max([1.0, nan, 3.0])` | **`nan`** — NaN is *not* dropped |
| `hl.nanmax([1.0, nan, 3.0])` | `3.0` |
| `hl.agg.max(x)` over a NaN row | **`3.0`** — ignores NaN |
| `hl.agg.mean(x)` over a NaN row | **`nan`** — does not |

Two traps stacked. First, missing and NaN are handled *oppositely* by
`hl.max`: missing is silently skipped, NaN poisons the result. Second,
`hl.agg.max` ignores NaN as well ("for back-compatibility reasons, in
contrast with `hl.max`" —
[agg.max docs](https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.max)),
so the scalar and the aggregator form of "the max" disagree on the same
column. A sanity-check `hl.agg.max` sitting next to a per-row `hl.max`
reports a plausible number while every per-row value is NaN. Identical
story for `min` / `nanmin`.

**Rule of thumb:** on float data that can carry NaN — any ratio with a
possibly-zero denominator, e.g. the EM outputs in the
[experimental-kernel entry](hail-idioms.md#hail-idioms) — reach for `hl.nanmax` /
`hl.nanmin` explicitly. Pass `filter_missing=False` when a missing input
should invalidate the result rather than be quietly skipped.

---

[← back to the index](README.md)
