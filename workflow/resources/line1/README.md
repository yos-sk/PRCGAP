# LINE-1 reference set for the `run_line1` annotation step

## Build this first

The three FASTAs are **not tracked** (`.gitignore`), so a fresh clone has only
this README and the build script. `run_line1` defaults to true and `line1_hap.sh`
reads all three directly, so the step fails until they exist:

```bash
bash workflow/resources/line1/scripts/build_line1_resources.sh
```

It needs `images/nanomonsv.sif` and nothing from the network. Run it once per
clone; the output is deterministic, so re-running it is also how you verify an
existing set.

| file | contents | source |
|---|---|---|
| `L1.3.fa` | L1.3, the canonical retrotransposition-competent full-length L1HS (6,059 bp) | GenBank **L19088** |
| `l1_3end.fa` | 67 `*_3end` LINE/L1 subunit models, human lineage | Dfam (CC0), via RepeatMasker's library in `images/nanomonsv.sif` |
| `l1_5end.fa` | 61 `*_5end` LINE/L1 subunit models, human lineage | same |

`L1.3.fa` is the detection query: it is 95.5-99% identical to L1PA2-L1PA5, so one
query finds every young subfamily.

The subunit models do two jobs. The `_3end` set decides the subfamily — that is
how RepeatMasker itself names a young L1, and scoring those 67 models reproduces
its labels at 97% on held-out loci. Including all 67 rather than only the five
wanted families is what rejects older elements: an L1PA6 wins on its own model
and is dropped instead of being forced into L1PA5. The `_5end` set sets the 5'
boundary, which the whole-length query cannot do because blastn does not reach
the 5'UTR of the older subfamilies.

Dfam has no full-length model for L1PA2-L1PA5 (22 L1PA entries, all `_3end` or
`_5end`), which is why the detection query is L1.3 and not a per-family
consensus. RepBase has full-length consensus sequences but needs a license.

The RepeatMasker library covers every species, so the build script keeps only
the clades on the path from the root to *Homo sapiens* and drops the 51
rodent-specific L1 subunit models. Note the library spells one of those clades
`Theria_mammals`, not `Theria`.

Validation and the parameter choices behind this step are in
`workspace/repeat_bench/README.md`, Experiment 8.
