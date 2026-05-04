# clocker v2.2.0 -> v2.2.4 — complete bundle

This is a **complete package** including all fixes since v2.2.0. Just
drop into your repo (overwriting), run `devtools::document()` and
`devtools::install()`, and push.

## Why a complete bundle this time

Looking at your last run output, the v2.2.3 incremental patch I sent
was clearly never applied — the install pulled commit `9d82560` which
doesn't include the EpiDISH CP capture or the PC-Clocks message
capture. To eliminate any ambiguity, here's the complete current
state of every file.

## What's fixed since v2.2.0

| Fix                                                | First in   |
| -------------------------------------------------- | ---------- |
| roxygen `[0, 1]` link warning                      | v2.2.1     |
| EpiMitClocks `dataETOC3.l not found` error         | v2.2.1     |
| EpiMitClocks `data set X not found` warnings       | v2.2.2     |
| EpiMitClocks `[1] "Number of represented..."`     | v2.2.2     |
| EpiDISH CP `1\n2\n...\nN` sample-number flood      | v2.2.3     |
| methylCIPHER PC-Clocks `Calculating PC Clocks now` | v2.2.3     |
| `epiTOC2: Using N of M CpGs` ignoring `verbose`    | **v2.2.4** |
| PC-Clocks Age/Female fallback alerts ignoring      | **v2.2.4** |
|   `verbose` (now: message under verbose, warning   |            |
|   otherwise — so the alert is always visible but   |            |
|   doesn't pollute stdout)                          |            |

## Verbose vs. quiet behavior

After this patch:

**`verbose = TRUE`**: full log including all upstream-package output
(captured and reformatted as indented `message()` lines). What you
saw working in your real-data run.

**`verbose = FALSE`**: completely silent stdout. Pheno fallback alerts
("PC-Clocks used Horvath2 estimate as Age", etc.) become standard R
warnings, visible only via `warnings()` or `last_warnings()`. Any
upstream-package print() / cat() output is captured and dropped.

## To apply

```r
setwd("/Users/bc/Documents/GitHub/clocker/")

# Delete old R/ files first to be safe (no stale files left behind)
unlink("R/*.R")

# Drop in all files from this bundle (R/, DESCRIPTION, NAMESPACE, etc.)
# Or just `cp -r clocker_v2.2.4/* .` from the unzipped bundle

devtools::document()       # regenerates man/
devtools::install()        # rebuild
git add -A
git commit -m "v2.2.4: silence upstream package noise; respect verbose flag"
git push

# Then on the install machine:
# detach if loaded
detach("package:clocker", unload = TRUE)
remove.packages("clocker")
pak::pak("brianchengithub/clocker")
```

## Verifying it worked

After the new install, this should produce ZERO bare numbered lines
or unindented "Calculating PC Clocks now" output, even on a 39-sample
real-data run with verbose = FALSE:

```r
library(clocker)
result2 <- clocker(beta_matrix, verbose = FALSE)
# (output should be empty)
warnings()  # may show the Age/Female fallback alerts as warnings
```

With verbose = TRUE you should see the well-indented log we discussed
in the previous round.
