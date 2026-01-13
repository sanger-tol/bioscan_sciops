# bioscan_sciops
Script for generation of bioscan sciops manifest given a list of plate IDs

DEPRECATED, active version now at https://github.com/sanger-tol/bioscan-ops/tree/master/analysis/bioscan_sciops

## Installation

This script requires python=3.10, tol-sdk, and pandas

## Usage

```
python bioscan_sciops.py -p PLATES [-o OUTFILE] [-l]
```

By default, specimen (LILYS) manifest is generated. Use `-l` to generate lysate (LBSN) manifest instead

## Conventions
- collection country set as UK
- collection date set as YYYY, filled as 2023 if missing
- G12 is set to CONTROL_POS
- any other well that is 'blank sample' is CONTROL_NEG_LYSATE
- only plates found in ToL portal are returned
- entries sorted by plates (input list) and wells (A1,B1,...,G12,H12)

