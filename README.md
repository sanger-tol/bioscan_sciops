# bioscan_sciops
Script for generation of bioscan sciops manifest given a list of plate IDs

## Conventions
- collection country set as UK
- collection date set as YYYY, filled as 2023 if missing
- G12 is CONTROL_POS
- any other well that is 'blank sample' is CONTROL_NEG_LYSATE
- only plates found in STS returned
- entries sorted by plates (input list) and wells (A1,B1,...,G12,H12)
- other columns contents taken from template manifest
