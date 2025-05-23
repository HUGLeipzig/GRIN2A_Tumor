gnomAD Data: 

- Downloaded all GRIN2A variants from gnomAD v.2.1.1 in a csv
- converted to Excel csv (Semicolon)
- filtered for Transcript ENST00000396573.2 (1566 out of 1573 variants
- combined columns "Filters - exomes" and "Filters - genomes" into "Filters - all" --> All variants passed either Exome or Genome Filters, so no variant had to be exluded --> Deleted Column again
- Added Column "VariantID" according to gnomAD standard (Chromosome-Position-Reference-Alternate
- Saved as 01-1_gnomAD_v2.1.1_ENST00000396573.2_2022_11_14_07_59_07.csv


Cosmic Data

- Downloaded all GrCH37 GRIN2A variants from cosmic (https://cancer.sanger.ac.uk/cosmic/download) --> COSMIC MUtation Data --> Filter by Gene: GRIN2A
- removed whitespace in headers
- Filtered column "GRCH" != null
- Saved as "02-2_V96_37_MUTANT.csv"


ClinVar Data

- Downloaded all GRIN2A variants from ClinVar in a csv Filtered for Variant length < 1kb, single gene
- Splitted the "Name" Column into Transcript, HGVSc and HGVSp
- Filtered the three larger deletions with the g.Code
