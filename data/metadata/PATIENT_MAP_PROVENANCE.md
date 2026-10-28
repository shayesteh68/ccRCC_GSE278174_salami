# Patient map provenance

The `patient_map.tsv` crosswalk was generated from Supplementary Information 3 for the GSE278174 study. Only the `Phenotype` worksheet columns `Library.ID` and `Myriad` were used. `Library.ID` was matched to `coldata.tsv` `sample_id`; the corresponding de-identified `Myriad` case code was recorded as `patient_id`. No other supplementary fields were used to create the crosswalk.

Source: Cotta et al., *Scientific Reports* (2026), DOI: 10.1038/s41598-026-35021-7; Supplementary Information 3 (`41598_2026_35021_MOESM3_ESM.xlsx`). The downloaded source checksum was verified against the checksum in the article's PMC full-text XML (MD5 `5343c3663e1615f709483399d0a73d36`).

Validation: all 42 `sample_id` values match `coldata.tsv` exactly, with no duplicate or missing mapping; there are 19 distinct patient codes, each with exactly one primary sample and at least one metastatic sample (19 primary, 23 metastatic in total).
