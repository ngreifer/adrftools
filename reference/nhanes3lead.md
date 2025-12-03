# Data from NHANES III on blood lead levels and cognitive outcomes

This is a subsample of the data from NHANES III containing observations
of 2521 adolescents that participated in survey, and in particular the
portions of survey that included assessing participants' blood lead
levels and scores on several cognitive tests.

## Usage

``` r
nhanes3lead
```

## Format

An object of class `data.frame` with 2521 rows and 14 columns.

## Details

The treatment is `logBLL`, the natural log of blood levels measured in
μg/dL. The covariates are `Age`, `Male`, `Race`, `PIR`, `Enough_Food`,
`Smoke_in_Home`, `Smoke_Pregnant`, and `NICU`. The outcomes are `Math`,
`Reading`, `Block`, and `Digit`. `MEC_wt` are sampling weights.
