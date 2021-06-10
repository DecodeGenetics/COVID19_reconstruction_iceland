## Reconstruction of a large-scale outbreak of SARS-CoV-2 infection in Iceland informs vaccination strategies
---
This repository contains the code for the creation and analysis of the model presented in the manuscript *Reconstruction of a large-scale outbreak of SARS-CoV-2 infection in Iceland informs vaccination strategies*. The model is created with the script `run_model.R`, which accepts the following command line arguments:
```
[--nr_iter N]   (int) the number of MCMC iterations to be run.
[--chain N]     (int) the chain number, if running more than one chain.
[--head N]      (int) only the first N data points will be modeled.
[--mu_init X]   (float) the initial value of mu, the mutation rate.
[--move_mu]     indicates whether mu should be a parameter or fixed.
[--normalize]   indicates whether truncated time distributions should be normalized.
[--imputation]  indicates whether low-coverage and missing haplotypes should be imputed.
``` 
The resulting model will be stored in `results/full_model`. The vaccination strategy simulations can be run via `run_vaccination_simulation.R`, which accepts the following command line arguments:
```
--model_name            (str) name of the outbreaker model.
--sim_nr                (int) the number of the batch of simulations (for parallelising).
-N                      (int) the number of simulations to perform.
--first_dose_efficacy   (float) efficacy of first dose to use in the simulation.
--second_dose_efficacy  (float) efficacy of second dose to use in the simulation.
--out_dir               (str) path to output directory.
[--init_size X]         (int) size of initial outbreak.
[--step_size X]         (float) percentage interval at which to simulate.
[--replay]              indicates whether to replay de facto vaccinations or start from zero.
[--subtrees]            indicates whether we should run the simulation on subtrees or the whole tree.
```

### Data spec
---
The dataset should be located in `data/dat.tsv` and contains the following fields:
```
Wave (int) denotes the wave of infections
t_id (str) unique identifier per person
sam_id (str) unique identifier per sample
seq (str) denotes whether sequenced via ONT or Illumina
mutations_full (str) semicolon-delimited list of all mutations, LOCUS-REF-ALT
mutations (str) semicolon-delimited list of all mutations, LOCUS
REF (str) semicolon-delimited list of all mutations, REF
ALT (str) semicolon-delimited list of all mutations, ALT
clade (str) denotes the clade
mutations_full_extra (str) semicolon-delimited list of mutations on top of the clade, LOCUS-REF-ALT
mutations_extra (str) semicolon-delimited list of mutations on top of the clade, LOCUS
REF_extra (str) semicolon-delimited list of mutations on top of the clade, REF
ALT_extra (str) semicolon-delimited list of mutations on top of the clade, ALT
ct (float) the ct-value of the sample
score (float) the ratio of called bases in the sequence
connected_to_revised (str) comma-delimited list of the t_ids of reported contacts
TYPE (str) the type of infection, one of ["border", "family", "healthcare", "none", "school", "social", "travel", "work"]
age_at_diag (int) age at diagnosis
birth_year (int) birth year
child (bool) denotes whether person is a child or not
date (str) diagnosis date, YYYY-MM-DD
date_isolation (str) date when person isolated, YYYY-MM-DD
date_exposure (str) date of reported exposure, YYYY-MM-DD
date_symptoms (str) date of onset of symptoms, YYYY-MM-DD
date_q (str) date at the start of quarantine, YYYY-MM-DD
date_registered (str) date when the person was first registered in the contact tracing system, YYYY-MM-DD
date_ref (str) min(date, date_symptoms)
shows_symptoms (bool) denotes whether person was symptomatic at diagnosis
q_type_binary (str) denotes whether person was in quarantine at diagnosis, one of ['outside_q', 'in_q']
q_type (str) denotes how long the person was in quarantine, one of ['outside_q', 'short_q', 'long_q']
household_id (int) which curated household the person belongs to 
```

The files `data/depth_miseq_dat.tsv` and `data/depth_ont_dat.tsv` should contain the loci per sample with a sequencing depth of 30x or lower.
```
seq (str) one of ["MISEQ", "ONT"]
sam_id (str) the unique sample id
Pos (int) locus
depth (int) sequencing depth
```

The file `data/vaccination_data.xlsx` consists of a header containing the columns `age`, `dose`, and one column for each date of vaccinations. Rows 1-9 contain the number of people who received the first dose on the corresponding date in the age groups "0-15", "16-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", and "90+". Rows 11-18 contain the number of people in the same age groups who received the second dose on the corresponding date. Rows 19-27 contain the number of people in the same age groups who received a single-dose vaccine on the corresponding date.

The file `data/age_groups.tsv` contains a header `age_group, second, first, none` and denotes how many people in the age groups detailed above had received two doses, one dose, and zero doses of a vaccine as of April 28, 2021.
