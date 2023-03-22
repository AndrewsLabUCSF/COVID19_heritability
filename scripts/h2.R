library(tidyverse)
library(glue)
library(stringr)

# Define column types for summary statistics
coltypes = cols(
  ID = col_character(),
  CHROM = col_double(),
  POS = col_double(),
  REF = col_character(),
  ALT = col_character(),
  AF = col_double(),
  TRAIT = col_character(),
  BETA = col_double(),
  SE = col_double(),
  Z = col_double(),
  P = col_double(),
  N = col_double(),
  OR = col_double(),
  OR_L95 = col_double(),
  OR_U95 = col_double(),
  DIR = col_character(),
  G1000_ID = col_character(),
  G1000_VARIANT = col_character(),
  DBSNP_ID = col_character(),
  DBSNP_VARIANT = col_character(),
  OLD_ID = col_character(),
  OLD_VARIANT = col_character()
)

# LD File paths
ld_path = "resources/GenomicSEM/eur_w_ld_chr/"
hm3_path = "resources/GenomicSEM/w_hm3.snplist"

# Sample Sizes
## Longcovid 
longcovid_samplesize.raw <- read_csv("resources/LongCovidHGI_DF4_N.csv") %>%
  magrittr::set_colnames(c("cohort", "cases_n1", "ctrls_n1", "cases_n2", "ctrls_n2", "cases_w1", "ctrls_w1", "cases_w2", "ctrls_w2"))

longcovid_samplesize <- bind_rows(
  select(longcovid_samplesize.raw, cohort, cases_n1, ctrls_n1) %>%
    rename(case = cases_n1, ctrl = ctrls_n1) %>%
    mutate(trait = 'LongCovid_N1'), 
  select(longcovid_samplesize.raw, cohort, cases_n2, ctrls_n2) %>%
    rename(case = cases_n2, ctrl = ctrls_n2) %>%
    mutate(trait = 'LongCovid_N2'), 
  select(longcovid_samplesize.raw, cohort, cases_w1, ctrls_w1) %>%
    rename(case = cases_w1, ctrl = ctrls_w1) %>%
    mutate(trait = 'LongCovid_W1'), 
  select(longcovid_samplesize.raw, cohort, cases_w2, ctrls_w2) %>%
    rename(case = cases_w2, ctrl = ctrls_w2) %>%
    mutate(trait = 'LongCovid_W2'), 
) %>%
  relocate(trait, .before = cohort) %>%
  rowwise() %>%
  mutate(
    n = sum(case, ctrl, na.rm = T), 
    pi = case / n,
    EffN = 4*pi*(1-pi)*n
  ) %>%
  ungroup() %>%
  group_by(trait) %>%
  summarise(
    n = sum(n), 
    case = sum(case, na.rm = T),
    ctrl = sum(ctrl, na.rm = T),
    EffN = sum(EffN, na.rm = T), 
    preval = case / n
  ) %>%
  mutate(code = c("longcovid2022N1v4", "longcovid2022N2v4", "longcovid2022W1v4", "longcovid2022W2v4")) %>%
  arrange(code)

## COVID-HGI 
samplesize.raw <- readxl::read_xlsx("resources/283874_file06.xlsx", sheet = 1, skip = 1, na = c("", "-")) %>%
  janitor::clean_names() 

covid_samplesize <- bind_rows(
    select(samplesize.raw, site_abbreviation, ancestry, cases_critical_illness_due_to_covid_19, controls_critical_illness_due_to_covid_19) %>%
      rename(case = cases_critical_illness_due_to_covid_19, ctrl = controls_critical_illness_due_to_covid_19) %>%
      mutate(trait = 'critical_illness'), 
    select(samplesize.raw, site_abbreviation, ancestry, cases_hospitalized_due_to_covid_19, controls_hospitalized_due_to_covid_19) %>%
      rename(case = cases_hospitalized_due_to_covid_19, ctrl = controls_hospitalized_due_to_covid_19) %>%
      mutate(trait = 'hospitalized'), 
    select(samplesize.raw, site_abbreviation, ancestry, cases_reported_infection, controls_reported_infection) %>%
      rename(case = cases_reported_infection, ctrl = controls_reported_infection) %>%
      mutate(trait = 'reported_infection')
  ) %>%
  relocate(trait, .before = site_abbreviation) %>%
  filter(ancestry == 'EUR') %>%
  rowwise() %>%
  mutate(
    n = sum(case, ctrl, na.rm = T), 
    pi = case / n,
    EffN = 4*pi*(1-pi)*n
  ) %>%
  ungroup() %>%
  group_by(trait) %>%
  summarise(
    n = sum(n), 
    case = sum(case, na.rm = T),
    ctrl = sum(ctrl, na.rm = T),
    EffN = sum(EffN, na.rm = T), 
    preval = case / n
  ) %>%
  mutate(code = c("COVID19_HGI_C2", "COVID19_HGI_B2", "COVID19_HGI_A2")) %>%
  arrange(code)

# Munge Summary Statistics
## COVID HGI
a2_path = "resources/COVID19_HGI_A2_ALL_eur_20220602.txt.GRCh37.vcf.gz.tsv.gz"
b2_path = "resources/COVID19_HGI_B2_ALL_eur_20220602.txt.GRCh37.vcf.gz.tsv.gz"
c2_path = "resources/COVID19_HGI_C2_ALL_eur_20220602.txt.GRCh37.vcf.gz.tsv.gz"

covid_trait_names <- c(c("COVID19_HGI_A2", "COVID19_HGI_B2", "COVID19_HGI_C2"))

ss <- map(c(a2_path, b2_path, c2_path), read_tsv,
          col_select = c(rsid, '#CHR', POS, REF, ALT, all_meta_AF, all_inv_var_meta_beta, all_inv_var_meta_sebeta, 
                         all_inv_var_meta_p, all_inv_var_meta_cases, all_inv_var_meta_controls, all_inv_var_meta_effective
            )
          )

ss_out <- map2(ss, covid_trait_names, function(x, y){
    x %>%
      rename(CHR = '#CHR', AF = all_meta_AF, beta = all_inv_var_meta_beta, se = all_inv_var_meta_sebeta, 
             p = all_inv_var_meta_p, cases = all_inv_var_meta_cases, ctrls = all_inv_var_meta_controls, 
             effn = all_inv_var_meta_effective) %>%
      mutate(trait = y, 
             total_n = cases + ctrls, 
             maf = ifelse(AF > 0.5, 1-AF, AF), 
             ) %>%
      filter(nchar(REF) == 1, nchar(ALT) == 1) 
    }
  )

map2(ss_out, covid_trait_names, function(x,y){
  write_tsv(x, glue("data/{y}_ldsc.tsv.gz"))
})

## LongCovid 
n1_path = "resources/longcovid2022N1v4.chrall.CPRA_b37.tsv.gz"
n2_path = "resources/longcovid2022N2v4.chrall.CPRA_b37.tsv.gz"
w1_path = "resources/longcovid2022W1v4.chrall.CPRA_b37.tsv.gz"
w2_path = "resources/longcovid2022W2v4.chrall.CPRA_b37.tsv.gz"

longcovid_trait_names <- c(c("longcovid2022N1v4", "longcovid2022N2v4", "longcovid2022W1v4", "longcovid2022W2v4"))


longcovid_ss <- map(c(n1_path, n2_path, w1_path, w2_path), read_tsv,
                    comment = "##", col_types = coltypes, 
                    col_select = c(DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, SE, Z, P, N, TRAIT)
                    )

longcovid_ss_out <- map(longcovid_ss, function(x){
  x %>%
    mutate(maf = ifelse(AF > 0.5, 1-AF, AF), 
    ) %>%
    rename(
      rsid = DBSNP_ID, CHR = CHROM, beta = BETA, se = SE, p = P, total_n = N,
    )
})

map2(longcovid_ss_out, longcovid_trait_names, function(x,y){
  write_tsv(x, glue("data/{y}_ldsc.tsv.gz"))
})

# LDSC Munge 
trait_names = c(covid_trait_names, longcovid_trait_names)
sample_size = bind_rows(covid_samplesize, longcovid_samplesize)

## Total Sample Size 
GenomicSEM::munge(
  files = glue("data/{trait_names}_ldsc.tsv.gz"), 
  hm3 = hm3_path, 
  trait.names = glue("{trait_names}_N"), 
  maf.filter = 0.01, 
  column.names = list(
    SNP='rsid', 
    MAF='maf', 
    A1='ALT',
    A2='REF', 
    effect='beta', 
    N='total_n'
  ), 
  overwrite=TRUE
)

## Effective Sample Size 
GenomicSEM::munge(
  files = glue("data/{trait_names}_ldsc.tsv.gz"), 
  hm3 = hm3_path, 
  trait.names = glue("{trait_names}_Neff"), 
  maf.filter = 0.01, 
  N = pull(sample_size, EffN),
  column.names = list(
    SNP='rsid', 
    MAF='maf', 
    A1='ALT',
    A2='REF', 
    effect='beta'
    # N='cohort_neff'
  ), 
  overwrite=TRUE
)

ldsc_output <- list.files()
ldsc_files <- ldsc_output[str_detect(ldsc_output, "gz|log")]

file.copy(from = ldsc_files, to = glue("data/{ldsc_files}"))
file.remove(ldsc_files)


# Heritability
h2_strings <- "Liability scale results for|Total Liability Scale h2|Heritability Results for trait|Mean Chi^2 across remaining SNPs|Lambda GC|Intercept|Ratio|Total Observed Scale h2|h2 Z"

## Total Sample Size
n_ldsc_res <- GenomicSEM::ldsc(traits = glue("data/{trait_names}_N.sumstats.gz"),
                               sample.prev = pull(sample_size, preval),
                               population.prev = pull(sample_size, preval),
                               ld = ld_path, 
                               wld = ld_path,
                               trait.names = trait_names, 
                               ldsc.log = "data/covid_hgi_prev_observed_N"
)


log_n <- read_lines("data/covid_hgi_prev_observed_N_ldsc.log")

res_n_out <- log_n %>%
  as_tibble() %>%
  filter(str_detect(value, h2_strings)) %>%
  filter(!str_detect(value, "Cross trait Intercept")) %>%
  separate(value, c('term', 'estimate'), sep = ": ") %>%
  mutate(trait = case_when(
    term == "Heritability Results for trait" ~ estimate, 
    term == "Liability scale results for" ~ estimate, 
    TRUE ~ NA_character_
  )) %>%
  fill(trait, .direction = 'down') %>%
  filter(estimate != trait) %>%
  mutate(
    trait = str_replace_all(trait, "data/|_N.sumstats.gz", "")
  ) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  left_join(sample_size, by = c("trait" = 'code')) %>%
  select(-trait.y, -case, -ctrl) %>%
  relocate(n, EffN, preval, .after = trait) 

write_csv(res_n_out, "results/h2_n.csv")


## Effective Sample Size 
### Total Sample Prevalence (Observed)
neff_ldsc_res <- GenomicSEM::ldsc(traits = glue("data/{trait_names}_Neff.sumstats.gz"),
                             sample.prev = rep(0.5, nrow(sample_size)),
                             population.prev = pull(sample_size, preval),
                             ld = ld_path, 
                             wld = ld_path,
                             trait.names = trait_names, 
                             ldsc.log = "data/covid_hgi_prev_observed_Neff"
)

### Range of population prevalences
pop_prev = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
map(pop_prev, function(pop_prev){
    preval = rep(pop_prev, length(trait_names))
    GenomicSEM::ldsc(traits = glue("data/{trait_names}_Neff.sumstats.gz"),
                   sample.prev = rep(0.5, length(trait_names)),
                   population.prev = preval,
                   ld = ld_path, 
                   wld = ld_path,
                   trait.names = trait_names, 
                   ldsc.log = glue("data/covid_hgi_prev_{pop_prev}_Neff")
  )
})



### Munge Results 
res_out <- map_dfr(c("observed", pop_prev), function(pop_prev){
  
  res <- read_lines(glue("data/covid_hgi_prev_{pop_prev}_Neff_ldsc.log"))
  
  res_out <- res %>%
    as_tibble() %>%
    filter(str_detect(value, h2_strings)) %>%
    filter(!str_detect(value, "Cross trait Intercept")) %>%
    separate(value, c('term', 'estimate'), sep = ": ") %>%
    mutate(trait = case_when(
      term == "Heritability Results for trait" ~ estimate, 
      term == "Liability scale results for" ~ estimate, 
      TRUE ~ NA_character_
    )) %>%
    fill(trait, .direction = 'down') %>%
    filter(estimate != trait) %>%
    mutate(
      trait = str_replace_all(trait, "data/|_Neff.sumstats.gz", "")
    ) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    left_join(sample_size, by = c("trait" = 'code')) %>%
    select(-trait.y, -case, -ctrl) %>%
    relocate(n, EffN, preval, .after = trait) %>%
    mutate(pop_prev = pop_prev)
  
  res_out
})

write_csv(res_out, "results/h2_neff.csv")




















