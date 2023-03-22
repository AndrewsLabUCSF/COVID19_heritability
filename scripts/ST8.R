library(tidyverse)

res_out <- read_csv("results/h2_neff.csv") 

res_out %>%
  filter(str_detect(trait, "COVID19_HGI")) %>%
  mutate(
    pop_prev = as.numeric(pop_prev), 
    pop_prev = ifelse(is.na(pop_prev), preval, pop_prev), 
    label = case_when(
      trait == "COVID19_HGI_A2" ~ "Reported Infection vs population controls", 
      trait == "COVID19_HGI_B2" ~ "Hospitalization vs population controls", 
      trait == "COVID19_HGI_C2" ~ "Critical Illness vs population controls", 
    ), 
  ) %>%
  arrange(trait, pop_prev) %>%
  select(Analysis = label, N = n, EffN, Prevalance = pop_prev, 'Lambda GC', 'LDSC intercept (SE)' = Intercept, 
         'Ratio (SE)' = Ratio, 'Total Observed Scale h2 (SE)' = 'Total Observed Scale h2', 
         'h2 Z', 'Total Liability Scale h2 (SE)' = 'Total Liability Scale h2') %>%
  write_csv(., "results/ST8.csv")

## Ranges
res_out %>%
  janitor::clean_names() %>%
  separate(total_liability_scale_h2, c('l_h2', "l_h2_se"), sep = "\\(") %>%
  mutate(
    observed = ifelse(pop_prev == "observed", TRUE, FALSE),
    pop_prev = as.numeric(pop_prev), 
    pop_prev = ifelse(is.na(pop_prev), preval, pop_prev), 
    l_h2 = str_trim(l_h2), 
    l_h2_se = str_replace(l_h2_se, "\\)", ""), 
    l_h2 = as.numeric(l_h2), 
    l_h2_se = as.numeric(l_h2_se), 
    l_h2_lci = l_h2 - (1.96 * l_h2_se),
    l_h2_uci = l_h2 + (1.96 * l_h2_se)
  ) %>%
  mutate(
    label = case_when(
      trait == "COVID19_HGI_A2" ~ "Reported Infection", 
      trait == "COVID19_HGI_B2" ~ "Hospitalization", 
      trait == "COVID19_HGI_C2" ~ "Critical Illness", 
      trait == "longcovid2022N1v4" ~ "Test Verified, COVID+",
      trait == "longcovid2022N2v4" ~ "Test Verified, Population",
      trait == "longcovid2022W1v4" ~ "Any COVID, COVID+",
      trait == "longcovid2022W2v4" ~ "Any COVID, Population",
    ), 
    label = fct_relevel(label, "Reported Infection", "Hospitalization", "Critical Illness")
  ) %>%
  filter(pop_prev %in% c(0.01, 0.9)) %>%
  select(label, pop_prev, l_h2) %>%
  mutate(l_h2 = signif(l_h2, digits = 2)) %>%
  pivot_wider(names_from = pop_prev, values_from = l_h2) %>%
  unite(l_h2, c(`0.01`, `0.9`), sep = "–")

