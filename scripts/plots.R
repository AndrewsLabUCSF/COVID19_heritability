library(tidyverse)

## Import data 
res_out <- read_csv("results/h2_neff.csv") 

## Munge data for plotting 
h2_dat_p <- res_out %>%
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
  )

## COVID-HGI
ggplot(h2_dat_p %>% filter(str_detect(trait, "COVID19_HGI")), aes(x = pop_prev, y = l_h2, color = observed)) + 
  facet_wrap(vars(label)) + 
  geom_point(size = 0.75) + 
  geom_linerange(aes(ymin = l_h2_lci, ymax = l_h2_uci), linewidth = 0.5) + 
  scale_x_continuous(
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
    labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) + 
  scale_color_manual(values = c("black", "red")) + 
  labs(
    x = "Population Prevelance", 
    y = expression(paste('Liability Scale h'^2)), 
    color = "Observed Sample Prevalence") + 
  theme_bw() + 
  theme(
    legend.position = "bottom", 
    strip.background = element_blank(), 
    text = element_text(size = 8), 
    panel.grid.minor = element_blank()
  )

ggsave('results/plots/covid_h2.png', height = 2, width = 6, units = "in")
ggsave('results/plots/covid_h2.pdf', height = 2, width = 6, units = "in")

## Longcovid
ggplot(h2_dat_p %>% filter(str_detect(trait, "longcovid")), aes(x = pop_prev, y = l_h2, color = observed)) + 
  facet_wrap(vars(label)) + 
  geom_point(size = 0.75) + 
  geom_linerange(aes(ymin = l_h2_lci, ymax = l_h2_uci), linewidth = 0.5) + 
  scale_x_continuous(
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
    labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) + 
  scale_color_manual(values = c("black", "red")) + 
  labs(x = "Population Prevelance", 
       y = expression(paste('Liability Scale h'^2)), 
       color = "Observed Sample Prevalence") + 
  theme_bw()  + 
  theme(
    legend.position = "bottom", 
    strip.background = element_blank(), 
    text = element_text(size = 8), 
    panel.grid.minor = element_blank()
  )


ggsave('results/plots/longcovid_h2.png', height = 4, width = 4, units = "in")



