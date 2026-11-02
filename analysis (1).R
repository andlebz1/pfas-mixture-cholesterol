# =============================================================
# PFAS mixtures and total cholesterol in US adults
# NHANES 2017-2018 | single-chemical models + quantile g-computation
# Author: Andleeb Zahra
# =============================================================

# ---- Step 1: Packages ----
# install.packages(c("tidyverse", "haven", "broom", "qgcomp"))  # run once
library(tidyverse)
library(haven)
library(broom)
library(qgcomp)

# ---- Step 2: Data ----
# Download these four files from the NHANES website into a folder called "data":
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PFAS_J.xpt
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/TCHOL_J.xpt
#   https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.xpt
pfas  <- read_xpt("data/PFAS_J.xpt")   # serum PFAS (n = 2,133)
demo  <- read_xpt("data/DEMO_J.xpt")   # demographics (n = 9,254)
tchol <- read_xpt("data/TCHOL_J.xpt")  # total cholesterol (n = 7,435)
bmx   <- read_xpt("data/BMX_J.xpt")    # body measures (n = 8,704)

# ---- Step 3: Clean and merge ----
# PFOA and PFOS are measured as linear and branched forms; we sum them.
# Me-PFOSA-AcOH and PFUnDA are often below the detection limit and are not used.
dat <- pfas %>%
  mutate(
    PFOA  = LBXNFOA + LBXBFOA,
    PFOS  = LBXNFOS + LBXMFOS,
    PFHxS = LBXPFHS,
    PFNA  = LBXPFNA,
    PFDeA = LBXPFDE
  ) %>%
  select(SEQN, PFOA, PFOS, PFHxS, PFNA, PFDeA) %>%
  inner_join(select(demo, SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR), by = "SEQN") %>%
  inner_join(select(tchol, SEQN, LBXTC), by = "SEQN") %>%
  inner_join(select(bmx, SEQN, BMXBMI), by = "SEQN") %>%
  rename(age = RIDAGEYR, sex = RIAGENDR, ethnicity = RIDRETH3,
         income_ratio = INDFMPIR, chol = LBXTC, bmi = BMXBMI) %>%
  filter(age >= 20) %>%
  drop_na() %>%
  mutate(
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    ethnicity = factor(ethnicity)
  )

nrow(dat)      # 1,384 adults
summary(dat)

pfas_vars <- c("PFOA", "PFOS", "PFHxS", "PFNA", "PFDeA")
round(cor(dat[pfas_vars]), 2)   # PFAS are moderately correlated (r = 0.16-0.69)

# ---- Step 4: Single-chemical models ----
# Linear regression of cholesterol on log(PFAS), adjusted for confounders,
# with Benjamini-Hochberg FDR correction across the five PFAS.
single <- map_dfr(pfas_vars, function(v) {
  f <- as.formula(paste0("chol ~ log(", v, ") + age + sex + ethnicity + bmi + income_ratio"))
  tidy(lm(f, data = dat), conf.int = TRUE) %>%
    filter(str_detect(term, v)) %>%
    mutate(chemical = v)
})
single$p_fdr <- p.adjust(single$p.value, method = "BH")
single %>% select(chemical, estimate, conf.low, conf.high, p.value, p_fdr)

# ---- Step 5: Mixture model (quantile g-computation) ----
fit <- qgcomp.noboot(
  chol ~ PFOA + PFOS + PFHxS + PFNA + PFDeA + age + sex + ethnicity + bmi + income_ratio,
  expnms = pfas_vars,
  data   = dat,
  family = gaussian(),
  q      = 4
)
fit   # psi1 = change in cholesterol when all PFAS increase by one quartile together

# ---- Step 6: Figures ----
dir.create("figures", showWarnings = FALSE)

# Figure 1: single-PFAS estimates, rescaled to a doubling of PFAS
single_plot <- single %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x * log(2)))

p1 <- ggplot(single_plot, aes(x = estimate, y = reorder(chemical, estimate))) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Change in total cholesterol (mg/dL) per doubling of PFAS",
       y = NULL,
       title = "Each PFAS and total cholesterol",
       subtitle = "NHANES 2017-2018 adults (n = 1,384)\nAdjusted for age, sex, ethnicity, income and BMI") +
  theme_minimal(base_size = 13)
ggsave("figures/single_pfas.png", p1, width = 8, height = 4.5, dpi = 150)

# Figure 2: qgcomp weights
png("figures/qgcomp_weights.png", width = 1000, height = 600)
plot(fit)
dev.off()
