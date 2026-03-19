library(tidyverse)
library(fixest)
library(ggplot2)
library(rddensity)
library(vtable)
library(rdrobust)

# read file from raw_data folder
load("./raw_data/DWI_Data.rdata")

# Basic checks 
str(dwi)
names(dwi)

# Summary 
summary(dwi)

# Missing values - Zero missing values 
colMeans(is.na(dwi)) * 100

# Correlation between the 2, should be closer to 1 - to check for measurement discrepancies
cor(dwi$bac1, dwi$bac2, use = "complete.obs")

# Feature engineering
c <- 0.08  # cutoff (legal limit)

dwi <- dwi %>%
  mutate(
    min_BAC = pmin(bac1, bac2, na.rm = TRUE), #(legal rule in WA to consider minimum BAC of both tests)
    centered_BAC = min_BAC - c,  # centered running variable to 0
    Treated = as.integer(min_BAC >= c)
  )

# EDA PLOTS 
# Running var: centered_BAC (min_BAC - 0.08)
# Cutoff: 0 (Running Variable is Centered)
# Treatment indicator: Treated (1 if min_BAC >= 0.08 else 0)
# Outcome: recidivism

# Create BAC bins (binned means plots)
# Choose bin width. 0.001 BAC
# The bin width is chosen for visualization purposes
# to ensure the smoothness of the outcome and covariates around the cutoff
bin_w <- 0.001

dwi_binned <- dwi %>%
  mutate(
    bac_bin = floor(centered_BAC / bin_w) * bin_w,
    bac_bin_mid = bac_bin + bin_w/2
  ) %>%
  group_by(bac_bin, bac_bin_mid) %>%
  summarize(
    n = n(),
    mean_recid = mean(recidivism, na.rm = TRUE),
    mean_treat = mean(Treated, na.rm = TRUE),
    mean_male  = mean(male, na.rm = TRUE),
    mean_white = mean(white, na.rm = TRUE),
    mean_age   = mean(aged, na.rm = TRUE),
    mean_acc   = mean(acc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

# Plot Outcome vs Running variable (binned means) with quadratic and linear fit
# to identify the jump and check for curvature for best fit functional form
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_recid)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Linear fit 
  geom_smooth(
    aes(group = bac_bin_mid < 0),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 1
  ) +
  
  # Quadratic fit
  geom_smooth(
    aes(group = bac_bin_mid < 0),
    method = "lm",
    formula = y ~ poly(x, 2, raw = TRUE),
    se = FALSE,
    color = "red",
    linetype = "dotted",
    linewidth = 1
  ) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean recidivism (within BAC bin)",
    title = "RDD Binned Means: Linear (Black) vs Quadratic (Red Dotted)"
  ) +
  theme_minimal()


# Plot Treatment probability vs running variable
# Checks whether treatment assignment jumps at the cutoff (sharp/fuzzy RD)
# to check if the treatment is based perfectly on being below or above cutoff
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_treat)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(aes(group = bac_bin_mid < 0), method = "lm", se = FALSE) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean Treated (within BAC bin)",
    title = "Probability of Harsher Punishment vs Centered BAC"
  ) +
  theme_minimal()

# Plot placebo test graphs (covariate smoothness around cutoff) - RDD assumption check
# If these jump at cutoff, RD assumptions are not satisfied

# Male
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_male)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(aes(group = bac_bin_mid < 0), method = "lm", se = FALSE) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean male (within BAC bin)",
    title = "Placebo Check: Male vs Centered BAC"
  ) +
  theme_minimal()

# White
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_white)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(aes(group = bac_bin_mid < 0), method = "lm", se = FALSE) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean white (within BAC bin)",
    title = "Placebo Check: White vs Centered BAC"
  ) +
  theme_minimal()

# Age
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_age)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(aes(group = bac_bin_mid < 0), method = "lm", se = FALSE) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean age (within BAC bin)",
    title = "Placebo Check: Age vs Centered BAC"
  ) +
  theme_minimal()

# Accident
ggplot(dwi_binned, aes(x = bac_bin_mid, y = mean_acc)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(aes(group = bac_bin_mid < 0), method = "lm", se = FALSE) +
  labs(
    x = "Centered BAC (min_BAC - 0.08)",
    y = "Mean accident (within BAC bin)",
    title = "Placebo Check: Accident vs Centered BAC"
  ) +
  theme_minimal()

# Manipulation/Bunching check: Density histogram + formal rddensity test 
h  <- 0.05
c0 <- 0

dwi_local <- dwi %>%
  filter(abs(centered_BAC) < h)

# Histogram (within bandwidth)
ggplot(dwi_local, aes(x = centered_BAC)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = c0, linetype = "dashed") +
  labs(
    x = "Centered BAC (BAC - 0.08)",
    y = "Count",
    title = paste0("Density Check (Histogram within ±", h, "):Centered BAC Around Cutoff")
  ) +
  theme_minimal()

# Formal density discontinuity test
dens_test <- dwi_local %>%
  pull(centered_BAC) %>%
  rddensity(c = c0)

summary(dens_test)

# Define the RD analysis sample and weights:
# - X is the centered running variable (BAC - 0.08), so the cutoff is at X = 0
# - Treated = 1 if BAC is at/above the legal limit (harsher punishment)
# - Restrict to a local bandwidth |X| < 0.05 to compare drivers near the cutoff
# - Apply triangular kernel weights so observations closest to the cutoff get more weight
# h = bandwidth
h <- 0.05

dwi_rd <- dwi %>%
  mutate(
    X = centered_BAC,
    Treated = as.integer(X >= 0),
    w_tri = pmax(0, 1 - abs(X)/h)
  ) %>%
  filter(abs(X) < h) 

# PLACEBO (BALANCE) TEST MODELS
# Goal: In a valid RD, people just below and just above the cutoff should look similar
# in "predetermined" characteristics (things the cutoff cannot cause).
# So we run the SAME RD regression, but we replace the outcome (recidivism) with covariates
# like male/age/race/accident/year. These should NOT jump at the cutoff.
# If they DO jump, it suggests sorting/manipulation or non-comparability around the cutoff,
# which threatens RD validity

male_placebo <- feols(
  male ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)

aged_placebo <- feols(
  aged ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)

white_placebo <- feols(
  white ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)

acc_placebo <- feols(
  acc ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)

year_placebo <- feols(
  year ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)

etable(male_placebo, aged_placebo, white_placebo, acc_placebo, year_placebo)

#Base model- No controls
m1 <- feols(
  recidivism ~ Treated + X + Treated:X,
  data = dwi_rd,
  vcov = "hetero"
)
etable(m1)

#Bandwidth Selection Models
m2 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .07), vcov = "hetero")
m3 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .06), vcov = "hetero")
m4 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .04), vcov = "hetero")
m5 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .03), vcov = "hetero")
m6 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .025), vcov = "hetero")
m7 <- feols(recidivism ~ Treated*X, data = dwi_rd %>% filter(abs(X) < .02), vcov = "hetero")

etable(m1,m2,m3,m4,m5,m6,m7, keep = "Treated")

# BANDWIDTH SELECTION / ROBUSTNESS MODELS
# Goal: RD estimates can change depending on the bandwidth (how close to the cutoff we look).
# We re-run the MAIN RD model multiple times using different windows around the cutoff.
# If the estimated treatment effect (Treated at cutoff) stays similar across bandwidths,
# results are more credible (not driven by arbitrary bandwidth choice).
mods <- list(
  `0.07`  = m2,
  `0.06`  = m3,
  `0.05`  = m1,
  `0.04`  = m4,
  `0.03`  = m5,
  `0.025` = m6,
  `0.02`  = m7
)

# Extract the cutoff jump
# Goal of this chunk:
#   1) Build a tidy table (sens) with, for each bandwidth h:
#        - h (the bandwidth)
#        - n (number of observations used by that model)
#        - jump (estimated treatment effect at the cutoff = coef on Treated)
#        - se (standard error of that estimate)
#        - ci_low / ci_high (95% confidence interval bounds)
#   2) Plot the estimate vs bandwidth to visually assess robustness/sensitivity.
sens <- bind_rows(lapply(names(mods), function(h) {
  m <- mods[[h]]
  data.frame(
    h = as.numeric(h),
    n = nobs(m),
    jump = coef(m)["Treated"],
    se = se(m)["Treated"]
  )
})) %>%
  mutate(
    ci_low  = jump - 1.96 * se,
    ci_high = jump + 1.96 * se
  ) %>%
  arrange(h)

# SENSITIVITY PLOT: RD estimate (jump) vs bandwidth-
# Goal: Visual robustness check.
# If the RD effect is credible, the "jump" should be reasonably stable across a range of h's.
# Typical pattern:
#   - As h decreases (more local): bias decreases but variance increases (wider CI)
#   - As h increases (less local): variance decreases but bias risk increases
# This plot shows that tradeoff and whether your conclusion depends on bandwidth choice
sens
ggplot(sens, aes(x = h, y = jump)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0) +
  theme_minimal() +
  labs(
    x = "Bandwidth (|X| ≤ h)  — smaller h = more local",
    y = "Estimated jump in recidivism at cutoff (coef on Treated)",
    title = "Sensitivity of RD estimate to bandwidth choice",
    subtitle = "Points are estimates; vertical bars are 95% confidence intervals"
  )

#Model with controls and chosen bandwidth 0.05 (robustness check) 
m8 <- feols(
  recidivism ~ Treated + X + Treated:X + white + male + aged | year,
  data = dwi_rd,
  vcov = "hetero"
)

etable(m1, m8)

#Model Linear with Triangular Kernel (Main Specification Model)
m9 <- feols(recidivism ~ Treated + X + Treated:X,
            data = dwi_rd,
            weights = ~ w_tri,
            vcov = "hetero")

etable(m1, m9)

# Model Quadratic RD with Triangular Kernel (robustness check)
m10 <- feols(
  recidivism ~ Treated + X + I(X^2) + Treated:X + Treated:I(X^2),
  data = dwi_rd,
  weights = ~ w_tri,     # <-- apply triangular weights
  vcov = "hetero"
)

# Compare linear vs quadratic
etable(
  m9, m10,
  headers = c("Linear", "Quadratic"),
  se.below = TRUE
)

#Joint wald test for quadratic terms
wald(m10, keep = "X\\^2")

# RDROBUST 
y <- dwi$recidivism
x <- dwi$centered_BAC

# Model RDROBUST default bandwidth (Robustness Check) 
m11 <- rdrobust(y = y, x = x, c = 0)
summary(m11)

# Model RDROBUST Manual bandwidth 0.05 (Robustness check)
m12 <- rdrobust(y = dwi$recidivism, x = dwi$centered_BAC, c = 0, h = 0.05)
summary(m12)

#Plot rdrobust default bandwidth
rdplot(dwi$recidivism, dwi$centered_BAC, c = 0)

#Plot rdrobust manual bandwidth 0.05
rdplot(dwi$recidivism, dwi$centered_BAC, c = 0, h=0.05)


