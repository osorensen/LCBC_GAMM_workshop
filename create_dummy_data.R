library(tidyverse)
library(gamm4)

load(file.path(Sys.getenv("MOAS_PATH"), "MOAS.RData"))

dat <- MOAS %>% 
  select(CrossProject_ID, Age, MRI_aseg_left_hippocampus,
         MRI_aseg_right_hippocampus, 
         MRI_aseg_estimated_total_intra_cranial_vol, Sex) %>% 
  rename(
    id = CrossProject_ID,
    age = Age,
    sex = Sex,
    icv = MRI_aseg_estimated_total_intra_cranial_vol  
  ) %>% 
  mutate(
    hippocampus = MRI_aseg_left_hippocampus + 
      MRI_aseg_right_hippocampus
  ) %>% 
  select(-starts_with("MRI")) %>% 
  na.omit() %>% 
  mutate(
    across(c(icv, age), list(z = ~ (. - mean(.)) / sd(.)))
  )

mod <- gamm(hippocampus ~ s(age_z, bs = "cr") + sex + icv_z,
            data = dat, random = list(id =~ 1))

varcomps <- VarCorr(mod$lme)[c("(Intercept)", "Residual"), "StdDev"] %>% 
  as.numeric()

icv_quantiles <- quantile(dat$icv, probs = seq(from = 0, to = 1, by = .001), names = FALSE)
age_quantiles <- dat %>% 
  group_by(id) %>% 
  filter(age == min(age)) %>% 
  distinct(age) %>% 
  pull(age) %>% 
  quantile(probs = seq(from = 0, to = 1, by = .001), names = FALSE)

# Create new data
n <- 2000
simdat <- tibble(
  id = seq_len(n),
  age = age_quantiles[sample(seq_along(age_quantiles), n, replace = TRUE)],
  sex = sample(c("Female", "Male"), size = n, replace = TRUE),
  time = map(
    sample(1:4, size = n, replace = TRUE, 
           prob = c(.5, .3, .15, .05)), 
    ~ cumsum(c(0, runif(.x - 1, min = .5, max = 5)))),
  random_intercept = rnorm(n, sd = varcomps[[1]]),
  icv = icv_quantiles[sample(seq_along(icv_quantiles), n, replace = TRUE)],
  icv_z = (icv - mean(dat$icv)) / sd(dat$icv),
  genotype = sample(c("A", "B", "C"), size = n, 
                    replace = TRUE, prob = c(.5, .4, .1)),
  education = runif(n, min = 10, max = 25)
) %>% 
  unnest(cols = time) %>% 
  mutate(
    age = age + time,
    age_z = (age - mean(dat$age)) / sd(dat$age),
    residual = rnorm(nrow(.), sd = varcomps[[2]])
  ) %>% 
  select(-time) %>% 
  mutate(
    hippocampus = as.numeric(predict(mod$gam, newdata = .) +
      random_intercept + residual) + 50 * (genotype == "A") -
      30 * (genotype == "C"),
  )

ggplot(simdat, aes(x = age, y = hippocampus, group = id)) +
  geom_point() +
  geom_line()

saveRDS(simdat, file = "data/simdat.rds")
