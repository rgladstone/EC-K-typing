
# Pneumo lineage vs. serotype effect
# Logistic mixed effects model with serotype as fixed effect and lineage as random effect

# Outcome is binary disease / carriage, 
# Build a model with two hierarchical levels: serotypes as fixed effects and lineage as random effect
# Take also the country into account (USA / South-Africa)

library(dplyr)
require(lme4)
library(lattice)
library(ggplot2)
library(optimx)
library(tidyr)
library(gplots)
library(circlize)
library(RColorBrewer)
library(tidyverse)

####### data
dat <- read.csv("~/Documents/Research/Rebecca_OR_logistic_hier_model/Abstract_ISPPD13/NT_update_Pneumo_lineage_invasiveness.csv", sep=";")


# 'Clinical_Manifest' is either carriage or disease
# The Serotype is 'In_Silico_Serotype'
# The lineage definition I would like to use is ‘Lineage-GPSC’ (ST’s sit in GPSCs, so it’s a slightly broader definition)

dat$Lineage.GPSC <- as.factor(dat$Lineage.GPSC)
dat$In_Silico_Serotype_NT_corrected <- as.factor(dat$In_Silico_Serotype_NT_corrected)
dat$Country <- as.factor(dat$Country)
dat$Clinical_Manifest <- as.factor(dat$Clinical_Manifest)
dat$invasive.bin <- ifelse(dat$Clinical_Manifest == "Carriage", 0, 1)

###### filter data

# In the logistic mixed model there is no need to filter lineages (random effects), but for the model convergence it is good to ensure
# that you have enough isolates per serotype (fixed effect), as well as avoid serotypes with 0 carriage or disease isolates

# filter such that smaller of the carriage/inf # of isolates is at least "thres"
# and the total number of isolates per serotype is at least "tot"

# these example values would lead to quite wide confidence intervals for some serotypes, so they could be stricker
thres <- 2
tot <- 15

st.inf <- as.data.frame(table(dat$In_Silico_Serotype_NT_corrected, dat$invasive.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(ST = Var1)

st.filt <- st.inf$ST

dat.filt <- dat %>%
  filter(In_Silico_Serotype_NT_corrected %in% st.filt)


# set NT as the reference category for the serotype variable (factor)

# unique(dat.filt$In_Silico_Serotype_NT_corrected)

dat.filt$In_Silico_Serotype.ref.NT <- factor(dat.filt$In_Silico_Serotype_NT_corrected,
                                             levels = c("NT", as.character(unique(dat.filt$In_Silico_Serotype_NT_corrected)[-which(unique(dat.filt$In_Silico_Serotype_NT_corrected) == "NT")])))


######### two-level hierarchical model

mm.random.lineage <- glmer(invasive.bin ~ In_Silico_Serotype.ref.NT + Country + (1 | Lineage.GPSC), 
                           data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                      optCtrl=list(maxfun=1e5)), nAGQ = 10)

summary(mm.random.lineage)
#vcov <- vcov(mm.random.lineage)

##### fixed effects for serotypes

# combine the fixed effect estimates + OR + CI + p-val in one data frame

se <- sqrt(diag(vcov(mm.random.lineage)))
tab <- cbind(Est = fixef(mm.random.lineage), LL = fixef(mm.random.lineage) - 1.96 * se, UL = fixef(mm.random.lineage) + 1.96 *se)
exp(round(tab, 3))

results <- as.data.frame(coef(summary(mm.random.lineage))) %>%
  mutate(parameter = row.names(coef(summary(mm.random.lineage)))) %>%
  mutate(short.parameter = gsub('In_Silico_Serotype.ref.NT', '', parameter)) %>%
  select(short.parameter, Estimate, "Std. Error", "Pr(>|z|)" )

ci <- as.data.frame(exp(round(tab, 3))) %>%
  mutate(parameter = row.names(exp(round(tab, 3)))) %>%
  mutate(short.parameter = gsub('In_Silico_Serotype.ref.NT', '', parameter)) %>%
  select(short.parameter, Est, LL, UL) %>%
  rename(OR = Est)

res.comb <- left_join(results, ci, by = "short.parameter") %>%
  select(short.parameter, Estimate, "Std. Error", OR, LL, UL, "Pr(>|z|)") %>%
  rename(p.value = "Pr(>|z|)", std.error = "Std. Error", Serotype = short.parameter)

res.comb$significance <- NA
res.comb$significance[res.comb$p.value < 0.1] <- "."
res.comb$significance[res.comb$p.value < 0.05] <- "*"
res.comb$significance[res.comb$p.value < 0.01] <- "**"
res.comb$significance[res.comb$p.value < 0.001] <- "***"
res.comb$significance[is.na(res.comb$significance)] <- ""

##### random effects for lineages

# plot of random effects for lineages (they are assumed to be normally ditributed around zero)
p1 <- lattice::dotplot(ranef(mm.random.lineage, which = "Lineage.GPSC", condVar = TRUE))
p1

# random effects in a data frame

ran <- ranef(mm.random.lineage, which = "Lineage.GPSC", condVar = TRUE)

rr1 <- ranef(mm.random.lineage)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "lineage"
dd.sort.sel <- dd.sort[,c(3:5)]



