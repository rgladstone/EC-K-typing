
# Capsular K-antigens and invasive potential of Escherichia coli

# Logistic mixed models with infection/carriage as outcome, k-antigen as fixed eefect and lineage as random effect run on different data:

# 1) The full dataset 2003-2017 with the updated K-designations

# + sensitivity analyses:

# 2) BSAC restricted to adults (>19 years old) only n=1550/2036
# 3) BSAC elderly adults only (>59 year old) n=1202/2036
# 4) BSAC with ST393 subsampled in 2003-2013 down to the same proportion as 2014-2017 of 0.02 n=2008/2036
# 5) Adjusting for the BSAC CTX-M prevalence of each lineage n=2036 (it is not possible to determine the CTX-M status of individual babybiome assemblies so a lineage based penalty is the only way to do it)



# Outcome is binary (bloodstream infection / carriage), 
# and these two levels are equal to the two different collections (BSAC = infection, BB = carriage)
# Build a model with two hierarchical levels, SC_clade -> K

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


######## 1) The full dataset 2003-2017 with the updated K-designations

# data
BSAC <- read.csv("BSAC_FULL_pp25_2003_2018_clades_020625.csv")
BB <- read.csv("BB_SCandK_020625.csv")

# correct one clade typo
BB$ST131.clade[which(BB$ST131.clade == "C")] <- "C0"

# use the CC name only for the lineage:
# - in BSAC: "pp_CCname"
# - in BB: "pp_CCname"

BSAC$ST131.clade[which(BSAC$ST131.clade == "PP195_CCnew")] <- ""

# full BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp_CCname, ST131.clade, Designation, collection, infection.bin, clinical_manifest) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id)

bsac <- BSAC %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))

dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)


#### filter capsules

# filter such that smaller of the carriage/inf # of isolates is at least thres
thres <- 5
tot <- 20

k.inf <- as.data.frame(table(dat$K, dat$infection.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(K = Var1)

k.filt <- k.inf$K

dat.filt <- dat %>%
  filter(K %in% k.filt)
unique(dat.filt$cc_clade)
sort(table(dat.filt$cc_clade), decreasing = T)

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))

# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + (1 | cc_clade), 
                            data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                       optCtrl=list(maxfun=1e5)), nAGQ = 10)

s <- summary(mm.random.SC_clade)


write.table(s$coefficients, "1_Full_data_set/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

# ORs with CIs for fixed effects
se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
exp(round(tab, 3))
write.table(exp(round(tab, 3)), "1_Full_data_set/log.mixed.K_OR.txt", quote = F)

fixed.k.est <- exp(round(tab, 3))

# random effects
ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "1_Full_data_set/random.effects.for.SC_clade.txt", quote = F, row.names = F)


################# Heatmap of the results

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove the first K from the name
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))
# remove NT
st.fixed.log.scale.no.nt <- st.fixed.log.scale[-which(st.fixed.log.scale$st == "Untypeable"),]
# order according to the estimate
st.fixed.log.scale.no.nt$st <- factor(st.fixed.log.scale.no.nt$st, levels=(st.fixed.log.scale.no.nt$st)[order(st.fixed.log.scale.no.nt$est, decreasing = T)])
st.fixed.log.scale.ord <- st.fixed.log.scale.no.nt[order(st.fixed.log.scale.no.nt$est),]

lin.ran.log.scale <- dd.sort.sel %>%
                    select(condval)

rownames(lin.ran.log.scale) <- dd.sort.sel[,1]
lin.ran.log.scale$lineage <- rownames(lin.ran.log.scale)
lin.ran.log.scale$condval <- as.numeric(lin.ran.log.scale$condval)
lin.ran.log.scale.ord <- lin.ran.log.scale %>% map_df(rev)

rownames(lin.ran.log.scale.ord) <- lin.ran.log.scale.ord$lineage
lin.ran.log.scale.ord$row <- as.factor(rep("row1", nrow(lin.ran.log.scale.ord)))
lin.ran.log.scale.ord$lineage <- factor(lin.ran.log.scale.ord$lineage, levels=(lin.ran.log.scale.ord$lineage)[order(lin.ran.log.scale.ord$condval, decreasing = T)])


# plot onløy top 40 largest lineages
lineage.top <- 40

# select only lineages that have capsule information

lineage.temp <- as.data.frame(sort(table(dat.filt$cc_clade), decreasing = T)) %>%
  rename(lineage = Var1)  %>%
  top_n(lineage.top)

lineage.filt <- lineage.temp$lineage

plot.dat <- dat.filt %>%
  filter(cc_clade %in% lineage.filt)

# combinations of serotype and lineage
temp <- data.frame(comb = unique(plot.dat$cc_clade_k))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}


temp$comb.est <- NA

# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")

st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)

for(i in 1:nrow(temp)){
  temp$comb.est[i] <- st.fixed.log.scale.ord$est[which(st.fixed.log.scale.ord$st == temp$serotype[i])] +
    lin.ran.log.scale$condval[which(lin.ran.log.scale$lineage == temp$lineage[i])]
}



# fill matrix
# all combinations of lineage and serotype

# lineage estimate data filtered
lin.ran.log.scale.filt <- lin.ran.log.scale %>%
  filter(lineage %in% lineage.filt)

all <- expand.grid(lineage = lin.ran.log.scale.filt$lineage, serotype = st.fixed.log.scale.ord$st)
heat.comb <- left_join(all, temp, by = c("lineage", "serotype"))
heat.comb$lineage <- factor(heat.comb$lineage, levels=(lin.ran.log.scale.filt$lineage)[order(lin.ran.log.scale.filt$condval, decreasing = F)])
heat.comb$serotype <- factor(heat.comb$serotype, levels=(st.fixed.log.scale.ord$st)[order(st.fixed.log.scale.ord$est, decreasing = T)])

# sts over all k

lin.over.st <- lin.ran.log.scale.filt
lin.over.st$serotype <- rep(" ", nrow(lin.over.st))
lin.over.st$lineage <- lin.ran.log.scale.filt$lineage
lin.over.st$comb <- NA
names(lin.over.st)[which(names(lin.over.st) == "condval")] <- "comb.est"
lin.over.st <- lin.over.st %>% 
  select(lineage, serotype, comb, comb.est)

# ks over all st
st.over.lin <- st.fixed.log.scale.ord
st.over.lin$lineage <- rep(" ", nrow(st.over.lin))
st.over.lin$comb <- NA
names(st.over.lin)[which(names(st.over.lin) == "est")] <- "comb.est"
st.over.lin <- st.over.lin %>% 
  select(lineage, st, comb, comb.est) %>%
  rename(serotype = st)

# combine
heatmap.all <- rbind(heat.comb, lin.over.st, st.over.lin)
heatmap.all$lineage <- factor(heatmap.all$lineage, levels=levels(heatmap.all$lineage))
heatmap.all$serotype <- factor(heatmap.all$serotype, levels = c(" ", levels(heatmap.all$serotype)[-length(levels(heatmap.all$serotype))]))


# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# heatmap

display.brewer.pal(n = 9, name = 'Oranges')
oranges <- brewer.pal(n = 9, name = 'Oranges')

heat <- ggplot(heatmap.all2, aes(x = serotype, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE",
                                   "#D1E5F0",
                                   "#FDDBC7","#F4A582", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603"),
                       na.value = "white",
                       values=scales:::rescale(c(-1.5, -0.5,
                                                 0,
                                                 1, 1.5, 2, 2.5, 3, 4, 5)))+
  
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")
  ) +
  labs(x="Lineage", y="In silico serotype", fill="Estimate")
heat
ggsave("1_Full_data_set/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)





######## 2) BSAC adults only (>19 years old)


# use the CC name only for the lineage:
# - in BSAC: "pp_CCname"
# - in BB: "pp_CCname"

BSAC$ST131.clade[which(BSAC$ST131.clade == "PP195_CCnew")] <- ""

# full BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp_CCname, ST131.clade, Designation,collection, infection.bin, clinical_manifest) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id)


# BSAC adults only
bsac <- BSAC %>%
  filter(age.subsets == ">=20-59" | age.subsets == ">=60+") %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))


dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)


#### filter capsules

# filter such that smaller of the carriage/inf # of isolates is at least thres
thres <- 5
tot <- 20

k.inf <- as.data.frame(table(dat$K, dat$infection.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(K = Var1)

k.filt <- k.inf$K

dat.filt <- dat %>%
  filter(K %in% k.filt)
unique(dat.filt$cc_clade)
sort(table(dat.filt$cc_clade), decreasing = T)

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))


# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + (1 | cc_clade), 
                            data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                       optCtrl=list(maxfun=1e5)), nAGQ = 10)

s <- summary(mm.random.SC_clade)

write.table(s$coefficients, "2_BSAC_adults_only/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
exp(round(tab, 3))
write.table(exp(round(tab, 3)), "2_BSAC_adults_only/log.mixed.K_OR.txt", quote = F)

fixed.k.est <- exp(round(tab, 3))

ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "2_BSAC_adults_only/random.effects.for.SC_clade.txt", quote = F, row.names = F)


################# heatmap of the results

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove the first K
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))
# remove NT
st.fixed.log.scale.no.nt <- st.fixed.log.scale[-which(st.fixed.log.scale$st == "Untypeable"),]
# order according to the estimate
st.fixed.log.scale.no.nt$st <- factor(st.fixed.log.scale.no.nt$st, levels=(st.fixed.log.scale.no.nt$st)[order(st.fixed.log.scale.no.nt$est, decreasing = T)])
st.fixed.log.scale.ord <- st.fixed.log.scale.no.nt[order(st.fixed.log.scale.no.nt$est),]

lin.ran.log.scale <- dd.sort.sel %>%
  select(condval)
rownames(lin.ran.log.scale) <- dd.sort.sel[,1]
lin.ran.log.scale$lineage <- rownames(lin.ran.log.scale)
lin.ran.log.scale$condval <- as.numeric(lin.ran.log.scale$condval)
lin.ran.log.scale.ord <- lin.ran.log.scale %>% map_df(rev)
rownames(lin.ran.log.scale.ord) <- lin.ran.log.scale.ord$lineage
lin.ran.log.scale.ord$row <- as.factor(rep("row1", nrow(lin.ran.log.scale.ord)))
lin.ran.log.scale.ord$lineage <- factor(lin.ran.log.scale.ord$lineage, levels=(lin.ran.log.scale.ord$lineage)[order(lin.ran.log.scale.ord$condval, decreasing = T)])


# plot only top 40 lineages
lineage.top <- 40

# select only lineages that have capsule information

lineage.temp <- as.data.frame(sort(table(dat.filt$cc_clade), decreasing = T)) %>%
  rename(lineage = Var1)  %>%
  top_n(lineage.top)
lineage.filt <- lineage.temp$lineage
plot.dat <- dat.filt %>%
  filter(cc_clade %in% lineage.filt)

# combinations of serotype and lineage
temp <- data.frame(comb = unique(plot.dat$cc_clade_k))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}


temp$comb.est <- NA


# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")

st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)


for(i in 1:nrow(temp)){
  temp$comb.est[i] <- st.fixed.log.scale.ord$est[which(st.fixed.log.scale.ord$st == temp$serotype[i])] +
    lin.ran.log.scale$condval[which(lin.ran.log.scale$lineage == temp$lineage[i])]
}


# fill matrix
# all combinations of lineage and serotype

# lineage estimate data filtered
lin.ran.log.scale.filt <- lin.ran.log.scale %>%
  filter(lineage %in% lineage.filt)

all <- expand.grid(lineage = lin.ran.log.scale.filt$lineage, serotype = st.fixed.log.scale.ord$st)
heat.comb <- left_join(all, temp, by = c("lineage", "serotype"))
heat.comb$lineage <- factor(heat.comb$lineage, levels=(lin.ran.log.scale.filt$lineage)[order(lin.ran.log.scale.filt$condval, decreasing = F)])
heat.comb$serotype <- factor(heat.comb$serotype, levels=(st.fixed.log.scale.ord$st)[order(st.fixed.log.scale.ord$est, decreasing = T)])

# sts over all k

lin.over.st <- lin.ran.log.scale.filt
lin.over.st$serotype <- rep(" ", nrow(lin.over.st))
lin.over.st$lineage <- lin.ran.log.scale.filt$lineage
lin.over.st$comb <- NA
names(lin.over.st)[which(names(lin.over.st) == "condval")] <- "comb.est"
lin.over.st <- lin.over.st %>% 
  select(lineage, serotype, comb, comb.est)

# ks over all st
st.over.lin <- st.fixed.log.scale.ord
st.over.lin$lineage <- rep(" ", nrow(st.over.lin))
st.over.lin$comb <- NA
names(st.over.lin)[which(names(st.over.lin) == "est")] <- "comb.est"
st.over.lin <- st.over.lin %>% 
  select(lineage, st, comb, comb.est) %>%
  rename(serotype = st)

# combine
heatmap.all <- rbind(heat.comb, lin.over.st, st.over.lin)
heatmap.all$lineage <- factor(heatmap.all$lineage, levels=levels(heatmap.all$lineage))
heatmap.all$serotype <- factor(heatmap.all$serotype, levels = c(" ", levels(heatmap.all$serotype)[-length(levels(heatmap.all$serotype))]))


# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# heatmap

display.brewer.pal(n = 9, name = 'Oranges')
oranges <- brewer.pal(n = 9, name = 'Oranges')

heat <- ggplot(heatmap.all2, aes(x = serotype, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE",
                                   "#D1E5F0",
                                   "#FDDBC7","#F4A582", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603"),
                       na.value = "white",
                       values=scales:::rescale(c(-1.5, -0.5,
                                                 0,
                                                 1, 1.5, 2, 2.5, 3, 4, 5)))+
  
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")
  ) +
  labs(x="Lineage", y="In silico serotype", fill="Estimate")
heat
ggsave("2_BSAC_adults_only/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)




######## 3) BSAC elderly adults only (>59 year old)



# full BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp_CCname, ST131.clade, Designation,collection, infection.bin, clinical_manifest) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id)


# BSAC elderly adults only
bsac <- BSAC %>%
  filter(age.subsets == ">=60+") %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))


dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)

#### filter capsules

# filter such that smaller of the carriage/inf # of isolates is at least thres

thres <- 5
tot <- 20

k.inf <- as.data.frame(table(dat$K, dat$infection.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(K = Var1)

k.filt <- k.inf$K

dat.filt <- dat %>%
  filter(K %in% k.filt)
unique(dat.filt$cc_clade)
sort(table(dat.filt$cc_clade), decreasing = T)

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))


# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + (1 | cc_clade), 
                            data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                       optCtrl=list(maxfun=1e5)), nAGQ = 10)

s <- summary(mm.random.SC_clade)


write.table(s$coefficients, "3_BSAC_elderly_adults_only/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
exp(round(tab, 3))
write.table(exp(round(tab, 3)), "3_BSAC_elderly_adults_only/log.mixed.K_OR.txt", quote = F)

fixed.k.est <- exp(round(tab, 3))

ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "3_BSAC_elderly_adults_only/random.effects.for.SC_clade.txt", quote = F, row.names = F)


################# heatmap of the results

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove the first K
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))
# remove NT
st.fixed.log.scale.no.nt <- st.fixed.log.scale[-which(st.fixed.log.scale$st == "Untypeable"),]
# order according to the estimate
st.fixed.log.scale.no.nt$st <- factor(st.fixed.log.scale.no.nt$st, levels=(st.fixed.log.scale.no.nt$st)[order(st.fixed.log.scale.no.nt$est, decreasing = T)])
st.fixed.log.scale.ord <- st.fixed.log.scale.no.nt[order(st.fixed.log.scale.no.nt$est),]

lin.ran.log.scale <- dd.sort.sel %>%
  select(condval)
rownames(lin.ran.log.scale) <- dd.sort.sel[,1]
lin.ran.log.scale$lineage <- rownames(lin.ran.log.scale)
lin.ran.log.scale$condval <- as.numeric(lin.ran.log.scale$condval)
lin.ran.log.scale.ord <- lin.ran.log.scale %>% map_df(rev)
rownames(lin.ran.log.scale.ord) <- lin.ran.log.scale.ord$lineage
lin.ran.log.scale.ord$row <- as.factor(rep("row1", nrow(lin.ran.log.scale.ord)))
lin.ran.log.scale.ord$lineage <- factor(lin.ran.log.scale.ord$lineage, levels=(lin.ran.log.scale.ord$lineage)[order(lin.ran.log.scale.ord$condval, decreasing = T)])


# plot top 40 lineages
lineage.top <- 40

# select only lineages that have capsule information

lineage.temp <- as.data.frame(sort(table(dat.filt$cc_clade), decreasing = T)) %>%
  rename(lineage = Var1)  %>%
  top_n(lineage.top)
lineage.filt <- lineage.temp$lineage
plot.dat <- dat.filt %>%
  filter(cc_clade %in% lineage.filt)

# combinations of serotype and lineage
temp <- data.frame(comb = unique(plot.dat$cc_clade_k))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}


temp$comb.est <- NA


# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")


st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)


for(i in 1:nrow(temp)){
  temp$comb.est[i] <- st.fixed.log.scale.ord$est[which(st.fixed.log.scale.ord$st == temp$serotype[i])] +
    lin.ran.log.scale$condval[which(lin.ran.log.scale$lineage == temp$lineage[i])]
}


# fill matrix
# all combinations of lineage and serotype

# lineage estimate data filtered
lin.ran.log.scale.filt <- lin.ran.log.scale %>%
  filter(lineage %in% lineage.filt)

all <- expand.grid(lineage = lin.ran.log.scale.filt$lineage, serotype = st.fixed.log.scale.ord$st)
heat.comb <- left_join(all, temp, by = c("lineage", "serotype"))
heat.comb$lineage <- factor(heat.comb$lineage, levels=(lin.ran.log.scale.filt$lineage)[order(lin.ran.log.scale.filt$condval, decreasing = F)])
heat.comb$serotype <- factor(heat.comb$serotype, levels=(st.fixed.log.scale.ord$st)[order(st.fixed.log.scale.ord$est, decreasing = T)])

# sts over all k

lin.over.st <- lin.ran.log.scale.filt
lin.over.st$serotype <- rep(" ", nrow(lin.over.st))
lin.over.st$lineage <- lin.ran.log.scale.filt$lineage
lin.over.st$comb <- NA
names(lin.over.st)[which(names(lin.over.st) == "condval")] <- "comb.est"
lin.over.st <- lin.over.st %>% 
  select(lineage, serotype, comb, comb.est)

# ks over all st
st.over.lin <- st.fixed.log.scale.ord
st.over.lin$lineage <- rep(" ", nrow(st.over.lin))
st.over.lin$comb <- NA
names(st.over.lin)[which(names(st.over.lin) == "est")] <- "comb.est"
st.over.lin <- st.over.lin %>% 
  select(lineage, st, comb, comb.est) %>%
  rename(serotype = st)

# combine
heatmap.all <- rbind(heat.comb, lin.over.st, st.over.lin)
heatmap.all$lineage <- factor(heatmap.all$lineage, levels=levels(heatmap.all$lineage))
heatmap.all$serotype <- factor(heatmap.all$serotype, levels = c(" ", levels(heatmap.all$serotype)[-length(levels(heatmap.all$serotype))]))


# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# heatmap

display.brewer.pal(n = 9, name = 'Oranges')
oranges <- brewer.pal(n = 9, name = 'Oranges')

heat <- ggplot(heatmap.all2, aes(x = serotype, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE",
                                   "#D1E5F0",
                                   "#FDDBC7","#F4A582", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603"),
                       na.value = "white",
                       values=scales:::rescale(c(-1.5, -0.5,
                                                 0,
                                                 1, 1.5, 2, 2.5, 3, 4, 5)))+
  
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")
  ) +
  labs(x="Lineage", y="In silico serotype", fill="Estimate")
heat
ggsave("3_BSAC_elderly_adults_only/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)



######## 4) BSAC ST393 subsampled


# full BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp_CCname, ST131.clade, Designation,collection, infection.bin, clinical_manifest) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id)

# BSAC: select only subsampled
bsac <- BSAC %>%
  filter(ST393.subsample != "Exclude") %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))


dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)


#### filter capsules

# filter such that smaller of the carriage/inf # of isolates is at least thres

thres <- 5
tot <- 20

k.inf <- as.data.frame(table(dat$K, dat$infection.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(K = Var1)

k.filt <- k.inf$K

dat.filt <- dat %>%
  filter(K %in% k.filt)
unique(dat.filt$cc_clade)
sort(table(dat.filt$cc_clade), decreasing = T)

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))


# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + (1 | cc_clade), 
                            data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                       optCtrl=list(maxfun=1e5)), nAGQ = 10)

s <- summary(mm.random.SC_clade)


write.table(s$coefficients, "4_BSAC_ST393_subsampled/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
exp(round(tab, 3))
write.table(exp(round(tab, 3)), "4_BSAC_ST393_subsampled/log.mixed.K_OR.txt", quote = F)

fixed.k.est <- exp(round(tab, 3))

ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "4_BSAC_ST393_subsampled/random.effects.for.SC_clade.txt", quote = F, row.names = F)


################# heatmap of the results, version without NT in the marginal and not intercept in the combined estimate

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove the first K
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))
# remove NT
st.fixed.log.scale.no.nt <- st.fixed.log.scale[-which(st.fixed.log.scale$st == "Untypeable"),]
# order according to the estimate
st.fixed.log.scale.no.nt$st <- factor(st.fixed.log.scale.no.nt$st, levels=(st.fixed.log.scale.no.nt$st)[order(st.fixed.log.scale.no.nt$est, decreasing = T)])
st.fixed.log.scale.ord <- st.fixed.log.scale.no.nt[order(st.fixed.log.scale.no.nt$est),]

lin.ran.log.scale <- dd.sort.sel %>%
  select(condval)
rownames(lin.ran.log.scale) <- dd.sort.sel[,1]
lin.ran.log.scale$lineage <- rownames(lin.ran.log.scale)
lin.ran.log.scale$condval <- as.numeric(lin.ran.log.scale$condval)
lin.ran.log.scale.ord <- lin.ran.log.scale %>% map_df(rev)
rownames(lin.ran.log.scale.ord) <- lin.ran.log.scale.ord$lineage
lin.ran.log.scale.ord$row <- as.factor(rep("row1", nrow(lin.ran.log.scale.ord)))
lin.ran.log.scale.ord$lineage <- factor(lin.ran.log.scale.ord$lineage, levels=(lin.ran.log.scale.ord$lineage)[order(lin.ran.log.scale.ord$condval, decreasing = T)])


# plot top 40 lineages
lineage.top <- 40

# select only lineages that have capsule information

lineage.temp <- as.data.frame(sort(table(dat.filt$cc_clade), decreasing = T)) %>%
  rename(lineage = Var1)  %>%
  top_n(lineage.top)
lineage.filt <- lineage.temp$lineage
plot.dat <- dat.filt %>%
  filter(cc_clade %in% lineage.filt)

# combinations of serotype and lineage
temp <- data.frame(comb = unique(plot.dat$cc_clade_k))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}


temp$comb.est <- NA


# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")


st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)


for(i in 1:nrow(temp)){
  temp$comb.est[i] <- st.fixed.log.scale.ord$est[which(st.fixed.log.scale.ord$st == temp$serotype[i])] +
    lin.ran.log.scale$condval[which(lin.ran.log.scale$lineage == temp$lineage[i])]
}



# fill matrix
# all combinations of lineage and serotype

# lineage estimate data filtered
lin.ran.log.scale.filt <- lin.ran.log.scale %>%
  filter(lineage %in% lineage.filt)

all <- expand.grid(lineage = lin.ran.log.scale.filt$lineage, serotype = st.fixed.log.scale.ord$st)
heat.comb <- left_join(all, temp, by = c("lineage", "serotype"))
heat.comb$lineage <- factor(heat.comb$lineage, levels=(lin.ran.log.scale.filt$lineage)[order(lin.ran.log.scale.filt$condval, decreasing = F)])
heat.comb$serotype <- factor(heat.comb$serotype, levels=(st.fixed.log.scale.ord$st)[order(st.fixed.log.scale.ord$est, decreasing = T)])

# sts over all k

lin.over.st <- lin.ran.log.scale.filt
lin.over.st$serotype <- rep(" ", nrow(lin.over.st))
lin.over.st$lineage <- lin.ran.log.scale.filt$lineage
lin.over.st$comb <- NA
names(lin.over.st)[which(names(lin.over.st) == "condval")] <- "comb.est"
lin.over.st <- lin.over.st %>% 
  select(lineage, serotype, comb, comb.est)

# ks over all st
st.over.lin <- st.fixed.log.scale.ord
st.over.lin$lineage <- rep(" ", nrow(st.over.lin))
st.over.lin$comb <- NA
names(st.over.lin)[which(names(st.over.lin) == "est")] <- "comb.est"
# remove country
st.over.lin <- st.over.lin %>% 
  select(lineage, st, comb, comb.est) %>%
  rename(serotype = st)

# combine
heatmap.all <- rbind(heat.comb, lin.over.st, st.over.lin)
heatmap.all$lineage <- factor(heatmap.all$lineage, levels=levels(heatmap.all$lineage))
heatmap.all$serotype <- factor(heatmap.all$serotype, levels = c(" ", levels(heatmap.all$serotype)[-length(levels(heatmap.all$serotype))]))


# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# heatmap

display.brewer.pal(n = 9, name = 'Oranges')
oranges <- brewer.pal(n = 9, name = 'Oranges')

heat <- ggplot(heatmap.all2, aes(x = serotype, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE",
                                   "#D1E5F0",
                                   "#FDDBC7","#F4A582", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603"),
                       na.value = "white",
                       values=scales:::rescale(c(-1.5, -0.5,
                                                 0,
                                                 1, 1.5, 2, 2.5, 3, 4, 5)))+
  
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")
  ) +
  labs(x="Lineage", y="In silico serotype", fill="Estimate")
heat
ggsave("/Users/maijupesonen/Documents/Research/Rebecca_OR_logistic_hier_model/E_coli/Review_June2025/4_BSAC_ST393_subsampled/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)




############ 5) Weighting for the resistance


# resistance probabilities
BSAC_prop_CTX.M <- read.csv2("BSAC_prop_CTX-M.csv")

ctxm <- BSAC_prop_CTX.M %>%
  mutate(inv.weight.ctxm = 1-Proportion) %>%
  mutate(weight.ctxm = Proportion) %>%
  select(Strain, weight.ctxm, inv.weight.ctxm) %>%
  rename(ppCC = Strain)

# use the CC name only for the lineage:
# - in BSAC: "pp_CCname"
# - in BB: "pp_CCname"


# full BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp2025_CC, pp_CCname, ST131.clade, Designation,collection, infection.bin, clinical_manifest) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id) %>%
  rename(ppCC = pp2025_CC)


bsac <- BSAC %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))


dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)

dat <- left_join(dat, ctxm, by = "ppCC") %>%
  filter(!is.na(weight.ctxm))


#### filter capsules

# filter such that smaller of the carriage/inf # of isolates is at least thres

thres <- 5
tot <- 20

k.inf <- as.data.frame(table(dat$K, dat$infection.bin)) %>%
  mutate(inf = ifelse(Var2 == 0, "carriage", "infection")) %>%
  select(Var1, inf, Freq) %>%
  pivot_wider(., names_from = inf, values_from = Freq) %>%
  rowwise() %>%
  mutate(min = min(carriage, infection), tot.n = sum(carriage, infection)) %>%
  filter(min > thres & tot.n > tot) %>%
  rename(K = Var1)

k.filt <- k.inf$K

dat.filt <- dat %>%
  filter(K %in% k.filt)
unique(dat.filt$cc_clade)
sort(table(dat.filt$cc_clade), decreasing = T)

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))


# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + weight.ctxm + (1 | cc_clade), 
                            data = dat.filt, family = binomial,
                            control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5)), nAGQ = 10)


s <- summary(mm.random.SC_clade)


write.table(s$coefficients, "5_Adjusting_CTXM_prevalence/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
exp(round(tab, 3))
write.table(exp(round(tab, 3))[-20,], "5_Adjusting_CTXM_prevalence/log.mixed.K_OR.txt", quote = F)

fixed.k.est <- exp(round(tab, 3))[-20,]

ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
  ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)
}

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "5_Adjusting_CTXM_prevalence/random.effects.for.SC_clade.txt", quote = F, row.names = F)


################# heatmap of the results

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[-20,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove the first K
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))
# remove NT
st.fixed.log.scale.no.nt <- st.fixed.log.scale[-which(st.fixed.log.scale$st == "Untypeable"),]
# order according to the estimate
st.fixed.log.scale.no.nt$st <- factor(st.fixed.log.scale.no.nt$st, levels=(st.fixed.log.scale.no.nt$st)[order(st.fixed.log.scale.no.nt$est, decreasing = T)])
st.fixed.log.scale.ord <- st.fixed.log.scale.no.nt[order(st.fixed.log.scale.no.nt$est),]

lin.ran.log.scale <- dd.sort.sel %>%
  select(condval)
rownames(lin.ran.log.scale) <- dd.sort.sel[,1]
lin.ran.log.scale$lineage <- rownames(lin.ran.log.scale)
lin.ran.log.scale$condval <- as.numeric(lin.ran.log.scale$condval)
lin.ran.log.scale.ord <- lin.ran.log.scale %>% map_df(rev)
rownames(lin.ran.log.scale.ord) <- lin.ran.log.scale.ord$lineage
lin.ran.log.scale.ord$row <- as.factor(rep("row1", nrow(lin.ran.log.scale.ord)))
lin.ran.log.scale.ord$lineage <- factor(lin.ran.log.scale.ord$lineage, levels=(lin.ran.log.scale.ord$lineage)[order(lin.ran.log.scale.ord$condval, decreasing = T)])


# plot top 40 lineages
lineage.top <- 40

# select only lineages that have capsule information

lineage.temp <- as.data.frame(sort(table(dat.filt$cc_clade), decreasing = T)) %>%
  rename(lineage = Var1)  %>%
  top_n(lineage.top)
lineage.filt <- lineage.temp$lineage
plot.dat <- dat.filt %>%
  filter(cc_clade %in% lineage.filt)

# combinations of serotype and lineage

temp <- data.frame(comb = unique(plot.dat$cc_clade_k))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"),
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("K54_K96", temp$comb[i], fixed = T) | grepl("K6_K24", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}


temp$comb.est <- NA


# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")


st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)


for(i in 1:nrow(temp)){
  temp$comb.est[i] <- st.fixed.log.scale.ord$est[which(st.fixed.log.scale.ord$st == temp$serotype[i])] +
    lin.ran.log.scale$condval[which(lin.ran.log.scale$lineage == temp$lineage[i])]
}


# fill matrix
# all combinations of lineage and serotype

# lineage estimate data filtered
lin.ran.log.scale.filt <- lin.ran.log.scale %>%
  filter(lineage %in% lineage.filt)

all <- expand.grid(lineage = lin.ran.log.scale.filt$lineage, serotype = st.fixed.log.scale.ord$st)
heat.comb <- left_join(all, temp, by = c("lineage", "serotype"))
heat.comb$lineage <- factor(heat.comb$lineage, levels=(lin.ran.log.scale.filt$lineage)[order(lin.ran.log.scale.filt$condval, decreasing = F)])
heat.comb$serotype <- factor(heat.comb$serotype, levels=(st.fixed.log.scale.ord$st)[order(st.fixed.log.scale.ord$est, decreasing = T)])

# sts over all k

lin.over.st <- lin.ran.log.scale.filt
lin.over.st$serotype <- rep(" ", nrow(lin.over.st))
lin.over.st$lineage <- lin.ran.log.scale.filt$lineage
lin.over.st$comb <- NA
names(lin.over.st)[which(names(lin.over.st) == "condval")] <- "comb.est"
lin.over.st <- lin.over.st %>% 
  select(lineage, serotype, comb, comb.est)

# ks over all st
st.over.lin <- st.fixed.log.scale.ord
st.over.lin$lineage <- rep(" ", nrow(st.over.lin))
st.over.lin$comb <- NA
names(st.over.lin)[which(names(st.over.lin) == "est")] <- "comb.est"

st.over.lin <- st.over.lin %>% 
  select(lineage, st, comb, comb.est) %>%
  rename(serotype = st)

# combine
heatmap.all <- rbind(heat.comb, lin.over.st, st.over.lin)
heatmap.all$lineage <- factor(heatmap.all$lineage, levels=levels(heatmap.all$lineage))
heatmap.all$serotype <- factor(heatmap.all$serotype, levels = c(" ", levels(heatmap.all$serotype)[-length(levels(heatmap.all$serotype))]))


# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# heatmap

display.brewer.pal(n = 9, name = 'Oranges')
oranges <- brewer.pal(n = 9, name = 'Oranges')

heat <- ggplot(heatmap.all2, aes(x = serotype, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  
  scale_fill_gradientn(colours = c("#2166AC", "#92C5DE",
                                   "#D1E5F0",
                                   "#FDDBC7","#F4A582", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603"),
                       na.value = "white",
                       values=scales:::rescale(c(-1.5, -0.5,
                                                 0,
                                                 1, 1.5, 2, 2.5, 3, 4, 5)))+
  
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")
  ) +
  labs(x="Lineage", y="In silico serotype", fill="Estimate")
heat
ggsave("5_Adjusting_CTXM_prevalence/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)





