

# Group 2 and 3 ABC-transporter dependant capsular K-loci contribute significantly to variation in the invasive potential of Escherichia coli

# R-code for the logistic mixed model to estimate marginal 
# and combined invasive potential of different K-loci and lineages.
# The model includes clinical manifest (bloodstream infection / carriage) as a binary outcome variable, 
# K-loci indicator variable as a fixed effect (constant for each isolate) and
# lineage indicator variable as a random effect. 
# Untypeable is set as a reference category for K-loci.


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
library(plyr)


######## Data

# with clades for ST131
BSAC <- read.csv("~/Data/BSAC_OKH_2003_2018_clades_071124.csv")
BB <- read.csv("~/Data/BB_SCandK_lost_K_061124.csv") 


# Outcome is binary (bloodstream infection / carriage), 
# and these two levels are equal to the two different collections (BSAC = infection, BB = carriage)

bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, Strain, Designation, O, H, ST131.clade, collection, infection.bin, clinical_manifest) %>%
  rename(lineage = Strain) %>%
  rename(K = Designation) %>%
  rename(lane = lane_id)

# fix one wrong K name in BB data
bb$K[which(bb$K == "KL3A2_K2ab")] <- "KL3_K2ab"

bsac <- BSAC %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection")

dat <- rbind(bb, bsac)
dat$ST131.clade[which(dat$ST131.clade == 0)] <- ""

dat <- dat %>% 
  mutate(ST = sapply(strsplit(dat$lineage,split='_', fixed=TRUE), "[[", 2)) %>%
  mutate(CC = gsub("ST", "CC", ST, fixed = TRUE)) %>%
  mutate(SC = sapply(strsplit(dat$lineage,split='_', fixed=TRUE), "[[", 1)) %>%
  mutate(SC_ST_clade = ifelse(ST131.clade != "", paste(SC, ST, ST131.clade, sep = "_"), paste(SC, ST, sep = "_"))) %>%
  mutate(SC_ST_clade_K = paste(SC_ST_clade, K, sep = "_"))

dat$SC_ST_clade <- as.factor(dat$SC_ST_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)


######## Filter capsules

# Filter such that smaller of the carriage/inf # of isolates is at least thres.
# No need to filter by lineage (random effect).

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

# proportions of carriage isolates per K
prop.carriage.k <- sort(table(dat.filt$clinical_manifest, dat.filt$K)["carriage",] / colSums(table(dat.filt$clinical_manifest, dat.filt$K)), decreasing = T)
# set the least invasive K as the reference: untypeable
dat.filt$K <- factor(dat.filt$K, levels = names(prop.carriage.k))

######## Model

# Two levels: SC_clade --> K

# two-level hierarchical model
mm.random.SC_clade <- glmer(infection.bin ~ K + (1 | SC_ST_clade), 
                            data = dat.filt, family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                                       optCtrl=list(maxfun=1e5)), nAGQ = 10)

s <- summary(mm.random.SC_clade)

write.table(s$coefficients, "~/Results_thres5_tot20/log.mixed.fixed.K.estimates.log.scale.txt", quote = F)

# OR + CI of the fixed effects
se <- sqrt(diag(vcov(mm.random.SC_clade)))
tab <- cbind(Est = fixef(mm.random.SC_clade), LL = fixef(mm.random.SC_clade) - 1.96 * se, UL = fixef(mm.random.SC_clade) + 1.96 *se)
fixed.k.est <- exp(round(tab, 3))
write.table(fixed.k.est, "~/Results_thres5_tot20/log.mixed.K_OR.txt", quote = F)

# random effect estimates
ran <- ranef(mm.random.SC_clade, which = "SC_ST_clade", condVar = TRUE)

# plot of the random effects
p1 <- lattice::dotplot(ran)
pdf("~/Results_thres5_tot20/SC_clade_Raneffs.stKcomb.pdf",  width = 12, height = 21) 
p1
dev.off()

# estimates with condsd
rr1 <- ranef(mm.random.SC_clade)
str(dd <- as.data.frame(rr1))
# if (require(ggplot2)) {
#   ggplot(dd, aes(y=grp,x=condval)) +
#     geom_point() + facet_wrap(~term,scales="free_x") +
#     geom_errorbarh(aes(xmin=condval -2*condsd,
#                        xmax=condval +2*condsd), height=0)
# }

dd.sort <- dd[order(dd$condval, decreasing = T),]
names(dd.sort)[3] <- "SC_ST_clade"
dd.sort.sel <- dd.sort[,c(3:5)]

write.table(dd.sort.sel, "~/Results_thres5_tot20/random.effects.for.SC_clade.txt", quote = F, row.names = F)


######## Heatmap plot

# data for heatmap
# estimates on log-scale

s <- summary(mm.random.SC_clade)

st.fixed.log.scale <- data.frame(est = s$coefficients[,1])
rownames(st.fixed.log.scale)[1] <- "Untypeable"
# remove first K from the names
rownames(st.fixed.log.scale) <- c(rownames(st.fixed.log.scale)[1], sub('.', '', rownames(st.fixed.log.scale)[-1]))

st.fixed.log.scale$est <- as.numeric(st.fixed.log.scale$est)
st.fixed.log.scale$st <- as.factor(rownames(st.fixed.log.scale))
st.fixed.log.scale$row <- as.factor(rep("row1", nrow(st.fixed.log.scale)))

# remove NT (reference)
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


# plot top 40 lineages (by size)
lineage.top <- 40

lineage.temp <- as.data.frame(sort(table(dat.filt$SC_ST_clade), decreasing = T)) %>%
              rename(lineage = Var1)  %>%
              top_n(lineage.top)
lineage.filt <- lineage.temp$lineage
plot.dat <- dat.filt %>%
            filter(SC_ST_clade %in% lineage.filt)

# combinations of serotype and lineage for picking up the right estimates for summing up the coefficients
temp <- data.frame(comb = unique(plot.dat$SC_ST_clade_K))

temp$lineage <- NA
for(i in 1:length(temp$comb)){
  temp$lineage[i] <- ifelse(grepl("like", temp$comb[i], fixed = T) | grepl("KL", temp$comb[i], fixed = T), 
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-2)), collapse = "_"), 
                            paste(c(sapply(strsplit(temp$comb[i], "_"), head, lengths(strsplit(temp$comb[i], "_"))-1)), collapse = "_"))
  temp$serotype[i] <- ifelse(grepl("like", temp$comb[i], fixed = T) | grepl("KL", temp$comb[i], fixed = T),
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 2)), collapse = "_"), 
                             paste(c(sapply(strsplit(temp$comb[i], "_"), tail, 1)), collapse = "_"))
}

temp$comb.est <- NA

# remove all rows with NT
temp <- temp %>%
  filter(serotype != "Untypeable")

# fix some of the K-names in temp file:
temp$serotype[which(temp$comb == "SC20_ST59_KL99")] <- "KL99"
temp$lineage[which(temp$comb == "SC20_ST59_KL99")] <- "SC20_ST59"
temp$serotype[which(temp$comb == "SC6_ST69_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC6_ST69_KL7")] <- "SC6_ST69"
temp$serotype[which(temp$comb == "SC6_ST69_KL99")] <- "KL99"
temp$lineage[which(temp$comb == "SC6_ST69_KL99")] <- "SC6_ST69"
temp$serotype[which(temp$comb == "SC15_ST648")] <- "KL7"
temp$lineage[which(temp$comb == "SC15_ST648")] <- "SC15_ST648"
temp$serotype[which(temp$comb == "SC3_ST73_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC3_ST73_KL7")] <- "SC3_ST73"
temp$serotype[which(temp$comb == "SC16_ST38_KL99")] <- "KL99"
temp$lineage[which(temp$comb == "SC16_ST38_KL99")] <- "SC16_ST38"
temp$serotype[which(temp$comb == "SC14_ST405_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC14_ST405_KL7")] <- "SC14_ST405"
temp$serotype[which(temp$comb == "SC22_ST398_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC22_ST398_KL7")] <- "SC22_ST398"
temp$serotype[which(temp$comb == "SC16_ST38_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC16_ST38_KL7")] <- "SC16_ST38"
temp$serotype[which(temp$comb == "SC15_ST648_KL7")] <- "KL7"
temp$lineage[which(temp$comb == "SC15_ST648_KL7")] <- "SC15_ST648"


st.fixed.log.scale.ord$st <- as.character(st.fixed.log.scale.ord$st)

# sum the lineage + k estimates (fixed + random effect)
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

# add margins
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

# keep only lineages that have capsule(s) that passed the filtering
# remove SC_STs with no capsule information
temp2 <- heatmap.all %>%
  group_by(lineage) %>%
  mutate(n.lineage = n()) %>%
  mutate(n.mis.comb = sum(is.na(comb))) %>%
  ungroup()

heatmap.all2 <- temp2 %>%
  filter(n.lineage != n.mis.comb | lineage == " ")


# Heatmap

# Mistake in data:
# change KL99 --> KL13
heatmap.all2$serotype2 <- revalue(heatmap.all2$serotype, c("KL99" = "KL13"))

heat <- ggplot(heatmap.all2, aes(x = serotype2, y = lineage, fill = comb.est)) +
  geom_tile(color = "grey") +
  scale_fill_gradient2(low = "darkgreen",
                       mid = "#FFFFCC",
                       high = "#FF0000", na.value = 'white') +
  coord_fixed()+
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
ggsave("~/Results_thres5_tot20/heatmap.top40_lineages.thres5_tot20.pdf", heat,  width = 10, height = 15)



