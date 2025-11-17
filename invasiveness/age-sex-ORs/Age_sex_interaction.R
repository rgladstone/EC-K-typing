
######## Age + sex -interaction

# data
BSAC <- read.csv("~/data/BSAC_FULL_pp25_2003_2018_clades_241025.csv")
BB <-  read.csv("~/data/BB_SCandK_241025.csv")

# correct one clade typo
BB$ST131.clade[which(BB$ST131.clade == "C")] <- "C0"
BSAC$ST131.clade[which(BSAC$ST131.clade == "PP195_CCnew")] <- ""

# BB
bb <- BB %>%
  mutate(collection = "BB1.BB2") %>%
  mutate(infection.bin = 0) %>%
  mutate(clinical_manifest = "carriage") %>%
  select(lane_id, pp_CCname, ST131.clade, Designation, collection, infection.bin, clinical_manifest, Sex, Age, Phylogroup) %>%
  rename(lineage = pp_CCname) %>%
  rename(K = Designation, sex = Sex, age = Age, phylogroup = Phylogroup) %>%
  rename(lane = lane_id)

# BSAC
bsac <- BSAC %>%
  mutate(collection = "BSAC") %>%
  mutate(infection.bin = 1) %>%
  mutate(clinical_manifest = "infection") %>%
  rename(lineage = pp_CCname) %>%
  select(lane, lineage, ST131.clade, K, collection, infection.bin, clinical_manifest, Sex, Age, Phylogroup) %>%
  rename(sex = Sex, age = Age, phylogroup = Phylogroup)

dat <- rbind(bb, bsac)

dat <- dat %>% 
  mutate(cc_clade = ifelse(ST131.clade != "", paste(lineage, ST131.clade, sep = "_"), lineage)) %>%
  mutate(cc_clade_k = paste(cc_clade, K, sep = "_"))

dat$cc_clade <- as.factor(dat$cc_clade)
dat$K <- as.factor(dat$K)
dat$clinical_manifest <- as.factor(dat$clinical_manifest)
dat$sex <- ifelse(dat$sex == "Missing", NA, dat$sex) 
dat$sex <- factor(dat$sex, levels = c("M", "F"), labels = c("Male", "Female"))

# age groups
dat <- dat %>%
  mutate(age.groups = case_when(age == "Missing" ~ NA,
                                age == "60-69 yrs" ~ ">=60",
                                age == "70-79 yrs" ~ ">=60",
                                age == "80+ yrs" ~ ">=60",
                                .default = "<60"))


##### only BSAC (BB is mostly <1y)

dat.bsac <- dat %>%
  filter(collection == "BSAC")

bsac.st.k.comb <- as.data.frame(sort(table(dat.bsac$cc_clade_k), decreasing = T))

# data frame to save results
# the ones with at least 50 isolates
results <- data.frame(cc_k = c("CC95_K1", "CC73_K5", "CC73_K2", "CC69_K52", "CC131_C2_K5", "CC1193_K5", "CC73_K7", "CC131_C2_K100"),
                      OR=NA, CI_lower = NA, CI_upper = NA, p.value=NA)

######### CC95_K1, n=233

# males
m.cc95.k1 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC95_K1") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc95.k1 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC95_K1") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc95.k1 <- rbind(f.cc95.k1, m.cc95.k1)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc95.k1$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC95_K1")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC95_K1")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC95_K1")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC95_K1")] <- ft$p.value



######### CC73_K5,  n=172

# males
m.cc73.k5 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC73_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc73.k5 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC73_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc73.k5 <- rbind(f.cc73.k5, m.cc73.k5)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc73.k5$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC73_K5")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC73_K5")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC73_K5")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC73_K5")] <- ft$p.value


######### CC73_K2, n=88

# males
m.cc73.k2 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC73_K2") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc73.k2 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC73_K2") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc73.k2 <- rbind(f.cc73.k2, m.cc73.k2)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc73.k2$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft


results$OR[which(results$cc_k == "CC73_K2")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC73_K2")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC73_K2")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC73_K2")] <- ft$p.value

######### CC69_K52, n=58

# males
m.cc69.k52 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC69_K52") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc69.k52 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC69_K52") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc69.k52 <- rbind(f.cc69.k52, m.cc69.k52)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc69.k52$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC69_K52")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC69_K52")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC69_K52")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC69_K52")] <- ft$p.value

######### CC131_C2_K5, n=54

# males
m.cc131c2.k5 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC131_C2_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc131c2.k5 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC131_C2_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc131c2.k5 <- rbind(f.cc131c2.k5, m.cc131c2.k5)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc131c2.k5$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC131_C2_K5")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC131_C2_K5")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC131_C2_K5")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC131_C2_K5")] <- ft$p.value


######### CC1193_K5, n=53

# males
m.cc1193.k5 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC1193_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc1193.k5 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC1193_K5") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc1193.k5 <- rbind(f.cc1193.k5, m.cc1193.k5)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc1193.k5$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC1193_K5")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC1193_K5")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC1193_K5")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC1193_K5")] <- ft$p.value


######### CC73_K7, n=52

# males
m.cc73.k7 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC73_K7") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc73.k7 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC73_K7") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc73.k7 <- rbind(f.cc73.k7, m.cc73.k7)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc73.k7$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC73_K7")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC73_K7")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC73_K7")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC73_K7")] <- ft$p.value


######## CC131_C2_K100, n=50

# males
m.cc131c2.k100 <- dat.bsac %>%
  filter(sex == "Male") %>%
  filter(cc_clade_k == "CC131_C2_K100") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "male")

# females
f.cc131c2.k100 <- dat.bsac %>%
  filter(sex == "Female") %>%
  filter(cc_clade_k == "CC131_C2_K100") %>%
  group_by(age.groups) %>%
  summarise(n.age.groups = n()) %>%
  mutate(sex = "female")

cc131c2.k100 <- rbind(f.cc131c2.k100, m.cc131c2.k100)

# odds ratio males/females (old/young)
# p-value by fisher test
tab <- matrix(cc131c2.k100$n.age.groups, nrow = 2, byrow = T,
              dimnames = list(sex = c("female", "male"),
                              age  = c("<60", ">=60")))
tab
ft <- fisher.test(tab)
ft

results$OR[which(results$cc_k == "CC131_C2_K100")] <- ft$estimate
results$CI_lower[which(results$cc_k == "CC131_C2_K100")] <- ft$conf.int[1]
results$CI_upper[which(results$cc_k == "CC131_C2_K100")] <- ft$conf.int[2]
results$p.value[which(results$cc_k == "CC131_C2_K100")] <- ft$p.value

