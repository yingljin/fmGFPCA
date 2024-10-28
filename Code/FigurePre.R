
# Produce figures used for presentation
#### Set up ####
library(here)
library(tidyverse)
library(plotly)
library(RColorBrewer)
# https://rviews.rstudio.com/2020/12/14/plotting-surfaces-with-r/


#### MS Lesion data example ####
# getwd()
# list.files(here("RawData/"))
# load(here("RawData/df_scans.RData"))
load(here("RawData/df_comb.RData"))
head(df_scan)

df_exp <- df_scan %>%
  filter(sub == "sub-01001" & site == "NIH" & z==165)
class(df_exp$t1_intensity1_std)

df_exp %>% group_by(z) %>% summarise(n=n()) %>% arrange(desc(n))

RColorBrewer::brewer.pal(2, "Set2")

Mat1 <- df_exp %>% select(x, y ,t1_intensity1) %>% 
  pivot_wider(names_from = y, values_from = t1_intensity1) %>% 
  select(-x) %>% as.matrix()
Mat2 <- df_exp %>% select(x, y ,t1_intensity2) %>% 
  pivot_wider(names_from = y, values_from = t1_intensity2) %>% 
  select(-x) %>% as.matrix()

# plot 

display.brewer.pal(2, "Accent") 

plot_ly(df_exp, colors = "Accent", size = 0.5) %>%
  add_markers(x=~x, y=~y, z=~t1_intensity1, color = "Scan 1") %>%
  add_markers(x=~x, y=~y, z=~t1_intensity2, color = "Scan 2") %>%
  layout(scene = list(xaxis = list(title = ''),
                      yaxis = list(title = ''),
                      zaxis = list(title = '')),
         title = list(text = "T1 intensity"),
         legend = list(orientation = "h",   
                       xanchor = "center",
                       x = 0.5))
  
add_mesh(df_exp,  x = ~x, y=~y, z = ~t1_intensity1, col = "Scan 1")
  
head(df_exp)


#### Examples of gmFPCA  ####
# simulated data
load(here("Data/SimData/SimData_J200.RData"))

library(lme4)

data <- lme4data <- data_list_allM_J200[[1]]
rm(data_list_allM_J200)

bin_w=10; L=4
t_vec <- unique(data$t) # measurement grid
J <- max(as.numeric(unique(data$visit))) # number of visits per subject
K <- length(t_vec) # number of observations per visit

# add a function index
# data$ind <- rep(1:K, J)

# Step 1
print("Step 1: bin data")
## create an index frame for bins
df_bin <- data.frame(t = t_vec) %>% 
  mutate(tind = 1:K) %>% 
  mutate(bin = ceiling(tind/bin_w))

## bin midpoints
bin_mid <- df_bin %>% group_by(bin) %>% summarise(bin_mid = median(t)) 

## data
data <- data %>% 
  left_join(df_bin, by = "t")

# Step 2
print("Step 2: local GLMM")
df_s2 <- data %>% 
  filter(complete.cases(.)) %>% ## make it accommodate missing values
  group_by(bin) %>%
  group_map(~{
    fit_local <- glmer(Y~1+(1|id)+(1|id:visit), data = .x, 
                       family = binomial, nAGQ=0)
    eta_hat <- predict(fit_local, type = "link")
    .x$eta_hat <- eta_hat
    return(data.frame(.x))
  }, .keep = T) %>% bind_rows()

cols <- c(brewer.pal(3, "Set2"), "#000000") # define a color palette 
names(cols) <- c(1:3, "True")

tend <- df_bin %>% group_by(bin) %>% 
  summarize(t_end = round(max(t), 3))

df_s2 %>% filter(id %in% 1:4 & bin %in% 1:5) %>% 
  group_by(bin) %>%
  mutate(t_mid = median(t)) %>% 
  mutate(t_max = round(max(t), 2)) %>% 
  mutate_at(c("eta_hat"), function(x){exp(x)/{1+exp(x)}}) %>% #head()
  ggplot()+
  geom_point(aes(x=t, y=Y, col = visit), size = 0.5)+
  geom_point(aes(x=t_mid, y=eta_hat, col = visit))+
  geom_line(aes(x=t_mid, y=eta_hat, col = visit))+
  geom_vline(xintercept = tend$t_end[1:5], linetype = "dashed", linewidth = 0.3)+
  facet_wrap(~id, nrow = 1)+
  scale_color_manual(values = cols)+
  labs(col = "Series")
ggsave(filename = here("Images/bin_exp.png"),
       width = 14, height = 4, bg = "white")

# step 3
print("Step 4: re-evaluation")
df_s3 <- df_s2 %>% select(id, visit, bin, eta_hat) %>% 
  distinct(.) %>% 
  pivot_wider(., names_from = bin, values_from = eta_hat) 

## number of knots
library(refund)
nknot <- min(nrow(bin_mid)-4, 30)
fit_mfpca <- mfpca.face(as.matrix(df_s3 %>% select(-id, -visit)),
                        id = df_s3$id, visit = df_s3$visit, 
                        argvals = as.vector(bin_mid$bin_mid),
                        npc = L, knots = nknot)

# step 4
print("Step 4: re-evaluation")
source(here("Code/Functions/ReEvalEfuncs.R"))
## eigenfunctions:
reeval_efunc_l1 <- reeval_efunctions(
  knots = nrow(bin_mid)/2,
  argvals_bin = bin_mid$bin_mid,
  argvals = t_vec,
  efunctions = fit_mfpca$efunctions$level1,
  npc = ncol(fit_mfpca$efunctions$level1)
)

reeval_efunc_l2 <- reeval_efunctions(
  knots = nrow(bin_mid)/2,
  argvals_bin = bin_mid$bin_mid,
  argvals = t_vec,
  efunctions = fit_mfpca$efunctions$level2,
  npc = ncol(fit_mfpca$efunctions$level2)
)
colnames(reeval_efunc_l1) <- paste0("PC", 1:L)
colnames(reeval_efunc_l2) <- paste0("PC", 1:L)
head(reeval_efunc_l1)

# plot of PC funcitons
cols <- c(brewer.pal(4, "Set2")) # define a color palette 

p1 <- reeval_efunc_l1 %>% data.frame() %>% 
  mutate(t = df_bin$t) %>% 
  pivot_longer(1:4) %>%
  ggplot()+
  geom_line(aes(x=t, y=value, col=name))+
  scale_color_manual(values = cols)+
  theme_minimal()+
  labs(title = "Subject level PC functions", x="t", y="", color = "")
  
p2 <- reeval_efunc_l2 %>% data.frame() %>% 
  mutate(t = df_bin$t) %>% 
  pivot_longer(1:4) %>%
  ggplot()+
  geom_line(aes(x=t, y=value, col=name))+
  scale_color_manual(values = cols)+
  theme_minimal()+
  labs(title = "Subject-series level PC functions", x="t", y="", color = "")

ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
ggsave(filename = here("Images/gmfpca_fpc.png"),
       width = 12, height = 4, bg = "white")
  