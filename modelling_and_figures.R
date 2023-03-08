#Script for all manuscript figures and modelling of primary objectives
#1) Days to Asymptomatic status
#2) Days to Recovery

#libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(rstan)
library(bayesplot)
library(rethinking)
library(tidybayes)
library(ggrepel)
library(RColorBrewer)
library(gtsummary)
library(gt)
library(rstanarm)
library(tidybayes.rethinking)
library(gganimate)

#for speed with mcmc
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#in-house functions
z_score <- function(x) { 
  
  out <- (x - mean(x))/sd(x)
  
}

#load relevant files
#for days to asymptomatic status
m1_df <- read.csv('m1_df.csv',stringsAsFactors = FALSE)

#for days to recovery
m2_df <- read.csv('m2_df.csv',stringsAsFactors = FALSE)

#for group contrasts on symptoms (Table 2 modelling)
symp_df <- read.csv('symp_df.csv',stringsAsFactors = FALSE)
#change column names
colnames(symp_df)[c(5:7)] <- c("Total Symptoms",'Symptom Severity',
                               '% of Normal')


#####DAYS TO ASYMPTOMATIC STATUS MODEL#####
#list data
dat <- list(days = m1_df$days_injury,
            reach_norm = m1_df$reach_norm,
            group = m1_df$randomization_group)

#create stan model call with rethinking package
m1 <- ulam(
  alist(
    days|reach_norm==1 ~ exponential(lambda),
    days|reach_norm==0 ~ custom(exponential_lccdf(!Y|lambda)),
    lambda <- 1/mu,
    log(mu) <- a[group],
    a[group] ~ normal(3.15,0.6)
  ),data = dat, chains = 4,iter = 20000, cores=4)

precis(m1,2) 
#extract posterior
pm1 <- data.frame(extract.samples(m1))

#create contrast
pm1$days1 <- exp(pm1$a.1)
pm1$days2 <- exp(pm1$a.2)
pm1$diff <- pm1$days1 - pm1$days2

#summary
precis(pm1,2,prob = 0.90)
mean(pm1$diff>0)
summary(pm1$diff)

#priors
prior_m1 <- data.frame(extract.prior(m1))
precis(prior_m1,2)
prior_m1$days1 <- exp(prior_m1$a.1)
prior_m1$days2 <- exp(prior_m1$a.2)
prior_m1$diff <- prior_m1$days1 - prior_m1$days2
precis(prior_m1,2,prob = 0.90)

#####DAYS TO RECOVERY MODEL#####
dat2 <- list(rtp_days = m2_df$days,
             reach_rtp =m2_df$reach_rtp,
             group = m2_df$randomization_group)

str(dat2)

#create stan model call with rethinking package
m2 <- ulam(
  alist(
    rtp_days|reach_rtp==1 ~ exponential(lambda),
    rtp_days|reach_rtp==0 ~ custom(exponential_lccdf(!Y|lambda)),
    lambda <- 1/mu,
    log(mu) <- a[group],
    a[group] ~ normal(3.15,0.6)
  ),data = dat2, chains = 4,iter = 20000, cores=4)

precis(m2,2) 

#extract posterior
pm2 <- data.frame(extract.samples(m2))

#create contrast
pm2$days1 <- exp(pm2$a.1)
pm2$days2 <- exp(pm2$a.2)
pm2$diff <- pm2$days1 - pm2$days2

#summary
precis(pm2,2, prob = 0.9)
mean(pm2$diff>0)

####PLOT PREPPING#######
#m1
#Counterfactual data for m1
days <- data.frame(rep(c(1:28),2))
group <- rep(c(1,2),each = 28)
fake_dat <- cbind.data.frame(days,group)
colnames(fake_dat) <- c("days",'randomization_group')

fake_m1 <- data.frame(link(m1,dat = fake_dat))
#separate lambda's for each day value, by group
df1 <- fake_m1[c(1)] #all the columns of lambda are the same
df1$group <- "UCEP"
df2 <- fake_m1[c(29)]
df2$group <- "SAEP"

#code to create draws of all possibilities of the CCDF for days 1:28 in Group 1
l1 <- lapply(1:28, function(i) list(x = exp(-df1$lambda.1*i),
                                    y = rep(i,length(df1$lambda.1))))
df1<- do.call(rbind.data.frame, l1)

#same in group 2
l2 <- lapply(1:28, function(i) list(x = exp(-df2$lambda.29*i),
                                    y = rep(i,length(df2$lambda.29))))
df2<- do.call(rbind.data.frame, l2)

df1$Treatment <- "UCEP"
df2$Treatment <- "SAEP"

#bind dfs
df_pprob <- rbind.data.frame(df1,df2)
colnames(df_pprob)
colnames(df_pprob) <- c("Proportion Remaining",'Days','Treatment')
df_pprob$`Proportion Remaining` <- df_pprob$`Proportion Remaining`*100

#for prior simulation
#create lambda
prior_m1$lambda1 <- 1/exp(prior_m1$a.1)
prior_m1$lambda2 <- 1/exp(prior_m1$a.2)

#function to create CCDF draws
lprior1 <- lapply(1:28, function(i) list(x = exp(-prior_m1$lambda1*i),
                                         y = rep(i,length(prior_m1$lambda1))))
lprior1_df<- do.call(rbind.data.frame, lprior1)
lprior2 <- lapply(1:28, function(i) list(x = exp(-prior_m1$lambda2*i),
                                         y = rep(i,length(prior_m1$lambda2))))
lprior2_df<- do.call(rbind.data.frame, lprior2)

#create group variable
lprior1_df$Treatment <- "UCEP"
lprior2_df$Treatment <- "SAEP"

#bind prior df2
df_prior_prob <- rbind.data.frame(lprior1_df,lprior2_df)
colnames(df_prior_prob)
colnames(df_prior_prob) <- c("Proportion Remaining",'Days','Treatment')
df_prior_prob$`Proportion Remaining` <- df_prior_prob$`Proportion Remaining`*100

#real data for m1
#split to give each treatment their own CCDF
m1_df1 <- m1_df[m1_df$randomization_group==1,]
m1_df2 <- m1_df[m1_df$randomization_group==2,]

#give each its own CCDF
#UCEP
m1_df1 <- m1_df1[order(m1_df1$days_injury),]
m1_df1$ccdf <- cumsum(m1_df1$reach_norm)
m1_df1$ccdf <- m1_df1$ccdf/19
m1_df1$ccdf <- 1-m1_df1$ccdf

#SAEP
m1_df2 <- m1_df2[order(m1_df2$days_injury),]
m1_df2$ccdf <- cumsum(m1_df2$reach_norm)
m1_df2$ccdf <- m1_df2$ccdf/19
m1_df2$ccdf <- 1-m1_df2$ccdf
m1_df3 <- rbind.data.frame(m1_df1,m1_df2)

#format
colnames(m1_df3)[c(2,3,5)] <- c("Treatment","Days","Proportion Remaining")
m1_df3$Treatment <- factor(m1_df3$Treatment)
m1_df3$Treatment <- ifelse(m1_df3$Treatment==1,"UCEP",'SAEP')
m1_df3$`Proportion Remaining` <- m1_df3$`Proportion Remaining`*100

#Counterfactual data for m2
days <- data.frame(rep(c(1:150),2))
group <- rep(c(1,2),each = 150)
fake_dat2 <- cbind.data.frame(days,group)
colnames(fake_dat2) <- c("Final_Days",'randomization_group')

fake_m2 <- data.frame(link(m2,dat = fake_dat2))
#separate lambda's for each day value, by group
df1 <- fake_m2[c(1)] #all the columns of lambda are the same
df1$group <- "UCEP"
df2 <- fake_m2[c(151)]
df2$group <- "SAEP"

# code  to create draws of all possibilities of the CCDF for days 1:28 in Group 1
l1 <- lapply(1:150, function(i) list(x = exp(-df1$lambda.1*i),
                                     y = rep(i,length(df1$lambda.1))))
df1<- do.call(rbind.data.frame, l1)

#same in group 2
l2 <- lapply(1:150, function(i) list(x = exp(-df2$lambda.151*i),
                                     y = rep(i,length(df2$lambda.151))))
df2<- do.call(rbind.data.frame, l2)

df1$Treatment <- "UCEP"
df2$Treatment <- "SAEP"

#bind dfs
df_pprob2 <- rbind.data.frame(df1,df2)
colnames(df_pprob2)
colnames(df_pprob2) <- c("Proportion Remaining",'Days','Treatment')
df_pprob2$`Proportion Remaining` <- df_pprob2$`Proportion Remaining`*100

#real data for m2
#split to give each treatment their own CCDF
#usual care
m2_df1 <- m2_df[m2_df$randomization_group==1,]
m2_df2 <- m2_df[m2_df$randomization_group==2,]

#Give each their own CCDF
#UCEP
m2_df1 <- m2_df1[order(m2_df1$days),]
m2_df1$ccdf <- cumsum(m2_df1$randomization_group)
m2_df1$ccdf <- m2_df1$ccdf/19
m2_df1$ccdf <- 1-m2_df1$ccdf

#SAEP
m2_df2 <- m2_df2[order(m2_df2$days),]
m2_df2$ccdf <- 1:19
m2_df2$ccdf <- m2_df2$ccdf/19
m2_df2$ccdf <- 1-m2_df2$ccdf
m2_df3 <- rbind.data.frame(m2_df1,m2_df2)

#format
colnames(m2_df3)[c(3,2,5)] <- c("Treatment","Days","Proportion Remaining")
m2_df3$Treatment <- factor(m2_df3$Treatment)
m2_df3$Treatment <- ifelse(m2_df3$Treatment==1,'UCEP','SAEP')
m2_df3$`Proportion Remaining` <- m2_df3$`Proportion Remaining`*100

##########PLOTS###########

#m1 plot of symptom resolution
theme_set(theme_tidybayes() + panel_border())
m1_df_plot_ppc = df_pprob%>%
  ggplot(aes(x = Days, y = `Proportion Remaining`, fill = Treatment)) +
  stat_lineribbon(aes(fill = Treatment), .width = c(.50,.75,.90),alpha = 1/5) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent")+
  xlab("Days to symptom resolution") +
  ylab('Proportion remaining (%)')+
  xlim(0,30)+ 
  theme(axis.title.x = element_text(face='bold'),
        axis.title.y = element_text(face='bold'),
        legend.title = element_text(face = 'bold'))+
  geom_point(data=m1_df3,aes(x = Days,
                             y = `Proportion Remaining`,
                             color = Treatment))
m1_df_plot_ppc

#histogram plot m1
#prep
colnames(pm1)
ucep <- cbind.data.frame(pm1$days1,rep('UCEP',length(pm1)))
saep <- cbind.data.frame(pm1$days2, rep('SAEP',length(pm1)))
diff <- cbind.data.frame(pm1$diff,rep('Difference (Contrast)',length(pm1)))
colnames(ucep) <- c("value",'condition')
colnames(saep) <- colnames(ucep)
colnames(diff) <- colnames(saep)
histdf <- rbind.data.frame(ucep,saep,diff)

hist1 <- histdf %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank(),
                                            axis.title.x = element_text(face='bold'),
                                            axis.text.y = element_text(face='bold'))+
  scale_fill_manual(values=c('orange','lightgreen','purple'))+
  xlab('Days to symptom resolution')
hist1 

#combine
fig2_combined <- plot_grid(hist1,m1_df_plot_ppc,labels="AUTO",ncol=2)

#save
ggsave('fig2.jpg',fig2_combined,dpi = 600, width =12)

#prior simulation plot for Fig 1 (and for both m1 and m2)
prior_sim_line <- df_prior_prob %>%
  ggplot(aes(x = Days, y = `Proportion Remaining`, Fill = Treatment)) +
  xlab('Days to symptom resolution') + ylab("Proportion remaining (%)")+
  stat_lineribbon(aes(fill=Treatment),.width = c(.50,.75,.90),alpha = 1/5) +
  scale_color_brewer(palette = "Dark2") 
prior_sim_line

#prior histograms
#prep
colnames(prior_m1)
ucep <- cbind.data.frame(prior_m1$days1,rep('UCEP',length(prior_m1$days1)))
saep <- cbind.data.frame(prior_m1$days2, rep('SAEP',length(prior_m1$days2)))
diff <- cbind.data.frame(prior_m1$diff,rep('Difference (Contrast)',length(prior_m1$diff)))
colnames(ucep) <- c("value",'condition')
colnames(saep) <- colnames(ucep)
colnames(diff) <- colnames(saep)
prior_hist_df <- rbind.data.frame(ucep,saep,diff)

prior_sim_hist <- prior_hist_df %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank())+
  xlab('Days to symptom resolution')
prior_sim_hist

#combine prior sims for supplementary Figure 2
prior_sims_plot <- plot_grid(prior_sim_hist,prior_sim_line,
                labels="AUTO",rel_widths = c(1, 1.3), ncol = 2, nrow = 1)

#save
ggsave('s2_fig.jpg',prior_sims_plot,width = 12)


#m2 plot for rtp
m2_df_plot_ppc = df_pprob2%>%
  ggplot(aes(x = Days, y = `Proportion Remaining`, fill = Treatment)) +
  stat_lineribbon(aes(fill = Treatment), .width = c(.50,.75,.90),alpha = 1/5) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent")+
  xlab("Days to medical clearance") +
  ylab('Proportion remaining (%)')+
  xlim(0,150)+ 
  theme(axis.title.x = element_text(face='bold'),
        axis.title.y = element_text(face='bold'),
        legend.title = element_text(face = 'bold'))+
  geom_point(data=m2_df3,aes(x = Days,
                                   y = `Proportion Remaining`,
                                   color = Treatment))
m2_df_plot_ppc

#prep histogram for m2
colnames(pm2)
ucep <- cbind.data.frame(pm2$days1,rep('Usual Care',length(pm2)))
saep <- cbind.data.frame(pm2$days2, rep('Exercise',length(pm2)))
diff <- cbind.data.frame(pm2$diff,rep('Difference (Contrast)',length(pm2)))
colnames(ucep) <- c("value",'condition')
colnames(saep) <- colnames(ucep)
colnames(diff) <- colnames(saep)
histdf2 <- rbind.data.frame(ucep,saep,diff)
precis(pm2)

#plot hist for m2
hist2 <- histdf2 %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank(),
                                            axis.title.x = element_text(face='bold'),
                                            axis.text.y = element_text(face='bold'))+
  scale_fill_manual(values=c('orange','lightgreen','purple'))+
  xlab('Days to medical clearance')
hist2

fig_m2_combined <- plot_grid(hist2,m2_df_plot_ppc,labels="AUTO",ncol=2)
fig_m2_combined

#save
ggsave('s1_fig.jpg',fig_m2_combined,dpi = 600, width =12)

#MCMC CHECKS 
#m1
m1_stan <- m1@stanfit
color_scheme_set("mix-blue-red")
m1_trace <-mcmc_trace(m1_stan, pars = c('a[1]','a[2]'),
                      facet_args = list(ncol = 1, strip.position = "left")) +
  xlab('Posterior Draw') + ylab('Coefficient Value')

#Gelman convergence checks
m1_rhat <- rhat(m1_stan)
mcmc_rhat(m1_rhat)

#m2
m2_stan <- m2@stanfit
color_scheme_set("mix-blue-red")

m2_trace <-mcmc_trace(m2_stan, pars = c('a[1]','a[2]'),
                      facet_args = list(ncol = 1, strip.position = "left"))
m2_trace

trace_plot <- plot_grid(m1_trace,m2_trace,labels = "AUTO",ncol=2)
trace_plot

#save
ggsave('sup_fig4.jpg',trace_plot,dpi=600,width = 12)

#Data for the results section from m1 (time to asymptomatic status)
#1-exp(-lambda*days) #cumulative prob dist.
ucep_day_28 <- 1-(exp(-(1/exp(pm1$a.1))*28)) 
saep_day_28 <- 1-(exp(-(1/exp(pm1$a.2))*28))
ucep_day_14 <- 1-(exp(-(1/exp(pm1$a.1))*14)) 
saep_day_14 <- 1-(exp(-(1/exp(pm1$a.2))*14) )
counter_factuals <- cbind.data.frame(ucep_day_14,saep_day_14,ucep_day_28,
                                     saep_day_28)
counter_factuals$diff_28 <- counter_factuals$saep_day_28 - counter_factuals$ucep_day_28
precis(counter_factuals,2,prob = 0.9)
mean(counter_factuals$diff_28>0)

###Fig 3###
#prep
dot_df <- symp_df

facet_labs <- list('0' = 'Enrolment',
                   '1' = 'Assessment 1',
                   '2' = 'Assessment 2',
                   '3' = 'Assessment 3',
                   '4' = 'Assessment 4')

facet_labeller <- function(variable,value){
  return(facet_labs[value])}

#make plot
dot_plot <- ggplot(data = dot_df,
                   aes(x = factor(randomization_group), 
                       y = `Symptom Severity`,fill = factor(randomization_group))) +
  geom_dotplot(binaxis='y', stackdir='center',position = 'dodge',binwidth = 2,
               dotsize = 1.4) + 
  geom_violin(trim = FALSE,alpha = 0.1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5)+
  facet_wrap(~factor(days_bin),ncol = 5,labeller = facet_labeller)+
  theme_bw()+
  theme(strip.text = element_text(face = 'bold',size = 12),
        strip.background = element_rect(fill = "white",color = 'white'),
        panel.spacing = unit(2, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face='bold'),
        axis.text.y = element_text(face='bold'),
        axis.title.y=element_text(face='bold'),
        legend.position = 'none')+
  scale_fill_manual(values = c('green','purple'))
dot_plot

#save
ggsave('Fig3.jpg',dot_plot,dpi = 600, width = 12)

###modelling for symptom comparisons (from Table 2)###
#SYMPTOM SEVERITY
#unpooled model of symptoms over time
symptoms_model <- symp_df[c(1,4,6,8)]
str(symptoms_model) 
symptoms_model$ids <- as.factor(symptoms_model$ids)
symptoms_model$ids<- as.integer(symptoms_model$ids)
#have to remove enrolment, they aren't part of a group yet, and the estimates 
#are biased as a result.
symptoms_model <- symptoms_model[!(symptoms_model$redcap_repeat_instance=='Enrolment'),]
symptoms_model$redcap_repeat_instance <-
  factor(symptoms_model$redcap_repeat_instance,levels = 
           c("Assessment 1 (MED = 6.5 Days)",
             'Assessment 2 (MED = 13 Days)','Assessment 3 (MED = 20 Days)',
             'Assessment 4 (MED = 28 Days)'))
symptoms_model$redcap_repeat_instance <-
  as.integer(symptoms_model$redcap_repeat_instance)
symptoms_model$randomization_group <- 
  ifelse(symptoms_model$randomization_group=='UCEP',1,2)

#checking
symptoms_model <- na.omit(symptoms_model)
by(symptoms_model$`Symptom Severity`[symptoms_model$randomization_group==2],
   symptoms_model$redcap_repeat_instance[symptoms_model$randomization_group==2],median)

#sqrt transform
dat_symp <- list(subj = symptoms_model$ids,
                 time = as.integer(symptoms_model$redcap_repeat_instance),
                 group = as.integer(symptoms_model$randomization_group),
                 symptoms = as.numeric(sqrt(symptoms_model$`Symptom Severity`)))
str(dat_symp)

#model symptom values
m_symp <- ulam(
  alist(
    symptoms ~ dnorm(mu,sig),
    mu <- a1[group]+a2[time],
    a1[group] ~ dnorm(2,1),
    a2[time] ~ dnorm(3,1),
    sig ~ dexp(1)
  ), data = dat_symp, iter = 20000,  chains = 4, log_lik = TRUE  )

precis(m_symp,2)

#extract posterior
pm_msymp <- data.frame(extract.samples(m_symp))

#recover values
pm_msymp$uc_w1 <- pm_msymp$a1.1+pm_msymp$a2.1
pm_msymp$uc_w2 <- pm_msymp$a1.1+pm_msymp$a2.2
pm_msymp$uc_w3 <- pm_msymp$a1.1+pm_msymp$a2.3
pm_msymp$uc_w4 <- pm_msymp$a1.1+pm_msymp$a2.4
pm_msymp$ex_w1 <- pm_msymp$a1.2+pm_msymp$a2.1
pm_msymp$ex_w2 <- pm_msymp$a1.2+pm_msymp$a2.2
pm_msymp$ex_w3 <- pm_msymp$a1.2+pm_msymp$a2.3
pm_msymp$ex_w4 <- pm_msymp$a1.2+pm_msymp$a2.4
pm_msymp[c(8:15)] <- data.frame(apply(pm_msymp[c(8:15)],2,function(x) x^2))
pm_msymp$diff_w1 <- pm_msymp$uc_w1-pm_msymp$ex_w1
pm_msymp$diff_w2 <- pm_msymp$uc_w2-pm_msymp$ex_w2
pm_msymp$diff_w3 <- pm_msymp$uc_w3-pm_msymp$ex_w3
pm_msymp$diff_w4 <- pm_msymp$uc_w4-pm_msymp$ex_w4
precis(pm_msymp,2,prob=0.9)
mean(pm_msymp$diff_w3>0)

#priors
prior_msymp <- data.frame(extract.prior(m_symp))
precis(prior_msymp,2,prob=0.9)
prior_msymp$uc_w1 <- prior_msymp$a1.1+prior_msymp$a2.1
prior_msymp$uc_w2 <- prior_msymp$a1.1+prior_msymp$a2.2
prior_msymp$uc_w3 <- prior_msymp$a1.1+prior_msymp$a2.3
prior_msymp$uc_w4 <- prior_msymp$a1.1+prior_msymp$a2.4
prior_msymp$ex_w1 <- prior_msymp$a1.2+prior_msymp$a2.1
prior_msymp$ex_w2 <- prior_msymp$a1.2+prior_msymp$a2.2
prior_msymp$ex_w3 <- prior_msymp$a1.2+prior_msymp$a2.3
prior_msymp$ex_w4 <- prior_msymp$a1.2+prior_msymp$a2.4
prior_msymp[c(8:15)] <- data.frame(apply(prior_msymp[c(8:15)],2,function(x) x^2))
prior_msymp$diff_w1 <- prior_msymp$uc_w1-prior_msymp$ex_w1
prior_msymp$diff_w2 <- prior_msymp$uc_w2-prior_msymp$ex_w2
prior_msymp$diff_w3 <- prior_msymp$uc_w3-prior_msymp$ex_w3
prior_msymp$diff_w4 <- prior_msymp$uc_w4-prior_msymp$ex_w4
precis(prior_msymp,2,prob=0.9)

#Pooled model to check whether there is information gain through a 
#hierarchical structure
m_symp_p <- ulam(
  alist(
    symptoms ~ dnorm(mu,sig),
    mu <- a1[group]+mu_time + z[time]*sig_time,
    a1[group] ~ dnorm(0,0.5),
    z[time] ~ dnorm(0,0.5),
    mu_time ~ dnorm(0,0.5),
    sig ~ dexp(1),
    sig_time ~ dexp(1),
    gq > vector[time]:a2 <<- mu_time + z*sig_time
  ), data = dat_symp, iter = 20000,  chains = 4, log_lik = TRUE  )

precis(m_symp_p,2,prob = 0.9)

#extract posterior
pm_pooled <- data.frame(extract.samples(m_symp_p))

#recover values
pm_pooled$uc_w1 <- pm_pooled$a1.1+pm_pooled$a2.1
pm_pooled$uc_w2 <- pm_pooled$a1.1+pm_pooled$a2.2
pm_pooled$uc_w3 <- pm_pooled$a1.1+pm_pooled$a2.3
pm_pooled$uc_w4 <- pm_pooled$a1.1+pm_pooled$a2.4
pm_pooled$ex_w1 <- pm_pooled$a1.2+pm_pooled$a2.1
pm_pooled$ex_w2 <- pm_pooled$a1.2+pm_pooled$a2.2
pm_pooled$ex_w3 <- pm_pooled$a1.2+pm_pooled$a2.3
pm_pooled$ex_w4 <- pm_pooled$a1.2+pm_pooled$a2.4
pm_pooled[c(14:21)] <- 
  data.frame(apply(pm_pooled[c(14:21)],2,function(x)  x^2))
pm_pooled$diff_w1 <- pm_pooled$uc_w1-pm_pooled$ex_w1
pm_pooled$diff_w2 <- pm_pooled$uc_w2-pm_pooled$ex_w2
pm_pooled$diff_w3 <- pm_pooled$uc_w3-pm_pooled$ex_w3
pm_pooled$diff_w4 <- pm_pooled$uc_w4-pm_pooled$ex_w4
precis(pm_pooled,2,prob=0.9)

compare(m_symp,m_symp_p) #pooled model offers no benefit

#TOTAL SYMPTOMS
#unpooled model of total symptoms over time
symptom_t_model <- symp_df[c(1,4,5,8)]
str(symptom_t_model) 
symptom_t_model$ids <- as.factor(symptom_t_model$ids)
symptom_t_model$ids<- as.integer(symptom_t_model$ids)
#have to remove enrolment, they aren't part of a group yet, and the estimates 
#are biased as a result.
symptom_t_model <- symptom_t_model[!(symptom_t_model$redcap_repeat_instance=='Enrolment'),]
symptom_t_model$redcap_repeat_instance <-
  factor(symptom_t_model$redcap_repeat_instance,levels = 
           c("Assessment 1 (MED = 6.5 Days)",
             'Assessment 2 (MED = 13 Days)','Assessment 3 (MED = 20 Days)',
             'Assessment 4 (MED = 28 Days)'))
symptom_t_model$redcap_repeat_instance <-
  as.integer(symptom_t_model$redcap_repeat_instance)
symptom_t_model$randomization_group <- 
  ifelse(symptom_t_model$randomization_group=='UCEP',1,2)

#checking
symptom_t_model <- na.omit(symptom_t_model)
by(symptom_t_model$`Total Symptoms`[symptom_t_model$randomization_group==1],
   symptom_t_model$redcap_repeat_instance[symptom_t_model$randomization_group==1],median)

#sqrt transform
dat_tot_symp <- list(subj = symptom_t_model$ids,
                     time = as.integer(symptom_t_model$redcap_repeat_instance),
                     group = as.integer(symptom_t_model$randomization_group),
                     tot_symptoms = as.numeric(sqrt(symptom_t_model$`Total Symptoms`)))
str(dat_tot_symp)

m_tot_symp <- ulam(
  alist(
    tot_symptoms ~ dnorm(mu,sig),
    mu <- a1[group]+a2[time],
    a1[group] ~ dnorm(1.25,1),
    a2[time] ~ dnorm(2.25,1),
    sig ~ dexp(1)
  ), data = dat_tot_symp, iter = 20000,  chains = 4, log_lik = TRUE  )

precis(m_tot_symp,2,prob = 0.9)
pm_m_tot_symp <- data.frame(extract.samples(m_tot_symp))

#recover values
pm_m_tot_symp$uc_w1 <- pm_m_tot_symp$a1.1+pm_m_tot_symp$a2.1
pm_m_tot_symp$uc_w2 <- pm_m_tot_symp$a1.1+pm_m_tot_symp$a2.2
pm_m_tot_symp$uc_w3 <- pm_m_tot_symp$a1.1+pm_m_tot_symp$a2.3
pm_m_tot_symp$uc_w4 <- pm_m_tot_symp$a1.1+pm_m_tot_symp$a2.4
pm_m_tot_symp$ex_w1 <- pm_m_tot_symp$a1.2+pm_m_tot_symp$a2.1
pm_m_tot_symp$ex_w2 <- pm_m_tot_symp$a1.2+pm_m_tot_symp$a2.2
pm_m_tot_symp$ex_w3 <- pm_m_tot_symp$a1.2+pm_m_tot_symp$a2.3
pm_m_tot_symp$ex_w4 <- pm_m_tot_symp$a1.2+pm_m_tot_symp$a2.4
pm_m_tot_symp[c(8:15)] <- data.frame(apply(pm_m_tot_symp[c(8:15)],2,function(x) x^2))
pm_m_tot_symp$diff_w1 <- pm_m_tot_symp$uc_w1-pm_m_tot_symp$ex_w1
pm_m_tot_symp$diff_w2 <- pm_m_tot_symp$uc_w2-pm_m_tot_symp$ex_w2
pm_m_tot_symp$diff_w3 <- pm_m_tot_symp$uc_w3-pm_m_tot_symp$ex_w3
pm_m_tot_symp$diff_w4 <- pm_m_tot_symp$uc_w4-pm_m_tot_symp$ex_w4
precis(pm_m_tot_symp,2,prob=0.9)
mean(pm_m_tot_symp$diff_w4>0)

#priors
prior_m_tot_symp <- data.frame(extract.prior(m_tot_symp))
precis(prior_m_tot_symp,2,prob=0.90)
prior_m_tot_symp$uc_w1 <- prior_m_tot_symp$a1.1+prior_m_tot_symp$a2.1
prior_m_tot_symp$uc_w2 <- prior_m_tot_symp$a1.1+prior_m_tot_symp$a2.2
prior_m_tot_symp$uc_w3 <- prior_m_tot_symp$a1.1+prior_m_tot_symp$a2.3
prior_m_tot_symp$uc_w4 <- prior_m_tot_symp$a1.1+prior_m_tot_symp$a2.4
prior_m_tot_symp$ex_w1 <- prior_m_tot_symp$a1.2+prior_m_tot_symp$a2.1
prior_m_tot_symp$ex_w2 <- prior_m_tot_symp$a1.2+prior_m_tot_symp$a2.2
prior_m_tot_symp$ex_w3 <- prior_m_tot_symp$a1.2+prior_m_tot_symp$a2.3
prior_m_tot_symp$ex_w4 <- prior_m_tot_symp$a1.2+prior_m_tot_symp$a2.4
prior_m_tot_symp[c(8:15)] <- data.frame(apply(prior_m_tot_symp[c(8:15)],2,function(x) x^2))
prior_m_tot_symp$diff_w1 <- prior_m_tot_symp$uc_w1-prior_m_tot_symp$ex_w1
prior_m_tot_symp$diff_w2 <- prior_m_tot_symp$uc_w2-prior_m_tot_symp$ex_w2
prior_m_tot_symp$diff_w3 <- prior_m_tot_symp$uc_w3-prior_m_tot_symp$ex_w3
prior_m_tot_symp$diff_w4 <- prior_m_tot_symp$uc_w4-prior_m_tot_symp$ex_w4
precis(prior_m_tot_symp,2,prob=0.90)

#FEEL
#unpooled model of % feel over time
feel_model <- symp_df[c(1,4,7,8)]
str(feel_model) 
feel_model$ids <- as.factor(feel_model$ids)
feel_model$ids<- as.integer(feel_model$ids)
#have to remove enrolment, they aren't part of a group yet, and the estimates 
#are biased as a result.
feel_model <- feel_model[!(feel_model$redcap_repeat_instance=='Enrolment'),]
feel_model$redcap_repeat_instance <-
  factor(feel_model$redcap_repeat_instance,levels = 
           c("Assessment 1 (MED = 6.5 Days)",
             'Assessment 2 (MED = 13 Days)','Assessment 3 (MED = 20 Days)',
             'Assessment 4 (MED = 28 Days)'))
feel_model$redcap_repeat_instance <-
  as.integer(feel_model$redcap_repeat_instance)
feel_model$randomization_group <- 
  ifelse(feel_model$randomization_group=='UCEP',1,2)

#checking
feel_model <- na.omit(feel_model)
by(feel_model$`% of Normal`[feel_model$randomization_group==2],
   feel_model$redcap_repeat_instance[feel_model$randomization_group==2],median)

#model with z transform
dat_feel <- list(subj = feel_model$ids,
                 time = as.integer(feel_model$redcap_repeat_instance),
                 group = as.integer(feel_model$randomization_group),
                 feel = as.numeric(z_score(feel_model$`% of Normal`)))
str(dat_feel)

m_feel<- ulam(
  alist(
    feel ~ dnorm(mu,sig),
    mu <- a1[group]+a2[time],
    a1[group] ~ dnorm(0,0.5),
    a2[time] ~ dnorm(0,0.5),
    sig ~ dexp(1)
  ), data = dat_feel, iter = 20000,  chains = 4, log_lik = TRUE  )

precis(m_feel,2,prob = 0.9)
pm_feel <- data.frame(extract.samples(m_feel))

#recover values
pm_feel$uc_w1 <- pm_feel$a1.1+pm_feel$a2.1
pm_feel$uc_w2 <- pm_feel$a1.1+pm_feel$a2.2
pm_feel$uc_w3 <- pm_feel$a1.1+pm_feel$a2.3
pm_feel$uc_w4 <- pm_feel$a1.1+pm_feel$a2.4
pm_feel$ex_w1 <- pm_feel$a1.2+pm_feel$a2.1
pm_feel$ex_w2 <- pm_feel$a1.2+pm_feel$a2.2
pm_feel$ex_w3 <- pm_feel$a1.2+pm_feel$a2.3
pm_feel$ex_w4 <- pm_feel$a1.2+pm_feel$a2.4
pm_feel[c(8:15)] <- data.frame(apply(pm_feel[c(8:15)],2,
                                     function(x) (x*sd(feel_model$`% of Normal`))+
                                       mean(feel_model$`% of Normal`)))
pm_feel$diff_w1 <- pm_feel$uc_w1-pm_feel$ex_w1
pm_feel$diff_w2 <- pm_feel$uc_w2-pm_feel$ex_w2
pm_feel$diff_w3 <- pm_feel$uc_w3-pm_feel$ex_w3
pm_feel$diff_w4 <- pm_feel$uc_w4-pm_feel$ex_w4
precis(pm_feel,2,prob=0.9)
mean(pm_feel$diff_w1<0)
mean(pm_feel$diff_w2<0)
mean(pm_feel$diff_w3<0)
mean(pm_feel$diff_w4<0)

####Prior predictive chescks found in the supplement.
#priors
prior_feel <- data.frame(extract.prior(m_feel))
precis(prior_feel,2,prob=0.9)
prior_feel$uc_w1 <- prior_feel$a1.1+prior_feel$a2.1
prior_feel$uc_w2 <- prior_feel$a1.1+prior_feel$a2.2
prior_feel$uc_w3 <- prior_feel$a1.1+prior_feel$a2.3
prior_feel$uc_w4 <- prior_feel$a1.1+prior_feel$a2.4
prior_feel$ex_w1 <- prior_feel$a1.2+prior_feel$a2.1
prior_feel$ex_w2 <- prior_feel$a1.2+prior_feel$a2.2
prior_feel$ex_w3 <- prior_feel$a1.2+prior_feel$a2.3
prior_feel$ex_w4 <- prior_feel$a1.2+prior_feel$a2.4
prior_feel[c(8:15)] <- data.frame(apply(prior_feel[c(8:15)],2,
                                        function(x) (x*sd(feel_model$`% of Normal`))+
                                          mean(feel_model$`% of Normal`)))
prior_feel$diff_w1 <- prior_feel$uc_w1-prior_feel$ex_w1
prior_feel$diff_w2 <- prior_feel$uc_w2-prior_feel$ex_w2
prior_feel$diff_w3 <- prior_feel$uc_w3-prior_feel$ex_w3
prior_feel$diff_w4 <- prior_feel$uc_w4-prior_feel$ex_w4
precis(prior_feel,2,prob=0.9)

#SYMPTOM SEVERITY
usual1 <- cbind.data.frame(prior_msymp$uc_w1,rep('UCEP A1',length(prior_msymp$uc_w1)))
usual2 <- cbind.data.frame(prior_msymp$uc_w2,rep('UCEP A2',length(prior_msymp$uc_w2)))
usual3 <- cbind.data.frame(prior_msymp$uc_w3,rep('UCEP A3',length(prior_msymp$uc_w3)))
usual4 <- cbind.data.frame(prior_msymp$uc_w4,rep('UCEP A4',length(prior_msymp$uc_w4)))
exercise1 <- 
  cbind.data.frame(prior_msymp$ex_w1,rep('SAEP A1',length(prior_msymp$ex_w1)))
exercise2 <- 
  cbind.data.frame(prior_msymp$ex_w2,rep('SAEP A2',length(prior_msymp$ex_w2)))
exercise3 <- 
  cbind.data.frame(prior_msymp$ex_w3,rep('SAEP A3',length(prior_msymp$ex_w3)))
exercise4 <- 
  cbind.data.frame(prior_msymp$ex_w4,rep('SAEP A4',length(prior_msymp$ex_w4)))
diff1 <- 
  cbind.data.frame(prior_msymp$diff_w1,rep('Contrast A1',length(prior_msymp$diff_w1)))
diff2 <- 
  cbind.data.frame(prior_msymp$diff_w2,rep('Contrast A2',length(prior_msymp$diff_w2)))
diff3 <- 
  cbind.data.frame(prior_msymp$diff_w3,rep('Contrast A3',length(prior_msymp$diff_w3)))
diff4 <- 
  cbind.data.frame(prior_msymp$diff_w4,rep('Contrast A4',length(prior_msymp$diff_w4)))

#rename
colnames(usual1) <- c("value",'condition')
colnames(usual2) <- c("value",'condition')
colnames(usual3) <- c("value",'condition')
colnames(usual4) <- c("value",'condition')
colnames(exercise1) <- c("value",'condition')
colnames(exercise2) <- c("value",'condition')
colnames(exercise3) <- c("value",'condition')
colnames(exercise4) <- c("value",'condition')
colnames(diff1) <- c("value",'condition')
colnames(diff2) <- c("value",'condition')
colnames(diff3) <- c("value",'condition')
colnames(diff4) <- c("value",'condition')

hist_symp_sev_t2 <- rbind.data.frame(usual1,usual2,usual3,usual4,
                                     exercise1,exercise2,exercise3,exercise4,
                                     diff1,diff2,diff3,diff4)
hist_symp_sev_t2$condition <- factor(hist_symp_sev_t2$condition,
                                     levels = c('Contrast A4','Contrast A3','Contrast A2','Contrast A1',
                                                'SAEP A4','SAEP A3','SAEP A2','SAEP A1',
                                                'UCEP A4','UCEP A3','UCEP A2','UCEP A1'))

hist_prior_symp_sev<- hist_symp_sev_t2 %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank())+
xlab('Symptom severity')
hist_prior_symp_sev

#TOTAL SYMPTOMS
colnames(prior_m_tot_symp)
usual1 <- cbind.data.frame(prior_m_tot_symp$uc_w1,rep('UCEP A1',length(prior_m_tot_symp$uc_w1)))
usual2 <- cbind.data.frame(prior_m_tot_symp$uc_w2,rep('UCEP A2',length(prior_m_tot_symp$uc_w2)))
usual3 <- cbind.data.frame(prior_m_tot_symp$uc_w3,rep('UCEP A3',length(prior_m_tot_symp$uc_w3)))
usual4 <- cbind.data.frame(prior_m_tot_symp$uc_w4,rep('UCEP A4',length(prior_m_tot_symp$uc_w4)))
exercise1 <- 
  cbind.data.frame(prior_m_tot_symp$ex_w1,rep('SAEP A1',length(prior_m_tot_symp$ex_w1)))
exercise2 <- 
  cbind.data.frame(prior_m_tot_symp$ex_w2,rep('SAEP A2',length(prior_m_tot_symp$ex_w2)))
exercise3 <- 
  cbind.data.frame(prior_m_tot_symp$ex_w3,rep('SAEP A3',length(prior_m_tot_symp$ex_w3)))
exercise4 <- 
  cbind.data.frame(prior_m_tot_symp$ex_w4,rep('SAEP A4',length(prior_m_tot_symp$ex_w4)))
diff1 <- 
  cbind.data.frame(prior_m_tot_symp$diff_w1,rep('Contrast A1',length(prior_m_tot_symp$diff_w1)))
diff2 <- 
  cbind.data.frame(prior_m_tot_symp$diff_w2,rep('Contrast A2',length(prior_m_tot_symp$diff_w2)))
diff3 <- 
  cbind.data.frame(prior_m_tot_symp$diff_w3,rep('Contrast A3',length(prior_m_tot_symp$diff_w3)))
diff4 <- 
  cbind.data.frame(prior_m_tot_symp$diff_w4,rep('Contrast A4',length(prior_m_tot_symp$diff_w4)))

#rename
colnames(usual1) <- c("value",'condition')
colnames(usual2) <- c("value",'condition')
colnames(usual3) <- c("value",'condition')
colnames(usual4) <- c("value",'condition')
colnames(exercise1) <- c("value",'condition')
colnames(exercise2) <- c("value",'condition')
colnames(exercise3) <- c("value",'condition')
colnames(exercise4) <- c("value",'condition')
colnames(diff1) <- c("value",'condition')
colnames(diff2) <- c("value",'condition')
colnames(diff3) <- c("value",'condition')
colnames(diff4) <- c("value",'condition')

hist_symp_tot_t2 <- rbind.data.frame(usual1,usual2,usual3,usual4,
                                     exercise1,exercise2,exercise3,exercise4,
                                     diff1,diff2,diff3,diff4)
hist_symp_tot_t2$condition <- factor(hist_symp_tot_t2$condition,
                                     levels = c('Contrast A4','Contrast A3','Contrast A2','Contrast A1',
                                                'SAEP A4','SAEP A3','SAEP A2','SAEP A1',
                                                'UCEP A4','UCEP A3','UCEP A2','UCEP A1'))

hist_prior_symp_tot<- hist_symp_tot_t2 %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank())+
  xlab('Total symptoms')
hist_prior_symp_tot

#FEEL
colnames(prior_feel)
usual1 <- cbind.data.frame(prior_feel$uc_w1,rep('UCEP A1',length(prior_feel$uc_w1)))
usual2 <- cbind.data.frame(prior_feel$uc_w2,rep('UCEP A2',length(prior_feel$uc_w2)))
usual3 <- cbind.data.frame(prior_feel$uc_w3,rep('UCEP A3',length(prior_feel$uc_w3)))
usual4 <- cbind.data.frame(prior_feel$uc_w4,rep('UCEP A4',length(prior_feel$uc_w4)))
exercise1 <- 
  cbind.data.frame(prior_feel$ex_w1,rep('SAEP A1',length(prior_feel$ex_w1)))
exercise2 <- 
  cbind.data.frame(prior_feel$ex_w2,rep('SAEP A2',length(prior_feel$ex_w2)))
exercise3 <- 
  cbind.data.frame(prior_feel$ex_w3,rep('SAEP A3',length(prior_feel$ex_w3)))
exercise4 <- 
  cbind.data.frame(prior_feel$ex_w4,rep('SAEP A4',length(prior_feel$ex_w4)))
diff1 <- 
  cbind.data.frame(prior_feel$diff_w1,rep('Contrast A1',length(prior_feel$diff_w1)))
diff2 <- 
  cbind.data.frame(prior_feel$diff_w2,rep('Contrast A2',length(prior_feel$diff_w2)))
diff3 <- 
  cbind.data.frame(prior_feel$diff_w3,rep('Contrast A3',length(prior_feel$diff_w3)))
diff4 <- 
  cbind.data.frame(prior_feel$diff_w4,rep('Contrast A4',length(prior_feel$diff_w4)))

#rename
colnames(usual1) <- c("value",'condition')
colnames(usual2) <- c("value",'condition')
colnames(usual3) <- c("value",'condition')
colnames(usual4) <- c("value",'condition')
colnames(exercise1) <- c("value",'condition')
colnames(exercise2) <- c("value",'condition')
colnames(exercise3) <- c("value",'condition')
colnames(exercise4) <- c("value",'condition')
colnames(diff1) <- c("value",'condition')
colnames(diff2) <- c("value",'condition')
colnames(diff3) <- c("value",'condition')
colnames(diff4) <- c("value",'condition')

hist_feel_t2 <- rbind.data.frame(usual1,usual2,usual3,usual4,
                                 exercise1,exercise2,exercise3,exercise4,
                                 diff1,diff2,diff3,diff4)
hist_feel_t2$condition <- factor(hist_feel_t2$condition,
                                 levels = c('Contrast A4','Contrast A3','Contrast A2','Contrast A1',
                                            'SAEP A4','SAEP A3','SAEP A2','SAEP A1',
                                            'UCEP A4','UCEP A3','UCEP A2','UCEP A1'))

hist_prior_feel<- hist_feel_t2 %>%
  ggplot(aes(x = value, y = condition,fill = condition)) +
  stat_halfeye(.width = c(.50, .90)) +theme(legend.position ='none',
                                            axis.title.y = element_blank())+
  xlab('% of Normal')
hist_prior_feel

#combine
t2_prior_sim <- plot_grid(hist_prior_symp_sev,hist_prior_symp_tot,
                          hist_prior_feel,labels = "AUTO",ncol=3)
t2_prior_sim

#save 
ggsave('s3_fig.jpg',t2_prior_sim,dpi=600,width = 10)







