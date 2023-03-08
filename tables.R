#script to create Table 2

#Note: Table 1 script is not available as it contains demographic information...
#that we are not permitted to share publicly. If needed, please contact the...
#corresponding author.

#Libraries required
library(gtsummary)
library(gt)
library(tidyverse)

#load files
t2_df <- read.csv("t2_df.csv",stringsAsFactors = FALSE)
#change colnames
colnames(t2_df) <- c('randomization_group','redcap_repeat_instance',
                     'Total Symptoms','Symptom Severity','% of Normal')

####TABLE 2####
#prep according to assessment
enrolment_df <- t2_df[t2_df$redcap_repeat_instance=="Enrolment",]
ass1_df <-  t2_df[t2_df$redcap_repeat_instance=="Assessment 1 (MED = 6.5 Days)",]
ass2_df <-  t2_df[t2_df$redcap_repeat_instance=="Assessment 2 (MED = 13 Days)",]
ass3_df <-  t2_df[t2_df$redcap_repeat_instance=="Assessment 3 (MED = 20 Days)",]
ass4_df <-  t2_df[t2_df$redcap_repeat_instance=="Assessment 4 (MED = 28 Days)",]

#gather tables for tbl_strata
#enrolment
enrolment_dfg <-gather(enrolment_df,`Symptom Evaluation Type`, 
                       Score,`Total Symptoms`:`% of Normal`)
enrolment_dfg $`Symptom Evaluation Type` <- 
  factor(enrolment_dfg$`Symptom Evaluation Type`,
  levels = c('Symptom Severity','Total Symptoms','% of Normal'))
enrolment_dfg $redcap_repeat_instance <- NULL
colnames(enrolment_dfg )[3] <-"Enrolment"

#Assessment 1
ass1_dfg <-gather(ass1_df,`Symptom Evaluation Type`, 
                  Score,`Total Symptoms`:`% of Normal`)
ass1_dfg $`Symptom Evaluation Type` <- factor(ass1_dfg$`Symptom Evaluation Type`,
                                              levels = c('Symptom Severity',
                                                         'Total Symptoms',
                                                         '% of Normal'))
ass1_dfg $redcap_repeat_instance <- NULL
colnames(ass1_dfg)[3] <-"Assessment 1 (MED = 6.5 days)"

#Assessment 2
ass2_dfg <-gather(ass2_df,`Symptom Evaluation Type`, 
                  Score,`Total Symptoms`:`% of Normal`)
ass2_dfg $`Symptom Evaluation Type` <- factor(ass2_dfg$`Symptom Evaluation Type`,
                                              levels = c('Symptom Severity',
                                                         'Total Symptoms',
                                                         '% of Normal'))
ass2_dfg $redcap_repeat_instance <- NULL
colnames(ass2_dfg)[3] <-"Assessment 2 (MED = 13 days)"

#Assessment 3
ass3_dfg <-gather(ass3_df,`Symptom Evaluation Type`, 
                  Score,`Total Symptoms`:`% of Normal`)
ass3_dfg $`Symptom Evaluation Type` <- factor(ass3_dfg$`Symptom Evaluation Type`,
                                              levels = c('Symptom Severity',
                                                         'Total Symptoms',
                                                         '% of Normal'))
ass3_dfg $redcap_repeat_instance <- NULL
colnames(ass3_dfg)[3] <-"Assessment 3 (MED = 20 days)"

#Assessment 4
ass4_dfg <-gather(ass4_df,`Symptom Evaluation Type`, 
                  Score,`Total Symptoms`:`% of Normal`)
ass4_dfg $`Symptom Evaluation Type` <- factor(ass4_dfg$`Symptom Evaluation Type`,
                                              levels = c('Symptom Severity',
                                                         'Total Symptoms',
                                                         '% of Normal'))
ass4_dfg $redcap_repeat_instance <- NULL
colnames(ass4_dfg)[3] <-"Assessment 4 (MED = 28 days)"


#Table Creation
theme_gtsummary_journal(journal = "nejm")
reset_gtsummary_theme()
enrolment_tbl <- 
  tbl_strata(enrolment_dfg,
             strata = `Symptom Evaluation Type` ,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = `randomization_group`, missing = "no",
                           type = c('Enrolment') ~ 'continuous',
                           digits = all_continuous() ~ 0)) %>% 
  bold_labels()%>%
  modify_caption("**Table 2. Participant Symptoms Across Assessments**") 

ass1_tbl <- 
  tbl_strata(ass1_dfg,
             strata = `Symptom Evaluation Type` ,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = `randomization_group`,missing = "no",
                           type = c('Assessment 1 (MED = 6.5 days)') ~ 'continuous',
                           digits = all_continuous() ~ 0)) %>% 
  bold_labels()


ass2_tbl <- 
  tbl_strata(ass2_dfg,
             strata = `Symptom Evaluation Type` ,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = `randomization_group`,missing = "no",
                           type = c('Assessment 2 (MED = 13 days)') ~ 'continuous',
                           digits = all_continuous() ~ 0)) %>% 
  bold_labels()


ass3_tbl <- 
  tbl_strata(ass3_dfg,
             strata = `Symptom Evaluation Type` ,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = `randomization_group`,missing = "no",
                           type = c('Assessment 3 (MED = 20 days)') ~ 'continuous',
                           digits = all_continuous() ~ 0)) %>% 
  bold_labels()


ass4_tbl <- 
  tbl_strata(ass4_dfg,
             strata = `Symptom Evaluation Type` ,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = `randomization_group`,missing = "no",
                           type = c('Assessment 4 (MED = 28 days)') ~ 'continuous',
                           digits = all_continuous() ~ 0)) %>% 
  bold_labels()


t2_stacked <- tbl_stack(list(enrolment_tbl,ass1_tbl,ass2_tbl,ass3_tbl,
                             ass4_tbl))%>%
  modify_header(update = list(label ~ '**Assessment**'))%>%
  as_gt() %>%
  tab_source_note(md('MED, median; SAEP, structured aerobic exercise protocol;
                     UCEP, usual care exercise protocol.'))%>%
  tab_source_note(md('Data presented as the Median (IQR), or n (%);
  IQR = interquartile range.'))%>%
  tab_source_note(md('Due to loss to follow-up and missing assessments, 
                     sample sizes differed across assessments: Enrolment, 
                     SAEP = 20 vs, UCEP = 19; assessment 1, SAEP = 19 vs. UCEP = 19; 
                     assessment 2, SAEP = 16 vs. UCEP = 19; assessment 3, 
                     SAEP = 16 vs. UCEP = 18; assessment 4, SAEP = 14 vs. UCEP = 18.'))%>%
  gt::gtsave(filename = "t2.html") 
