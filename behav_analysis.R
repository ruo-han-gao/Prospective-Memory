library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)

wd = "~/Documents/RPTU/SS_2025/LB_Neuroscience/data_analysis/behavioral"
setwd(wd)

######################################################
#####    Functions for triggers and accuracy    ######
######################################################

get_trigger = function(df){
  df %>%
    mutate(trigger = ifelse(Valence == "Neutral" & Arousal == "Neutral", Stimulus + 20,
                            ifelse(Valence == "Positive" & Arousal == "High", Stimulus + 40,
                                   ifelse(Valence == "Positive" & Arousal == "Low", Stimulus + 60,
                                          ifelse(Valence == "Negative" & Arousal == "High", Stimulus + 80,
                                                 ifelse(Valence == "Negative" & Arousal == "Low", Stimulus + 100, NA))))))
}

get_accuracy <- function(df) {
  df %>%
    mutate(
      correct_response = case_when(
        Trial %in% c(1, 2) ~ NA_real_,
        trigger %in% PM.cues ~ 2,
        participant %% 2 == 1 & Stimulus == LastStimuli ~ 1,
        participant %% 2 == 1 & Stimulus != LastStimuli ~ 3,
        participant %% 2 == 0 & Stimulus == LastStimuli ~ 3,
        participant %% 2 == 0 & Stimulus != LastStimuli ~ 1,
        TRUE ~ NA_real_
      ),
      
      accuracy = case_when(
        !(Trial %in% c(1, 2)) & RT == 0 ~ 0,
        Response == correct_response ~ 1,
        Response != correct_response ~ 0,
        TRUE ~ NA_real_
      )
    )
}

# PM triggers
PM.cues = c(33, 34, 53, 54, 73, 74, 93, 94, 113, 114)


######################################################
#####     Prepare data from main PM blocks      ######
######################################################

# get behavioral response files
pm.mainfiles <- list.files(path = "behav_data", pattern = "^MainExperiment_[0-9]+\\.txt$", full.names = TRUE)

# combine the data across participants
pm.data <- do.call(rbind, lapply(pm.mainfiles, function(file) {
  participant_num <- sub(".*MainExperiment_([0-9]+)\\.txt$", "\\1", file)
  
  df <- read.table(file, header = TRUE, sep = "\t")
  
  #add participant number for random effect analysis
  df$participant <- as.integer(participant_num)
  return(df)
}))

# compute trigger and accuracy
pm.data = get_trigger(pm.data)
pm.data = get_accuracy(pm.data)


# check accuracy
acc.check = pm.data %>%
  group_by(participant) %>%
  summarise(
    accuracy_percent = round(mean(accuracy, na.rm = TRUE) * 100, 2)
  )

# Unusual low accuracy found in participant 40, 41, and 46. 

# exclude those with accuracy lower than 70%
acc.lower.70 <- acc.check %>%
  filter(accuracy_percent < 70) %>%
  pull(participant)

data <- pm.data %>%
  filter(!participant %in% acc.lower.70)

length(unique(data$participant)) # 21 participants remain


######################################################
#####                  Accuracy                 ######
######################################################

# general accuracy
mean(data$accuracy, na.rm = TRUE) #0.8913

# accuracy per person
acc.per.person = data %>%
  group_by(participant) %>%
  summarise(
    accuracy_percent = round(mean(accuracy, na.rm = TRUE) * 100, 2)
  )
range(acc.per.person$accuracy_percent)  # 77.35-94.12

# accuracy per condition
acc.per.cond <- data %>%
  group_by(Valence) %>%
  summarise(
    accuracy_percent = mean(accuracy, na.rm = TRUE) * 100,
    n = n(),
    .groups = "drop"
  )
acc.per.cond$Valence # "Negative" "Neutral"  "Positive"
round(acc.per.cond$accuracy_percent,2) # 89.01 86.62 90.51

# set contrast: dummy coding use "Neutral" as reference
data$Valence = as.factor(data$Valence) # "Negative", "Neutral", "Positive"
contrasts(data$Valence) = contr.treatment(3, base = 2)
colnames(contrasts(data$Valence)) = c("+Neg vs. N", "+Pos vs. N")
contrasts(data$Valence)

# compare accuracy across conditions
mod.acc <- glmer(accuracy ~ Valence + (1|participant), data = data, family = binomial())
summary(mod.acc)

exp(fixef(mod.acc))
# Neg vs. N: 25% higher odds of a correct response
# Pos vs. N: 48% higher odds of a correct response


######################################################
#####             Reaction times                ######
######################################################

# Trim data
data.trim <- data %>%
  filter(accuracy == 1 & RT != 0) %>%
  group_by(participant) %>%
  filter(between(RT,
                 mean(RT, na.rm = TRUE) - 2.5 * sd(RT),
                 mean(RT, na.rm = TRUE) + 2.5 * sd(RT))) %>%
  ungroup()

hist(data.trim$RT)

# 14.84% of data were trimmed out
round((1-length(data.trim$RT)/length(data$RT))*100,2)

# Plot RT
ggplot(data.trim, aes(x = Valence, y = RT, group = 1)) +  # group=1 ensures line connects points
  stat_summary(fun = mean, geom = "point", size = 3, color = "darkblue") +
  stat_summary(fun = mean, geom = "line", size = 0.5, color = "darkblue") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "darkblue") +
  labs(
    title = " ",
    x = "Valence",
    y = "Reaction time (ms)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


# Set contrasts: dummy coding use "Neutral" as reference
data.trim$Valence = as.factor(data.trim$Valence) # "Negative", "Neutral", "Positive"
contrasts(data.trim$Valence) = contr.treatment(3, base = 2)
colnames(contrasts(data.trim$Valence)) = c("+Neg vs. N", "+Pos vs. N")
contrasts(data.trim$Valence)

# compare RTs across conditions
mod.rt = lmer(RT ~ Valence + (Valence| participant), data = data.trim)

# check residuals
qqnorm(resid(mod.rt))
qqline(resid(mod.rt))

# model summary
summary(mod.rt)

