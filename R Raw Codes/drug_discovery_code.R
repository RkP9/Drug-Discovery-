library(tidyverse)
library(dplyr)
library(ggplot2)
library(corrplot)


dd <- read.csv("/Users/rahulpurnapatre/Downloads/drug_discovery.csv")
glimpse(dd)
str(dd)

print(colSums(is.na(dd))) ##This shows the number of missing values (60 in hydrophobicity)

dd$logp[is.na(dd$logp)] <- median(dd$logp, na.rm = TRUE)
dd$hydrophobicity[is.na(dd$hydrophobicity)] <- median(dd$hydrophobicity, na.rm = TRUE)
dd$polar_surface_area[is.na(dd$polar_surface_area)] <- median(dd$polar_surface_area, na.rm = TRUE)

print(colSums(is.na(dd)))

summary(dd[,c("molecular_weight","logp","binding_affinity")]) #why only[]?#


ggplot(dd, aes(y = molecular_weight)) +
  geom_boxplot(fill = "skyblue", outlier.color = "red", outlier.shape = 18, outlier.size = 2) +
  labs(title = "Distribution of Molecular Weight", y = "Molecular Weight (g/mol)") +
  theme_minimal()

ggplot(dd, aes(y=logp))+
  geom_boxplot(fill="lightgreen",outlier.color = "red",outlier.shape = 18,outlier.size = 2)+
  labs(title = "Distribution of LogP",y="LogP")+
  theme_minimal()

ggplot(dd, aes(y = binding_affinity)) +
  geom_boxplot(fill = "salmon", outlier.color = "red", outlier.shape = 18, outlier.size = 2) +
  labs(title = "Distribution of Binding Affinity", y = "Binding Affinity (pKi)") +
  theme_minimal()


summary_stats <- dd %>% 
  summarise(
    Variable= c("Molecular Affinity","LogP","Binding Affinity"),
    Mean= c(mean(molecular_weight),mean(logp),mean(binding_affinity)),
    Median= c(median(molecular_weight),median(logp),median(binding_affinity)),
    SD = c(sd(molecular_weight), sd(logp), sd(binding_affinity)),
    Min = c(min(molecular_weight), min(logp), min(binding_affinity)),
    Max = c(max(molecular_weight), max(logp), max(binding_affinity))
  )
print(summary_stats, row.names = FALSE)


ggplot(dd, aes(x = molecular_weight)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean(molecular_weight)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Molecular Weight", x = "Molecular Weight (g/mol)", y = "Frequency") +
  theme_minimal()


ggplot(dd, aes(x = logp)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  geom_vline(aes(xintercept = mean(logp)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of LogP", x = "LogP", y = "Frequency") +
  theme_minimal()


ggplot(dd, aes(x = binding_affinity)) +
  geom_histogram(bins = 30, fill = "salmon", color = "black") +
  geom_vline(aes(xintercept = mean(binding_affinity)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Binding Affinity", x = "Binding Affinity (pKi)", y = "Frequency") +
  theme_minimal()


dd$active_status <- factor(dd$active,
                           levels = c(0,1),
                           labels = c("Inactive","Active")) #EXPLAIN?#


activity_grouped_summary <- dd %>% 
  group_by(active_status) %>% 
  summarize(
    Count = n(),
    Mean_MW = mean(molecular_weight),
    Median_MW = median(molecular_weight),
    Mean_LogP = mean(logp),
    Median_LogP = median(logp),
    Mean_BA= mean(binding_affinity),
    Median_BA= median(binding_affinity)
  )
print(activity_grouped_summary)

ggplot(dd, aes(x = active_status, y = molecular_weight, fill = active_status)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Molecular Weight vs. Compound Activity",
       x = "Compound Activity",
       y = "Molecular Weight (g/mol)") +
  theme_minimal() +
  scale_fill_manual(values = c("Inactive" = "red", "Active" = "blue")) +
  theme(legend.position = "none") # Hide legend as colors are self-explanatory

ggplot(dd, aes(x = active_status, y = logp, fill = active_status)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "LogP vs. Compound Activity",
       x = "Compound Activity",
       y = "LogP") +
  theme_minimal() +
  scale_fill_manual(values = c("Inactive" = "red", "Active" = "blue")) +
  theme(legend.position = "none")

ggplot(dd,aes(x=active_status, y=binding_affinity, fill=active_status))+
  geom_boxplot(alpha=0.8)+
  labs(title= "Binding Affinity vs Compound Activity",
       x="Compound Activity",
       y="Binding Affinity")+
  theme_minimal()+
  scale_fill_manual(values = c("Inactive" = "red", "Active"= "blue"))+
  theme(legend.position= "none" )


numeric_dd <- dd %>% 
  select(molecular_weight,logp,binding_affinity)


correlation_matrix <- cor(numeric_dd)
print(round(correlation_matrix,2))



corrplot(correlation_matrix,
         method="color",
         type="upper",
         order="hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         diag = FALSE)

ggplot(dd, aes(x = logp, y = binding_affinity)) +
  geom_point(aes(color = active_status), alpha = 0.6) + 
  geom_smooth(method = "lm", color = "black", se = FALSE) 
  scale_color_manual(values = c("Inactive" = "red", "Active" = "blue")) +
  labs(title = "Binding Affinity vs. LogP",
       x = "LogP",
       y = "Binding Affinity (pKi)",
       color = "Activity Status") +
  theme_minimal()

logp_ttest <- t.test(logp ~ active_status, data = dd)
print(logp_ttest)


protein_activity_table <- table(dd$protein_id, dd$active_status)


chi_test <- chisq.test(protein_activity_table)

print(chi_test)
