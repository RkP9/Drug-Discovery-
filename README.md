# Drug-Discovery-
ID relationship between molecular properties 

## Objective
The primary objective of this project was to analyze a dataset of 2000 virtual compounds to identify key relationships between molecular properties (e.g., LogP, molecular weight), target proteins, and compound activity. The analysis employed descriptive statistics, group comparisons, correlation analysis, and formal statistical testing to uncover significant trends.

## Key Findings & Insights
The analysis yielded several clear and statistically significant findings that point to the properties driving compound activity.

### Lipophilicity (LogP) is the Dominant Factor for Activity 
The most significant discovery is the strong link between a compound's lipophilicity (LogP) and its biological activity. Active compounds are significantly more lipophilic than their inactive counterparts, with a mean LogP of 4.81 versus 2.90 for inactive compounds. A t-test confirmed this difference is highly statistically significant (p < 2.2e-16). This is further supported by a strong, positive correlation (r = 0.60) between LogP and binding affinity.

### Molecular Weight Does Not Influence Activity 
In contrast to LogP, a compound's molecular weight showed no meaningful relationship with its activity. The mean molecular weight was nearly identical for both active and inactive groups, and the correlation with binding affinity was negligible (r = -0.01). This suggests that compound size is not a critical design parameter for activity against the targets in this dataset.

### Activity is Independent of Specific Protein Targets 
A Chi-Squared test was performed to check for an association between the specific protein ID and compound activity. The result was not statistically significant (p = 0.83), indicating that the proportion of active compounds is consistent across all tested proteins. This implies that the importance of high LogP is a general requirement for activity in this chemical space, rather than being a phenomenon tied to a specific protein.

## Conclusion
This exploratory and statistical analysis robustly demonstrates that lipophilicity is the most critical molecular property associated with compound activity in this dataset.

For a real-world drug discovery project based on these findings, the primary recommendation would be to prioritize the design and synthesis of compounds with higher LogP values. This strategy would significantly increase the probability of identifying potent hits against the target class. Molecular weight can be considered a secondary factor, used for optimizing other properties without being a primary driver for activity itself.

