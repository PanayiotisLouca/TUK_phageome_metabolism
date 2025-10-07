# Characterising the gut phageome and its links with host metabolism in the general population

*Panayiotis Louca, Mohammadali Khan Mirzaei, Afroditi Kouraki, Erfan Khamespanah, Yu Lin, Xue Peng, Robert Pope, Alessia Visconti, Francesco Asnicar, Daniel Kirk, Ricardo Costeira, Nicola Segata, Mario Falchi, Jordana T. Bell, Tim D. Spector, Lindsey A. Edwards, Li Deng, Ana M. Valdes, Cristina Menni*

---

  Analysis scripts used in our comprehensive study of the human gut phageome in the general populations and its links to metabolic health (DOI: To be updated). 
    
> *If you use this code, please cite:*
> > Louca, P. et al. (2025). *Characterising the gut phageome and its links with host metabolism in the general population*. [DOI: To be updated]


---

  ## 📂 Repository Structure
```
├── 1.de_novo_Assembly_and_Profiling
│   └── 1.1.Preprocessing_Host_Decontamination.txt
│   └── 1.2.Viral_Genome_Assembly_and_Detection.txt
│   └── 1.3.Clustring_Viral_Contigs_Mapping_Reads.txt
│   └── 1.4.Read_Count_Extraction_TPM_Calculation.txt
│   └── 1.5.Taxa_Host_Rep_cycle_AMG_Prediction.txt
├── 2.Phage_Heritability
│   └── 2.1.Phage_alpha_div_heritability_SCRIPT.R
│   └── 2.2.Phage_beta_div_heritability_SCRIPT.R
│   └── 2.3.Individual_phage_heritability_SCRIPT_CREATE.R
├── 3.Phage_Taxa_Network_Analysis
│   ├── 3.1.Phage_species_ggraph_network_POSITIVE_SCRIPT.R
│   ├── 3.2.Phage_species_ggraph_network_NEGATIVE_SCRIPT.R
├── 4.Phage_Bacterial_Metabolite_Associations_and_Network
│   ├── 4.1.Phage_bile_acids_CREATE_ARRAY_SCRIPT.R
│   ├── 4.2.Phage_SCFA_CREATE_ARRAY_SCRIPT.R
│   └── 4.3.Phage_SCFA_bile_acid_network_SCRIPT.R
├── 5.Phage_Metabolic_Health_Associations
│   ├── 5.1.Phages_metabolic_health_Maaslin2_SCRIPT_CREATE.R
└── 6.Phage_Diet_Associations
    └── 6.1.Phages_diet_CREATE_ARRAY_SCRIPT.R
```

  ## ℹ️ Repository Information 

  - **1. de novo Assembly and Phageome Profiling**
      - `1.1.Preprocessing_Host_Decontamination.txt`: Pipeline for preprocessing paired-end sequencing reads.
       - `1.2.Viral_Genome_Assembly_and_Detection.txt`: Pipeline for viral genome assembly and detection.
       - `1.3.Clustring_Viral_Contigs_Mapping_Reads.txt`: Pipeline for clustering viral contigs and mapping reads.
       - `1.4.Read_Count_Extraction_TPM_Calculation.txt`: Pipeline for Extracting Read Counts from SAM/BAM and Calculating TPM.
       - `1.4.Read_Count_Extraction_TPM_Calculation.txt`: Pipeline for Extracting Read Counts from SAM/BAM and Calculating TPM.
       - `1.5.Taxa_Host_Rep_cycle_AMG_Prediction.txt`: Pipeline for Taxonomic, host, replication cycle, and Auxiliary Metabolic Gene Assignment.
       
  - **2. Phage Heritability**
      - `2.1.Phage_alpha_div_heritability_SCRIPT.R`: Uses the `mets` package to estimate heritability of alpha diversity of the gut phageome (Shannon index and observed taxonomic richness).
      - `2.2.Phage_beta_div_heritability_SCRIPT.R`: As above for principal components of phageome beta-diversity (Bray-Curtis index).
      - `2.3.Individual_phage_heritability_SCRIPT_CREATE.R`: As above for individual viral contigs.

  - **3. Phage Taxa Network Analysis**
      - `3.1.Phage_species_ggraph_network_POSITIVE_SCRIPT.R`: Constructs network for positive correlations (Figure 3A).
      - `3.2.Phage_species_ggraph_network_NEGATIVE_SCRIPT.R`: Builds network for negative correlations between phages and taxa (Figure 3B).

  - **4. Phage Bacterial Metabolite Associations and Network**
      - `4.1.Phage_bile_acids_CREATE_ARRAY_SCRIPT.R`: Assesses associations between viral contigs and bile acids.
      - `4.2.Phage_SCFA_CREATE_ARRAY_SCRIPT.R`: As above, for short-chain fatty acids.
      - `4.3.Phage_SCFA_bile_acid_network_SCRIPT.R`: Build correlation network between viral contigs and short-chain fatty acids & bile acids.

  - **5. Phage Metabolic Health Associations**
      - `5.1.Phages_metabolic_health_Maaslin2_SCRIPT_CREATE.R`: Tests associations between viral contigs & metabolic health parameters (Triglycerides, TyG index, BMI, & glucose) using linear mixed effect modelling.
      - `5.2.pathways_MaAslin2_SCRIPT.R`: Models associations between neoplasia and MetaCyc pathways.

  - **6. Phage Diet Associations**
      - `6.1.Phages_diet_CREATE_ARRAY_SCRIPT.R`: MTests associations between viral contigs & alpha diveristy of the gut phageome and dietary nutrient intakes and the healthy eating index using linear mixed effect models.
