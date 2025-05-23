Readme
================

# About

This repository contains data, code and plots for our paper on somatic
GRIN2A variants (doi pending)

# Code

## Download files

### gnomAD

- Downloaded all *GRIN2A* variants from gnomAD v.2.1.1 in a csv
- converted to Excel csv (Semicolon)
- filtered for Transcript ENST00000396573.2 (1566 out of 1573 variants)
- combined columns “Filters - exomes” and “Filters - genomes” into
  “Filters - all” –\> All variants passed either Exome or Genome
  Filters, so no variant had to be exluded –\> Deleted Column again
- Added Column “VariantID” according to gnomAD standard
  (Chromosome-Position-Reference-Alternate)
- Saved as 01-1_gnomAD_v2.1.1_ENST00000396573.2_2022_11_14_07_59_07.csv

### Cosmic

- Downloaded the GRCh37 VCF Files (coding and non-coding mutations) from
  <https://cancer.sanger.ac.uk/cosmic/download>
- Indexed the VCF with `tabix -p vcf CosmicCodingMuts.vcf.gz` in WSL2
- Extracted the Region with GRIN2A based on Emsembl Coordinates with
  `tabix CosmicCodingMuts.vcf.gz 16:9,847,261-10,276,785 > 02_Cosmic_GRIN2A.vcf`
- Downloaded the Metadata (for sample information etc) from COSMIC
  Mutation Data –\> Filter by Gene: GRIN2A
- removed whitespace in headers
- Filtered column “GRCH” != null
- Saved as “02-1_V96_37_MUTANT.csv”

### ClinVar

- Use the ClinVar vcf and tsv from the ftp site

## Load libraries

``` r
library(pandoc)
library(readr)
library(dplyr)
library(vcfR)
library(stringr)
library(fuzzyjoin)
library(tidyr)
library(ggplot2)
library(ggvenn)
library(cowplot)
library(tibble)
```

## Read and prepare cosmic, gnomAD and clinVar

First read the data files. The .vcf files have to be edited a little bit
to extract the metadata. Filter the files for *GRIN2A* and/or the right
transcript. We also filter the **Cosmic** dataset for **confirmed
somatic** variants. In **ClinVar** we **exclude benign and somatic**
variants.

``` r
gnomAD = read_csv2("00_Data/01_gnomAD_Data/01-1_gnomAD_v2.1.1_ENST00000396573.2_2022_11_14_07_59_07.csv")
cosmic_meta = read_csv2("00_Data/02_Cosmic_Data/02-1_V96_37_MUTANT.csv")
cosmic_vcf = read.vcfR("00_Data/02_Cosmic_Data/02_Cosmic_GRIN2A.vcf", verbose = T, limit = 6e+09)
#clinvar_vcf = read.vcfR("00_Data/03_ClinVar_Data/03_clinvar.vcf.gz", verbose = T, limit = 6e+09)
#save(clinvar_vcf, file = "00_Data/00_RData/clinvar_vcf.RData")
load(file = "00_Data/00_RData/clinvar_vcf.RData")


### Edit ClinVar Data ###
clinvar_fix = as_tibble(getFIX(clinvar_vcf)) %>% 
  mutate(Key = seq(1:nrow(.))) %>% 
  mutate(Variant = paste(CHROM, POS, REF, ALT, sep = "-"))

#clinvar_meta = extract_info_tidy(clinvar_vcf, info_fields = NULL,
#                                             info_types = TRUE, info_sep = ";")
#save(clinvar_meta, file = "00_Data/00_RData/clinvar_meta.RData")
load(file = "00_Data/00_RData/clinvar_meta.RData")

clinvar = clinvar_fix %>% 
  left_join(clinvar_meta, by = c("Key" = "Key")) %>% 
  filter(grepl("GRIN2A", GENEINFO)) %>% 
  filter(ORIGIN != 2) %>% # somatic variants
  filter(CLNSIG == "Pathogenic" | CLNSIG == "Likely_pathogenic" | CLNSIG == "Pathogenic/Likely_pathogenic")

### Edit Cosmic Data ###
cosmic_vcf_fix = as_tibble(getFIX(cosmic_vcf)) %>% 
  mutate(Key = seq(1:nrow(.))) %>% 
  mutate(Variant = paste(CHROM, POS, REF, ALT, sep = "-"))

cosmic_vcf_meta = extract_info_tidy(cosmic_vcf, info_fields = NULL,
                                             info_types = TRUE, info_sep = ";")

cosmic_vcf_combined = cosmic_vcf_fix %>% 
  left_join(cosmic_vcf_meta, by = c("Key" = "Key")) %>% 
  filter(GENE == "GRIN2A")

cosmic = cosmic_vcf_combined %>% 
  left_join(cosmic_meta, by = c("LEGACY_ID" = "LEGACY_MUTATION_ID"), na_matches = "never") %>% 
  filter(MUTATION_SOMATIC_STATUS == "Confirmed somatic variant")

gnomAD_variants = unique(sort(gnomAD$VariantID))
cosmic_variants = unique(sort(cosmic$Variant))
clinvar_variants = unique(sort(clinvar$Variant))
```

## Build variant table and generate vcf for annotation

Now we have variants per dataset, but all have different forms of
annotation.

So we **generate a dummy vcf** and **annotate with Ensembl VEP**.

``` r
variantvector = c(clinvar_variants, 
                  cosmic_variants, 
                  gnomAD_variants) %>% 
  sort() %>% 
  unique()

Variants = tibble(Variant = variantvector)

Variants_vcf = Variants %>% 
  select(Variant) %>% 
  separate(Variant, into = c("CHROM", "POS", "REF", "ALT"), sep = "-") %>% 
  mutate(Dummy1 = ".") %>% 
  mutate(Dummy2 = ".") %>% 
  mutate(Dummy3 = ".") %>% 
  mutate(Dummy4 = ".") %>% 
  relocate(CHROM, POS, Dummy1, REF, ALT, Dummy2, Dummy3, Dummy4)
  

write_tsv(Variants_vcf, "01_Ensembl_annotation/01_ExampleData/Test_incl_intronic.vcf")
```

***Paste this data to the online Ensembl VEP and download the results*,
filter for transcript NM_001134407.3**.

<https://grch37.ensembl.org/Homo_sapiens/Tools/VEP>

Parameters: RefSeq Transcripts, include HGVS. Save in folder
01_Ensembl_annotation/02_EnsemblVEP_output.

## Edit VEP output

The output files from the VEP need to be edited a little bit

``` r
# Get the Columns from the output and put them in a vector
# The output will be split into columns accordingly
Columns = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|REFSEQ_MATCH|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|HGVS_OFFSET|AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
Columns2 = str_replace_all(Columns, "\\|", ", ")
Columns3 = str_split(Columns2, ", ")[[1]]

Variants_VEP = read.vcfR("01_Ensembl_annotation/02_EnsemblVEP_output/F4DbdVwaWKSOtBzG.Feature_is_NM_001134407.3.vcf_", 
                         verbose = T, limit = 6e+09)

# read the VEP file and separate the info column
Variants_VEP_tidy = (vcfR2tidy(Variants_VEP, info_only = T))$fix %>% 
  separate(CSQ, into = Columns3, sep = "\\|", extra = "drop", fill = "right", convert = T) %>% 
  mutate(ID = paste(CHROM, POS, REF, ALT, sep = "-")) %>% # add the ID for joining
  mutate(HGVSp = str_replace_all(HGVSp, "%3D", "=")) %>% # tidy up synonymous variants
  # join the datasets
  left_join(y = (clinvar %>% select(Variant) %>% mutate(ClinVar = T)),
            by = c("ID" = "Variant")) %>% 
  left_join(y = (gnomAD %>% select(VariantID) %>% mutate(gnomAD = T)), 
            by = c("ID" = "VariantID")) %>% 
  left_join(y = (cosmic %>% select(Variant) %>% mutate(Cosmic = T)), 
            by = c("ID" = "Variant")) %>% 
  distinct() %>% 
  # extract the protein position for plotting
  mutate(Protein_position2 = str_remove_all(Protein_position, "-[0-9].*")) %>% 
  mutate(Protein_position2 = if_else(Protein_position2 == "", NA_character_, Protein_position2)) %>% 
  filter(IMPACT == "MODERATE" | IMPACT == "HIGH") %>% # keep variants with high & moderate impact
  mutate(Protein_position2 = as.numeric(Protein_position2)) %>% 
  mutate(Consequence = if_else(str_detect(pattern = "intron", string = Consequence), 
                               "splice_variant_non_canonical", 
                               Consequence)) %>% # restructure consequence
  mutate(Consequence = if_else(str_detect(pattern = "splice_region_variant", 
                                          string = Consequence), 
                               "missense_variant", Consequence)) %>% 
  mutate(Null_variant = case_when(
         Consequence == "stop_lost" ~ T,
         Consequence == "stop_gained" ~ T,
         Consequence == "frameshift_variant" ~ T,
         Consequence == "start_lost" ~ T,
         Consequence == "splice_acceptor_variant" ~ T,
         Consequence == "splice_donor_variant" ~ T)) %>% 
  separate(HGVSc, into = c("Transcript_c", "c_code"), sep = ":") %>% 
  separate(HGVSp, into = c("Transcript_p", "p_code"), sep = ":")
```

## Read domains

Get the domain data which will be used for plotting

``` r
domains = read_tsv("00_Data/04_Domain_Structure/GRIN2A_domains.txt") %>% 
  mutate(mid = (End-Start)/2+Start)
domains$Domain = factor(domains$Domain, levels = domains$Domain)
```

## Plot Data

Now we have the data ready and can do some plots

### Variant density

We need to know where the variants are located. Do we see a hotspot? Do
somatic variants cluster differently?

First we write a plot function and then plot the data for: 1. All
variants 2. missense only 3. Null variants only

``` r
### write plot function ###
plot_variant_density = function(input_df, label_prefix = ""){
  g = ggplot(input_df, aes(x = Protein_position2)) + 
  geom_density(data = (input_df %>% select(Protein_position2, gnomAD) %>% filter(!is.na(gnomAD))), 
               aes(y = -after_stat(density), fill = "gnomAD"), alpha = 0.5) + 
  geom_density(data = (input_df %>% select(Protein_position2, ClinVar) %>% filter(!is.na(ClinVar))), 
               aes(fill = "ClinVar"), alpha = 0.5) + 
  geom_density(data = (input_df %>% select(Protein_position2, Cosmic) %>% filter(!is.na(Cosmic))), 
               aes(fill = "Cosmic"), alpha = 0.5) + 
  geom_segment(data = domains, 
               aes(x = Start, xend = End, color = Domain), 
               y = 0, yend = 0, 
               inherit.aes = F, size = 5, show.legend = F) + 
  geom_text(data = domains, 
            aes(x = mid, y = 0, label = Domain), 
            color = "white", fontface = "bold") + 
  labs(fill = "Dataset", y = paste(label_prefix, "Variants (Density)"), x = "Amino Acid Position") +
  scale_fill_manual(values = c("indianred4", "lightblue3", "darkslategrey")) + 
  scale_color_viridis_d(option = "E", end = 0.9) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.25))) + 
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12))
  
  return(g)
}

### missense variants only ###
Variants_VEP_tidy_missense = Variants_VEP_tidy %>% 
  filter(Consequence == "missense_variant")

### Null variants only ###
Variants_VEP_tidy_null = Variants_VEP_tidy %>% 
  filter(Null_variant == T)

### plot data ###
g_density = plot_variant_density(Variants_VEP_tidy)
g_density_missense = plot_variant_density(Variants_VEP_tidy_missense, label_prefix = "Missense")
g_density_null = plot_variant_density(Variants_VEP_tidy_null, label_prefix = "Null")
```

Now lets have a look at the plots:

#### All variants

![](Readme_files/figure-gfm/Print%20all%20variants-1.png)<!-- -->

#### Missense variants

![](Readme_files/figure-gfm/Print%20missense%20variants-1.png)<!-- -->

#### Null variants

![](Readme_files/figure-gfm/Print%20null%20variants-1.png)<!-- -->

### Venn diagramm

With a Venn diagramm we can visualise which variants are present in
which dataset

``` r
x = list(gnomAD = (Variants_VEP_tidy %>% filter(!is.na(gnomAD)))$ID, 
         Cosmic = (Variants_VEP_tidy %>% filter(!is.na(Cosmic)))$ID, 
         ClinVar = (Variants_VEP_tidy %>% filter(!is.na(ClinVar)))$ID)

v = ggvenn(x, fill_color = c("darkslategrey", "lightblue3", "indianred4")) + 
  theme(plot.background = element_rect(fill = "white"))
```

![](Readme_files/figure-gfm/Print%20Venn%20diagramm-1.png)<!-- -->

### Domain data

We can also incorporate the domain data. This way we can determine if
specific variant types are more frequent in different domains.

#### Prepare domain dataset

First we need to prepare the dataset by assigning the correct domain to
every variant

``` r
### join the domain data ###
# with the function fuzzy_join we can assign the correct domain based on the variant position
Variants_VEP_tidy_domains = Variants_VEP_tidy %>% 
  fuzzy_left_join(domains, by = c("Protein_position2" = "Start", "Protein_position2" = "End"), 
             match_fun = list(`>=`, `<=`)) %>% 
  select(-Start, -End, -mid)

### melt dataset into long format ###
Variants_VEP_tidy_domains_melted = Variants_VEP_tidy_domains %>%
  select(Protein_position2, gnomAD, Cosmic, ClinVar, Consequence, Domain) %>%
  gather("Dataset", "TrueFalse", gnomAD:ClinVar) %>%
  filter(TrueFalse == T) %>%
  select(-TrueFalse) %>% 
  mutate(Consequence = str_replace_all(Consequence, "_", " ")) %>% 
  mutate(Domain = str_replace_na(Domain, replacement = "unspecified"))

Variants_VEP_tidy_domains_melted$Domain = factor(Variants_VEP_tidy_domains_melted$Domain, 
                                                       levels = c(as.vector(domains$Domain), "unspecified"))

### count number of variant type per domain and dataset ###
Variants_VEP_tidy_domains_melted = Variants_VEP_tidy_domains_melted %>% 
  group_by(Domain, Dataset) %>% 
  mutate(Number_of_Variants = n()) %>% 
  ungroup()
```

#### Plot domain diagramms

This dataset can now be used to plot the variant types per domain. We
will use pie charts and split them by domain and dataset, and annotate
with the variant consequence.

``` r
g_domain = ggplot(Variants_VEP_tidy_domains_melted, aes(fill = Consequence, x = "")) + 
  geom_bar(position = "fill") + 
  facet_grid(cols = vars(Dataset), 
             rows = vars(Domain), 
             switch = "y") + 
  coord_polar(theta = "y") + 
  theme_minimal() + 
  scale_y_continuous(labels = NULL, breaks = NULL, position = "bottom") +
  labs(x = "Domain", y = "Dataset", fill = "Variant Consequence") + 
  scale_fill_viridis_d() + 
  theme(strip.text = element_text(size = 12), 
        panel.grid = element_blank(), 
        plot.background = element_rect(fill = "white"), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12)) + 
  geom_text(x = 0.5, y = 0, aes(label = Number_of_Variants), 
            vjust = 0.5, hjust = 0.5, color = "white")
```

![](Readme_files/figure-gfm/Print%20domain%20diagramms-1.png)<!-- -->

## Statistics

Here we validate our observations with statistical analyses. In
particular, we test for a uniform variant distribution with the
Kolmogorov Smirnov Test

``` r
##### Kolmogorov Smirnov Test for uniform distribution #####

plot_ks = function(input_vector){
  x = input_vector
  dd = data.frame(x)
  ed = ecdf(dd$x)
  maxdiffidx <- which.max(abs(ed(dd$x)-punif(dd$x,0,1469)))
  maxdiffat <- dd$x[maxdiffidx]
  p<-ggplot(aes(x),data=dd)+
    stat_ecdf()+
    theme_bw()+
    stat_function(fun=punif,args=list(0,1469))
  return(p)
}

# per Base
# 16:9,847,261-10,276,785

test = Variants_VEP_tidy %>% 
  select(ID, POS, CDS_position, ClinVar, gnomAD, Cosmic, Protein_position2) %>% 
  mutate(CDS_position = as.numeric(str_extract(CDS_position, "[0-9]+")))

test_cosmic_vector = (test %>% 
  filter(Cosmic == T))$CDS_position
test_clinvar_vector = (test %>% 
  filter(ClinVar == T))$CDS_position
test_gnomad_vector = (test %>% 
  filter(gnomAD == T))$CDS_position

#ks.test(test_cosmic_vector, "punif", 1,4407)
plot_ks = function(input_vector, title){
  x = input_vector
  dd = data.frame(x)
  ed = ecdf(dd$x)
  ks_output = ks.test(input_vector, "punif", 1,4407)
  p_value = ks_output$p.value
  D = round(ks_output$statistic[["D"]], digits = 5)
  maxdiffidx <- which.max(abs(ed(dd$x)-punif(dd$x,1,4407)))
  maxdiffat <- dd$x[maxdiffidx]
  p<-ggplot(aes(x),data=dd)+
    stat_ecdf()+
    theme_minimal()+
    stat_function(fun=punif,args=list(1,4407)) + 
    labs(title = title) + 
    geom_label(x = 3500, y = 0.2, label = paste0("D = ", D)) +
    geom_label(x = 3500, y = 0.1, label = paste0("p = ", p_value))
  return(p)
}

cosmic_plot = plot_ks(input_vector = test_cosmic_vector, title = "Cosmic")
clinvar_plot = plot_ks(input_vector = test_clinvar_vector, title = "ClinVar")
gnomad_plot = plot_ks(input_vector = test_gnomad_vector, title = "gnomAD")

combi = plot_grid(cosmic_plot, 
                  clinvar_plot, 
                  gnomad_plot, 
                  labels = "AUTO", 
                  ncol = 3) + 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
```

![](Readme_files/figure-gfm/Print%20Statistics%20diagramms-1.png)<!-- -->
