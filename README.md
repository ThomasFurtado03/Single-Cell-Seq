# Single-Cell-Seq

### INTRODUCTION

Influenza viruses are a major cause of acute respiratory infection, affecting a large proportion of the global population each year and contributing largely to morbidity and mortality (Dumm et al., 2019). Infection occurs primarily in the respiratory epithelium, where viral replication, immune cell infiltration, and inflammatory signaling all disrupt tissue structure and function. The airway epithelium serves as a first line of defense against inhaled pathogens, and its interaction with immune cells plays a critical role in determining disease progression and recovery. In addition to epithelial responses, adaptive immune cells such as B cells are essential for viral removal, and long term immunity via production of antibodies (Mathew et al., 2021; MacLean et al., 2022).

A central component of the early immune response to influenza infection is the activation of the interferon signaling pathways. Interferon-stimulated genes, including members of the IFITM family, are rapidly induced and function to inhibit viral entry and replication (Brass et al., 2009). These antiviral mechanisms are essential for controlling infection, as disruption of IFITM function has been shown to increase disease severity and viral strength. (Everitt et al., 2012). Additional interferon-responsive genes, such as *Ifi27l2a*, are also upregulated during early infection and are associated with immune cell activation and infiltration, particularly in macrophages and lymphocytes (Tantawy et al., 2014). The timing and magnitude of these responses are critical, as early interferon signaling alters downstream adaptive immune processes, including antigen presentation and lymphocyte activation (Rowe et al., 2024).

Despite extensive knowledge of influenza immunology, the complexity of cellular responses within infected tissues is not completely understood. Traditional bulk RNA sequencing approaches average gene expression across large populations of cells. This obscures the heterogeneity of individual cell types and states. In contrast, single-cell RNA sequencing allows for resolution of transcriptional differences at the level of individual cells, allowing for the identificiation of distinct cell populations and changes in gene expression over time. Recent studies have demonstrated that both epithelial and immune cells exhibit highly heterogenous responses to viral infection within the respiratory tract, which signifies the importance of cellular context in shaping immune outcomes (Wang et al., 2025).

However, single-cell approaches also present certain limitations. These methods can be computationally intensive and sensitive to parameter selection, including clustering resolution and dimensionality reduction techniques, which can influence the identification or interpretation of cell populations. In addition, differential expression analysis at the single-cell level may be affected by changes in cell composition, which may make it challenging to distinguish between true transcriptional regulation and population abundance (Tantawy et al., 2014). These trade-offs must be considered while interpreting results from single-cell datasets.

The objective of this study is to characterize transcriptional changes associated with influenza infection, at single-cell resolution. This will be done with a particular focus on identifying cell populations that exhibit dynamic responses across infection time points. This study will be using clustering, differential expression analysis, and functional enrichment. Through these methods, this study aims to identify distinct cell populations within the respiratory tract, determine how gene expression changes across infection stages, and relate these changes to known biological processes regarding antiviral defense and immune activation.

### METHODS 

Single-cell RNA sequencing data were provided as a pre-processed Seurat object containing gene expression cunts and associated metadata, including tissue origin (*organ_custom*), infection time point (*time*), and bio sample identifiers (*biosample_id*). Due to incomplete entries in the *mouse_id* field, this variable was not sued for downstream analysis. Additional quality control filtering was not performed, as the dataset was provided in a pre-processed format.

All analyses were conducted in R (version 4.5.1) using the Seurat package (version 5.4.0) (Stuart et al., 2019). Highly variable genes were identified using the FindVariableFeatures() function with default parameters. To remedy and reduce computational burden, scaling was restricted to the identified variable using ScaleData(). Principal component analysis (PCA) was performed using RunPCA() on the variable gene set to reduce dimensionality.

A shared nearest neighbor graph was constructed using FindNeighbors() with the first 12 principal components, selected based on inspection of the elbow plot. Clustering was performed using FindClusters() with a resolution parameter of 0.3 to identify transcriptionally distinct cell populations. Uniform Manifold Approximation (UMAP) was applied using RunUMAP() with the same principal components to visualize the data in two dimensions (McInnes et al., 2018).

Cluster composition across experimental conditions was assessed by tabulating cluster membership across infection time points, Cluster 5 was selected for downstream analysis due to its differential distribution across time.

Differential expression analysis was performed within Cluster 5 using Seurat's FindMarkers() function, which compared cells from D14 and D05 time points. Parameters included a minimum expression threshold (min.pct = 0.1) and a log fold-change threshold (logfc.threshold = 0.25). Ribosomal protein genes (Rps/Rpl) were later filtered form downstream interpretation due to their ubiquitous expression, and lack of biological specificity.

Gene ontology enrichment analysis was performed using the clusterProfiler package (version 4.18.4) (Yu et al., 2012). Differentially expressed genes were analyzed using the enrichGO() function with biological process ontology (ont = "BP"), with gene annotations obtained from the *org.Mm.eg.db* database (Bioconductor). Enrichment results were visualized using dot plots.

Marker gene identification for cluster annotation was performed using a down sampled subset of cells to reduce computational strain. Cells from Cluster 5 and a representative subset of other clusters were compared using FindMarkers() with positive markers only. The resulting gene lists were filtered to remove ribosomal and mitochondrial genes prior to interpretation.


### RESULTS

![Figure 1](Figures/UMAP1.png)

**Figure 1. UMAP visualization of clustered single-cell transcriptomic data.** UMAP projection of all cells showing 23 transcriptionally distinct clusters identified using Seurat clustering. Each point represents a single cell, colored by cluster identity. The clear separation of clusters indicates distinct gene expression profiles across cell populations.
<br>

![Figure 2](Figures/UMAP_organcustom.png)

**Figure 2. Distribution of cells across tissue types.** UMAP projection colored by tissue of origin: respiratory mucosa (RM) , olfactory mucosa (OM) , and lateral nasal gland (LNG). Cells from different tissues are distributed across multiple clusters, with some regions showing tissue-specific enrichment, indicating both shared and tissue-specific cellular populations.
<br>

![Figure 3](Figures/UMAP_time.png)

**Figure 3 Distribution of cells across infection time points.** UMAP projection colored by infection stage (Naive, D02, D05, D08, D14). Cells from different time points are broadly distributed across clusters, with localized enrichment in specific regions, suggesting dynamic transcriptional changes during infection progression.

<br><br>

Single-cell RNA sequencing data were analyzed to characterize cellular differences and infection-associated changes in transcription across nasal tissues in mice. Dimensionality reduction and clustering identified 23 distinct cell populations, visualized using UMAP (Figure 1). Clusters were well separated, indicating distinct transcriptional profiles across cell populations. Overlay of metadata revealed that cells from different tissues (RM, OM, LNG) and different time points were broadly distributed across clusters. Although, certain regions showed enrichment for specific conditions, suggesting that there are some condition-specific cellular states observed. (Figures 2-3).
<br>

To identify populations associated with the dynamics of infection, cluster composition across time points was examined. Cluster 5 displayed a notable shift in abundance between the early (D05) and late (D14) infection stages, with a significant amount of transcription at the early stage. This pattern suggested that Cluster 5 may be involved in infection-related responses and was therefore selected for downstream analysis.
<br>
### Table 1. Differentially expressed genes in Cluster 5 (D14 vs D05)

| Gene    | p_val   | avg_log2FC | pct.1 | pct.2 | p_val_adj |
|---------|--------|------------|-------|-------|-----------|
| Obp1a   | 2.74E-98 | -310.216  | 0.278 | 0.731 | 6.88E-94 |
| Ifitm2  | 5.87E-40 | -3.012    | 0.117 | 0.341 | 1.48E-35 |
| Ifi213  | 3.27E-37 | -4.889    | 0.025 | 0.150 | 8.22E-33 |
| Ifi27l2a| 1.73E-35 | -4.236    | 0.391 | 0.648 | 4.35E-31 |
| Lef1    | 7.60E-35 | -7.226    | 0.019 | 0.128 | 1.91E-30 |
| Il7r    | 1.53E-34 | -1.839    | 0.027 | 0.150 | 3.84E-30 |
| Map1b   | 1.08E-33 | -10.133   | 0.850 | 0.585 | 2.73E-29 |
| H2-Aa   | 2.42E-33 | 6.331     | 0.910 | 0.762 | 6.07E-29 |
| Vpreb1  | 5.86E-33 | -14.464   | 0.012 | 0.101 | 1.47E-28 |
| Pafah1b3| 4.85E-32 | -6.104    | 0.122 | 0.327 | 1.22E-27 |

[View full Table 1 (CSV)](DE_results.csv)

<br>

![Figure 4](Figures/Cluster5_selectgenes.png)

**Figure 4. Expression of select differentially expressed genes in Cluster 5.** Feature plots showing expression of representative genes (*Ifitm2*, *Ifi27l2a*, *Il7r*, and *Obp1a*) within Cluster 5. Expression is localized to subsets of cells, indicating transcriptional heterogeneity within this population and supporting differential expression results across infection stages.
<br><br>

Differential expression analysis within Cluster 5 comparing D14 and D05 revealed notable transcriptional changes (Table 1). Several genes associated with antiviral responses were seen to be significantly enriched at the earlier time points, including interferon-stimulated genes such as *Iftim2* and *Ifi27l2a*. Additional immune-related genes including *Il7r* and *Lef1*, were also differentially expressed, indicating involvement of lymphocyte-associated processes. Visualization of selected genes using FeaturePlot demonstrated that their expression was localized to subsets of cells with Cluster 5, indicating heterogeneity in this population (Figure 4).

<br>

![Figure 5](Figures/ORA_dotplot.png)

**Figure 5. Gene ontology enrichment analysis of differentially expressed genes.** Dot plot showing enriched biological processes for genes differentially expressed between D05 and D14 within Cluster 5. Dot size represents gene count and color indicates adjusted p-value. Enriched pathways include antiviral response, innate immune activation, lymphocyte differentiation, and antigen processing, consistent with an immune response to infection.
<br><br>

Functional enrichment analysis of differentially expressed genes further supported these observations. Gene ontology analysis revealed significant enrichment of pathway related to antiviral response, innate immune activation, lymphocyte differentiation, and antigen processing and presentation (Figure 5). These results strongly indicate that cluster 5 is involved in coordinated immune responses during infection.

<br>

### Table 2. Top 15 marker genes for Cluster 5 with proportional expression values between time points

| Gene    | log2FC | pct_D14 | pct_D05 |
|---------|--------|---------|---------|
| Igkc    | 166.40 | 0.92    | 0.02    |
| Iglc1   | 133.67 | 0.36    | 0.00    |
| Iglc3   | 99.63  | 0.61    | 0.01    |
| Iglc2   | 74.13  | 0.65    | 0.00    |
| Dnajc7  | 45.47  | 0.41    | 0.40    |
| Ighm    | 42.22  | 0.92    | 0.11    |
| Vpreb3  | 38.07  | 0.52    | 0.03    |
| Ms4a1   | 31.45  | 0.68    | 0.00    |
| Fau     | 27.97  | 1.00    | 1.00    |
| Jund    | 23.64  | 0.99    | 0.92    |
| Ighd    | 20.36  | 0.56    | 0.00    |
| Cd79b   | 17.89  | 0.87    | 0.03    |
| Smarca4 | 17.61  | 0.21    | 0.32    |
| Cd37    | 17.47  | 0.79    | 0.08    |
| Fcmr    | 16.73  | 0.53    | 0.00    |

[View full Table 2 (CSV)](cluster5_top15_markers.csv)
<br>

To further elucidate the identity of this population, marker gene analysis was performed. Cluster 5 displayed strong expression of canonical B-cell markers. These included immunoglobulin genes (*Igkc*, *Ighm*, *Ighd*) and B-cell receptor components (Cd79b, Ms4a1). These markers were highly enriched relative to other cell populations, confirming that Cluster 5 represents a B-cell population. A summary of the top marker genes is provided in Table 2.


### DISCUSSION

The results of this study demonstrate that distinct immune cell populations, in this case B-cells, undergo dynamic transcriptional changes across the course of influenza infection. The identification of Cluster 5 as a B-cell population is strongly supported by the enrichment of known B-cell markers, including immunoglobulin genes and B-cell receptor components. These findings are consistent with previous single-cell studies showing that influenza infection drives the emergence of antigen-specific B-cell populations with distinct transcriptional profiles that evolve over time (Mathew et al., 2021). Furthermore, lung-resident memory B-cells have been shown to rapidly respond to infection by migrating to site of viral exposure and differentiating into antibody-secreting cells, highlighting their role in localized immune defense (MacLean et al., 2022).

A key observation in this study is the temporal shift in gene expression within Cluster 5. Early time points were characterized by strong upregulation of interferon-stimulated genes (such as *Iftm2* and *Ifi27l2a*), indicating activation of innate immunity pathways. IFITM family proteins are well-established mediators of intrinsic antiviral defense, which function to restrict viral entry and replication across multiple viruses, including influenza (Brass et al., 2009). In vivo studies further demonstrate that disruption of IFITM function leads increased viral effects and severity, which emphasizes the importance of these pathways (Everitt et al., 2012). The presence of these genes in early-stage cells suggests that B cells (or closely related immune populations) participate in antiviral defense beyond just their classic antibody-producing roles.

The expression pattern of *Ifi27l2a* further supports a temporally regulated immune response. This gene has been shown to be strongly induced during early influenza infection, particularly in macrophages and lymphocytes, accompanied by a depletion in expression by the later stages (Tantawy et al., 2014). This aligns closely with the observed differential expression between D05 and D14, suggesting that the transcriptional differences identified in this study actually reflect progression through stages of infection rather than static cell-type differences. Notably, previous work has shown that increased expression of such genes may reflect changes in immune cell composition rather than solely transcriptional regulation within individual cells (Tantawy et al., 2014).

The broader interferon response observed in this study is consistent with the established models of influenza infection in the respiratory tract specifically. Early interferon signaling plays a crucial role in activating downstream immune pathways, including antigen presentation and adaptive immune activation (Rowe et al., 2024). Differences in the magnitude and timing of interferon responses have been shown to significantly influence disease progression and immune outcomes, For example, influenza A virus induces stronger and earlier interferon and inflammatory responses compared to the influenza B virus, which can delay adaptive immune activation (Rowe et al., 2024). These findings support the interpretation that early antiviral gene expression in this dataset reflects am active innate immune response which later shapes into adaptive immunity.

In addition to immune cell dynamics, the role of epithelial cells in shaping the immune response must also be taken into consideration. The nasal and upper respiratory epithelium represents the primary stie of infection and is known to exhibit stored antiviral responses at the single-cell level (Wang et al., 2025). These tissue-specific responses ay contribute to the distribution of cells observed across clusters and explain why cells from multiple tissues are dispersed within transcriptionally defined populations. Epithelial cells can also survive viral infection and undergo long-term transcriptional reprogramming, which may influence immune signalling and tissue recovery (Dumm et al., 2019). These findings highlight the importance of considering both immune and non-immune cell contributions when looking at transcriptional changes.

The enrichment of pathways related to antigen processing, lymphocyte activation, and immune signaling suggests that Cluster 5 plays an active role in coordinating immune responses during infection. The transition from strong interferon-driven responses at earlier time points to more adaptive signature later on is consistent with the known progression of influenza infection. Which is described by innate immunity initiating viral control, and adaptive immunity mediates clearance and long-term protection. The presence of B-cell markers alongside interferon-stimulated genes indicates that these processes may occur simultaneously within specific cell populations or reflect interactions between closely associated immune cells.


