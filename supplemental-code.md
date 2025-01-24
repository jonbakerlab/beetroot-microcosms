# daniela-AHS-beetroot-project Supplemental Code
# Preliminary Steps
From the nanopore sequencing, separate feature tables (OTU tables) for each flow cell run were generated. To merge the feature tables together, the `feature-table merge` function in the qiime2 conda environment was used.

## Convert .txt files to .qza files
Before merging the feature tables, first they must be converted the .txt tables into compatible .qza files:

```sh
#!/bin/bash
# first convert featureTable into a .biom format:
conda activate qiime2-amplicon-2024.5
biom convert -i featureTable.txt -o featureTable.biom --to-hdf5

# import the featureTable into qiime to get a .qza file:
qiime tools import \
 --type 'FeatureTable[Frequency]' \
 --input-path merged-featureTable.biom \
 --output-path merged-featureTable.qza    
```

## Merge the feature tables
With the .qza featureTables, the following command can be run to merge the featureTables together:

```sh
#!/bin/bash
# merge feature tables 1 - N together:
qiime feature-table merge \
    --i-tables featureTable1.qza featureTable2.qza ...featureTableN.qza \
    --p-overlap-method sum \
    --o-merged-table merged-featureTable.qza
```

## Alpha rarefaction curve
Now that all of the featureTables have been combined into `merged-featureTable.qza`, the following code generates an alpha rarefaction curve visual with the qiime2 environment's `diversity alpha-rarefaction` command. This is helpful to visualize and ensure the all the taxa present has been observed in each sample.

```sh
#!/bin/bash
qiime diversity alpha-rarefaction \
    --i-table merged-featureTable.qza \
    --p-max-depth 60000 \
    --p-steps 40 \
    --m-metadata-file metadata.txt \
    --o-visualization alpha-rarefaction.qzv
```
## Core metrics generation
In addition, the following code uses the `diversity core-metrics` function in qiime2's environment to generate a variety of diversity metric artifacts, a starting point for downstream analysis. Of particular importance is the `rarefied_table.qza` artifact. This is a normalized version of the `merged-featureTable.qza` file dictated by the parameter `--p-sampling-depth`. This rarefies the abundances present in each sample to the number inputted (in this case 4400 reads) so samples with a higher number of reads don't skew the data away from those with a lower number of reads.

**NOTE: `gemelli` is a plugin downloaded and used in downstream analyses (RCPA and CTF) however, downloading the gemelli plugin downgrades the scikit-bio package to v0.5.9 from the base qiime scikit-bio v0.6.0. This downgraded version of scikit-bio is incompatible with the `diversity core-metrics` command specifically.**

```sh
#!/bin/bash
qiime diversity core-metrics \
    --i-table ./feature-tables/merged-featureTable.qza \
    --p-sampling-depth 4400 \
    --m-metadata-file metadata.txt \
    --output-dir core-metrics-results
```
`diversity core-metrics` exports the following artifact files:

* `core-metrics-results/bray_curtis_pcoa_results.qza`
* `core-metrics-results/shannon_vector.qza`
* `core-metrics-results/rarefied_table.qza`
* `core-metrics-results/jaccard_pcoa_results.qza`
* `core-metrics-results/observed_features_vector.qza`
* `core-metrics-results/jaccard_distance_matrix.qza`
* `core-metrics-results/evenness_vector.qza`
* `core-metrics-results/bray_curtis_distance_matrix.qza`

As well as the following visual files:
* `core-metrics-results/jaccard_emperor.qzv`
* `core-metrics-results/bray_curtis_emperor.qzv`

## Phylogenetic Tree Generation
A phylogenetic tree artifact is required to for the Weighted and Unweighted Unifrac metrics, therefore, before starting, a properly formatted `.fasta` file must be extracted from the 16S nanopore sequencing data:

### Downloading and formatting the RefSeq 16S and 18S database
The EPI2ME 16S workflow by default uses the RefSeq 16S/18S database.  This can be downloaded by going to https://www.ncbi.nlm.nih.gov/refseq/targetedloci/ and clicking on '16S RefSeq Nucleotide sequence records'. Then click 'Send to:' and select 'Complete Record' 'File' Format: FASTA Sort By: Organism Name and Create File. The created file will have a bunch of unecessary newlines, which can be removed with the following uses the seqkit toolkit (available as a conda environment):

```sh
seqkit seq -w 0 sequence.fasta > sequence-edit.fasta
```

### Extract 16S sequences from the database
First, we will want to get our list of features that we want to extract 16S sequences of. You can just copy the contents of the 'species' column of the feature table if the contents are in this format "Streptococcus mutans". You will need to remove the 'Unknown' row from this, as well as remove it from the feature table in the downstream phylogenetic analyses. Copy the list of taxa names into a .txt file, we can call it patterns.txt. You can use the following to then extract 16S sequences of interest from the database.

```sh
while IFS= read -r pattern; do     
    grep -m 1 -A 1 "$pattern" sequence-edit.fasta >> output.fasta
done < patterns.txt
```

Now you can simplify the FASTA deflines:

```sh
sed -E 's/^>([^ ]+) ([^ ]+) ([^ ]+).*/>\2_\3/' output.fasta > output-simple.fasta
```
### Import to QIIME and generate your trees
With the properly formatted .fasta file, it can now be imported into qiime:

```sh
#!/bin/bash
# import .fasta file:
qiime tools import \
    --type FeatureData[Sequence] \
    --input-path output-simple.fasta \
    --output-path 16S-seqs.qza

# generate our trees:
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 16S-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```
### Filtering the featureTable and rooted phylogenetic tree
Through some trouble shooting, it was found that there were feature IDs in the feature tables that did not have a corresponding leaf in the rooted tree file. Fortunately, these missing species were relatively low in abundance. These species were filtered out based on a frequency count of reads less than 200. This was done using the `feature-table filter-feature` function: 

```sh
#!/bin/bash
# filter out features with < 200 reads:
qiime feature-table filter-features \
    --i-table ./feature-tables/rarefied_table.qza \
    --p-min-frequency 200 \
    --o-filtered-table ./feature-tables/filtered-rarefied-featureTable.qza
```

With the new `filtered-featureTable.qza` file, the `rooted-tree.qza` was also filtered with `phylogeny filter-tree` based on the remaining species present in the feature table:

```sh
#!/bin/bash
# filter tree based on remaining features:
qiime phylogeny filter-tree \
  --i-tree ./16s-phylogeny/rooted-tree.qza \
  --i-table ./feature_tables/filtered-rarefied-featureTable.qza \
  --o-filtered-tree ./16s-phylogeny/filtered-rarefied-tree.qza
```

# Alpha Diversity Analysis
## Vectors of alpha diversity
Before performing the alpha diversity analyses, vector files for all the methods of alpha diversity analysis must be created. The shannon, chao1, and simpson metrics are all of interest, however, the vector for the shannon metric was created and exported using the `core-metrics` function. To create vectors for chao1 and simpson, the `diversity alpha` function was used:

```sh
#!/bin/bash
# simpson metric:
qiime diversity alpha \
    --i-table ./feature-tables/rarefied_table.qza \
    --p-metric simpson \
    --o-alpha-diversity ./core-metrics-results/simpson_vector.qza

# chao1 metric
qiime diversity alpha \
    --i-table ./feature-tables/rarefied_table.qza \
    --p-metric chao1 \
    --o-alpha-diversity ./core-metrics-results/chao1_vector.qza 
```
### Vectors with phylogenetic trees
For vectors associated with a phylogenetic tree, its a similar process using the qiime `diversity alpha-phylogeny` function:

```sh
#!/bin/bash
# faith_pd metric:
qiime diversity alpha-phylogenetic \
    --i-table ./feature-tables/filtered-rarefied-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-rarefied-tree.qza \
    --p-metric faith_pd \
    --o-alpha-diversity ./core-metrics-results/faith_pd_vector.qza
```

## Alpha diversity group significance
With the alpha diversity vector artifacts, a group significance analysis using the `diversity alpha-group-significance` function can be run for each of the metrics, resepectively:

```sh
#!/bin/bash
# shannon metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./core-metrics-results/shannon_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha-diversity/shannon-group-significance.qzv

# simpson metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./core-metrics-results/simpson_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha-diversity/simpson-group-significance.qzv

# chao1 metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./core-metrics-results/chao1_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha-diversity/chao1-group-significance.qzv

# faith_pd metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha-diversity/faith-group-significance.qzv
```

# Beta Diversity Analysis
## Bray Curtis emperor plot qnd permanova
Using the `bray_curtis_distance_matrix.qza` artifact previously exported from the `diversity core-metrics` function, a permanova analysis can be conducted:

```sh
#!/bin/bash
# generate beta diversity permanova
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results//bray_curtis_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./beta_diversity/bray_curtis_permanova.qzv
```

##  Jaccard emperor plot and permanova
An identical process was done for the Jaccard metric by substituting the `jaccard` argument for `--p-metric` in the `diversity beta-group-significance` function:
```sh
#!/bin/bash
# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./beta_diversity/jaccard_permanova.qzv
```

## Weighted and Unweighted Unifrac emperor plots and permanova
### Weighted Unifrac emperor plot and permanova
Since the `filtered-rarefied-featureTable.qza` and `filtered-rarefied-tree.qza` have the same corresponding feature ID's, a similar sequence of codes as the Bray Curtis and Jaccard beta diversity analyses can be used for the weighted unifrac analysis using `diversity beta-phylogenetic`:

```sh
#!/bin/bash
# distance matrix:
qiime diversity beta-phylogenetic \
    --i-table ./feature-tables/filtered-rarefied-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-rarefied-tree.qza \
    --p-metric weighted_unifrac \
    --o-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza

# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization  ./beta-diversity/weighted_unifrac_permanova.qzv

# pcoa:
qiime diversity pcoa \
    --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --o-pcoa ./core-metrics-results/weighted_unifrac_pcoa.qza

# emperor plot:
qiime emperor plot \
    --i-pcoa ./core-metrics-results/weighted_unifrac_pcoa.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./beta-diversity/weighted_unifrac_emperor.qzv
```
### Unweighted Unifrac emperor plot and permanova
And similarly for the Unweighted Unifrac analysis by using the `--p-metric unweighted_unifrac` argument:

```sh
#!/bin/bash
# distance matrix:
qiime diversity beta-phylogenetic \
    --i-table ./feature-tables/filtered-rarefied-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-rarefied-tree.qza \
    --p-metric unweighted_unifrac \
    --o-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza

# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization  ./beta-diversity/unweighted_unifrac_permanova.qzv

# pcoa:
qiime diversity pcoa \
    --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --o-pcoa ./core-metrics-results/unweighted_unifrac_pcoa.qza

# emperor plot:
qiime emperor plot \
    --i-pcoa ./core-metrics-results/unweighted_unifrac_pcoa.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./beta-diversity/unweighted_unifrac_emperor.qzv
```

Finally, it is worth noting that the function `diversity core-metrics-phylogenetic` is would provide the results from `diversity core-metrics` in addition to the Weighted and Unweighted Unifrac artifacts. However, due to the order in which these were conducted (creating the rarefied feature table) the phylogenetic analyses were done independently.

# Gemelli Plugin Analyses
The `gemelli` plugin provides other methods for visualizing significant differences among groups: Robust Aitchison PCA (RPCA) and Compositional Tensor Factorization (CTF).

## Robust Aitchison PCA (RPCA)

RPCA is another metric for investigating beta diversity which is done through use of the `gemelli` plugin. A separate qiime enivronment was created specifically for the `gemelli` due to the aforementioned version differences. 

**Note: For both the RPCA and phylogenetic RPCA (as well as the CTF analyses below) use an in-house normalization method for a feature table. Therefore, in the following sections, the overall `merged-featureTable.qza` was used.**

```sh
#!/bin/bash
conda activate qiime2-gemelli_env
```
This process is similar to the beta diversity sequence of codes, just using the `gemelli rpca` function in the `gemelli` plugin:

```sh
#!/bin/bash
# use gemelli plugin to create ordination and distance matrix:
qiime gemelli rpca \
    --i-table merged-featureTable.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot ./core-metrics-results/ordination.qza \
    --o-distance-matrix ./core-metrics-results/rcpa_distance.qza
```
The exported files-

* `ordination.qza`
* `rcpa_distance.qza`

-provide the pcoa and distance matrix artifacts, respectively, to use for the `emperor biplot` and `diversity beta-grou-significance` functions in the base qiime environment:

```sh
# generate rpca emperor biplot:
qiime emperor biplot \
    --i-biplot ./core-metrics-results/ordination.qza \
    --m-sample-metadata-file metadata.txt \
    --o-visualization ./beta_diversity/rpca_biplot.qzv \
    --p-number-of-features 8

# rpca permanova analysis:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/rcpa_distance.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./beta_diversity/rcpa_permanova.qzv
```
### Interlude: Creating a taxonomy metadata file
In order to do a phylogenetic RCPA analysis, a taxonomy table is required as one of the input files. In order to create this, the raw data from the sequencer contains taxonomic information for each level, i.e., kingdom, order, phylum, genus, and species. taking each level individually, they were copied, pasted, and adjusted into the right format required by the `FeatureData[Taxonomy]` file type. For example, kingdom would be 'k__bacteria; '. Each of the levels were then concatonated into a single string associated with the species and converted into a .tsv file. Now, the taxonomy file is ready to be imported as an artifact into qiime:

```sh
#!/bin/bash
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path filtered-taxonomy.tsv \
    --output-path filtered-taxonomy.qza
```

To view that the filtered-taxonomy.qza was imported correctly, a .qzv can be created using `metadata tabulate` for visual verification:

```sh
#!/bin/bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
**NOTE: Since the filtered featureTable and filtered phylogenetic tree files are being used, the taxonomy must *only* contain the species present in both.**

### Phylogenetic RPCA analysis
Prior to starting this analysis, a new filtered phylogenetic tree must be created. This is because, since the `merged-featureTable.qza` is being inputted, a phylogenetic tree must be filtered based on this specific feature table rather than the `rarefied_table.qza`. Luckily, this process is the exact same: filter `merged-featureTable.qza` and filter `rooted-tree.qza`.

```sh
#!/bin/bash
# filter out features with < 200 reads:
qiime feature-table filter-features \
    --i-table ./feature-tables/merged-featureTable.qza \
    --p-min-frequency 200 \
    --o-filtered-table ./feature-tables/filtered-merged-featureTable.qza

# filter tree based on remaining features:
qiime phylogeny filter-tree \
  --i-tree ./16s-phylogeny/rooted-tree.qza \
  --i-table ./feature_tables/filtered-merged-featureTable.qza \
  --o-filtered-tree ./16s-phylogeny/filtered-merged-tree.qza
```

With the newly generated filtered artifacts and the previously formatted taxonomy artifact, all the required files are created and can be inputted as arguments in the `gemelli phylogenetic-rcpa-with-taxonomy` function:

```sh
#!/bin/bash
# create all the metrics needed for downstream analysis:
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table ./feature-tables/filtered-merged-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-merged-tree.qza \
    --m-taxonomy-file ./16s-phylogeny/taxonomy/filtered-taxonomy.qza \
    --o-biplot ./core-metrics-results/phylo-ordination.qza \
    --o-distance-matrix ./core-metrics-results/phylo-distance.qza \
    --o-counts-by-node-tree ./core-metrics-results/phylo-tree.qza \
    --o-counts-by-node ./core-metrics-results/phylo-table.qza \
    --o-t2t-taxonomy ./core-metrics-results/phylo-taxa.qza
```
This function outputs the following files:
* `phylo-ordination.qza`
* `phylo-distance.qza`
* `phylo-tree.qza`
* `phylo-table.qza`
* `phylo-taxa.qza`

These files are required to generate an emperor plot and a permanova, as previously done. However, instead of an emperor plot, with the phylogenetic and taxonomic data, an empress plot can be generated. `empress` is an additional plugin that displays a phylolegenetic tree data and an emperor biplot side-by-side in the `phylo-empress.qzv` output file; this can be done through the `empress community-plot` function as seen below:

```sh
#!/bin/bash
# generate an empress community plot:
qiime empress community-plot \
    --i-tree ./core-metrics-results/phylo-tree.qza \
    --i-feature-table ./core-metrics-results/phylo-table.qza \
    --i-pcoa ./core-metrics-results/phylo-ordination.qza \
    --m-sample-metadata-file ./metadata.txt \
    --m-feature-metadata-file ./core-metrics-results/phylo-taxa.qza \
    --p-filter-missing-features \
    --o-visualization ./beta_diversity/phylo-empress.qzv
```

The permanova can be created in the base qiime environment, as done previously with the RCPA:

```sh
#!/bin/bash
# phylogeny/taxonomy permanova analysis:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/phylo-distance.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column 'group' \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./beta_diversity/phylo-permanova.qzv
```

## Compositional Tensor Factorization (CTF)
As mentioned in the main document, CTF helps aleviate some of the subject-to-subject oral microbiome variability by comparing repeated measures of the same subject. Much like the RCPA analysis, the CTF analysis must be done in the `gemelli` environment.  In addition, like RCPA, the CTF workflow uses its own in-house normalizion methods, therefore, the feature table before rarefaction (`merged-featureTable.qza`) was used. This applies to both the CTF and phylogenetic CTF analyses.

```sh
#!/bin/bash
qiime gemelli ctf \
    --i-table merged-featureTable.qza \
    --m-sample-metadata-file metadata.txt \
    --p-individual-id-column subject \
    --p-state-column state \
    --output-dir ctf-results
```

The output directory contains the following files:
* `subject_biplot.qza`
* `state_biplot.qza`
* `distance_matrix.qza`
* `state_subject_ordination.qza`
* `state_feature_ordination.qza`

Using some of these output files, an emperor biplot can be generated using `emperor biplot` in the base qiime environment:

```sh
#!/bin/bash
qiime emperor biplot \
    --i-biplot ./ctf-results/subject_biplot.qza \
    --m-sample-metadata-file metadata.txt \
    --p-number-of-features 5 \
    --o-visualization subject_biplot.qvz
```

## Phylogenetic CTF analysis
Similarly to the phylogenetic RPCA analysis, the `filtered-taxonomy.qza` taxonomic metadata and `filtered-merged-tree.qza` phylogeny tree files are required.

```sh
#!/bin/bash
qiime gemelli phylogenetic-ctf-with-taxonomy \
    --i-table ./feature-tables/filtered-merged-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-merged-tree.qza \
    --m-sample-metadata-file metadata.txt \
    --m-taxonomy-file ./16s-phylogeny/taxonomy/filtered-taxonomy.qza \
    --p-state-column state \
    --p-individual-id-column subject \
    --output-dir ./phylo-ctf-results
```
The output directory contains the following files:
* `subject_biplot.qza`
* `state_biplot.qza`
* `distance_matrix.qza`
* `state_subject_ordination.qza`
* `state_feature_ordination.qza`
* `counts_by_node_tree.qza`
* `counts_by_node.qza`
* `t2t_taxonomy.qza`
* `subject_table.qza`

With the files above, an empress plot can be generated using `empress community-plot` since the data contains phylogenetic and taxonomic data:

```sh
# empress plot:
qiime empress community-plot \
    --i-tree ./phylo-ctf-results/counts_by_node_tree.qza \
    --i-feature-table ./phylo-ctf-results/subject_table.qza \
    --i-pcoa ./phylo-ctf-results/subject_biplot.qza \
    --m-sample-metadata-file metadata.txt \
    --m-feature-metadata-file ./phylo-ctf-results/t2t_taxonomy.qza \
    --p-filter-missing-features \
    --o-visualization ./phylo-ctf-results/phylo-subject-empress.qzv
```
# Subject Analysis
In addition to the alpha and beta diversity analysis on the different conditions, i.e. smoker/nonsmoker and saliva/growth, alpha and beta diversity analyses were conducted based on subjects response to the initial survey to tease out any correlations between subject metadata and oral microbiome. These include gender, age, coffee consumption, sugar consumption, alcohol consumption, and pH. All except gender were recorded in the metadata on a numerical scale, therefore, a different functions in qiime2 were used to test for alpha and beta diversity significance.

For these analyses, only the saliva samples were considered. This is due to the fact that those were the original samples provided by the subjects and no selectivity had occured during a growth process, as is the case for the microcosms. Hence, using the rarefied feature table as a basis, `saliva-featureTable.qza` and `subject-metadata.qza` were used as new feature tables and metadata tables, respectively, to exclude all microcosm data.

## Numerical alpha diversity
Much like with categorical metadata (done previously), an evenness vector can be generated using `qiime diversity alpha` for the Shannon, Simpson, and Chao1 metrics:

```sh
#!/bin/bash
qiime diversity alpha \
    --i-table ./feature-tables/saliva-featureTable.qza \
    --p-metric simpson \
    --o-alpha-diversity ./subject-analysis/metrics/simpson-saliva-vector.qza

qiime diversity alpha \
    --i-table ./feature-tables/saliva-featureTable.qza \
    --p-metric shannon \
    --o-alpha-diversity ./subject-analysis/metrics/shannon-saliva-vector.qza

qiime diversity alpha \
    --i-table ./feature-tables/saliva-featureTable.qza \
    --p-metric chao1 \
    --o-alpha-diversity ./subject-analysis/metrics/chao1-saliva-vector.qza
 ```

With the vector artifacts, rather than using `qiime diversity alpha-group-significance` for categorical metadata, numerical metadata used the `qiime alpha-correlation` function. Again, this was done below for each metric: Shannon, Simpson, Chao1.

```sh
#!/bin/bash
# simpson metric:
qiime diversity alpha-correlation \
    --i-alpha-diversity ./subject-analysis/metrics/simpson-saliva-vector.qza \
    --m-metadata-file  subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/simpson-numeric-group-significance.qzv

# shannon metric:
qiime diversity alpha-correlation \
    --i-alpha-diversity ./subject-analysis/metrics/shannon-saliva-vector.qza \
    --m-metadata-file  subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/shannon-numeric-group-significance.qzv

# chao1 metric:
qiime diversity alpha-correlation \
    --i-alpha-diversity ./subject-analysis/metrics/chao1-saliva-vector.qza \
    --m-metadata-file  subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/chao1-numeric-group-significance.qzv
```
The output files `-numeric-group-significance.qzv` is a visual artifact that also provides statistical significance information.

## Numerical beta diversity
Similarly to categorical metadata, a distance matrix artifact must be generated using `qiime diversity pcoa` for numerical data. As previously done, the distance matrix for each metric: Bray Curtis and Jaccard.

```sh
#!/bin/bash
# for Bray Curtis
qiime diversity beta \
    --i-table ./feature-tables/saliva-featureTable.qza \
    --p-metric braycurtis \
    --o-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza

# for Jaccard
qiime diversity beta \
    --i-table ./feature-tables/saliva-featureTable.qza \
    --p-metric jaccard \
    --o-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza
```

With the distance matrix artifacts, instead of using `qiime diversity beta-group-significance` for categorical metadata, numerical metadata significance is shown through a mantel plot generated via `qiime diversity beta-correlation`. This function is applied to both Bray Curtis and Jaccard as well as to each of the metadata groups listed above.

```sh
#!/bin/bash
# for Bray Curtis
# mantle plot for numerical 'pH' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column pH \
    --o-metadata-distance-matrix ./subject-analysis/metrics/bc-pH-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/bc-pH-mantel.qzv

# mantel plot for numerical 'age' metadata 
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column age \
    --o-metadata-distance-matrix ./subject-analysis/metrics/bc-age-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/bc-age-mantel.qzv

# mantel plot for numerical 'sugar' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column sugar \
    --o-metadata-distance-matrix ./subject-analysis/metrics/bc-sugar-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/bc-sugar-mantel.qzv

# mantel plot for numerical 'coffee' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column coffee \
    --o-metadata-distance-matrix ./subject-analysis/metrics/bc-coffee-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/bc-coffee-mantel.qzv

# mantel plot for numerical 'alcohol' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column alcohol \
    --o-metadata-distance-matrix ./subject-analysis/metrics/bc-alcohol-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/bc-alcohol-mantel.qzv

# for Jaccard
# mantle plot for numerical 'pH' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column pH \
    --o-metadata-distance-matrix ./subject-analysis/metrics/jaccard-pH-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/jaccard-pH-mantel.qzv

# mantel plot for numerical 'age' metadata 
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column age \
    --o-metadata-distance-matrix ./subject-analysis/metrics/jaccard-age-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/jaccard-age-mantel.qzv

# mantel plot for numerical 'sugar' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column sugar \
    --o-metadata-distance-matrix ./subject-analysis/metrics/jaccard-sugar-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/jaccard-sugar-mantel.qzv

# mantel plot for numerical 'coffee' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column coffee \
    --o-metadata-distance-matrix ./subject-analysis/metrics/jaccard-coffee-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/jaccard-coffee-mantel.qzv

# mantel plot for numerical 'alcohol' metadata
qiime diversity beta-correlation \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column alcohol \
    --o-metadata-distance-matrix ./subject-analysis/metrics/jaccard-alcohol-distance-matrix.qza \
    --o-mantel-scatter-visualization ./subject-analysis/beta-diversity/jaccard-alcohol-mantel.qzv
```
Here, the `-mantel.qzv` output file provides the mantel plot as well as statistical significance information.

## Gender alpha and beta diversity
Due to the categorical nature of gender (male/female), the original method of alpha and beta diversity analyses can be employed. 

For alpha diversity:
```sh 
#!/bin/bash
# for shannon metric...
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./subject-analysis/metrics/shannon-saliva-vector.qza \
    --m-metadata-file subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/shannon-categorical-group-significance.qzv

# for simpson metric...
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./subject-analysis/metrics/simpson-saliva-vector.qza \
    --m-metadata-file subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/simpson-categorical-group-significance.qzv

# for chao1 metric...
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./subject-analysis/metrics/chao1-saliva-vector.qza \
    --m-metadata-file subject-metadata.txt \
    --o-visualization ./subject-analysis/alpha-diversity/chao1-categorical-group-significance.qzv

```
And for beta diversity:
```sh
#!/bin/bash
# bray curtis permanova
qiime diversity beta-group-significance \
    --i-distance-matrix ./subject-analysis/metrics/bc-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column gender \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./subject-analysis/beta-diversity/bc-gender-permanova.qzv

# jaccard permanoava
qiime diversity beta-group-significance \
    --i-distance-matrix ./subject-analysis/metrics/jaccard-saliva-distance-matrix.qza \
    --m-metadata-file subject-metadata.txt \
    --m-metadata-column gender \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./subject-analysis/beta-diversity/jaccard-gender-permanova.qzv
```

# Differential Abundance Analysis
The differential abundance analysis was executed using two programs: ANCOM-BC which is integrated in the base package of qiime2, and songbird, in which the plugin version for the qiime2 environment, was used. Due to the inconsistent nature of various differential abundance analysis methods, using two methods and looking for corresponding results gives greater confidence in the outcome.

**Note: Like the RPCA and CTF analyses, these differential abundance tests filter the feature tables in the process. Therefore, the `merged-featureTable.qza` artifact is used as the feature table input.**

## ANCOM-BC
By inputting the feature table and metadata files into the `qiime composition ancombc` function, a differential abundance artifact is generated. This can be used for downstream visualizations. It is important to note a few arguments in this function:

* `--p-formula`: the categorical metadata column for the differential analysis. In this case it is `group` which divides the samples into their respective conditions: 'saliva-control', 'saliva-smoker', 'control-growth', 'control-beetroot', 'smoker-growth', and 'smoker-beetroot'.
* `--p-reference-levels`: the level in which all other conditions in the specified formula above will be compared to. The default selects the reference level in alphabetical order. Here it is specified as `group::saliva-control` meaning the the 'saliva-control' group will be the reference.
```sh
#!/bin/bash
qiime composition ancombc \
    --i-table merged-featureTable.qza \
    --m-metadata-file metadata.txt \
    --p-formula 'group' \
    --p-reference-levels group::saliva-control \
    --o-differentials ./differential_abundance/differential.qza
```
### Differential abundance bar plot and table
With the `differential.qza` outputted from the function above, this can now be inputted into `qiime composition da-barplot` and `qiime composition tabulate` to generate a barplot and table visual, respectively.
```sh
#!/bin/bash
# differential abundance barplot
qiime composition da-barplot \
    --i-data ./differential_abundance/differential.qza \
    --o-visualization ./differential_abundance/da-barplot.qzv

# differential abundance table
qiime composition tabulate \
    --i-data ./differential_abundance/differential.qza \
    --o-visualization ./differential_abundance/da-table.qzv
```
This process was repeated using each group condition as the reference level.

## songbird
Songbird can be used either as a standalone environment or, as previously mentioned, a plugin for the qiime2 enviornment. The latter was used for this analysis. The songbird plugin is only compatible with qiime2 v2019.7 through v2020.6 so qiime2 v2020.6 was used. 

Songbird operates on a machine learning algorithm and hence input parameters were set and modified to iteratively improve the models fit until certain features were observable in the output model: an exponential decay followed by a plateau, the model line is "below" the null model baseline, and there is a positive Q^2 value. Because of the repetitiveness of this, the script below was created to expidite the process.

Finally, it is important to note that songbird doesn't allow for selection of which conditions to compare. Therefore, to bypass this issue, and allow a one-to-one comparison with ANCOM-BC, separate feature tables and metadata files were created for each comparison desired. For example, to compare just the nonsmokers beetroot to smokers beetroot samples, a new feature table and metadata file was filtered to contain just the information for those two conditions.
```sh
#!/bin/bash

# Runs Songbird through QIIME 2 on an OTU table, generating the null model
# for the Q^2 score and visualization. 
# Arguments:
#   1: OTU table as a .qza file
#   2: Path to QIIME 2 environment with Songbird install (e.g. "qiime2-2020.6")
#   3: Output iteration number as integer

# Songbird arguments
metadata_file="/home/exacloud/gscratch/BakerLab/daniela-AHS-project/exacloud/metadata.txt"
formula='group'
differential_prior='0.5'
batch_size='2'
test_examples='2'
summary_interval='1'
epochs='100000'

filename=$(basename $1 .qza)
working_dir=$(dirname $(realpath $1))
dir_songbird=$(echo "${working_dir}/${filename}-songbird-${3}")
dir_model=$(echo "${dir_songbird}/model/")
dir_nullmodel=$(echo "${dir_songbird}/nullmodel/")

mkdir $dir_songbird
mkdir $dir_model
mkdir $dir_nullmodel

file_model_differentials="${dir_model}/${filename}-model-differentials.qza"
file_model_regressionstats="${dir_model}/${filename}-model-regressionstats.qza"
file_model_regressionbiplot="${dir_model}/${filename}-model-regressionbiplot.qza"

file_nullmodel_differentials="${dir_nullmodel}/${filename}-nullmodel-differentials.qza"
file_nullmodel_regressionstats="${dir_nullmodel}/${filename}-nullmodel-regressionstats.qza"
file_nullmodel_regressionbiplot="${dir_nullmodel}/${filename}-nullmodel-regressionbiplot.qza"

file_regression_visualization="${dir_songbird}/${filename}-songbird-regression.qzv"

dir_model_export="${dir_model}/export/"
dir_nullmodel_export="${dir_nullmodel}/export/"

# Songbird
conda run -p $2 qiime songbird multinomial \
    --i-table $1 \
    --m-metadata-file $metadata_file \
    --p-formula $formula \
    --p-differential-prior $differential_prior \
    --p-summary-interval $summary_interval \
    --p-epochs $epochs \
    --p-num-random-test-examples $test_examples \
    --p-batch-size $batch_size \
    --verbose \
    --o-differentials $file_model_differentials \
    --o-regression-stats $file_model_regressionstats \
    --o-regression-biplot $file_model_regressionbiplot

# Null model
conda run -p $2 qiime songbird multinomial \
    --i-table $1 \
    --m-metadata-file $metadata_file \
    --p-formula "1" \
    --p-differential-prior $differential_prior \
    --p-summary-interval $summary_interval \
    --p-epochs $epochs \
    --p-num-random-test-examples $test_examples \
    --p-batch-size $batch_size \
    --verbose \
    --o-differentials $file_nullmodel_differentials \
    --o-regression-stats $file_nullmodel_regressionstats \
    --o-regression-biplot $file_nullmodel_regressionbiplot

# Visualize model
conda run -p $2 qiime songbird summarize-paired \
    --i-regression-stats $file_model_regressionstats \
    --i-baseline-stats $file_nullmodel_regressionstats \
    --o-visualization $file_regression_visualization

# Export model data
for file in $dir_model/*; do
    conda run -p $2 qiime tools export \
        --input-path $file \
        --output-path $dir_model_export
done

# Export null model data
for file in $dir_nullmodel/*; do
    conda run -p $2 qiime tools export \
        --input-path $file \
        --output-path $dir_nullmodel_export
done
```
This process was repeated for each combination of conditions to compare.