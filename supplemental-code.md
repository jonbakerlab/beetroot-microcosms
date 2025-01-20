# daniela-AHS-beetroot-project Supplemental Code
From the nanopore sequencing, we were able to generate separate feature tables (OTU tables) for each flow cell run. To merge the feature tables together, the `feature-table merge` function in the qiime2 conda environment was used.

## Convert .txt files to .qza files
Before merging the feature tables, first they must be converted the .txt tables into compatible .qza files:

```sh
#!/bin/bash
# first convert featureTable into a .biom format:
conda activate qiime2-metagenome-2024.5
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
## qiime2 Emperor Plots and Rarefaction Curve
Now that all of the featureTables have been combined into `merged-featureTable.qza`, the following code uses the `diversity core-metrics` function in qiime2's environment to generate a variety of diversity metric artifacts and visuals used in downstream analyses:

**NOTE: `gemelli` is a plugin downloaded and used in downstream analyses (RCPA and CTF) however, downloading the gemelli plugin downgrades the scikit-bio package to v0.5.9 from the base qiime scikit-bio v0.6.0. This downgraded version of scikit-bio is incompatible with the `diversity core-metrics` command specifically.**

```sh
#!/bin/bash
qiime diversity core-metrics \
    --i-table ./feature_tables/merged-featureTable.qza \
    --p-sampling-depth 4500 \
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

The following code generates an alpha rarefaction curve visual with the qiime2 environment's `diversity alpha-rarefaction` command. This is helpful to visualize and ensure the all the taxa present has been observed in each sample.

```sh
#!/bin/bash
qiime diversity alpha-rarefaction \
    --i-table merged-featureTable.qza \
    --p-max-depth 60000 \
    --p-steps 40 \
    --m-metadata-file metadata.txt \
    --o-visualization alpha-rarefaction.qzv
```

# Beta Diversity Analysis
## Bray Curtis emperor plot qnd permanova
Using the `bray_curtis_distance_matrix.qza` artifact previously exported from the `diversity core-metrics` function, a permanova analysis can be conducted. Alternatively, the following sequence of codes can be run independent of `diversity core-metrics` to obtain the necessary artifacts for an emperor plot and permanova visualization:

```sh
#!/bin/bash
# generate a braycurtis distance matrix
qiime diversity beta \
 --i-table merged-featureTable.qza \
 --p-metric braycurtis \
 --o-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza

# generate a braycurtis pcoa
qiime diversity pcoa \
 --i-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza \
 --o-pcoa ./core-metrics-results/bray_curtis_pcoa.qza

# generate the emperor plot
qiime emperor plot \
  --i-pcoa ./core-metrics-results/bray_curtis_pcoa.qza \
  --m-metadata-file metadata.txt \
  --o-visualization ./beta_diversity/bray_curtis_emperor.qzv

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
An identical process was done for the Jaccard metric by substituting the `jaccard` argument for `--p-metric` in the `diversity beta` function:
```sh
#!/bin/bash
# distance matrix:
qiime diversity beta \
 --i-table merged-featureTable.qza \
 --p-metric jaccard \
 --o-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza

# pcoa:
qiime diversity pcoa \
 --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
 --o-pcoa ./core-metrics-results/jaccard_pcoa.qza

# emperor plot:
qiime emperor plot \
  --i-pcoa ./core-metrics-results/jaccard_pcoa.qza \
  --m-metadata-file metadata.txt \
  --o-visualization ./beta_diversity/jaccard_emperor.qzv

# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization ./beta_diversity/jaccard_permanova.qzv
```

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
Now we can import the FASTA file we created into QIIME2.  Keep in mind you'll need to add underscores to your feature names in your feature table to match the Genus_species format we just created to make the tree (QIIME2 won't let us import the fasta file with spaces in the deflines):

With the properly formatted .fasta file, it can now be imported it into qiime:

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

## Weighted and Unweighted Unifrac emperor plots and permanova
### Filtering the featureTable and rooted phylogenetic tree
Through some trouble shooting, it was found that there were feature IDs in the feature tables that did not have a corresponding leaf in the rooted tree file. Fortunately, these missing species were relatively low in abundance. These species were filtered out based on a frequency count of reads less than 200. This was done using the `feature-table filter-feature` function: 

```sh
#!/bin/bash
# filter out features with < 200 reads:
qiime feature-table filter-features \
    --i-table ./feature-tables/merged-featureTables.qza \
    --p-min-frequency 200 \
    --o-filtered-table ./feature-tables/filtered-featureTable.qza
```

With the new `filtered-featureTable.qza` file, the `rooted-tree.qza` was also filtered with `phylogeny filter-tree` based on the remaining species present in the feature table:

```sh
#!/bin/bash
# filter tree based on remaining features:
qiime phylogeny filter-tree \
  --i-tree ./16s-phylogeny/rooted-tree.qza \
  --i-table ./feature_tables/filtered-featureTable.qza \
  --o-filtered-tree ./16s-phylogeny/filtered-tree.qza
```

### Weighted Unifrac emperor plot and permanova
Now that the `filtered-featureTable.qza` and `filtered-tree.qza` have the same corresponding feature ID's, a similar sequence of codes as the Bray Curtis and Jaccard beta diversity analyses can be used for the weighted unifrac analysis using `diversity beta-phylogenetic`:

```sh
#!/bin/bash
# distance matrix:
qiime diversity beta-phylogenetic \
    --i-table ./feature_tables/filtered-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-tree.qza \
    --p-metric weighted_unifrac \
    --o-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza

# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization  ./beta_diversity/weighted_unifrac_permanova.qzv

# pcoa:
qiime diversity pcoa \
    --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --o-pcoa ./core-metrics-results/weighted_unifrac_pcoa.qza

# emperor plot:
qiime emperor plot \
    --i-pcoa ./core-metrics-results/weighted_unifrac_pcoa.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./beta_diversity/weighted_unifrac_emperor.qzv
```
### Unweighted Unifrac emperor plot and permanova
And similarly for the Unweighted Unifrac analysis by using the `--p-metric unweighted_unifrac` argument:

```sh
#!/bin/bash
# distance matrix:
qiime diversity beta-phylogenetic \
    --i-table ./feature_tables/filtered-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-tree.qza \
    --p-metric unweighted_unifrac \
    --o-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza

# permanova:
qiime diversity beta-group-significance \
    --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --p-method permanova \
    --o-visualization  ./beta_diversity/unweighted_unifrac_permanova.qzv

# pcoa:
qiime diversity pcoa \
    --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --o-pcoa ./core-metrics-results/unweighted_unifrac_pcoa.qza

# emperor plot:
qiime emperor plot \
    --i-pcoa ./core-metrics-results/unweighted_unifrac_pcoa.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./beta_diversity/unweighted_unifrac_emperor.qzv
```
Finally, it is worth noting that the function `diversity core-metrics-phylogenetic` is would provide the results from `diversity core-metrics` in addition to the Weighted and Unweighted Unifrac artifacts. However, due to the order in which these were conducted, the phylogenetic analyses were done independently.

# Alpha Diversity Analysis
## Vectors of alpha diversity
Before performing the alpha diversity permanova, vector files for all the methods of alpha diversity analysis must be created. The shannon, chao1, and simpson metrics are all of interest, however, the vector for the shannon metric was created and exported using the `core-metrics` function. To create vectors for chao1 and simpson, `diversity alpha` function was used:

```sh
#!/bin/bash
# simpson metric:
qiime diversity alpha \
    --i-table ./feature_tables/merged-featureTable.qza \
    --p-metric simpson \
    --o-alpha-diversity ./core-metrics-results/simpson_vector.qza

# chao1 metric
qiime diversity alpha \
    --i-table ./feature_tables/merged-featureTable.qza \
    --p-metric chao1 \
    --o-alpha-diversity ./core-metrics-results/chao1_vector.qza 

# shannon metric:
qiime diversity alpha \
    --i-table ./feature_tables/merged-featureTable.qza \
    --p-metric shannon \
    --o-alpha-diversity ./core-metrics-results/shannon_vector.qza
```
### Vectors with phylogenetic trees
For vectors associated with a phylogenetic tree, its a similar process using the qiime `diversity alpha-phylogeny` function:

```sh
#!/bin/bash
# faith_pd metric:
qiime diversity alpha-phylogenetic \
    --i-table ./feature_tables/filtered-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-tree.qza \
    --p-metric faith_pd \
    --o-alpha-diversity ./core-metrics-results/faith_pd_vector.qza
```
## Alpha diversity group significance
With the alpha diversity vector artifacts, a group significance analysis using the `diversity alpha-group-significance` function can be run for each of the metrics, resepectively:

```sh
#!/bin/bash
# shannon metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./alpha_diversity/shannon_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha_diversity/shannon-group-significance.qzv

# simpson metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./alpha_diversity/simpson_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha_diversity/simpson-group-significance.qzv

# chao1 metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./alpha_diversity/chao1_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha_diversity/chao1-group-significance.qzv

# faith_pd metric:
qiime diversity alpha-group-significance \
    --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization ./alpha_diversity/faith-group-significance.qzv
```

# Gemelli Plugin Analysis
The `gemelli` plugin provides other methods for visualizing significant differences among groups: Robust Aitchison PCA (RPCA) and Compositional Tensor Factorization (CTF).

## Robust Aitchison PCA (RPCA)
**Insert something about RPCA**

... For this section of analysis, a separate qiime enivronment was created specifically for the gemelli plugin due to the aforementioned version differences:
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
With the newly formatted taxonomy artifact, all the required files are created and can be inputted as arguments in the `gemelli phylogenetic-rcpa-with-taxonomy` function:

```sh
#!/bin/bash
# create all the metrics needed for downstream analysis:
qiime gemelli phylogenetic-rpca-with-taxonomy \
    --i-table ./feature_tables/filtered-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-tree.qza \
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

The permanova can be created in the base qiime environment, as done previously with the RCPA, Jaccard, and Bray Curtis metrics:

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
**Insert something about CTF**

...Much like the RCPA analysis, the CTF analysis must be done in the `gemelli` environment.  Using the `gemelli ctf` function... **something about what the subject and state columns are**
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
Similarly to the phylogenetic RPCA analysis, the taxonomic metadata is required. This essentially combines the Phylogenetic RPCA section and the CTF analysis, i.e., uses the state and subject columns.
```sh
#!/bin/bash
qiime gemelli phylogenetic-ctf-with-taxonomy \
    --i-table ./feature_tables/filtered-featureTable.qza \
    --i-phylogeny ./16s-phylogeny/filtered-tree.qza \
    --m-sample-metadata-file metadata.txt \
    --m-taxonomy-file ./16s-phylogeny/taxonomy/filtered-taxonomy.qza \
    --p-state-column 'group' \
    --p-individual-id-column sample \
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

# Differential Abundance Analysis
## ANCOMBC
```sh
#!/bin/bash
qiime composition ancombc \
    --i-table merged-featureTable.qza \
    --m-metadata-file metadata.txt \
    --p-formula 'group' \
    --p-reference-levels group::saliva-control \
    --o-differentials ./differential_abundance/differential.qza
```
### differential abundance bar plots
```sh
#!/bin/bash
qiime composition da-barplot \
    --i-data ./differential_abundance/differential.qza \
    --o-visualization ./differential_abundance/da-barplot.qzv
```

### differential abundance table
```sh
#!/bin/bash
qiime composition tabulate \
    --i-data ./differential_abundance/differential.qza \
    --o-visualization ./differential_abundance/da-table.qzv
```

## songbird

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