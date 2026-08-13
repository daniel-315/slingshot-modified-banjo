# Slingshot + Monocle 3 + Modified Banjo Gene Regulatory Network Analysis

This repository contains the code, input files, configuration files, and results used to infer dynamic gene regulatory networks from Slingshot pseudotime-ordered single-cell gene-expression data using a modified version of Banjo (Bayesian Network Inference with Java Objects).

Two primary analyses are included:

1. A network inferred from the 10 most densely connected genes in the supplied Gold Standard.
2. A network inferred from the 15 most densely connected genes in the supplied Gold Standard.

The repository preserves all top-scoring non-identical networks returned by Banjo for both analyses rather than saving only the single highest-scoring network.

---

# 1. Project Overview

Single-cell experiments measure many individual cells at a limited number of experimental time points. Cells collected at the same experimental time do not necessarily represent exactly the same biological stage.

For this project, Slingshot was used to infer a continuous pseudotime ordering from gene-expression patterns. The cells were then reordered according to this inferred trajectory.

The pseudotime-ordered expression matrix was supplied to Modified Banjo. Banjo was used as a dynamic Bayesian network (DBN) inference method to identify lagged statistical relationships among genes.

The primary analysis used:

- Original number of cells: 960
- Cell excluded from the primary Slingshot analysis: Cell 591
- Cells remaining: 959
- Measured genes in the original RNA dataset: 45
- Pseudotime method: Slingshot
- Network-inference method: Modified Banjo
- Markov lag: 1
- Usable lag-1 transitions: 958
- Discretization policy: i2
- Maximum parents per gene: 5
- Mandatory identity lag: 1
- Maximum requested non-identical networks: 5
- Banjo search time for the reported analyses: 5 minutes

Two subsets of the original genes were analyzed:

- Top 10 Gold-Standard-density genes
- Top 15 Gold-Standard-density genes

---

# 2. Repository Structure

The main repository structure is:

~~~
slingshot-modified-banjo/
|
|-- README.md
|-- compile_banjo.sh
|-- run_top10.sh
|-- run_top15.sh
|
|-- modified_banjo/
|   |-- source/
|   |   `-- edu/duke/cs/banjo/
|   `-- build/
|
|-- configs/
|   |-- banjo_top10.txt
|   `-- banjo_top15.txt
|
|-- data/
|   `-- banjo_inputs/
|       |-- top10_slingshot_input.txt
|       `-- top15_slingshot_input.txt
|
|-- gene_sets/
|   |-- top10_genes.txt
|   `-- top15_genes.txt
|
|-- scripts/
|   `-- gene_selection/
|       `-- make_banjo_subset.R
|
|-- docs/
|   `-- gold_standard_gene_ranking.csv
|
`-- results/
    |-- top10/
    |   |-- banjo_full_report.txt
    |   |-- banjo_input_top10.txt
    |   |-- banjo_settings_top10.txt
    |   |-- top10_genes.txt
    |   `-- all_networks/
    |
    `-- top15/
        |-- banjo_full_report.txt
        |-- banjo_input_top15.txt
        |-- banjo_settings_top15.txt
        |-- top15_genes.txt
        `-- all_networks/
~~~

The folders `results/top10/` and `results/top15/` contain the original reported Banjo results.

When the provided run scripts are executed again, newly generated results are written separately under `results/reproduced/`.

---

# 3. Why These Specific Genes Were Selected

The Top-10 and Top-15 gene subsets were not chosen randomly.

The goal was to focus network inference on genes that were relatively densely connected according to the supplied Gold Standard regulatory network.

The professor's guidance was to identify genes with more positive regulatory connections in the Gold Standard because selecting densely connected genes could make network inference more informative than using an arbitrary subset.

## 3.1 Gold Standard Relationships

The supplied Gold Standard contains regulatory relationships between genes.

For the gene-density ranking used in this project:

- A positive relationship is represented by `1`.
- A negative relationship is represented by `-1`.
- Only positive relationships were counted for the density score.
- Only relationships where both genes were among the 45 genes measured in the RNA dataset were included.

Relationships involving genes outside the measured 45-gene dataset were not used for ranking.

## 3.2 Connectivity Score

For each measured gene, two positive-connectivity quantities were counted.

Positive outgoing degree:

    Number of positive Gold Standard edges from the gene
    to another one of the 45 measured genes.

Positive incoming degree:

    Number of positive Gold Standard edges from another
    measured gene into the gene.

The total connectivity score was then calculated as:

    Total positive connectivity
    =
    positive outgoing degree
    +
    positive incoming degree

Both incoming and outgoing relationships were included because a gene may be highly connected either as a regulator, as a regulatory target, or both.

Genes were ranked from highest to lowest total positive connectivity.

If genes had equal total positive connectivity, positive outgoing degree was used as the first tie-breaking criterion.

---

# 4. Top-10 Gene Set

The 10 highest-ranked genes were:

| Rank | Gene | Positive Outgoing | Positive Incoming | Total Connectivity |
|---:|---|---:|---:|---:|
| 1 | PPARG | 7 | 28 | 35 |
| 2 | NFATC2 | 9 | 23 | 32 |
| 3 | BCL6 | 5 | 27 | 32 |
| 4 | MAFB | 5 | 23 | 28 |
| 5 | PRDM1 | 0 | 23 | 23 |
| 6 | SNAI3 | 2 | 19 | 21 |
| 7 | EGR2 | 11 | 9 | 20 |
| 8 | ETS1 | 3 | 16 | 19 |
| 9 | STAT1 | 9 | 8 | 17 |
| 10 | FOS | 8 | 9 | 17 |

Therefore, the Top-10 gene set is:

    PPARG
    NFATC2
    BCL6
    MAFB
    PRDM1
    SNAI3
    EGR2
    ETS1
    STAT1
    FOS

The same list is stored in:

    gene_sets/top10_genes.txt

---

# 5. Top-15 Gene Set

The Top-15 analysis contains all 10 genes above plus five additional highly connected genes.

| Rank | Gene | Positive Outgoing | Positive Incoming | Total Connectivity |
|---:|---|---:|---:|---:|
| 11 | NFYC | 13 | 3 | 16 |
| 12 | IRF8 | 2 | 14 | 16 |
| 13 | TCF3 | 14 | 1 | 15 |
| 14 | PPARD | 6 | 8 | 14 |
| 15 | SNAI1 | 10 | 3 | 13 |

The complete Top-15 set is:

    PPARG
    NFATC2
    BCL6
    MAFB
    PRDM1
    SNAI3
    EGR2
    ETS1
    STAT1
    FOS
    NFYC
    IRF8
    TCF3
    PPARD
    SNAI1

The same list is stored in:

    gene_sets/top15_genes.txt

The Top-15 set is therefore an extension of the Top-10 set rather than an independently selected group.

## 5.1 Rank-15 Tie

SNAI1 and RREB1 both had total positive connectivity values of 13.

Their positive outgoing degrees differed:

- SNAI1: 10 positive outgoing relationships
- RREB1: 9 positive outgoing relationships

Because positive outgoing degree was used as the first tie-breaking criterion, SNAI1 was ranked above RREB1 and selected as gene 15.

---

# 6. Slingshot Pseudotime

Slingshot is a trajectory-inference method designed to estimate progression through biological states from single-cell data.

Pseudotime is not literal elapsed clock time. It is an inferred ordering that represents progression along a trajectory based on similarities and changes in gene-expression profiles.

For the primary Slingshot workflow used in this project:

1. The original RNA dataset contained 960 cells.
2. The dataset contained 45 measured genes plus the experimental-time variable.
3. Cell 591 was excluded from the primary analysis as an extreme multigene/PCA outlier.
4. PCA was recalculated using the remaining 959 cells.
5. Expression values were centered and scaled during PCA.
6. The first 10 principal components were used for the clustering step.
7. K-means clustering was evaluated.
8. Two clusters were selected for the primary trajectory.
9. Slingshot was used to infer the trajectory.
10. Cells were assigned pseudotime values.
11. Cells were sorted from lower to higher pseudotime.
12. The original expression values were reordered according to this pseudotime ordering.
13. The reordered expression matrix was supplied to Banjo.

The Slingshot pseudotime ordering showed a positive relationship with the original experimental-time variable, supporting the use of pseudotime as an alternative ordering of the cells.

---

# 7. Banjo Input Files

The repository already contains the pseudotime-ordered Banjo input files:

    data/banjo_inputs/top10_slingshot_input.txt

and:

    data/banjo_inputs/top15_slingshot_input.txt

Therefore, someone who only wants to reproduce the Banjo network inference does not need to rerun the entire Slingshot analysis first.

Both files contain:

- 959 ordered observations

The Top-10 file contains:

- 10 genes

The Top-15 file contains:

- 15 genes

Because a first-order dynamic Bayesian network uses transitions between adjacent observations:

    959 observations - 1 = 958 lag-1 transitions

Banjo therefore reports 959 observations in the file but 958 observations used for learning the dynamic Bayesian network.

This is expected behavior.

---

# 8. What Is Banjo?

Banjo stands for:

**Bayesian Network Inference with Java Objects**

Banjo was developed at Duke University for Bayesian-network and dynamic-Bayesian-network inference.

In this project, Banjo is used to learn relationships among genes using the order supplied by Slingshot pseudotime.

The analysis uses a Markov lag of 1.

Conceptually, a relationship such as:

    ETS1(t) -> PPARG(t+1)

means that the state of ETS1 at one ordered observation contributes to the statistical model for PPARG at the following ordered observation.

An inferred Banjo edge represents a statistical dependency in the learned dynamic Bayesian network.

It should not automatically be interpreted as experimental proof of direct biological causation.

---

# 9. Banjo Settings Used

The major settings used in these analyses include:

    minMarkovLag = 1
    maxMarkovLag = 1
    dbnMandatoryIdentityLags = 1
    maxParentCount = 5
    discretizationPolicy = i2
    nBestNetworks = 5
    bestNetworksAre = nonidentical
    maxTime = 5 m
    computeInfluenceScores = yes
    createDotOutput = yes

For Top 10:

    variableCount = 10
    observationCount = 959

For Top 15:

    variableCount = 15
    observationCount = 959

The portable settings files are:

    configs/banjo_top10.txt
    configs/banjo_top15.txt

---

# 10. Mandatory Identity Edges

The setting:

    dbnMandatoryIdentityLags = 1

requires each gene to retain a lag-1 relationship with itself.

Conceptually:

    Gene(t) -> Gene(t+1)

These identity relationships are different from inferred cross-gene relationships.

A network can therefore be valid even if it contains no additional cross-gene edges beyond the mandatory identity relationships.

---

# 11. What Was Modified in Banjo?

The source code in this repository is based on Banjo 2.2.0.

The project-specific modification is located in:

    modified_banjo/source/edu/duke/cs/banjo/utility/PostProcessor.java

The original post-processing code assumed that the DOT representation of every inferred network contained a section containing non-identity edges.

During this project, Banjo sometimes inferred a network containing only mandatory identity/self-lag relationships and no additional cross-gene relationships.

In that situation, the original post-processing logic could attempt to access a DOT-output section that did not exist.

The original behavior effectively assumed:

    A non-identity edge section exists
    -> parse that section

The modified behavior is:

    Check whether the expected DOT section exists
    -> check whether it contains a "->" edge
    -> parse it only if such an edge exists

The relevant logic now safely initializes an empty edge list when no cross-gene edges are present.

This prevents post-processing from failing on networks containing no non-identity edges.

Importantly, this modification affects post-processing and graph generation. It does not change Banjo's underlying network-search or BDe scoring algorithm.

---

# 12. Why Some Network Images Can Be Empty

The DOT visualizations focus on cross-gene relationships.

Because mandatory identity-lag relationships are treated separately, an inferred network containing only identity relationships can produce an essentially empty cross-gene DOT graph.

This does not necessarily indicate a failed Banjo run.

The Modified Banjo patch allows such a network to be written safely as an empty graph rather than causing the Java post-processing code to fail.

---

# 13. Requirements

The easiest environment for reproducing the analyses is Ubuntu Linux or Ubuntu through Windows Subsystem for Linux (WSL).

For the already prepared Banjo inputs, the main requirements are:

- Java Development Kit
- Graphviz
- Git, if cloning from GitHub

R is only necessary if the user wants to perform additional data processing or regenerate gene subsets.

---

# 14. Installing the Required Software on Ubuntu or WSL

Open an Ubuntu terminal.

Run:

    sudo apt update

Then install Java, Graphviz, and Git:

    sudo apt install -y openjdk-17-jdk graphviz git

Verify Java:

    java -version

Verify the Java compiler:

    javac -version

Verify Graphviz:

    dot -V

If these commands print version information, the main required software is installed.

---

# 15. Downloading the Repository

Using Git:

    git clone https://github.com/daniel-315/slingshot-modified-banjo.git

Then enter the repository:

    cd slingshot-modified-banjo

Replace `https://github.com/daniel-315/slingshot-modified-banjo.git` with the actual private GitHub repository URL.

Alternatively, the repository can be downloaded as a ZIP file from GitHub and extracted manually.

---

# 16. Compiling Modified Banjo

The repository includes:

    compile_banjo.sh

This script automatically finds and compiles the Java source files.

Run:

    ./compile_banjo.sh

A successful build should end with output similar to:

    Compilation successful.
    Compiled class files:
    138

Modern versions of Java may display warnings that Java source/target version 8 is obsolete.

These are compiler warnings rather than compilation failures.

The script uses:

    --release 8

to compile the older Banjo Java code with Java 8 compatibility.

---

# 17. Running the Top-10 Analysis

From the repository root, run:

    ./run_top10.sh

The script:

1. Checks that Java is installed.
2. Checks that Graphviz is installed.
3. Checks whether Banjo has been compiled.
4. Automatically compiles Banjo if necessary.
5. Loads the Top-10 input.
6. Loads the Top-10 settings.
7. Runs the Banjo network search.
8. Requests up to five non-identical networks.
9. Runs for approximately five minutes.
10. Stores newly generated files separately from the original reported results.

New results are written to:

    results/reproduced/top10/

The original reported results remain under:

    results/top10/

---

# 18. Running the Top-15 Analysis

Run:

    ./run_top15.sh

The same process is used with:

- 15 genes
- 959 observations
- 958 usable lag-1 transitions
- Up to 5 non-identical networks
- Approximately 5 minutes of Banjo search time

Newly reproduced results are written to:

    results/reproduced/top15/

The original reported results remain under:

    results/top15/

---

# 19. Understanding the Banjo Output

Banjo generates several output types.

## 19.1 dynamic.report File

A file with a name similar to:

    dynamic.report.YYYY.MM.DD.HH.MM.SS.txt

is the primary Banjo report.

It contains information such as:

- Input file
- Number of observations
- Number of variables
- Discretization settings
- Markov lag
- Maximum parent count
- Network scores
- Network structures
- Search statistics
- Influence scores
- DOT graph output

This report is the best place to verify exactly how a run was performed.

## 19.2 DOT Files

Files ending in:

    .dot

contain Graphviz descriptions of inferred cross-gene network relationships.

## 19.3 PNG Files

Files ending in:

    .png

are graphical network images produced from the DOT files using Graphviz.

These are usually the easiest network outputs to inspect visually.

## 19.4 SMV Files

Files with names such as:

    SMVGRAPH0.smv

contain additional network representations produced during Banjo post-processing.

They are preserved so that the complete original output is available.

---

# 20. Original Top-10 Results

The original Top-10 Banjo analysis requested:

    nBestNetworks = 5

with:

    bestNetworksAre = nonidentical

However, Banjo returned only two non-identical top-scoring networks.

This is valid behavior.

The setting requests up to five non-identical networks; it does not guarantee that five distinct networks will be returned.

The Top-10 results therefore contain two network outputs.

They are preserved under:

    results/top10/all_networks/

The full Banjo report is:

    results/top10/banjo_full_report.txt

The highest-scoring Top-10 network had a Banjo score of approximately:

    -2262.2199

One cross-gene relationship appearing in the highest-scoring Top-10 network was:

    ETS1 -> PPARG

with a reported influence score of approximately:

    +0.0637

The second returned Top-10 network contained no additional cross-gene relationships beyond the identity structure.

---

# 21. Original Top-15 Results

The Top-15 Banjo analysis successfully returned all five requested non-identical networks.

The five network scores were approximately:

    Network 1: -3889.6251
    Network 2: -3889.8552
    Network 3: -3891.3869
    Network 4: -3891.6365
    Network 5: -3892.2044

The full set of Top-15 outputs is stored under:

    results/top15/all_networks/

The full Banjo report is:

    results/top15/banjo_full_report.txt

All five returned network results were preserved rather than keeping only Network 1.

---

# 22. Reproducibility Test

A short smoke test was performed after reorganizing the project into this repository.

For this test:

- Modified Banjo was compiled directly from the source stored in this repository.
- 138 Java class files were produced.
- The portable Top-10 configuration was used.
- The repository correctly found the Top-10 Slingshot input.
- Banjo recognized 959 observations.
- Banjo recognized 10 variables.
- Banjo used 958 observations for lag-1 DBN learning.
- The temporary search time was reduced to 10 seconds.
- Banjo successfully produced report, SMV, DOT, and PNG outputs.
- A network containing no non-identity edges was processed successfully.

This confirms that the repository does not depend on the original local Banjo directory in order to compile and run the included analysis.

---

# 23. Original Results Versus Reproduced Results

The repository deliberately separates the exact original results from new runs.

Original results:

    results/top10/
    results/top15/

Results from a new execution:

    results/reproduced/top10/
    results/reproduced/top15/

Users should not overwrite the original reported result directories when testing the program.

---

# 24. Random Search and Reproducibility

Banjo's search procedure uses a starting random seed.

Because of this, running the same five-minute analysis again may not always produce exactly the same sequence of searched networks or identical final results unless the same random seed is explicitly controlled.

For that reason, the exact outputs produced for this project have been preserved in:

    results/top10/
    results/top15/

The reproduced-output directories are intended for newly executed runs.

---

# 25. Interpretation of Results

The inferred networks represent statistical relationships learned by a dynamic Bayesian network.

Several cautions are important.

First, an inferred edge is not automatically proof of a direct biological regulatory interaction.

Second, pseudotime is an inferred ordering rather than experimentally measured continuous time.

Third, the signs and magnitudes of Banjo influence scores should not automatically be interpreted as direct biological activation or inhibition without additional validation.

Finally, the Gold Standard was used to select the densely connected gene subsets.

Because gene selection itself used Gold Standard information, evaluation against the same Gold Standard is not a completely independent or unbiased benchmark.

The Top-10 and Top-15 analyses are therefore best understood as targeted, connectivity-enriched network-inference experiments.

---

# 26. Troubleshooting

## Java command not found

Install the Java Development Kit:

    sudo apt update
    sudo apt install -y openjdk-17-jdk

Then verify:

    java -version
    javac -version

## Graphviz or dot command not found

Install Graphviz:

    sudo apt install -y graphviz

Verify:

    dot -V

## Permission denied when running a script

Run:

    chmod +x compile_banjo.sh run_top10.sh run_top15.sh

Then try the script again.

## Banjo has not been compiled

Run:

    ./compile_banjo.sh

The Top-10 and Top-15 run scripts also attempt to compile Banjo automatically if the main compiled Banjo class cannot be found.

## Only two Top-10 networks are present

This is expected.

The original Top-10 analysis requested up to five non-identical networks but returned two.

The two networks in the results directory are the complete set returned by that analysis.

## A network graph appears empty

This can happen when Banjo returns a network without additional cross-gene edges beyond the required identity-lag structure.

The modification to `PostProcessor.java` allows this type of network to be processed without causing a DOT post-processing error.

---

## Using WSL or a Command Prompt Terminal

A lab member independently tested this repository successfully using WSL launched through a Windows command-prompt workflow.

When running the provided shell scripts from the repository directory, WSL/Linux may require `./` in front of the script name.

For example, use:

    ./compile_banjo.sh
    ./run_top10.sh
    ./run_top15.sh

rather than:

    compile_banjo.sh
    run_top10.sh
    run_top15.sh

The `./` prefix tells the Linux shell to execute the script located in the current directory.

This applies to the included shell scripts. The Java executable itself is still called normally as:

    java

If a command such as `run_top10.sh` is reported as "command not found" even though the file is present, first try:

    ./run_top10.sh

Also make sure you are inside the repository directory:

    cd slingshot-modified-banjo

---

# 27. Quick Start for a New User

For someone using a new Ubuntu or WSL environment, the shortest procedure is:

    sudo apt update
    sudo apt install -y openjdk-17-jdk graphviz git

Then clone the repository:

    git clone https://github.com/daniel-315/slingshot-modified-banjo.git
    cd slingshot-modified-banjo

Compile Banjo:

    ./compile_banjo.sh

Run Top 10:

    ./run_top10.sh

Run Top 15:

    ./run_top15.sh

The reproduced networks will then appear under:

    results/reproduced/top10/
    results/reproduced/top15/

---

# 28. Software Used

The workflow uses:

- Slingshot for trajectory and pseudotime inference
- Monocle 3 for an independent pseudotime and trajectory comparison
- R for single-cell data processing and gene subsetting
- Banjo 2.2.0 for dynamic Bayesian network inference
- Java for compiling and running Banjo
- Graphviz for rendering DOT network files as PNG images
- Ubuntu/WSL as the primary command-line environment

---

# 29. Banjo License Notice

Banjo is software developed at Duke University.

The original Banjo source files included in this project retain their original copyright and licensing notices.

The source identifies Banjo as licensed from Duke University and contains copyright notices for Alexander J. Hartemink.

The modification included in this repository does not imply ownership of the original Banjo software.

Anyone redistributing or using the Banjo source should review and comply with the original Banjo licensing terms.

Because of these licensing considerations and because the repository contains research data and derived inputs, this project should remain private unless appropriate permission is obtained for public redistribution.

---

# 30. Analysis Summary

The workflow represented in this repository can be summarized as:

    Original single-cell RNA expression
                    |
                    v
             Slingshot analysis
                    |
                    v
          Pseudotime cell ordering
                    |
                    v
              959 cells
                    |
              +-----+-----+
              |           |
              v           v
          Top 10       Top 15
         dense genes  dense genes
              |           |
              v           v
        Modified Banjo Modified Banjo
              |           |
              v           v
        All returned   All returned
          networks       networks

The Top-10 and Top-15 subsets were chosen according to a documented positive Gold Standard connectivity criterion instead of random selection.

All original Banjo results are preserved so that subsequent reruns do not overwrite the networks used for the reported project results.

---

# 31. Monocle 3 Comparison Analysis

A second pseudotime method was added to determine how sensitive the inferred
dynamic Bayesian networks are to the pseudotime ordering method.

The comparison was designed so that the pseudotime method is the primary
difference between the analyses.

The four network analyses are:

1. Slingshot Top-10
2. Monocle 3 Top-10
3. Slingshot Top-15
4. Monocle 3 Top-15

The same cells, expression matrix, gene-selection criterion, selected genes,
Banjo discretization policy, Markov lag, search algorithm, maximum parent
count, scoring method, and five-minute search duration were used wherever
applicable.

---

# 32. Monocle 3 Pseudotime Method

The Monocle comparison used:

- R 4.3.3
- Monocle 3 version 1.3.1
- 959 cells
- 45 measured genes
- Cell 591 excluded, matching the primary Slingshot analysis
- 10 PCA dimensions
- scaled PCA
- UMAP dimensionality reduction
- Monocle clustering
- principal graph learning
- programmatic trajectory rooting using the earliest experimental time point

The original expression matrix contains continuous non-integer expression
values rather than raw count values.

For that reason, the Monocle preprocessing step used:

    norm_method = "none"

This avoids applying an additional count-style normalization to expression
values that were already processed.

Scaling remained enabled for PCA:

    scaling = TRUE

The Monocle trajectory was generated with:

    preprocess_cds(...)
    reduce_dimension(...)
    cluster_cells(...)
    learn_graph(...)
    order_cells(...)

To obtain one pseudotime coordinate across the same 959 cells used by the
Slingshot/Banjo workflow, the graph was learned with:

    use_partition = FALSE

The trajectory root was not manually selected.

Instead, the principal graph node containing the greatest number of cells
from the earliest experimental time point, h = 0, was selected
programmatically.

The selected root principal node was:

    Y_13

There were 120 cells at h = 0.

Monocle assigned finite pseudotime to all 959 analyzed cells.

The Monocle pseudotime range was:

    0 to 15.09055

The Spearman correlation between Monocle pseudotime and experimental time was:

    rho = 0.6367

For comparison, the primary Slingshot ordering had a Spearman correlation of
approximately:

    rho = 0.718

Both pseudotime methods therefore produced positive agreement with the known
experimental-time progression, although Slingshot showed the stronger
monotonic association in this dataset.

---

# 33. Monocle 3 v1.3.1 Partition Compatibility Fix

During clustering with Monocle 3 version 1.3.1, the analysis encountered an
igraph error because the partition adjacency matrix contained NaN values.

This occurred when the partition calculation encountered a zero
between-cluster-edge case.

A small runtime compatibility fix was added to:

    scripts/pseudotime/10_run_monocle3.R

The fix converts NaN entries in the relevant partition matrices to zero
before the graph is passed to igraph.

This follows the behavior used by newer upstream Monocle 3 source code for
the same zero-total-edge condition.

The compatibility fix does not modify:

- expression values
- selected cells
- selected genes
- PCA dimensionality
- UMAP coordinates
- Banjo settings
- Banjo scoring
- Banjo network search

It only prevents the older Monocle 3 partition implementation from failing
during graph construction.

After the fix, Monocle produced two clusters:

    Cluster 1: 861 cells
    Cluster 2: 98 cells

All 959 cells subsequently received finite pseudotime values.

---

# 34. Monocle Files Included in This Repository

The Monocle analysis scripts are:

    scripts/pseudotime/10_run_monocle3.R
    scripts/pseudotime/11_make_monocle3_banjo_inputs.R

The exported pseudotime ordering is:

    data/pseudotime/monocle3_order.csv

The trajectory figures are:

    results/monocle3/trajectory/monocle3_trajectory_pseudotime.png
    results/monocle3/trajectory/monocle3_trajectory_experimental_time.png

The Banjo input matrices are:

    data/banjo_inputs/banjo_monocle3_top10.txt
    data/banjo_inputs/banjo_monocle3_top15.txt

The portable Banjo settings files are:

    configs/banjo_monocle3_top10.txt
    configs/banjo_monocle3_top15.txt

The original five-minute Banjo results are preserved under:

    results/monocle3/top10/
    results/monocle3/top15/

---

# 35. Controlled Top-10 and Top-15 Comparison

The Monocle analysis uses exactly the same selected gene sets as the
Slingshot analysis.

## Top-10 genes

    PPARG
    NFATC2
    BCL6
    MAFB
    PRDM1
    SNAI3
    EGR2
    ETS1
    STAT1
    FOS

## Top-15 genes

The Top-15 set contains the same Top-10 genes plus:

    NFYC
    IRF8
    TCF3
    PPARD
    SNAI1

This is important because changing both the pseudotime method and the gene
set would make the resulting networks difficult to interpret.

Here, the selected genes remain fixed while the cell ordering changes.

---

# 36. Banjo Settings for the Monocle Comparison

The Monocle-ordered input matrices used the same primary Banjo settings as
the Slingshot analyses:

    observationCount = 959
    minMarkovLag = 1
    maxMarkovLag = 1
    nBestNetworks = 5
    maxTime = 5 m
    computeInfluenceScores = yes

Additional shared settings include:

    DBN mandatory identity lag = 1
    maximum parent count = 5
    discretization policy = i2
    searcher = SearcherGreedy
    evaluator = EvaluatorBDe

Because a lag of one is used, Banjo learns relationships from one position
in the pseudotime-ordered sequence to the following position.

---

# 37. Top-10 Results: Slingshot Versus Monocle 3

## Slingshot Top-10

The original Slingshot Top-10 search returned two non-identical networks.

Best score:

    -2262.2199

The best network contained one additional cross-gene lag-1 relationship:

    ETS1 -> PPARG    +0.0637

The second returned network contained no additional cross-gene edges beyond
the mandatory identity-lag structure.

## Monocle 3 Top-10

The Monocle Top-10 run returned all five requested non-identical networks.

Scores:

    Network 1: -2254.2060
    Network 2: -2254.2629
    Network 3: -2254.4318
    Network 4: -2254.7464
    Network 5: -2255.3654

Banjo examined:

    55,967,761 networks

in:

    5.0 minutes

with:

    17,274 restarts

The best Monocle Top-10 network contained these additional cross-gene
lag-1 relationships:

    ETS1  -> PPARG    +0.0066
    ETS1  -> BCL6     +0.4659
    PRDM1 -> BCL6     +0.5064
    ETS1  -> MAFB     +0.1014
    PPARG -> SNAI3    -0.0278
    STAT1 -> EGR2     +0.0070
    PRDM1 -> ETS1     +0.4147
    ETS1  -> STAT1    -0.2147

The exact directed edge shared by the two best Top-10 networks is:

    ETS1 -> PPARG

The much denser Monocle Top-10 result demonstrates that the inferred network
can be sensitive to the pseudotime ordering supplied to the dynamic Bayesian
network procedure.

---

# 38. Top-15 Results: Slingshot Versus Monocle 3

## Slingshot Top-15

The Slingshot Top-15 analysis returned five non-identical networks.

Scores:

    Network 1: -3889.6251
    Network 2: -3889.8552
    Network 3: -3891.3869
    Network 4: -3891.6365
    Network 5: -3892.2044

The best Slingshot Top-15 network contained the following additional
cross-gene lag-1 relationships:

    ETS1  -> PPARG    +0.0637
    IRF8  -> BCL6     +0.2900
    PPARG -> BCL6     +0.1870
    IRF8  -> MAFB     +0.0667
    SNAI1 -> SNAI3    +0.5171
    IRF8  -> EGR2     approximately -0.0008
    IRF8  -> ETS1     +0.2901
    PPARG -> ETS1     +0.4083
    ETS1  -> NFYC     +0.1532
    ETS1  -> IRF8     -0.2595
    PPARG -> TCF3     -0.1029
    PRDM1 -> SNAI1    +0.5137

## Monocle 3 Top-15

The Monocle Top-15 analysis also returned five non-identical networks.

Scores:

    Network 1: -3863.1410
    Network 2: -3863.2844
    Network 3: -3864.0546
    Network 4: -3864.8857
    Network 5: -3865.7966

Banjo examined:

    62,620,951 networks

in:

    5.0 minutes

with:

    14,909 restarts

The best Monocle Top-15 network contained:

    ETS1  -> PPARG    +0.0066
    PPARD -> BCL6     +0.2185
    ETS1  -> MAFB     +0.1014
    IRF8  -> SNAI3    +0.2315
    STAT1 -> EGR2     +0.0070
    IRF8  -> ETS1     +0.0059
    ETS1  -> STAT1    -0.2147
    ETS1  -> NFYC     +0.1889
    ETS1  -> IRF8     -0.2991
    IRF8  -> PPARD    -0.0839
    PPARG -> SNAI1    -0.1958
    BCL6  -> SNAI1    +0.1969

The following directed edges occur in both best Top-15 networks:

    ETS1 -> PPARG
    IRF8 -> ETS1
    ETS1 -> NFYC
    ETS1 -> IRF8

The direction and sign also agree for these shared edges:

    ETS1 -> PPARG    positive
    IRF8 -> ETS1     positive
    ETS1 -> NFYC     positive
    ETS1 -> IRF8     negative

Other inferred edges depend strongly on the pseudotime ordering.

For example:

Slingshot inferred:

    IRF8  -> BCL6
    PPARG -> BCL6

whereas Monocle inferred:

    PPARD -> BCL6

Slingshot inferred:

    IRF8 -> MAFB

whereas Monocle inferred:

    ETS1 -> MAFB

These differences provide direct evidence that downstream dynamic Bayesian
network inference is sensitive to the pseudotime trajectory supplied as the
temporal ordering.

---

# 39. Interpretation of the Pseudotime Comparison

The comparison supports two conclusions.

First, both Slingshot and Monocle 3 produced pseudotime coordinates that
positively follow the known experimental-time progression.

The correlations were:

    Slingshot: rho approximately 0.718
    Monocle 3: rho = 0.6367

Second, the inferred Banjo networks were not identical.

Several relationships were stable across pseudotime methods, particularly
in the 15-gene analysis, while other edges changed substantially.

The shared edges can be viewed as comparatively robust to the choice of
pseudotime method, whereas method-specific edges should be interpreted more
cautiously.

The absolute Banjo BDe scores should not by themselves be interpreted as
proof that one pseudotime method is biologically superior to the other.
Within each Banjo search, the score ranks candidate networks for that
particular ordered dataset.

The purpose of this comparison is therefore network sensitivity and
agreement, not simply selecting the pseudotime method with the numerically
highest Banjo score.

The Top-10 and Top-15 genes were selected using connectivity information from
the supplied Gold Standard. Consequently, the same Gold Standard should not
be treated as a completely independent validation set for the selected-gene
analysis without acknowledging this selection dependence.

---

# 40. Reproducing the Monocle Banjo Runs

First compile Modified Banjo if necessary:

    ./compile_banjo.sh

Run the Monocle Top-10 network search:

    ./run_monocle3_top10.sh

Run the Monocle Top-15 network search:

    ./run_monocle3_top15.sh

Reproduced outputs are written to:

    results/reproduced/monocle3/top10/
    results/reproduced/monocle3/top15/

These directories are separate from the preserved original Monocle results:

    results/monocle3/top10/
    results/monocle3/top15/

This prevents a reproducibility test from overwriting the network files used
for the reported comparison.

To regenerate the Monocle trajectory and Banjo input ordering from the
original expression matrix, use:

    Rscript scripts/pseudotime/10_run_monocle3.R
    Rscript scripts/pseudotime/11_make_monocle3_banjo_inputs.R

These R scripts expect the original `rna.csv` expression matrix to be
available in the working directory.

---

# 41. Updated Comparison Workflow

The expanded workflow is:

    Original single-cell RNA expression
                    |
          +---------+---------+
          |                   |
          v                   v
      Slingshot           Monocle 3
          |                   |
          v                   v
     pseudotime          pseudotime
      ordering            ordering
          |                   |
          +---------+---------+
                    |
              same 959 cells
                    |
              +-----+-----+
              |           |
              v           v
           Top 10       Top 15
              |           |
              v           v
        Modified Banjo Modified Banjo
              |           |
              v           v
          compare networks across
           pseudotime methods

This design keeps the selected cells and genes fixed while changing the
pseudotime ordering method, allowing the downstream Banjo networks to be
compared directly.
