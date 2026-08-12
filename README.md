# Slingshot + Modified Banjo Gene Regulatory Network Analysis

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

    git clone <REPOSITORY-URL>

Then enter the repository:

    cd slingshot-modified-banjo

Replace `<REPOSITORY-URL>` with the actual private GitHub repository URL.

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

# 27. Quick Start for a New User

For someone using a new Ubuntu or WSL environment, the shortest procedure is:

    sudo apt update
    sudo apt install -y openjdk-17-jdk graphviz git

Then clone the repository:

    git clone <REPOSITORY-URL>
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
