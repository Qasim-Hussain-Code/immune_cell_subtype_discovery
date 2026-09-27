Machine Learning for Biology | Day 65
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

The rules, fixed before any code runs.

The data: the ten bead-enriched populations released by 10x Genomics alongside Zheng and colleagues, Nature Communications, 2017. Every cell the original pipeline called a cell, pooled into one table. Nothing removed on top by thresholds of my choosing.

The answer key goes in a locked drawer. At the very first step the labels move to a separate file, every cell gets a random code in place of its name, and the order is shuffled. No script that groups cells may open that file. It opens once, after the groups are final and committed to the repository.

The input is gene expression and nothing else. The settings are borrowed, not tuned: the processing recipe from the 2017 paper, which the Scanpy library ships under that paper's name. The 1,000 most variable genes, squeezed to 50 axes. No correction for sequencing batch, because here batch and label are the same thing. Tomorrow explains.

The number of groups comes from the data. Anything from 2 to 20 is allowed. The winner is the number where cells sit most clearly inside their own group rather than the next one over, a measure called the silhouette. It will never be set to ten because the key says ten. One extra run, labelled as borrowing from the key, forces ten groups, so a wrong count can be told apart from a wrong grouping.

The score: the adjusted Rand index. 0 is what random grouping scores on average, 1 is perfect agreement. It comes with a 95% confidence interval, and the 0 gets checked by shuffling the key 1,000 times. Scored twice, as declared on Day 63: against the ten labels and against the six lineages.

The baseline to beat: the same model, fed three numbers per cell that describe the measurement rather than the cell. Molecules captured, genes detected, and the share of molecules from mitochondrial genes. If gene expression cannot clearly beat that, the groups describe the machine.

Names before answers. Before the key opens, every group gets a name from its marker genes, twice: once by a fixed rule written down now, once by me. Both get committed, then scored.

Every setting, down to the random seed, goes into the repository before the analysis starts, along with a pass or fail test for each prediction from Day 63.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #Bioinformatics #Immunology #SingleCellRNAse
