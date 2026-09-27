Machine Learning for Biology | Day 63
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

A correction first. Day 62 said the ten populations were sorted by flow cytometry. They were not. They were pulled out of one donor's blood with antibody-coated beads, and flow cytometry only came in afterwards, to check how pure each tube was. Other published work describes this dataset as FACS-sorted too, which is how the slip travels.

The correction matters because it shows where the answer key came from. Every label here was made with an antibody against a protein on the cell surface: CD4, CD8, CD14, CD19, CD25, CD34, CD45RA, CD45RO, CD56.

The algorithm will never see a protein. It sees RNA, counted from the tail end of each message. Proteins and messages usually agree. Here are three places they do not.

One label is close to invisible. Naive and memory T cells were separated using CD45RA versus CD45RO. Those are not two genes. They are two versions of one gene, PTPRC, made by cutting different pieces out of the same message, and those pieces sit near the front, far from the tail this chemistry reads. A 2024 study tested exactly this on tail-end data from blood cells: even with very deep sequencing, fewer than one cell in ten showed either version.

Some labels overlap by definition. The helper T cells were selected on CD4 alone. Naive, memory and regulatory T cells all carry CD4, so the helper tube holds cells matching three other labels. The two CD8 populations nest the same way, with a twist: some of the most experienced killer T cells switch CD45RA back on, so "naive" here is not purely naive.

One label is mostly something else. The CD34+ tube was 45% pure by the dataset's own flow check.

This is the Chapter 7 version of Day 48. There, the gene did not follow the organism. Here, the label does not fully follow the RNA.

So, a prediction, written down before any code runs. Not a blind one: others have clustered these cells before. B cells and monocytes will come back largely intact. NK cells will partly blur into cytotoxic T cells, because both run the same killing programme. The six T cell labels will not separate cleanly, and the helper label cannot separate from the populations nested inside it. If it does, something other than biology is doing the separating. The CD34+ label will mostly hold together, as it did in a 2020 re-analysis. That is odd for a tube that is 45% pure, and odd results get checked against the run first.

One decision, made now rather than after a disappointing score. Results get scored against the ten labels, as promised, and against six broad lineages: B cells, monocytes, NK cells, CD34+ progenitors, CD4 T cells and CD8 T cells. Both are declared today, so neither can be chosen later to flatter the result.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #Bioinformatics #Immunology #SingleCellRNAseq
