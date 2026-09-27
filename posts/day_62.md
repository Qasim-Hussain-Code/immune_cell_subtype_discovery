Machine Learning for Biology | Day 62
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

Draw a tube of blood, spin it down, and the thin white layer in the middle holds T cells, B cells, NK cells, monocytes and dendritic cells, tangled together. 

Under a microscope, several of them look nearly identical.

Every chapter in this series so far has told the model the correct answer during training.

Logistic regression knew which genomes carried resistance. 

The network in Chapter 6 knew which peptides provoked a response. 

Day 2 of this series raised a different possibility and then left it alone for six chapters: biologists already sort cells into types without any of that supervision.

This chapter is where that gets tested properly.

The question: given only a cell's gene expression, nothing else, can an algorithm that has never seen a cell type label recover the immune cell populations a lab already knows are there?

The check is what makes this answerable rather than just interesting. 

The cells used here were sorted into ten populations by an entirely different method from transcriptomics. 

Those labels never touch the clustering. They exist only to be compared against afterward.

One limitation, named now. Clustering finds structure, not names. It can group T cells together without knowing to call them T cells, and it can just as easily split one real population into two or merge two distinct ones, especially where biology itself is gradual rather than discrete.

Whatever comes out of this still needs a human, and marker genes, to say what each group actually is.

hashtag#MachineLearningForBiology hashtag#MachineLearning hashtag#ArtificialIntelligence hashtag#Bioinformatics hashtag#Immunology hashtag#SingleCellRNAseq
