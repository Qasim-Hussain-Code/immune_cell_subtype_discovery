Machine Learning for Biology | Day 64
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

The model for the final chapter is k-means. Here is how it works, and why it got the job.

Each cell starts as counts for thousands of genes. Many of those genes rise and fall together, so they can be squeezed into 50 combined axes that keep the main patterns. That squeeze is called principal component analysis. Now picture every cell as a point in a space with those 50 directions. Similar cells sit close together.

k-means drops k markers into that space. Every cell joins its nearest marker. Each marker moves to the middle of the cells that joined it. Cells switch to whichever marker is now nearest. Repeat until nothing moves.

That is the whole algorithm. No answer key anywhere. The one instruction it gets is k, the number of groups.

Chapter 5's support vector machine searched for the widest gap between two outcomes it had been told about. k-means has no outcomes. It decides where the gaps are on its own.

It is not a toy choice for this data. The 2017 paper that produced these cells clustered with k-means, on 50 principal components, and for its 68,000-cell blood sample it set k to 10 after reading an error curve. Reading a curve like that is a judgement call made after seeing the data. This chapter replaces it with a rule fixed in advance, published tomorrow.

Two reasons it closes the series. It is the plainest version of the idea, the unsupervised counterpart to Chapter 1's logistic regression. And it drags the most important decision into the open. How many groups are there? Most methods answer that quietly, inside a setting. k-means will not start until someone answers it.

One weakness shows up at the start line: where the markers begin can change where they finish. So every run starts from 50 different random positions and keeps the tightest result.

k-means is not what single-cell labs usually reach for. Seurat and Scanpy, two of the most widely used toolkits, default to graph-based clustering: link each cell to its nearest neighbours, then find communities more tightly linked inside than out. The version used here is called Leiden. It needs no k, only a resolution dial that does the same job less visibly. Leiden runs alongside, its dial set by the same rule that picks k.

If the textbook method loses to the field's default, that gets reported too.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #Bioinformatics #Immunology #SingleCellRNAseq
