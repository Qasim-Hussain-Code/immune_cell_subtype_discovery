Machine Learning for Biology | Day 68
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

The answer key is open. Here is how the computer's groups compare with the lab's labels.

I scored every grouping against two keys: the ten labels, one per tube, and the six broad families from Day 63, which lump related tubes together (the four CD4 T cell tubes become one family).

The score runs from 0 (random grouping) to 1 (perfect agreement).

Every grouping beat random easily. I shuffled the labels 1,000 times, and not one shuffle scored above 0.003.

That was the warm-up. The real test was the "camera settings" baseline from Day 66: the same sorting method, k-means, but given only three numbers about how well each cell was measured, and nothing about the cell itself.

The rule, fixed in advance, was strict: gene expression would only win if the whole 95% interval for its lead over the baseline sat above zero.

First, k-means had to decide for itself how many groups there were. It chose two. At two groups, gene expression lost.

Ten labels: gene expression 0.01, camera settings 0.04
Six families: gene expression 0.05, camera settings 0.13

The reason is almost boring. Gene expression split off a small group, mostly monocytes. The camera settings split off the cells that gave the most RNA. The CD34+ tube gave far more RNA per cell than any other, so 89% of it landed in that group. Against the labels, the camera-settings split happens to score better.

Then I told k-means to make ten groups, the same number as the lab's labels. This time gene expression won clearly.

Ten labels: gene expression 0.44, camera settings 0.08
Six families: gene expression 0.55, camera settings 0.09

Leiden, the second sorting method from Day 66, was also left to pick its own number of groups. It settled on 23 and beat two-group k-means on both keys: 0.58 on the ten labels, 0.42 on the six families.

Against ten-group k-means, it was a split decision. Leiden scored higher on the ten labels (0.58 vs 0.44), and k-means scored higher on the six families (0.55 vs 0.42). That comparison was not planned in advance, so it stays a description, not a verdict.

So the headline fits in two sentences. When k-means has to decide how many groups exist, it gets the count wrong, and the camera settings beat it. When it is told the count, it recovers a large part of the lab's structure from RNA alone.

Day 66 called this "count wrong but sort right". That is what happened.

Code and data pipeline: https://lnkd.in/gvq7DXwV

hashtag#MachineLearningForBiology hashtag#MachineLearning hashtag#ArtificialIntelligence hashtag#Bioinformatics hashtag#Immunology hashtag#SingleCellRNAseq
