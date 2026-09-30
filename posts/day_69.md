Machine Learning for Biology | Day 69
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

On Day 63 I wrote down six predictions before any code ran. They were scored on the run told to make ten groups, because they are about ten labels.

Four passed. Two failed.

1. B cells come back intact. Pass. 99.9% of them landed in one group, and that group was 96% B cells.

2. Monocytes come back intact. Pass. 93% of them landed in one group.

3. NK cells partly blur into cytotoxic T cells. Fail. The prediction needed an overlap of at least 0.10. It came out at 0.05. Sharing a killing programme was not enough to mix them: 98% of NK cells sat in one group that was 96% NK cells.

4. The CD34+ tube mostly holds together. Fail, and the most interesting miss. The tube did not scatter into other cell types. It split into three groups of its own, each more than 98% CD34+ cells. But its biggest single group held only 30% of it, far short of the 80% needed. This was the one prediction shaped by the 2020 re-analysis, where most of these cells stayed in one cluster.

5. The six T cell labels do not separate cleanly. Pass. Not one of them got a group holding 80% of its cells at 80% purity.

6. The helper tube cannot be pulled apart from the naive T cells inside it. Pass, and this is the one that matters most. Helper against naive separated with a score of 0.05. Naive against memory, two genuinely different kinds of T cell, separated at 0.49. If the tubes, rather than the cells, had been driving the groups, that order would have flipped. It did not.

Two honest notes. These were not blind predictions: others had clustered these cells before, which made B cells and monocytes safe calls. And the two misses are the reason predictions get written down first: without Day 63 on record, it would be easy to claim afterwards that these results were expected.

hashtag#MachineLearningForBiology hashtag#MachineLearning hashtag#ArtificialIntelligence hashtag#Bioinformatics hashtag#Immunology hashtag#SingleCellRNAseq
