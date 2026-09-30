Machine Learning for Biology | Day 67
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

Before the answer key opened, everything below was frozen and committed to the repository.

94,655 blood cells went in. Each one was just a list of gene counts. No names attached.

The first job was to decide how many groups to make. The rule from Day 65 tried every number from 2 to 20 and kept the one where cells sat most clearly inside their own group.

It picked 2.

Not ten. Not six. Two. One group of 3,521 cells full of monocyte genes. One group of 91,134 cells holding everything else.

Day 66 warned about exactly this: in blood, the cleanest split is a coarse one. The clarity score for two groups was 0.63. For ten groups it was 0.21.

So the extra run, the one told to make exactly ten groups, now matters a lot. On their marker genes, its groups look like real cell families: B cells, NK cells, monocytes, several kinds of T cells, and three separate groups the naming rule called CD34+ progenitors.

Leiden, the kind of method most labs use, settled on 23 groups.

Then the 20 reruns from different random starting points. Leiden gave much the same answer each time. k-means did not. At two groups, 12 of the 20 reruns found a different split, and a slightly tighter one, than the run fixed in advance. The rules say the pre-committed run is the one that gets scored, so it is. The other answer gets reported next to it.

Finally, every group was named from its genes, twice. A fixed rule called the big group "CD8 T cells". The second naming, which Day 65 said I would do myself. It worked from the same marker genes, before the key opened, and called the big group "unresolved": a mix of several cell types. That change is logged in the repository as a deviation from the plan.

All of it was locked in before anyone looked at a label.

Tomorrow, the key opens.

Code and data pipeline: https://lnkd.in/gvq7DXwV

hashtag#MachineLearningForBiology hashtag#MachineLearning hashtag#ArtificialIntelligence hashtag#Bioinformatics hashtag#Immunology hashtag#SingleCellRNAseq
