Machine Learning for Biology | Day 66
Chapter 7: Immune cell subtype discovery from single-cell RNA-seq

Here are the ways this chapter can go wrong. 

1. It might sort by tube, not by cell.
Imagine ten piles of holiday photos, each taken with a different camera. Ask someone to sort them by who is in them, and they might sort them by camera instead, because one camera's photos are brighter and another's are blurrier.

The same risk exists here. Each cell type arrived in its own tube and was measured separately, some far more thoroughly than others.
Check: I give the computer only the "camera settings", three numbers about how well each cell was measured and nothing about the cell itself. The real data has to clearly beat that. If it cannot, the groups are about the camera, not the people.

Second check: one tube, the helper T cells, was collected with a broad net, so it already holds many cells from another tube, the naive T cells. Those two should be hard to separate. If the computer separates them more easily than two tubes that really are different, it is sorting by tube.

2. It might count wrong but sort right.
The computer decides how many groups there are by looking for the cleanest splits. In blood, the cleanest splits are the big ones, so it may find three or four groups where the lab defined ten.

Check: I publish its score for every group count from 2 to 20, plus one extra run where I tell it to make ten.

4. It never says "I am not sure".
Every cell must go into some group, even cells caught halfway through changing from one type to another. And where the computer starts is random, which can change where it ends.

Check: I run it 20 times from different starting points, and compare it with a second sorting method, Leiden, that many labs use.

5. A strange group is not a discovery.
A group that matches no label might be a new cell type. It is just as likely to be stray cells, two cells stuck together, or dying cells.

Check: any such group is called "unresolved" and its genes are listed. No new cell types will be claimed.

The labels stay locked away until every group is final.

#MachineLearningForBiology #MachineLearning #ArtificialIntelligence #Bioinformatics #Immunology #SingleCellRNAseq
