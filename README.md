# ALP
## Author: Newton Campbell
Repository for the ALP (A*, Landmarks, and Polygon Inequality) algorithm 

Point-to-point shortest path distance queries are a core operation in graph analytics. However, preprocessing algorithms that speed up these queries rely on large data structures for reference. In this paper, we discuss the computational challenge introduced by these data structures when using landmark-based preprocessing algorithms on large graphs. We introduce a new heuristic for the A* algorithm that references a data structure of size θ(|L| 2 + |V|), where L represents a set of strategically chosen landmark vertices and V the set of vertices in the graph. This heuristic's benefits are permitted by an approach for computing lower bounds based on generalized polygon inequalities. In this approach, each landmark stores the distances between the landmark and vertices within its graph partition. The heuristic is experimentally compared with a previous landmark heuristic in a fixed-memory environment, as an analog to an embedded system. The new heuristic demonstrates a reduction in overall computational time and memory requirements in this environment.

Further descriptions can be found in the following publications:

Landmark routing for large graphs in fixed-memory environments
N Campbell, MJ Laszlo, S Mukherjee
[2016 IEEE High Performance Extreme Computing Conference (HPEC), 1-7		2016](https://ieeexplore.ieee.org/abstract/document/7761581/)

Algorithmic Foundations of Heuristic Search using Higher-Order Polygon Inequalities
NH Campbell Jr
[Nova Southeastern University		2016](https://nsuworks.nova.edu/cgi/viewcontent.cgi?article=1372&context=gscis_etd&httpsredir=1&referer=/)

Using quadrilaterals to compute the shortest path
N Campbell Jr
[arXiv preprint arXiv:1603.00963	1	2016](https://arxiv.org/abs/1603.00963)

Computing shortest paths using A*, landmarks, and polygon inequalities
N Campbell Jr
[arXiv preprint arXiv:1603.01607](https://arxiv.org/abs/1603.01607)
