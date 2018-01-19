# smDBSCAN

smDBSCAN.m: applies the DBSCAN principle of density-reachability to group subsequent observations of the same immobile single-molecule over time.

DBSCAN_fwd_rev_cluster.m: performs pure forward and reverse-looking clustering using the DBSCAN principle followed by fusion via maximum core point overlap

DBSCAN_pot_link.m: finds points that potentially belong into one group

DBSCAN_hard_thresh.m: reduces the list of potential links to those that reside within the expected time window around each point

DBSCAN_point_score.m: calculates the density score for each point

DBSCAN_group.m: implements point clustering based on the principle of DBSCAN (Density-Based Spatial Clustering of Applications with Noise). DBSCAN is a data clustering algorithm proposed by Martin Ester, Hans-Peter Kriegel, Jörg Sander and Xiaowei Xu in 1996.

show_me_my_cluster.m: visualizes the results of the clustering process
