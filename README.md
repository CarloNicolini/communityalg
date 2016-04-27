# communityalg
Algorithms and functions in Matlab for community detection in networks. 
Expands BrainConnectivity toolbox.

# Matlab/Octave algorithms and functions

- **ami.m** Returns the adjusted mutual information between two membership vectors.
- **association_score.m** Returns the association score between pairs of communities specified by the graph and membership.
- **asymptotic_modularity.m** Compute asymptotic modularity of a graph with respect to a membership vector.
- **asymptotic_modularity_sum.m**  TODO
- **asymptotic_surprise.m** Compute asymptotic surprise of a graph with respect to a membership vector.
- **clique.m** Generate an adjacency matrix of a clique graph with `n` nodes.
- **cluster_similarity.m** Compare two membership vectors.
- **clustering_entropy.m** Compute the clustering entropy of an agreement matrix as in "Gfeller, Newman, 2006".
- **community_robustness_weighted.m** TODO
- **community_size2memb.m** Convert an array where every element is the size of a clique to the correspondig membership vector.
- **comm_mat.m** Returns the block matrix of a graph and its community structure as membership.
- **compute_surprise.m** Compute the surprise given surprise paramenters.
- **consensus_clustering.m** TODO
- **consensus_clustering_weighted.m** TODO
- **consensus_entropy.m** TODO
- **consensus_robustness.m** TODO
- **correlation_louvain.m** Adaptation of the BCT `community_louvain` method for correlation matrices, as described in MacMahon,2015.
- **count_comm.m** Plot the histogram with the community size given a membership vector.
- **cycle_graph.m** Generate the adjacency matrix of a cycle graph with n nodes.
- **effcommplot.m** TODO
- **generate_agreement.m** Generate the agreement matrix for a given community detection method.
- **generate_agreement_weighted.m** Generate the weighted agreement matrix for a given community detection method.
- **generate_connected_components.m** TODO
- **graph_JS_similarity.m** Compute the quantum Jensen Shannon divergence between two adjacency matrices.
- **graph_laplacian.m** Compute the graph combinatorial Laplacian matrix `L=D-A`.
- **group2membership.m** Convert a cell of arrays representing the nodes in the communities to a membership vector.
- **image_to_network.m** Convert a gray index image to its corresponding adjacency graph.
- **imagesctxt.m** Show a matrix like `imagesc` but with text values of elements displayed on the pixels.
- **isoctave.m** Returns true if using Octave, false if using Matlab.
- **jensen_shannon_sim.m** Returns the Jensen-Shannon symmetrized information theoretic distance between two graph Laplacians.
- **k_regular.m** Generate a `k`-regular graph, a graph where the degree of every vertex is `k`.
- **KL.m** Returns the binary Kullback-Leibler divergence between Bernoulli distribution `p` and `q`.
- **kullback_leibler_sim.m** Returns the Kullback-Leibler divergence between two graph Laplacians.
- **logHyperProbability.m** Compute the logarithm of the hypergeometric probability in base 10.
- **membership2groups.m** Convert a membership vector to a cell of arrays of nodes in every community.
- **membership_agreement.m** DEPRECATE
- **membership_similarity.m** TODO
- **method_best.m** Functio handle to the community detection method that returns the best value over a set of repetitions.
- **modularity.m** Returns the modularity of a graph with respect to a membership vector.
- **nearcorr.m** Returns the nearest correlation matrix of a square matrix. Implementation by Nick Higham.
- **nearestSPD.m** Returns the positive definite matrix of a square matrix. Implementation by Nick Higham.
- **norm_conf_mat.m** DEPRECATE
- **number_connected_components.m** Returns the number of connected components of a graph.
- **number_of_edges.m** Returns the number of edges of a binary or weighted graph.
- **paco.m** Function handle to the MEX implementation of PACO.
- **partition_params.m** Returns the partition parameters for use with `compute_surprise`.
- **quantum_density.m** Returns the quantum density of a graph, `Brauenstein et al.` "Ann. of Combinatorics, 10, no 3 (2006), 291-317."
- **reindex_membership.m** Transform a membership vector to have community indices sorted by community size from `1` to `|C|`
- **reorder_membership.m** Linearize a membership vector to have continuous indices of communities from `1` to `|C|`
- **ring_of_cliques.m** Returns a  network ring of cliques, with given number of cliques and clique size and its membership.
- **ring_of_custom_cliques.m** Returns a network ring of cliques, with given size of cliques specified as input and its membership.
- **rmtdecompose.m** Returns the Random Matrix Theory decomposition of a correlation matrix.
- **robustness_configuration_interp_und.m** 
- **robustness_edge_weight_und.m**
- **run_cluster_similarity.m**
- **rwalkent.m** Returns the random walk entropy of a graph as in `Estrada et al.` Walk entropies in graphs, "Linear Algebra and its Applications 443 (2014) 235–244"
- **significance.m** Returns the significance of a graph partitioning, `Traag (2013)`.
- **smi.m** Returns the standardized mutual information of two membership vectors.
- **sort_group_by_size.m**  Sort community groups by size.
- **star_of_custom_cliques.m** Returns a star of cliques, every clique is connected to all other cliques with one edge.
- **surprise.m** Returns the surprise of a graph partitioning.
- **threshold_by_giant_component.m** Returns the threshold over which the graph has more than one connected component.
- **threshold_by_num_edges.m** Returns a graph thresholded to have a specific number of edges.
- **vonneumann_entropy.m** Returns the VonNeumann quantum entropy of graph `Brauenstein et al.` "Ann. of Combinatorics, 10, no 3 (2006), 291-317."
- **write_brainet.m** Write a graph with coordinates of nodes and membership to Brainet format.
- **write_brainet_community.m** TODO
- **writetoEdgesList.m** TODO
- **writetoPAJ_labels.m**
- **writetoPAJ_labels_coords.m**

# Python algorithms
- **find_intersections.py** Find the intersections of edges in a bipartite graph (TODO)

# Datasets
This repository contains the 638 areas template used by Crossley (2013). It consists of two files:

1. `template.nii` is the NIFTI file describing the template, MNI space
2. `template_638.txt` is the description of anatomical a


Additionally the file `template_638_coords_abbr.txt`  is organized as follows

| Node |  X | Y | Z | Label | Abbr|
|------|----|---|---|-------|-----|


and differently from `template_638.txt` the NodeID starts from 0 to 637 (for indexing with Python).

