#pragma once

#include "tree.h"

/*
 * this implementation is based on the following description:
 * Our new initialization method put global information to use through adequate subgraphs
 * and thus meshes well with the refinement phase of the scoring procedure, which
 * uses adequate subgraphs to compute medians. We initialize internal nodes progressively:
 * in order for an internal node to be a candidate for initialization, two of its three
 * neighbors must be already initialized (or leaves). The third, while typically not be initialized,
 * is not devoid of information, so our method summarizes the data available in
 * the third subtree (rooted at the uninitialized neighbor) into a set of weighted adjacencies.
 * Thus information used in initializing a node consists of two 0-1 sets of adjacencies
 * from the two initialized neighbors and one weighted set of adjacencies from the third
 * neighbor. A suitable choice of the node to be initialized is thus a generalized version
 * of a median, one that takes into account the weighted nature of adjacencies in the third
 * node.
 *
 * we consider only
 * the directive nodes in the perspective, that is, initialized nodes or leaves connected to
 * the node of interest via a path of uninitialized nodes. In Fig. 1, nodes a, b, and c are
 * initialized and connected to node 1 through paths of uninitialized nodes and so are
 * directive nodes.
 **/

void
run_initialization();
void
init_directive_node();
