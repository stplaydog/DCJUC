#pragma once

#include "graph.h"
#include "list.h"

//the problem with bnb is that because of the ambiguity, it's hard to compute the bound based on the original way.
//but we can compute the bound using the following way:
//Upperbound: continue search in one direction to get a local result, and this local result is the upperbound.
//Lowerbound: connect the edge with pairs and the other by some other ways, and using the formular (d_12+d_13+d_23)/2 as lower bound.

void bnb(pg g);
