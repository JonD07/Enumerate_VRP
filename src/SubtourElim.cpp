/*
 * SubtoureLim.cpp
 *
 *  Created on: Jan 25, 2022
 *      Author: jonathan
 */

#include "SubtourElim.h"


Subtour::Subtour() {
}

Subtour::~Subtour() {
}

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.  Result is returned in 'tour', and length is
// returned in 'tourlenP'.
void Subtour::findsubtour(int n, double** sol, int* tourlenP, int* tour) {
	// Function variables
	int bestind, bestlen;
	int i, node, len, start;
	bestlen = n+1;
	bestind = -1;
	start = 0;
	node = 0;

	// Track which stops have been seen so far
	bool* seen = new bool[n];
	for(i = 0; i < n; i++) {
		seen[i] = false;
	}

	// Find all sub-tours, record which one is the best
	while(start < n) {
		// Find a next node that we haven't seen yet
		for(node = 0; node < n; node++) {
			if (!seen[node])
				break;
		}

		// If we have seen every node, stop algorithm
		if(node == n) {
			break;
		}

		for(len = 0; len < n; len++) {
			// Record this stop and mark it as seen
			tour[start+len] = node;
			seen[node] = true;

			// Find next stop on tour
			for(i = 0; i < n; i++) {
				if(sol[node][i] > 0.5 && !seen[i]) {
					node = i;
					break;
				}
			}

			// Check to see if we have traversed the entire tour
			if(i == n) {
				// Traversed this tour, (maybe) update best tour length found
				len++;
				if(len < bestlen) {
					bestlen = len;
					bestind = start;
				}

				// Update start marker
				start += len;
				break;
			}
		}
	}

	// Update tour array so that best sub-tour is at the beginning
	for(i = 0; i < bestlen; i++) {
		tour[i] = tour[bestind+i];
	}
	*tourlenP = bestlen;

	delete[] seen;
}


SubtourElim::SubtourElim(GRBVar** xvars, int xn) {
	vars = xvars;
	n = xn;

}

SubtourElim::~SubtourElim() {
	// TODO Auto-generated destructor stub
}


void SubtourElim::callback() {
	try {
		if(where == GRB_CB_MIPSOL) {
			// Found an integer feasible solution - does it visit every node?
			double **x = new double*[n];
			int *tour = new int[n];
			int i, j, len;

			// Get current solution from model
			for(i = 0; i < n; i++) {
				x[i] = getSolution(vars[i], n);	// Allocates memory!
			}

			// Check length of smallest sub-tour in solution
			sb.findsubtour(n, x, &len, tour);

			// If smallest sub-tour does not contain all stops, create a lazy constraint (cut off this node in the solver)
			if(len < n) {
				// Add sub-tour elimination constraint
				GRBLinExpr expr = 0;
				for (i = 0; i < len; i++) {
					for (j = i+1; j < len; j++) {
						expr += vars[tour[i]][tour[j]];
					}
				}
				addLazy(expr <= len-1);
			}

			// Cleanup memory
			for(i = 0; i < n; i++) {
				delete[] x[i];
			}
			delete[] x;
			delete[] tour;
		}
	}
	catch(GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch(...) {
		cout << "Error during callback" << endl;
	}
}
