/*
 * SubtoureLim.h
 *
 *  Created on: Jan 25, 2022
 *      Author: jonathan
 *
 *  Subtour elimination callback.  Whenever a feasible solution is found, find
 *  the smallest subtour, and add a subtour elimination constraint if the tour
 *  doesn't visit every node.
 */

#ifndef SRC_SUBTOURELIM_H_
#define SRC_SUBTOURELIM_H_

#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

class Subtour {
public:
	Subtour();
	virtual ~Subtour();

	void findsubtour(int n, double** sol, int* tourlenP, int* tour);
};

class SubtourElim : public GRBCallback {
public:
	SubtourElim(GRBVar** xvars, int xn);
	virtual ~SubtourElim();

	GRBVar** vars;
	int n;

protected:
	Subtour sb;

	void callback();
};

#endif /* SRC_SUBTOURELIM_H_ */
