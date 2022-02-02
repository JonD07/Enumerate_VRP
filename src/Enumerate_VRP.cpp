/* Copyright 2021, Gurobi Optimization, LLC */

/* Solve a traveling salesman problem on a randomly generated set of
   points using lazy constraints.   The base MIP model only includes
   'degree-2' constraints, requiring each node to have exactly
   two incident edges.  Solutions to this model may contain subtours -
   tours that don't visit every node.  The lazy constraint callback
   adds new constraints to cut them off. */

#include "gurobi_c++.h"
#include "SubtourElim.h"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>

using namespace std;

// To change the number of stops, edit NUM_STOPS and update the arrays below
#define NUM_STOPS	20
#define NUM_VHCLE	15
#define DEPOT_INDEX	NUM_STOPS

/*
 * Coordinates of the stops. By convention, the last entry is the depot.
 */
const double STOP_COORDS[NUM_STOPS + 1][2] = {
	{236, -54},
	{231, -17},
	{317, -32},
	{33, -126},
	{44, -22},
	{182, -26},
	{173, 25},
	{19, 184},
	{221, -39},
	{227, 170},
	{-25, -16},
	{273, 43},
	{368, -59},
	{-10, -195},
	{260, -109},
	{9, -37},
	{334, 176},
	{291, 117},
	{95, 123},
	{146, 120},
	{0.0, 0.0}	// Designated depot
};

/*
 * Demand at each stop. It was found manually that any feasible tour can have at
 * most 5 stops (not including the depot) due to truck capacity. If the demands
 * change then one must also re-evaluate how long the longest feasible tour is
 * and modify main() as needed.
 */
const double STOP_DEMAND[NUM_STOPS] = {
	0.25,
	0.33,
	0.39,
	0.40,
	0.27,
	0.70,
	0.28,
	0.43,
	0.50,
	0.22,
	0.21,
	0.68,
	0.16,
	0.19,
	0.22,
	0.38,
	0.26,
	0.29,
	0.17,
	0.31
};

string itos(int i) {
	stringstream s;
	s << i;
	return s.str();
}

double distance(double* x, double* y, int i, int j);
double distance(double x1, double y1, double x2, double y2);

// Euclidean distance between points 'i' and 'j'.
double distance(double* x, double* y, int i, int j) {
	double dx = x[i]-x[j];
	double dy = y[i]-y[j];

	return sqrt(dx*dx+dy*dy);
}

// Euclidean distance between points (x1, y1) and (x2, y2)
double distance(double x1, double y1, double x2, double y2) {
	double dx = x1-x2;
	double dy = y1-y2;

	return sqrt(dx*dx+dy*dy);
}

bool in_array(int l, int arr[], int arr_size) {
	for(int i = 0; i < arr_size; i++) {
		if(arr[i] == l) {
			return true;
		}
	}

	return false;
}


//vector<bool*> tours;
//vector<int*> tsp_sol;
//vector<double> tour_dists;

void solveTSP(int stops_to_visit[], int n, vector<bool*>& tours, vector<int*>& tsp_sol, vector<double>& tour_dists) {
	// Program variables
	GRBEnv *env = NULL;
	GRBVar **vars = NULL;
	int i;

	// Create an N x N array for decision variables
	vars = new GRBVar*[n];
	for(i = 0; i < n; i++) {
		vars[i] = new GRBVar[n];
	}

	try {
		int j;

		// Create Gurobi environment
		env = new GRBEnv();
		GRBModel model = GRBModel(*env);

		// Must set LazyConstraints parameter when using lazy constraints
		model.set(GRB_IntParam_LazyConstraints, 1);

		// Create the binary decision variables (link to model and store in NxN array)
		for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++) { // STOP_COORDS
				vars[i][j] =
						model.addVar(
								0.0,
								1.0,
								distance(STOP_COORDS[stops_to_visit[i]][0],
										STOP_COORDS[stops_to_visit[i]][1],
										STOP_COORDS[stops_to_visit[j]][0],
										STOP_COORDS[stops_to_visit[j]][1]),
								GRB_BINARY,
								"x_"+itos(i)+"_"+itos(j));
				vars[j][i] = vars[i][j];
			}
		}

		// Degree-2 constraints, each stop must have an edge in and an edge out
		for (i = 0; i < n; i++) {
			GRBLinExpr expr = 0;
			for (j = 0; j < n; j++)
				expr += vars[i][j];
			model.addConstr(expr == 2, "deg2_"+itos(i));
		}

		// Forbid edge from node back to itself
		for (i = 0; i < n; i++) {
			vars[i][i].set(GRB_DoubleAttr_UB, 0);
		}

		// Set callback function
		SubtourElim cb = SubtourElim(vars, n);
		model.setCallback(&cb);

		// Optimize model
		model.optimize();

		// Extract solution
		if (model.get(GRB_IntAttr_SolCount) > 0) {
			double **sol = new double*[n];
			for (i = 0; i < n; i++)
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

			int* tour = new int[n];
			int len;

			Subtour sb;
			sb.findsubtour(n, sol, &len, tour);
			assert(len == n);

			cout << endl << "Tour: ";
			for (i = 0; i < len; i++) {
				cout << stops_to_visit[tour[i]] + 1 << " ";
			}
			cout << endl;
			cout << "Length: " << model.get(GRB_DoubleAttr_ObjVal) << endl << endl;

			// Add this to our solution set
			bool* bn_sol = new bool[NUM_STOPS];
			for(int i = 0; i < NUM_STOPS; i++) {
				if(in_array(i, stops_to_visit, n)) {
					bn_sol[i] = 1;
				}
				else {
					bn_sol[i] = 0;
				}
			}
			tours.push_back(bn_sol);
			int* tspsol = new int[NUM_STOPS + 1];
			for(int i = 0; i <= NUM_STOPS; i++) {
				if(i < len) {
					tspsol[i] = stops_to_visit[tour[i]];
				}
				else {
					tspsol[i] = -1;
				}
			}
			tsp_sol.push_back(tspsol);
			tour_dists.push_back(model.get(GRB_DoubleAttr_ObjVal));

			for (i = 0; i < n; i++)
				delete[] sol[i];
			delete[] sol;
			delete[] tour;
		}

	}
	catch(GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch(...) {
		cout << "Error during optimization" << endl;
	}

	// Clean-up dynamic memory
	for (i = 0; i < n; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete env;
}

int main(int argc, char *argv[]) {
	int comb_count = 0, route_count = 0;
	vector<bool*> tours;
	vector<int*> tsp_sol;
	vector<double> tour_dists;

	// Try every combination of stops. We assume that no more that 5 stops can be
	// on a single route due to capacity constraints.
	for(int a = 0; a < NUM_STOPS; a++, comb_count++) {
		// A single stop is feasible, print distances
		{
			route_count++;
			printf("%d\n", a);
			double dist = distance(STOP_COORDS[a][0], STOP_COORDS[a][1], STOP_COORDS[DEPOT_INDEX][0], STOP_COORDS[DEPOT_INDEX][1]) * 2;
			cout << endl << "Tour: " << a+1 << " " << DEPOT_INDEX+1 << endl;
			cout << "Length: " << dist << endl << endl;

			// Add this to our solution set
			bool* sol = new bool[NUM_STOPS];
			for(int i = 0; i < NUM_STOPS; i++) {
				sol[i] = 0;
			}
			sol[a] = 1;
			tours.push_back(sol);
			int* tspsol = new int[NUM_STOPS + 1];
			for(int i = 0; i <= NUM_STOPS; i++) {
				tspsol[i] = -1;
			}
			tspsol[0] = NUM_STOPS;
			tspsol[1] = a;
			tsp_sol.push_back(tspsol);
			tour_dists.push_back(dist);
		}

		for(int b = a+1; b < NUM_STOPS; b++, comb_count++) {
			// Determine if this combo is feasible
			if((STOP_DEMAND[a] + STOP_DEMAND[b]) <= 1.0) {
				route_count++;
				printf("%d %d\n", a, b);

				// "Solve" TSP
				double dist = distance(STOP_COORDS[a][0], STOP_COORDS[a][1], STOP_COORDS[DEPOT_INDEX][0], STOP_COORDS[DEPOT_INDEX][1]) +
						distance(STOP_COORDS[a][0], STOP_COORDS[a][1], STOP_COORDS[b][0], STOP_COORDS[b][1]) +
						distance(STOP_COORDS[b][0], STOP_COORDS[b][1], STOP_COORDS[DEPOT_INDEX][0], STOP_COORDS[DEPOT_INDEX][1]);
				cout << endl << "Tour: " << a+1 << " " << b+1 << " " << DEPOT_INDEX+1 << endl;
				cout << "Length: " << dist << endl << endl;

				// Add this to our solution set
				bool* sol = new bool[NUM_STOPS];
				for(int i = 0; i < NUM_STOPS; i++) {
					sol[i] = 0;
				}
				sol[a] = 1;
				sol[b] = 1;
				tours.push_back(sol);
				int* tspsol = new int[NUM_STOPS + 1];
				for(int i = 0; i <= NUM_STOPS; i++) {
					tspsol[i] = -1;
				}
				tspsol[0] = NUM_STOPS;
				tspsol[1] = a;
				tspsol[2] = b;
				tsp_sol.push_back(tspsol);
				tour_dists.push_back(dist);
			}

			for(int c = b+1; c < NUM_STOPS; c++, comb_count++) {
				// Determine if this combo is feasible
				if((STOP_DEMAND[a] + STOP_DEMAND[b] + STOP_DEMAND[c]) <= 1.0) {
					route_count++;
					printf("%d %d %d\n", a, b, c);

					// Solve TSP
					int stops_to_visit[] = {a, b, c, DEPOT_INDEX};
					solveTSP(stops_to_visit, 4, tours, tsp_sol, tour_dists);
				}

				for(int d = c+1; d < NUM_STOPS; d++, comb_count++) {
					// Determine if this combo is feasible
					if((STOP_DEMAND[a] + STOP_DEMAND[b] + STOP_DEMAND[c] + STOP_DEMAND[d]) <= 1.0) {
						route_count++;
						printf("%d %d %d %d\n", a, b, c, d);

						// Solve TSP
						int stops_to_visit[] = {a, b, c, d, DEPOT_INDEX};
						solveTSP(stops_to_visit, 5, tours, tsp_sol, tour_dists);
					}

					for(int e = d+1; e < NUM_STOPS; e++, comb_count++) {
						// Determine if this combo is feasible
						if((STOP_DEMAND[a] + STOP_DEMAND[b] + STOP_DEMAND[c] + STOP_DEMAND[d]+ STOP_DEMAND[e]) <= 1.0) {
							// Route is feasible
							route_count++;
							printf("%d %d %d %d %d\n", a, b, c, d, e);

							// Solve TSP
							int stops_to_visit[] = {a, b, c, d, e, DEPOT_INDEX};
							solveTSP(stops_to_visit, 6, tours, tsp_sol, tour_dists);
						}
					}
				}
			}
		}
	}

	if(tours.size() == tsp_sol.size()) {
		// Print results to file
		FILE * pOutputFile;
		pOutputFile = fopen("enum_vrp.dat", "w");

		// Report the number of routes found
		fprintf(pOutputFile, "# Parameters\nparam N := %d;\n", route_count);
		// Set the number of stops for this problem
		fprintf(pOutputFile, "param M := %d;\n", NUM_STOPS);
		// Set the number of vehicles that can be used
		fprintf(pOutputFile, "param k := %d;\n", NUM_VHCLE);

		// Write out the routes found, and their TSP solutions
		fprintf(pOutputFile, "param r :\n");
		for(int i = 0; i < NUM_STOPS; i++) {
			fprintf(pOutputFile, "\t%d", i+1);
		}
		fprintf(pOutputFile, " :=\n");
		for(int t = 0; t < (int)tours.size(); t++) {
			fprintf(pOutputFile, "%d", t+1);
			for(int i = 0; i < NUM_STOPS; i++) {
				fprintf(pOutputFile, "\t%d", tours[t][i]);
			}
			fprintf(pOutputFile, "\n#");
			for(int i = 0; i < NUM_STOPS; i++) {
				if(tsp_sol[t][i] == -1) {
					break;
				}
				else {
					fprintf(pOutputFile, " %d", tsp_sol[t][i] + 1);
				}
			}
			fprintf(pOutputFile, "\n");
		}
		fprintf(pOutputFile, ";\n");

		// Write out costs of each tour
		fprintf(pOutputFile, "param d :=\n");
		for(int t = 0; t < (int)tour_dists.size(); t++) {
			fprintf(pOutputFile, "%d %f\n", t+1, tour_dists[t]);
		}
		fprintf(pOutputFile, ";\n");

		// Done writing, close file
		fclose(pOutputFile);
	}

	for(auto t : tours) {
		delete[] t;
	}

	for(auto arr : tsp_sol) {
		delete[] arr;
	}

	// Report statistics
	printf("\nTried %d combinations and found %d feasible tours\n", comb_count, route_count);

	return 0;
}
