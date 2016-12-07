// teste_algorithm_delaunay_voronoi_vs.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "..\..\source\voronoi_delaunayTriangulation.h"


int main()
{

	std::vector<std::vector<double>> matrixPoint {  {0.9501, 0.4057},
													{0.2311, 0.9355},
													{0.6068, 0.9169},
													{0.4860, 0.4103},
													{0.8913, 0.8936},
													{0.7621, 0.0579},
													{0.4565, 0.3529},
													{0.0185, 0.8132},
													{0.8214, 0.0099},
													{0.4447, 0.1389},
													{0.6154, 0.2028},
													{0.7919, 0.1987},
													{0.9218, 0.6038},
													{0.7382, 0.2722},
													{0.1763, 0.1988} };
	

	alg::DelaunayTriangulation dl = alg::generateDelaunayTriangulation(matrixPoint);

    return 0;
}

