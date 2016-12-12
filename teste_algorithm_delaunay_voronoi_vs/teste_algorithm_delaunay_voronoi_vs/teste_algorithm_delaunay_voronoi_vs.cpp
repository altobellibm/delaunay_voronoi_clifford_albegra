// teste_algorithm_delaunay_voronoi_vs.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "..\..\source\voronoi_delaunayTriangulation.h"


int main()
{

	//Multivector<double > n = 1.0*e(1) + 1.0*e(2) + 1.0*e(3) + +1.0*e(4);

	std::vector<std::vector<double>> matrixPoint { { -1.5,3.2 },
													{ 1.8,3.3 },
													{ -3.7,1.5 },
													{ -1.5,1.3 },
													{ 0.8,1.2 },
													{ 3.3,1.5 },
													{ -4.0,-1.0},
													{ -2.3,-0.7 },
													{ 0,-0.5 },
													{ 2.0,-1.5 },
													{ 3.7,-0.8 },
													{ -3.5,-2.9 },
													{ -0.9,-3.9 },
													{ 2.0,-3.5 },
													{ 3.5,-2.25 } };
	

	alg::DelaunayTriangulation dl = alg::generateDelaunayTriangulation2D_matlab(matrixPoint);

	alg::Voronoi vr = alg::generateVoronoi(dl);

	for (int i = 0; i < vr.m_region.size(); i++)
		for (int j = 0; j < vr.m_region[i].size(); j++) {
			std::cout << vr.m_map_id_simplex_pontoFinito[vr.m_region[i][j]].getCoefBladeBase(1) << ", ";
			std::cout << vr.m_map_id_simplex_pontoFinito[vr.m_region[i][j]].getCoefBladeBase(2) << std::endl;
		}

    return 0;
}

