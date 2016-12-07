
#include "Multivector.h"
#include <set>
#include "MatlabEng.h"

namespace alg {

	class DelaunayVertex {
	public:
		DelaunayVertex() {}
		DelaunayVertex(const Multivector<double> &pt) : m_pt(pt) {}
		DelaunayVertex(const DelaunayVertex &V) : m_pt(V.m_pt), m_triangles(V.m_triangles) {}

		DelaunayVertex &operator=(const DelaunayVertex &V) {
			if (&V != this) {
				m_pt = V.m_pt;
				m_triangles = V.m_triangles;
			}
			return *this;
		}


		Multivector<double> m_pt;
		std::set<int> m_triangles;
	};


	class DelaunaySimplex {
	public:
		DelaunaySimplex() {
		}

		DelaunaySimplex(std::vector<int> vtxIdx) {
			m_vertices = vtxIdx;
		}

		DelaunaySimplex(const DelaunaySimplex &T) {
			m_vertices = T.m_vertices;
			m_neighborTriangles = T.m_neighborTriangles;
		
		}

		DelaunaySimplex& operator=(const DelaunaySimplex &T) {
			if (this != &T) {
				m_vertices = T.m_vertices;
				
				m_neighborTriangles = T.m_neighborTriangles;
				
			}
			return *this;
		}

		std::vector<int> m_vertices;

		std::vector<int> m_neighborTriangles;

	};


	class DelaunayTriangulation {
	public:
		DelaunayTriangulation() {}
		DelaunayTriangulation(const DelaunayTriangulation &DT) :
			m_vertices(DT.m_vertices),
			m_triangles(DT.m_triangles) {}

		DelaunayTriangulation& operator=(const DelaunayTriangulation &DT) {
			if (this != &DT) {
				m_vertices = DT.m_vertices;
				m_triangles = DT.m_triangles;
			}
		}

		std::vector<DelaunayVertex> m_vertices;

		std::vector<DelaunaySimplex> m_triangles;
	};

	void calculateConvexHull(DelaunayTriangulation& dl, const int dim, const int No, const int Nii) {
		
		//calcula convex hull chama o matlab para realizar o cálculo

		CMatlabEng matlab;

		matlab.Open(NULL);
		matlab.SetVisible(FALSE);

		// assume c gets initialized, all of it 
		mxArray *mxc;
		mxc = mxCreateDoubleMatrix(dim + 1, dl.m_vertices.size(), mxREAL);
		double * ptrIn = (double *)mxGetPr(mxc);

		int i = 0;
		for (auto vertex : dl.m_vertices) {
			int j = 0;
			for (; j < dim; j++)
				ptrIn[i *(dim + 1) + j] = vertex.m_pt.getElementBase(j + 1);
			ptrIn[i *(dim + 1) + j] = vertex.m_pt.getElementBase(Nii);
			i++;
		}

		matlab.PutVariable("T", mxc);
		matlab.EvalString("K = convhulln(T');");
		matlab.EvalString("K = K';");
		matlab.EvalString("[dim1 dim2] = size(K);");

		mxArray* dim1 = matlab.GetVariable("dim1");
		double* a = (double *)mxGetData(dim1);
		mxArray* dim2 = matlab.GetVariable("dim2");
		double* b = (double *)mxGetData(dim2);

		mxArray* convex_hull = matlab.GetVariable("K");
		double * ptr = (double *)mxGetData(convex_hull);

		for (int i = 0; i < (int)b[0]; i++) { // para cada triângulo 
			dl.m_triangles.push_back(DelaunaySimplex());
			int indexTriang = dl.m_triangles.size() - 1;
			for (int j = 0; j < (int)a[0]; j++) {
				dl.m_triangles[indexTriang].m_vertices.push_back(ptr[i*(int)a[0] + j] - 1);
				dl.m_vertices[(ptr[i*(int)a[0] + j] - 1)].m_triangles.insert(indexTriang);
			}
		}
		
		mxDestroyArray(mxc);
		mxDestroyArray(convex_hull);
		mxDestroyArray(dim1);
		mxDestroyArray(dim2);

	}

	template <typename T>
	DelaunayTriangulation generateDelaunayTriangulation(const std::vector<std::vector<T>>& matrix_point) {
		
		DelaunayTriangulation dl;

		if (matrix_point.empty())
			return dl;

		const int qtdPoint = matrix_point.size();
		const int dim = matrix_point.begin()->size();

		const int No = matrix_point.begin()->size() + 1;
		const int Nii = matrix_point.begin()->size() + 2;

		for (size_t i = 0; i < matrix_point.size(); i++) {
			Multivector<double> P = (e(No)*1.0);
			for (size_t j = 0; j < matrix_point[i].size(); j++)
				P = P + e(j + 1)*matrix_point[i][j] + e(Nii)*((matrix_point[i][j] * matrix_point[i][j])*0.5);

			dl.m_vertices.push_back(DelaunayVertex(P));
		}

		calculateConvexHull(dl, dim, No, Nii);

		for (size_t i = 0; i < dl.m_triangles.size(); i++) {
			DelaunaySimplex &T = dl.m_triangles[i];

			for (size_t j = 0; j <  T.m_vertices.size(); j++) {

				const DelaunayVertex &V1 = dl.m_vertices[T.m_vertices[j]];
				const DelaunayVertex &V2 = dl.m_vertices[T.m_vertices[(j + 1) % 3]];

				// todo: use set union of STL
				for (std::set<int>::const_iterator I = V1.m_triangles.begin();
					I != V1.m_triangles.end(); I++) {
					if (*I == i) continue; // must not be this triangle!
					if (V2.m_triangles.find(*I) != V2.m_triangles.end()) {
						T.m_neighborTriangles.push_back(*I);
					}
				}
			}
		}


		return dl;
	}

	void generateVoronoi(const DelaunayTriangulation& _delaunay) {

	}

}