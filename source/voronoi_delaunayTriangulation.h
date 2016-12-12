
#include "Multivector.h"
#include <set>
#include "MatlabEng.h"
#include <algorithm>
#include <iterator>

namespace alg {

	class DelaunayVertex {
	public:
		DelaunayVertex() {}
		DelaunayVertex(const Multivector<double> &pt) : m_pt(pt) {}
		DelaunayVertex(const DelaunayVertex &V) : m_pt(V.m_pt), m_simplex(V.m_simplex) {}

		DelaunayVertex &operator=(const DelaunayVertex &V) {
			if (&V != this) {
				m_pt = V.m_pt;
				m_simplex = V.m_simplex;
			}
			return *this;
		}
		
		Multivector<double> m_pt;
		std::set<int> m_simplex;
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

		bool isBorder()const {
			for (int i : m_neighborTriangles)
				if (i == -1)
					return true;
			return false;
		}

		int getNextNeighbor(const std::set<int>& _neighbor) const{
			for (int index : m_neighborTriangles)
				if (_neighbor.find(index) != _neighbor.end())
					return index;
			return -1; // não achou vizinho
		}

		std::vector<int> m_vertices;

		std::vector<int> m_neighborTriangles;

	};


	class DelaunayTriangulation {
	public:
		DelaunayTriangulation() {}

		DelaunayTriangulation(Orthogonal<int> _ort) : m_metrica(_ort){
		}

		friend DelaunayTriangulation generateDelaunayTriangulation2D_matlab(const std::vector<std::vector<double>>& matrix_point);


		inline Multivector<double> No() const {
			return (0.5*e(m_metrica.getSize() - 1) + e(m_metrica.getSize()));
		}

		inline Multivector<double> Nii() const{
			return (1.0*e(m_metrica.getSize()) - e(m_metrica.getSize() - 1));
		}

		Multivector<double> hiperesfera(int _indexSimplex) const{
			Multivector<double> C = m_vertices[m_simplex[_indexSimplex].m_vertices[0]].m_pt;
			for (size_t i = 1; i < m_simplex[_indexSimplex].m_vertices.size(); i++)
				C = C^m_vertices[m_simplex[_indexSimplex].m_vertices[i]].m_pt;
			
			return C;
		}

		int selectInicialSimplex(const std::set<int>& _neighbors) const{
			for (const auto& index : _neighbors)
				if (m_simplex[index].isBorder())
					return index;
			return *_neighbors.begin();
		}

	public:
		std::vector<DelaunayVertex> m_vertices;

		std::vector<DelaunaySimplex> m_simplex;

		Orthogonal<int> m_metrica;

	};

	class Voronoi {
	  public:
		Voronoi(Multivector<double> Nii) {
			m_map_id_simplex_pontoFinito[-1] = Nii;
		}

		std::map<int ,Multivector<double>> m_map_id_simplex_pontoFinito;

		std::vector<std::vector<int>> m_region;
	};


	void calculateDelaunayMatlab2D(DelaunayTriangulation& dl, const int dim) {
		
		//calcular delaunay pelo matlab

		CMatlabEng matlab;

		matlab.Open(NULL);
		matlab.SetVisible(FALSE);

		// matlab 
		mxArray *mxc;
		mxc = mxCreateDoubleMatrix(dim, dl.m_vertices.size(), mxREAL);
		double * ptrIn = (double *)mxGetPr(mxc);

		int i = 0;
		for (auto vertex : dl.m_vertices) {
			int j = 0;
			for (; j < dim; j++) {
				ptrIn[i*(dim)+j] = vertex.m_pt.getCoefBladeBase(j + 1);
			}
			i++;
		}
			
		matlab.PutVariable("T", mxc);
		matlab.EvalString("K = delaunayTriangulation(T');");
		matlab.EvalString("K = K.ConnectivityList';");
		matlab.EvalString("[dim1 dim2] = size(K);");

		mxArray* dim1 = matlab.GetVariable("dim1");
		double* a = (double *)mxGetData(dim1);
		mxArray* dim2 = matlab.GetVariable("dim2");
		double* b = (double *)mxGetData(dim2);

		std::cout << "Qtd simplex" << (int)b[0] << std::endl;
		std::cout << "vertices" << (int)a[0] << std::endl;

		mxArray* convex_hull = matlab.GetVariable("K");
		double * ptr = (double *)mxGetData(convex_hull);

		for (int i = 0; i < (int)b[0]; i++) { // para cada triângulo 
			dl.m_simplex.push_back(DelaunaySimplex());
			int indexTriang = dl.m_simplex.size() - 1;
			for (int j = 0; j < (int)a[0]; j++) {
				dl.m_simplex[indexTriang].m_vertices.push_back(ptr[i*(int)a[0] + j] - 1);
				dl.m_vertices[(ptr[i*(int)a[0] + j] - 1)].m_simplex.insert(indexTriang);
			
			}
		}
		
		mxDestroyArray(mxc);
		mxDestroyArray(convex_hull);
		mxDestroyArray(dim1);
		mxDestroyArray(dim2);

	}

	DelaunayTriangulation generateDelaunayTriangulation2D_matlab(const std::vector<std::vector<double>>& matrix_point) {
		
		if (matrix_point.empty())
			return DelaunayTriangulation();

		std::vector<int> metric_values;
		for (size_t i = 0; i < matrix_point[0].size(); i++)
			metric_values.push_back(1);
		metric_values.push_back(1);
		metric_values.push_back(-1);

		Orthogonal<int> metric(metric_values);
		DelaunayTriangulation dl(metric);

		const int qtdPoint = matrix_point.size();
		const int dim = matrix_point[0].size();

		for (size_t i = 0; i < matrix_point.size(); i++) {
			Multivector<double> P = dl.No();
			double valueNii = 0.0;
			for (size_t j = 0; j < matrix_point[i].size(); j++) {
				P = P + e(j + 1)*matrix_point[i][j];
				valueNii += ((matrix_point[i][j] * matrix_point[i][j])*0.5);
			}
			P = P + valueNii*dl.Nii();
			dl.m_vertices.push_back(DelaunayVertex(P));
		}

		calculateDelaunayMatlab2D(dl, dim);
		
		// calcular vizinhança
		for (size_t i = 0; i < dl.m_simplex.size(); i++) {
			DelaunaySimplex &T = dl.m_simplex[i];

			for (size_t j = 0; j <  T.m_vertices.size(); j++) {
				std::set<int> simplex_neighbor = dl.m_vertices[T.m_vertices[j]].m_simplex; 
				
				for (size_t aux = 1; aux < T.m_vertices.size() - 1; aux++) {
					const DelaunayVertex &V2 = dl.m_vertices[T.m_vertices[(j + aux) % T.m_vertices.size()]];  // (j + aux ) mod qtd vertices , necessário para pecorrer o vetor de forma circular. 
					std::vector<int> intersection;
					std::set_intersection(simplex_neighbor.begin(), simplex_neighbor.end(), V2.m_simplex.begin(), V2.m_simplex.end(), std::back_inserter(intersection));
					simplex_neighbor.clear();
					simplex_neighbor.insert(intersection.begin(), intersection.end());
				}
				simplex_neighbor.erase(i);

				if (!simplex_neighbor.empty())
					T.m_neighborTriangles.push_back(*(simplex_neighbor.begin()));
				else
					T.m_neighborTriangles.push_back(-1); // não possui vizinho
			}
		}
		return dl;
	}

	Voronoi generateVoronoi(const DelaunayTriangulation& _delaunay) {
		
		Voronoi  v(_delaunay.Nii()); // blade direção 

		//calcula o ponto finito, ou seja, o centro de todas as hiperesferas
		for (size_t i = 0; i < _delaunay.m_simplex.size(); i++){
			
				const DelaunaySimplex &T = _delaunay.m_simplex[i];
				Multivector<double> hiperEsfera = _delaunay.hiperesfera(i);
				Multivector<double> pontoFinito = -0.5*(GP(GP(hiperEsfera, _delaunay.Nii(), _delaunay.m_metrica), hiperEsfera, _delaunay.m_metrica));
				Multivector<double> divisor = LCont(_delaunay.Nii(), hiperEsfera, _delaunay.m_metrica);
				Multivector<double> divisorFinal = GP(divisor, divisor, _delaunay.m_metrica);

				v.m_map_id_simplex_pontoFinito[i] = pontoFinito * (1.0 / divisorFinal.getCoefBladeBase(0));
		}

		// para cada vertices referência montar o polígono do voronoi usando os pontos finitos já calculados
		for (size_t j = 0; j < _delaunay.m_vertices.size(); j++) {

			const DelaunayVertex& vertex = _delaunay.m_vertices[j];
			std::set<int> neighbor = vertex.m_simplex;

			std::vector<int> region;
			
			int next_index = _delaunay.selectInicialSimplex(neighbor);
			region.push_back(next_index);
			neighbor.erase(next_index);
			
			while (!neighbor.empty()) {
				next_index = _delaunay.m_simplex[next_index].getNextNeighbor(neighbor);
				region.push_back(next_index); 
				neighbor.erase(next_index);
			}

			v.m_region.push_back(region);
		}
		return v;
	}

}