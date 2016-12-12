
#include <map>
#include <type_traits>
#include <bitset>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <utility>

#include "Utils.h"
#include "Metric.h"


//Biblioteca de álgebra geométrica 
//Toda implementação foi baseado no livro : <Esperando a publicação do livro do Leandro et.al para fornecer o link do livro>


enum class TYPE { VERSOR, BLADE, NO_STRUCTURE_OF_INTEREST};

template <typename T>
class Multivector{
public:
	Multivector(){}
	
    template<typename T1>
    friend Multivector<T1> e(unsigned int i);
	
    friend unsigned int take_grade(unsigned int masc); 

    template<typename T1>
    friend int take_grade(const Multivector<T1>& A); 

    template<typename T1, typename T2 >
    friend TYPE GetType(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension);

    template<typename T1>
    friend Multivector<T1> grade_extraction(const Multivector<T1>& M, unsigned int grade);

    template<typename T1>
    friend std::ostream& operator<<(std::ostream& out, const Multivector<T1>& A);

    friend Multivector<int> Pseudoscalar(unsigned int _dimension); 

    template<typename T1>
    friend bool IsZero(const Multivector<T1>& A);

	
	//equação 2.6.60
    template<typename T1>
    friend Multivector<T1> Involution(const Multivector<T1>& A);

	//equação 2.4.16
    template<typename T1>
    friend Multivector<T1> Reverse(const Multivector<T1>& A);

    // blades

	//algoritmo 5.1 do livro
	template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> operator+(const Multivector<T1>& A, const Multivector<T2>& B); //ok

	template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> operator-(const Multivector<T1>& A, const Multivector<T2>& B); //ok

	template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivector<T2>& A); // ok

	template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> operator*(const Multivector<T1>& A, const T2& l); // ok

	//algoritmo 5.1 do livro
	template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> operator^(const Multivector<T1>& A, const Multivector<T2>& B); // ok

    template <typename T1, typename T2>
    friend bool operator==(const Multivector<T1>& A, const Multivector<T2>& B);

	//algoritmo 5.5
    template <typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> RP(const Multivector<T1>& A, const Multivector<T2>& B, int dimention); //ok

    //geometric operations

	//algoritmo 5.6
    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort); // ok

	//subcaso do algoritmo 5.4
	template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> LCont(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort); //ok

	//subcaso do algoritmo 5.4
    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> RCont(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort); //ok
    
    //subcaso do algoritmo 5.4
    template<typename T1, typename T2, typename T3>
    friend typename std::common_type<T1, T2, T3>::type SCP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort); //ok


    template<typename T1, typename T2>
    friend typename std::common_type<T1, T2>::type SQR_Norm_Reverse(const Multivector<T1>& A, const Orthogonal<T2>& ort);

	//equação 2.4.21
    template<typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> INV(const Multivector<T1>& A, const Orthogonal<T2>& ort);


    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> IGP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

    //versores
	//algoritmo 5.8
	template<typename T1, typename T2, typename T3>
	friend Multivector<typename std::common_type<T1, T2, T3>::type> DeltaP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

	//definido logo abixo do algoritmo 5.9
    template<typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> NormalizeBlade(const Multivector<T1>& A, const Orthogonal<T2>& ort);

	//algoritmo 2.4.20
    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type> OrtoProjection(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

	//algoritmo 2.4.23
	template<typename T1, typename T2, typename T3>
	friend Multivector<typename std::common_type<T1, T2, T3>::type> OrtoRejection(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort);

	//algoritmo 5.9
	template<typename T1, typename T2>
	friend std::pair< typename std::common_type<T1, T2>::type, std::vector< Multivector<typename std::common_type<T1, T2>::type> > > factorization(const Multivector<T1>& A, const Orthogonal<T2>& ort);

	//algoritmo 5.10
    template<typename T1, typename T2, typename T3>
    friend std::pair<typename std::common_type<T1, T2, T3>::type, typename std::common_type<T1, T2, T3>::type > MeetJoinBlade(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& _orto, unsigned int dimension);

	//equação 2.7.63
    template<typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> Dual(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension);

	//equação 2.7.66
    template<typename T1, typename T2>
    friend Multivector<typename std::common_type<T1, T2>::type> UnDual(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension);

	T getCoefBladeBase(int index) {
		auto it = m.find(index);
		if (it != m.end())
			return it->second;

		return 0;
	}

private:
    template<typename T1, typename T2, typename T3>
    friend Multivector<typename std::common_type<T1, T2, T3>::type>  GP(const Orthogonal<T3>& ort, unsigned int masc_1, T1 cof_1, unsigned int masc_2, T2 cof_2);

	//algoritmo 5.4
	friend int canonical_order(unsigned int masc1, unsigned int masc2);
	
private:
    typedef std::map<unsigned int, T> MAP;
	MAP m;
};

template<typename T1 = int>
Multivector<T1> e(unsigned int i){

    Multivector<T1> a;
	if (i != 0)
		a.m[1 << (i - 1)] = 1;
	else
		a.m[i] = 1;

	return a;
}

template<typename T1>
std::ostream& operator<<(std::ostream& out, const Multivector<T1>& A){

    std::string strMultivector;

    for(auto it = A.m.begin(); it != A.m.end(); it++){

        it->second > 0 ? strMultivector+="+" : strMultivector;

        strMultivector+= std::to_string(it->second) + "*";

        unsigned int masc = it->first;

        if(masc == 0)
            strMultivector+= "e(0)";
        else
        {
            int indice = 1;
            unsigned int mascr = 1;
            strMultivector+= "( ";

            while (masc != 0){

                if ( (masc & mascr) != 0) 
                      strMultivector+= "e(" + std::to_string(indice) + ") ^ ";
                
                indice++;
                masc = masc >> 1;
            }

			strMultivector = strMultivector.substr(0, strMultivector.size() - 2);
            strMultivector+= " ) ";
        }
    }

    strMultivector.size() == 0 ? strMultivector : strMultivector.substr(0, strMultivector.size()-1);

	if (A.m.size() == 0)
		strMultivector = std::to_string(0);

    out <<  strMultivector;
    return out;
}

int canonical_order(unsigned int masc1, unsigned int masc2){

    unsigned int trocas = 0;
    masc1 = masc1 >> 1;
    while (masc1 != 0){
        trocas += Utils::HammingWeight(masc1 & masc2);
        masc1 = masc1 >> 1;
    }

    if ((trocas & 1) == 0)
        return +1;
    else
        return -1;
}

template<typename T>
T metric_factor(unsigned int masc, const Orthogonal<T>& orth){
	
	int indice = 0;
	T product_value = 1;
    unsigned int mascr = 1;

	do {

		if ((masc & mascr) != 0)
			product_value *= orth.eval(indice);

		masc = masc >> 1;
		indice++;

	} while (masc != 0);
		
	return product_value;
}

unsigned int take_grade(unsigned int masc){
    return Utils::HammingWeight(masc);
}

template<typename T1>
int take_grade(const Multivector<T1>& A){

    if(A.m.empty())
        return -1;

    auto it = A.m.begin();
    unsigned int grade = take_grade(it->first);
    it++;

    for (; it != A.m.end(); it++)
        if (take_grade(it->first) != grade)
            return -1;

    return grade;
}

template<typename T1, typename T2>
TYPE GetType(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension){

    if(GP(Involution(A),INV(A,ort),ort) == GP(INV(A,ort),Involution(A),ort) ){

        for(unsigned int i = 1; i <= dimension; i++){
            Multivector<typename std::common_type<T1, T2>::type> R = PG(PG(Involution(A),e(1),ort), Reverse(A),ort);
            if( (R.m.find(0) != R.m.end()))
                return TYPE::NO_STRUCTURE_OF_INTEREST;

            if(R.m[0] != (typename std::common_type<T1, T2>::type) 1 )
                return TYPE::NO_STRUCTURE_OF_INTEREST;
            }

        if(take_grade(A) == -1)
            return TYPE::VERSOR;
        else
            return TYPE::BLADE;
    }
    else
        return TYPE::NO_STRUCTURE_OF_INTEREST;
}

template<typename T1>
Multivector<T1> grade_extraction(const Multivector<T1>& mul, unsigned int grade){
    Multivector<T1> R;

	for (auto it = mul.m.begin(); it != mul.m.end(); it++) 
        if (take_grade(it->first) == grade)
			R.m[it->first] = it->second;

    return R;
}


Multivector<int> Pseudoscalar(unsigned int _dimension){
    Multivector<int> e;
    e.m[(~(unsigned int)0) >> (32 - _dimension) ] = 1;
    return e;
}

template<typename T1>
bool IsZero(const Multivector<T1>& A){
    return A.m.size() == 0;
}


 template <typename T1, typename T2>
 Multivector<typename std::common_type<T1, T2>::type> operator+(const Multivector<T1>& A, const Multivector<T2>& B){
	 Multivector<typename std::common_type<T1, T2>::type> C;

	 auto itA = A.m.begin();
	 auto itB = B.m.begin();

	 while (itA != A.m.end() && itB != B.m.end()){
		 if (itA->first < itB->first){
			 C.m[itA->first] = itA->second;
			 itA++;
		 }
		 else
			 if (itA->first > itB->first){
			 C.m[itB->first] = itB->second;
			 itB++;
			 }
			 else{
				 C.m[itB->first] = itA->second + itB->second;
				 itA++;
				 itB++;
			 }
	 }

	 while (itA != A.m.end()){
		 C.m[itA->first] = itA->second;
		 itA++;
	 }

	 while (itB != B.m.end()){
		 C.m[itB->first] = itB->second;
		 itB++;
	 }

	 return C;
 }

 template <typename T1, typename T2>
 Multivector<typename std::common_type<T1, T2>::type> operator-(const Multivector<T1>& A, const Multivector<T2>& B){
	 return A + uminus(B);
 }

template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator*(const T1& l, const Multivector<T2>& A){
	
	Multivector<typename std::common_type<T1, T2>::type> C;

    auto itA = A.m.begin();
	while (itA != A.m.end()){
        C.m[itA->first] = itA->second *l;
		itA++;
	}

	return C;
 }
 
template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator*(const Multivector<T1>& A, const T2& l){

    return l*A;
}

template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> operator^(const Multivector<T1>& A, const Multivector<T2>& B){

	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
        for (auto itB = B.m.begin(); itB != B.m.end(); itB++)
            if ((itA->first & itB->first) == 0){

				Multivector<typename std::common_type<T1, T2>::type> D;
				D.m[(itA->first | itB->first)] = canonical_order(itA->first, itB->first)*itA->second*itB->second;
				C = C + D;
			}

	return C;
}

template <typename T1, typename T2>
bool operator==(const Multivector<T1>& A, const Multivector<T2>& B){

    if ( A.m.size() == B.m.size() ){
        for (auto itA = A.m.begin(); itA != A.m.end(); itA++){
            auto itB = B.m.find(itA->first);
            if ((itB == B.m.end()) && ((typename std::common_type<T1, T2>::type) itA->second != (typename std::common_type<T1, T2>::type) itB->second))
                return false;
       }
        return true;
    }
    return false;
}


template <typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> RP(const Multivector<T1>& A, const Multivector<T2>& B, int dimention){

	Multivector<typename std::common_type<T1, T2>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
				Multivector<typename std::common_type<T1, T2>::type> D;

                unsigned int mascr = itA->first & itB->first;
				if ((take_grade(itA->first) + take_grade(itB->first) - take_grade(mascr)) == dimention){					
					D.m[mascr] = canonical_order(itA->first^mascr, itB->first^mascr)*itA->second*itB->second;
				}
				C = C + D;
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> GP(const Orthogonal<T3>& ort, unsigned int masc_1, T1 cof_1, unsigned int masc_2, T2 cof_2){

   Multivector<typename std::common_type<T1, T2, T3>::type> R;
   R.m[masc_1^masc_2] = canonical_order(masc_1, masc_2)*metric_factor(masc_1 & masc_2, ort)*cof_1*cof_2;
   return R;

}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> GP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){

     Multivector<typename std::common_type<T1, T2, T3>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++){
            const Multivector<typename std::common_type<T1, T2, T3>::type>& D = GP(ort,itA->first,itA->second,itB->first,itB->second);
            C.m[D.m.begin()->first] += D.m.begin()->second;
		}

	return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> LCont(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
	Multivector<typename std::common_type<T1, T2, T3>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++) {
			const Multivector<typename std::common_type<T1, T2, T3>::type>& D = GP(ort, itA->first, itA->second, itB->first, itB->second);
			if (take_grade(D.m.begin()->first) == (take_grade(itB->first) - take_grade(itA->first)) )
				C.m[D.m.begin()->first] += D.m.begin()->second;

		}

    return C;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> RCont(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
	Multivector<typename std::common_type<T1, T2, T3>::type> C;

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++) {
			const Multivector<typename std::common_type<T1, T2, T3>::type>& D = GP(ort, itA->first, itA->second, itB->first, itB->second);
			if( take_grade(D.m.begin()->first)  == (take_grade(itA->first) - take_grade(itB->first)) )
				C.m[D.m.begin()->first] += D.m.begin()->second;
		}

    return C;
}

template<typename T1, typename T2, typename T3>
typename std::common_type<T1, T2, T3>::type SCP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
 
    typename std::common_type<T1, T2, T3>::type coef(0);

	for (auto itA = A.m.begin(); itA != A.m.end(); itA++)
		for (auto itB = B.m.begin(); itB != B.m.end(); itB++) {
			const Multivector<typename std::common_type<T1, T2, T3>::type>& D = GP(ort, itA->first, itA->second, itB->first, itB->second);
			if (take_grade(D.m.begin()->first) ==  0)
                coef += D.m.begin()->second;
		}

    return coef;
}


template<typename T1>
Multivector<T1> Reverse(const Multivector<T1>& A){

	Multivector<T> multivector_r;
	for (auto it = A.m.begin(); it != A.m.end(); ++it) {
		int r = take_grade(it->first);
		multivector_r.M[it->first] = ((-1) ^ ((r*(r - 1)) >> 1)) * it->second;
	}
	return multivector_r;

}

template<typename T1>
Multivector<T1> Involution(const Multivector<T1>& A){
   
	Multivector<T> invol = A;
	for (auto it = invol.M.begin(); it != invol.M.end(); ++it) {
		invol.M[it->first] = ((-1) ^ (take_grade(A.m.begin()->first))) * it->second;
	}
	return invol;

}

template<typename T1, typename T2>
typename std::common_type<T1, T2>::type SQR_Norm_Reverse(const Multivector<T1>& A, const Orthogonal<T2>& ort){
    return SCP(A, Reverse(A), ort);
}

template<typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> INV(const Multivector<T1>& A, const Orthogonal<T2>& ort){
  
    return Reverse(A)*(1.0/(SQR_Norm_Reverse(A,ort)));
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> IGP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
  
    Multivector< typename std::common_type<T2,T3>::type > B_I;
    B_I.m[0] = INV(B, ort);
    return GP(A,B_I, ort);

}

template <typename T>
Multivector<T> uminus(const Multivector<T>& A){
    return -1*A;
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> DeltaP(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort) {

    Multivector<typename std::common_type<T1, T2, T3>::type> C = GP(A,B,ort);

    // get greater grade
    unsigned int greater_grade = 0;
    for(auto it = C.m.begin(); it != C.m.end(); it++){
        auto grade = take_grade(it->first);
        if(greater_grade < grade)
            greater_grade = grade;
    }

    return grade_extraction(C,greater_grade);
}

template<typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> NormalizeBlade(const Multivector<T1>& A, const Orthogonal<T2>& ort){
    return A*(1.0/sqrt(SQR_Norm_Reverse(A,ort)));
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> OrtoProjection(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort){
    return LCont(LCont(A,INV(B,ort),ort),B, ort);
}

template<typename T1, typename T2, typename T3>
Multivector<typename std::common_type<T1, T2, T3>::type> OrtoRejection(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort) {
	return LCont(LCont(A, B, ort), INV(B, ort), ort);
}

template<typename T1, typename T2>
std::pair< typename std::common_type<T1, T2>::type, std::vector< Multivector<typename std::common_type<T1, T2>::type> > > factorization(const Multivector<T1>& A, const Orthogonal<T2>& ort){

    std::vector< Multivector<typename std::common_type<T1, T2>::type> > factor_k;

    unsigned int masc;
    T1 greater_coef(0);

    // get masc from blade how has greater coeffcient.
    for(auto it = A.m.begin(); it != A.m.end(); it++)
        if(greater_coef < it->second){
            greater_coef = it->second;
            masc = it->first;
        }

    //T1& escalar, std::vector<Multivector<T1>>& _factor
    auto& scale = SQR_Norm_Reverse(A,ort);

    Multivector<T1> temp = NormalizeBlade;

    unsigned int grade = take_grade(masc);

   unsigned int mascr = 1;
   unsigned int indice = 0;

   for(unsigned int indice = 0; indice < (grade -1); indice++){

       if ((masc & mascr) != 0){
            Multivector<T1> ej;
            ej.m[1 << indice] = (T1)1.0;
            auto proj = OrtoProject(ej,temp, ort);

            factor_k.push_back(NormalizeBlade(proj,ort));
            temp = LCont(INV(factor_k.back(),ort),temp,ort);
        }
        masc = masc >> 1;
    }

    factor_k.push_back(NormalizeBlade(temp,ort));

    return std::make_pair(scale,factor_k);
}

template<typename T1, typename T2, typename T3>
std::pair<typename std::common_type<T1, T2, T3>::type, typename std::common_type<T1, T2, T3>::type > MeetJoinBlade(const Multivector<T1>& A, const Multivector<T2>& B, const Orthogonal<T3>& ort, unsigned int dimension){

    Multivector<typename std::common_type<T1, T2, T3>::type> Ar = A;
    Multivector<typename std::common_type<T1, T2, T3>::type> Br = B;

    unsigned int r = take_grade(Ar);
    unsigned int s = take_grade(Br);

    if (r > s){
        auto temp = Ar;
        Ar = Br;
        Br = temp;
    }

    auto delta = DeltaP(Ar,Br,ort);
    unsigned int t = 0.5*(r + s - take_grade(delta));

    auto& pair_scale_factor_k = factorization(Dual(delta,dimension),ort);

    auto INTER = 1*e(0);
    auto MEET = INV(Pseudoscalar(dimension),ort);

    const auto factors = pair_scale_factor_k.second;
    for(const auto factor_j : factors){

        auto PROJ = OrtoProjection(factor_j,Ar,ort);
        if(!IsZero(PROJ)){
           INTER = INTER ^ PROJ;
           //if(take_grade(*(INTER.m.begin())) == t){
               MEET = RCont(Ar,INV(INTER,ort),ort)^Br;
               continue;
           //}
        }

        auto& REJE = OrtoRejection(factor_j,Ar);
        if(!IsZero(REJE)){
            MEET = LCont(REJE,MEET,ort);
            //if take_grade(*(MEET.m.begin())) == (r + s - t){
                INTER = UnDual((Dual(Br,ort,dimension)^Dual(Ar,ort,dimension)), ort, dimension );
                continue;
           // }
        }
    }

    if (r > s){
        MEET = std::pow(-1,(r-t)*(s-t))*MEET;
        INTER = std::pow(-1,(r-t)*(s-t))*INTER;
    }

    return std::make_pair(MEET,INTER);
}

template<typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> Dual(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension){
    return LCont(A,INV(Pseudoscalar(dimension),ort),ort);
}

template<typename T1, typename T2>
Multivector<typename std::common_type<T1, T2>::type> UnDual(const Multivector<T1>& A, const Orthogonal<T2>& ort, unsigned int dimension){
    return LCont(A,Pseudoscalar(dimension),ort);
}

