#ifndef METRIC_H
#define METRIC_H

#include <vector>

template<typename T>
class Metric{

	virtual T eval(int i, int j) const = 0;

};

template<typename T>
class Orthogonal : public Metric<T> {

public:

	Orthogonal(){
	}

    Orthogonal(std::vector<T> a ){
		m = a;
    }

    T eval(int i, int) const {
        return m[i];
    }

    T eval(int i) const {
        return eval(i,i);
    }

	size_t getSize() const{
		return m.size();
	}

protected:
	void set_value(T value){
		m.push_back(value);
	}

private:
    std::vector<T> m;
};

class Orthonormal : public Orthogonal<int> {

public:
	Orthonormal(int dimension){
		for (int i=0; i < dimension; i++)
			set_value(1);
	}
};

#endif 
