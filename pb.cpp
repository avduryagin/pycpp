#include "pybind11\pybind11.h"
#include "constructor_stats.h"
#include "pybind11\numpy.h"
#include "pybind11_tests.h"
#include "Python.h"
#include <iostream>
#include <stdarg.h>

namespace py = pybind11;


int add(int i=0, int j=0)
{
	return(i + j);
}

struct Pet
{
	Pet(const std::string &name,const int age):name(name),age(age){}
	Pet() :name("Sobol"), age(5) {}
	void setName(std:: string &name)
	{
		this->name = name;
	}
	const std::string& getName()
	{
		return name;
	}
	void setAge(int& age) 
	{ 	
		this->age=age; 	

	}
	const int& getAge() { return age; }
	void increaseAge(int& a)
	{
		age = add(age, a);
	}

	std::string name;
	int age;
	
};

void print_dict(const py::dict& dict)
{
	for(auto k:dict)
	{
		std::cout << "key=" << std::string(py::str(k.first)) << "; value=" << std::string(py::str(k.second)) << "\n";
	}
}



//template <typename T, int ExtraFlags>
//class expanded_array:public py::array_t<T,ExtraFlags>{};
/*class expanded_array :public py::array_t<T, ExtraFlags>
{

public:
	//expanded_array():py::array_t<T, ExtraFlags>() {}
	ssize_t shape(ssize_t dim) const
	{
		return py::array_t<T, ExtraFlags>::shape(dim) + 1;
	}

	template <typename... Ix>
	const T& at(Ix ... args) const {
		//ssize_t i = i_ - 1;
		//ssize_t j = j_ - 1;
		//if (i < 0 | j < 0) return zero;		
		return py::array_t<T,ExtraFlags>::at(args...);



	}

};*/
template <typename T, int ExtraFlags>
class expanded_array
{
private:
	py::array_t<T, ExtraFlags> array;	
	
public:
	expanded_array(py::array_t<T, ExtraFlags>& array_) { array = array_; }
	py::ssize_t ndim() const 
	{ 
		return array.ndim();
	
	}
	py::ssize_t shape(int i) const
	{
		return array.shape(i) + 1;
	}
	const T& at(py::ssize_t i_,py::ssize_t j_) const
	{
		py::ssize_t i, j;
		i = i_-1;
		j = j_-1;
		if (i < 0 || j < 0) { return 0; }
		return array.at(i, j);
	}



};
void print_array(py::array_t<double, py::array::c_style | py::array::forcecast> array_)
{
	expanded_array<double, py::array::c_style | py::array::forcecast> array= expanded_array<double, py::array::c_style | py::array::forcecast>(array_);


	py::ssize_t ndim = array.ndim();
	if(ndim==2)
	{
		py::ssize_t n = array.shape(0);
		py::ssize_t m = array.shape(1);
		

		for(int i=0;i<n;i++)
		{
			for (int j = 0; j < m; j++)
			{
				std::cout << array.at(i,j)<<" ";

		    }
			std::cout << std::endl;
		}
	}

}

template <typename T>
class top_level
{
	T level;

public:
	top_level() { level = (T)0; }
	top_level(T val) { level = val; }
	T get_level() const
	{
		return level;
	}

	~top_level() { }


};
template <typename T>
class sub_level:top_level<T>
{public:
	sub_level():top_level<T>(){}
	sub_level(T val):top_level<T>(val){}
	T get_level()
	{
		return top_level<T>::get_level() + 1;
	}

};

template <typename T>
void cfill(T* x, int n, T val)
{
	for(int i=0;i<n;i++)
	{
		x[i] = val;
	}
}
double assignment(py::array_t<double, py::array::c_style | py::array::forcecast> array_){

	expanded_array<double, py::array::c_style | py::array::forcecast> a = expanded_array<double, py::array::c_style | py::array::forcecast>(array_);
	int k = (int)a.shape(0);
	int n = (int)a.shape(1);

	double* u = (double*)calloc(n + k + n, sizeof(double));
	double* v = u + k;
	double* minv = v + n;
	int* p = (int*)calloc(n + n + k, sizeof(int));
	int* way = p + n;
	int* c = way + n;
	bool* used = (bool*)calloc(n, sizeof(bool));
	double inf = std::numeric_limits<double>::infinity();

	for (int i = 1; i < k; i++)
	{
		p[0] = i;
		int j0 = 0;

		cfill<double>(minv, n, inf);
		cfill<bool>(used, n, false);
		do
		{
			used[j0] = true;
			int i0 = p[j0];
			double delta = inf;
			int j1 = 0;
			for (int j = 1; j < n; j++)
			{
				if (!used[j])
				{
					double cur = a.at(i0,j) - u[i0] - v[j];
					//double cur = a[i0][j] - u[i0] - v[j];
					if (cur < minv[j])
					{
						minv[j] = cur;
						way[j] = j0;
					}
					if (minv[j] < delta)
					{
						delta = minv[j];
						j1 = j;
					}
				}

			}
			for (int j = 0; j < n; j++)
			{
				if (used[j])
				{
					u[p[j]] = u[p[j]] + delta;
					v[j] = v[j] - delta;
				}
				else
				{
					minv[j] = minv[j] - delta;
				}
			}
			j0 = j1;

		} while (p[j0] != 0);
		do
		{
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;

		} while (j0 > 0);


	}
	for (int i = 0; i < n; i++)
	{
		c[p[i] - 1] = i - 1;
		//std::cout << c[i] << "  ";

	}


	/*for (int i = 0; i < (k - 1); i++)
	{
		//c[p[i] - 1] = i - 1;
		std::cout <<i<<"->"<< c[i] <<'\n';

	}*/
	//std::cout << "sum=" << v[0];
	double sum = v[0];

	//*sum = v[0];
	delete[] u;
	//delete p;
	delete[] used;
	return sum;
	//if (expand) { delete[] a; }
	//return c;
}

PYBIND11_MODULE(examples,m)
{
	m.doc()="my firs example";
	m.def("sum", &add, "returns sum of two integers",py::arg("i")=0, py::arg("j")=0);
	m.def("print_dict", &print_dict);
	m.def("print_array", &print_array);
	m.def("assignment", &assignment);
	
	
	py::class_<Pet>(m, "Pet")
		.def("__repr__", [](const Pet& a) {return "struct Pet; name=" + a.name + "; age= " + std::to_string(a.age); })
		.def_readonly("age",&Pet::age)
		.def_readonly("name", &Pet::name)
		.def(py::init<const std::string &,const int &>())
		.def(py::init())
		.def("setname", &Pet::setName)
		.def("getname", &Pet::getName)
		.def("setage", &Pet::setAge)
		.def("getage", &Pet::getAge)
		.def("increase", &Pet::increaseAge);

	py::class_<top_level<int>>(m, "Top")
		.def(py::init())
		.def("__repr__", [](const top_level<int>& t) { return "class level=" + std::to_string(t.get_level()); })
		.def("level", &top_level<int>::get_level);

	py::class_<sub_level<int>>(m, "Next")
		.def(py::init())
		.def("__repr__", [](sub_level<int>& t) { return "class level=" + std::to_string(t.get_level()); })
		.def("level", &sub_level<int>::get_level);




}