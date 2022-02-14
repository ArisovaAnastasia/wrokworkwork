#pragma once
#include <vector>

template <typename type_data>
class data3d {
	//std::vector<type_data> data;
	type_data* data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	/*void Create(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		data.resize((n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2));

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}*/
	data3d(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		data = new type_data[(n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2)];

	this->n = n; this->m = m; this->k = k;
	this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
}

	inline type_data& operator()(int i, int j, int l) {
		/*
		if (i<0 || j<0 || l<0 || i>=(n + 2 * delta_x + 2) || j>= (m + 2 * delta_y + 2) || l>= (k + 2 * delta_z + 2)){
			std::cout << "ERROR" << std::endl;
		}*/
		return data[i * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2) + j * (k + 2 * delta_z + 2) + l];
	}

	int Get_Nx() {
		return n;
	}
	int Get_Ny() {
		return m;
	}
	int Get_Nz() {
		return k;
	}
	int Get_deltaX() {
		return delta_x;
	}
	int Get_deltaY() {
		return delta_y;
	}
	int Get_deltaZ() {
		return delta_z;
	}
	~data3d() {
		delete[] data;
	}
};

template <typename type_data>
struct Fields {
	data3d<type_data> Ex, Ey, Ez, Bx, By, Bz;
	/*Fields(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Ex.Create(n, m, k, delta_x, delta_y, delta_z);
		Ey.Create(n, m, k, delta_x, delta_y, delta_z);
		Ez.Create(n, m, k, delta_x, delta_y, delta_z);

		Bx.Create(n, m, k, delta_x, delta_y, delta_z);
		By.Create(n, m, k, delta_x, delta_y, delta_z);
		Bz.Create(n, m, k, delta_x, delta_y, delta_z);
	}*/

	Fields(int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Ex(n, m, k, delta_x, delta_y, delta_z),
		Ey(n, m, k, delta_x, delta_y, delta_z),
		Ez(n, m, k, delta_x, delta_y, delta_z),
		Bx(n, m, k, delta_x, delta_y, delta_z),
		By(n, m, k, delta_x, delta_y, delta_z),
		Bz(n, m, k, delta_x, delta_y, delta_z) {}
};

template <typename type_data>
struct SplitFields {
	data3d<type_data> Exy, Exz, Eyx, Eyz, Ezx, Ezy,
						Bxy, Bxz, Byx, Byz, Bzx, Bzy;
	/*SplitFields(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Exy.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy.Create(n, m, k, delta_x, delta_y, delta_z);
	}*/
	SplitFields(int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Exy(n, m, k, delta_x, delta_y, delta_z),
		Exz(n, m, k, delta_x, delta_y, delta_z),
		Eyx(n, m, k, delta_x, delta_y, delta_z),
		Eyz(n, m, k, delta_x, delta_y, delta_z),
		Ezx(n, m, k, delta_x, delta_y, delta_z),
		Ezy(n, m, k, delta_x, delta_y, delta_z),

		Bxy(n, m, k, delta_x, delta_y, delta_z),
		Bxz(n, m, k, delta_x, delta_y, delta_z),
		Byx(n, m, k, delta_x, delta_y, delta_z),
		Byz(n, m, k, delta_x, delta_y, delta_z),
		Bzx(n, m, k, delta_x, delta_y, delta_z),
		Bzy(n, m, k, delta_x, delta_y, delta_z) {}
};

template <typename type_data>
struct Sigma {
	data3d<type_data> sigmaEx, sigmaEy, sigmaEz, sigmaBx, sigmaBy, sigmaBz;
	/*Sigma(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		sigmaEx.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaEy.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaEz.Create(n, m, k, delta_x, delta_y, delta_z);

		sigmaBx.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaBy.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaBz.Create(n, m, k, delta_x, delta_y, delta_z);
	}*/

	Sigma(int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		sigmaEx(n, m, k, delta_x, delta_y, delta_z),
		sigmaEy(n, m, k, delta_x, delta_y, delta_z),
		sigmaEz(n, m, k, delta_x, delta_y, delta_z),
		sigmaBx(n, m, k, delta_x, delta_y, delta_z),
		sigmaBy(n, m, k, delta_x, delta_y, delta_z),
		sigmaBz(n, m, k, delta_x, delta_y, delta_z) {}
};
template <typename type_data>
struct Coefficient {
	data3d<type_data> Exy1, Exz1, Eyx1, Eyz1, Ezx1, Ezy1, Bxy1, Bxz1, Byx1, Byz1, Bzx1, Bzy1,
		Exy2, Exz2, Eyx2, Eyz2, Ezx2, Ezy2, Bxy2, Bxz2, Byx2, Byz2, Bzx2, Bzy2;
	/*Coefficient(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Exy1.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy1.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy1.Create(n, m, k, delta_x, delta_y, delta_z);

		Exy2.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy2.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy2.Create(n, m, k, delta_x, delta_y, delta_z);
	}*/

	Coefficient(int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Exy1(n, m, k, delta_x, delta_y, delta_z),
		Exz1(n, m, k, delta_x, delta_y, delta_z),
		Eyx1(n, m, k, delta_x, delta_y, delta_z),
		Eyz1(n, m, k, delta_x, delta_y, delta_z),
		Ezx1(n, m, k, delta_x, delta_y, delta_z),
		Ezy1(n, m, k, delta_x, delta_y, delta_z),

		Bxy1(n, m, k, delta_x, delta_y, delta_z),
		Bxz1(n, m, k, delta_x, delta_y, delta_z),
		Byx1(n, m, k, delta_x, delta_y, delta_z),
		Byz1(n, m, k, delta_x, delta_y, delta_z),
		Bzx1(n, m, k, delta_x, delta_y, delta_z),
		Bzy1(n, m, k, delta_x, delta_y, delta_z),

		Exy2(n, m, k, delta_x, delta_y, delta_z),
		Exz2(n, m, k, delta_x, delta_y, delta_z),
		Eyx2(n, m, k, delta_x, delta_y, delta_z),
		Eyz2(n, m, k, delta_x, delta_y, delta_z),
		Ezx2(n, m, k, delta_x, delta_y, delta_z),
		Ezy2(n, m, k, delta_x, delta_y, delta_z),

		Bxy2(n, m, k, delta_x, delta_y, delta_z),
		Bxz2(n, m, k, delta_x, delta_y, delta_z),
		Byx2(n, m, k, delta_x, delta_y, delta_z),
		Byz2(n, m, k, delta_x, delta_y, delta_z),
		Bzx2(n, m, k, delta_x, delta_y, delta_z),
		Bzy2(n, m, k, delta_x, delta_y, delta_z) {}
};