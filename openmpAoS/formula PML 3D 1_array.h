#pragma once
#define _USE_MATH_DEFINES

#include "type data.h"

using namespace std;

template <typename type_data>
class data3d {
	//type_data* data;
	vector<type_data> data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	data3d(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		//this->data = (type_data*)malloc((n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2)*sizeof(type_data));
		data.resize((n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2));

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}
	inline type_data& operator()(int i, int j, int l) {
		return data[i * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2) + j * (k + 2 * delta_z + 2) + l];
	}
};

template <class ftypePML>
ftypePML distanceH(int N, int delta, int i)
{
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i);

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N) - (ftypePML)0.5;

	return dist / (ftypePML)delta;
}

template <class ftypePML>
ftypePML distanceE(int N, int delta, int i)
{ //���� ���������� ������������ PML-����
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i) - (ftypePML)0.5;

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N);

	return dist / (ftypePML)delta;
}


template <class ftype>
void Update_electric_field_three_dimen(data3d<Component<ftype>>& cube,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;

	tEx = (cube(i, j + 1, k).Bz - cube(i, j, k).Bz) * (dt_y)
		-(cube(i, j, k + 1).By - cube(i, j, k).By) * (dt_z)+cube(i, j, k).Ex;

	tEy = (cube(i, j, k + 1).Bx - cube(i, j, k).Bx) * (dt_z)
		-(cube(i + 1, j, k).Bz - cube(i, j, k).Bz) * (dt_x)+cube(i, j, k).Ey;

	tEz = (cube(i + 1, j, k).By - cube(i, j, k).By) * (dt_x)
		-(cube(i, j + 1, k).Bx - cube(i, j, k).Bx) * (dt_y)+cube(i, j, k).Ez;

	cube(i, j, k).Ex = tEx;
	cube(i, j, k).Ey = tEy;
	cube(i, j, k).Ez = tEz;
}

template <class ftype>
void Update_magnetic_field_three_dimen(data3d<Component<ftype>>& cube,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tBx, tBy, tBz;

	tBx = -(cube(i, j, k).Ez - cube(i, j - 1, k).Ez) * (dt_y)
		+(cube(i, j, k).Ey - cube(i, j, k - 1).Ey) * (dt_z)+cube(i, j, k).Bx;

	tBy = -(cube(i, j, k).Ex - cube(i, j, k - 1).Ex) * (dt_z)
		+(cube(i, j, k).Ez - cube(i - 1, j, k).Ez) * (dt_x)+cube(i, j, k).By;

	tBz = -(cube(i, j, k).Ey - cube(i - 1, j, k).Ey) * (dt_x)
		+(cube(i, j, k).Ex - cube(i, j - 1, k).Ex) * (dt_y)+cube(i, j, k).Bz;

	cube(i, j, k).Bx = tBx;
	cube(i, j, k).By = tBy;
	cube(i, j, k).Bz = tBz;
}

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//c1*(A + b*c2)

	tExy = cube_split(i, j, k).Exy + ((ftypePML)cube(i, j + 1, k).Bz - (ftypePML)cube(i, j, k).Bz) * Coeff(i, j, k).Exy2 * (_1dy);

	tExz = cube_split(i, j, k).Exz - ((ftypePML)cube(i, j, k + 1).By - (ftypePML)cube(i, j, k).By) * Coeff(i, j, k).Exz2 * (_1dz);

	tEyx = cube_split(i, j, k).Eyx - ((ftypePML)cube(i + 1, j, k).Bz - (ftypePML)cube(i, j, k).Bz) * Coeff(i, j, k).Eyx2 * (_1dx);

	tEyz = cube_split(i, j, k).Eyz + ((ftypePML)cube(i, j, k + 1).Bx - (ftypePML)cube(i, j, k).Bx) * Coeff(i, j, k).Eyz2 * (_1dz);

	tEzx = cube_split(i, j, k).Ezx + ((ftypePML)cube(i + 1, j, k).By - (ftypePML)cube(i, j, k).By) * Coeff(i, j, k).Ezx2 * (_1dx);

	tEzy = cube_split(i, j, k).Ezy - ((ftypePML)cube(i, j + 1, k).Bx - (ftypePML)cube(i, j, k).Bx) * Coeff(i, j, k).Ezy2 * (_1dy);

	tExy = tExy * Coeff(i, j, k).Exy1;
	tExz = tExz * Coeff(i, j, k).Exz1;
	tEyx = tEyx * Coeff(i, j, k).Eyx1;
	tEyz = tEyz * Coeff(i, j, k).Eyz1;
	tEzx = tEzx * Coeff(i, j, k).Ezx1;
	tEzy = tEzy * Coeff(i, j, k).Ezy1;

	cube_split(i, j, k).Exy = tExy;
	cube_split(i, j, k).Exz = tExz;
	cube_split(i, j, k).Eyx = tEyx;
	cube_split(i, j, k).Eyz = tEyz;
	cube_split(i, j, k).Ezx = tEzx;
	cube_split(i, j, k).Ezy = tEzy;

	cube(i, j, k).Ex = (ftype)(tExy + tExz);
	cube(i, j, k).Ey = (ftype)(tEyx + tEyz);
	cube(i, j, k).Ez = (ftype)(tEzx + tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split(i, j, k).Bxy - ((ftypePML)cube(i, j, k).Ez - (ftypePML)cube(i, j - 1, k).Ez) * Coeff(i, j, k).Bxy2 * (_1dy);

	tBxz = cube_split(i, j, k).Bxz + ((ftypePML)cube(i, j, k).Ey - (ftypePML)cube(i, j, k - 1).Ey) * Coeff(i, j, k).Bxz2 * (_1dz);

	tByx = cube_split(i, j, k).Byx + ((ftypePML)cube(i, j, k).Ez - (ftypePML)cube(i - 1, j, k).Ez) * Coeff(i, j, k).Byx2 * (_1dx);

	tByz = cube_split(i, j, k).Byz - ((ftypePML)cube(i, j, k).Ex - (ftypePML)cube(i, j, k - 1).Ex) * Coeff(i, j, k).Byz2 * (_1dz);

	tBzx = cube_split(i, j, k).Bzx - ((ftypePML)cube(i, j, k).Ey - (ftypePML)cube(i - 1, j, k).Ey) * Coeff(i, j, k).Bzx2 * (_1dx);

	tBzy = cube_split(i, j, k).Bzy + ((ftypePML)cube(i, j, k).Ex - (ftypePML)cube(i, j - 1, k).Ex) * Coeff(i, j, k).Bzy2 * (_1dy);

	tBxy = tBxy * Coeff(i, j, k).Bxy1;
	tBxz = tBxz * Coeff(i, j, k).Bxz1;
	tByx = tByx * Coeff(i, j, k).Byx1;
	tByz = tByz * Coeff(i, j, k).Byz1;
	tBzx = tBzx * Coeff(i, j, k).Bzx1;
	tBzy = tBzy * Coeff(i, j, k).Bzy1;

	cube_split(i, j, k).Bxy = tBxy;
	cube_split(i, j, k).Bxz = tBxz;
	cube_split(i, j, k).Byx = tByx;
	cube_split(i, j, k).Byz = tByz;
	cube_split(i, j, k).Bzx = tBzx;
	cube_split(i, j, k).Bzy = tBzy;

	cube(i, j, k).Bx = (ftype)(tBxy + tBxz);
	cube(i, j, k).By = (ftype)(tByx + tByz);
	cube(i, j, k).Bz = (ftype)(tBzx + tBzy);
}

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML_old(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//c1*(A + b*c2)

	tExy = cube_split(i, j, k).Exy * Coeff(i, j, k).Exy1 + ((ftypePML)cube(i, j + 1, k).Bz - (ftypePML)cube(i, j, k).Bz) * Coeff(i, j, k).Exy2 * (_1dy);

	tExz = cube_split(i, j, k).Exz * Coeff(i, j, k).Exz1 - ((ftypePML)cube(i, j, k + 1).By - (ftypePML)cube(i, j, k).By) * Coeff(i, j, k).Exz2 * (_1dz);

	tEyx = cube_split(i, j, k).Eyx * Coeff(i, j, k).Eyx1 - ((ftypePML)cube(i + 1, j, k).Bz - (ftypePML)cube(i, j, k).Bz) * Coeff(i, j, k).Eyx2 * (_1dx);

	tEyz = cube_split(i, j, k).Eyz * Coeff(i, j, k).Eyz1 + ((ftypePML)cube(i, j, k + 1).Bx - (ftypePML)cube(i, j, k).Bx) * Coeff(i, j, k).Eyz2 * (_1dz);

	tEzx = cube_split(i, j, k).Ezx * Coeff(i, j, k).Ezx1 + ((ftypePML)cube(i + 1, j, k).By - (ftypePML)cube(i, j, k).By) * Coeff(i, j, k).Ezx2 * (_1dx);

	tEzy = cube_split(i, j, k).Ezy * Coeff(i, j, k).Ezy1 - ((ftypePML)cube(i, j + 1, k).Bx - (ftypePML)cube(i, j, k).Bx) * Coeff(i, j, k).Ezy2 * (_1dy);

	cube_split(i, j, k).Exy = tExy;
	cube_split(i, j, k).Exz = tExz;
	cube_split(i, j, k).Eyx = tEyx;
	cube_split(i, j, k).Eyz = tEyz;
	cube_split(i, j, k).Ezx = tEzx;
	cube_split(i, j, k).Ezy = tEzy;

	cube(i, j, k).Ex = (ftype)(tExy + tExz);
	cube(i, j, k).Ey = (ftype)(tEyx + tEyz);
	cube(i, j, k).Ez = (ftype)(tEzx + tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML_old(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split(i, j, k).Bxy * Coeff(i, j, k).Bxy1 - ((ftypePML)cube(i, j, k).Ez - (ftypePML)cube(i, j - 1, k).Ez) * Coeff(i, j, k).Bxy2 * (_1dy);

	tBxz = cube_split(i, j, k).Bxz * Coeff(i, j, k).Bxz1 + ((ftypePML)cube(i, j, k).Ey - (ftypePML)cube(i, j, k - 1).Ey) * Coeff(i, j, k).Bxz2 * (_1dz);

	tByx = cube_split(i, j, k).Byx * Coeff(i, j, k).Byx1 + ((ftypePML)cube(i, j, k).Ez - (ftypePML)cube(i - 1, j, k).Ez) * Coeff(i, j, k).Byx2 * (_1dx);

	tByz = cube_split(i, j, k).Byz * Coeff(i, j, k).Byz1 - ((ftypePML)cube(i, j, k).Ex - (ftypePML)cube(i, j, k - 1).Ex) * Coeff(i, j, k).Byz2 * (_1dz);

	tBzx = cube_split(i, j, k).Bzx * Coeff(i, j, k).Bzx1 - ((ftypePML)cube(i, j, k).Ey - (ftypePML)cube(i - 1, j, k).Ey) * Coeff(i, j, k).Bzx2 * (_1dx);

	tBzy = cube_split(i, j, k).Bzy * Coeff(i, j, k).Bzy1 + ((ftypePML)cube(i, j, k).Ex - (ftypePML)cube(i, j - 1, k).Ex) * Coeff(i, j, k).Bzy2 * (_1dy);

	cube_split(i, j, k).Bxy = tBxy;
	cube_split(i, j, k).Bxz = tBxz;
	cube_split(i, j, k).Byx = tByx;
	cube_split(i, j, k).Byz = tByz;
	cube_split(i, j, k).Bzx = tBzx;
	cube_split(i, j, k).Bzy = tBzy;

	cube(i, j, k).Bx = (ftype)(tBxy + tBxz);
	cube(i, j, k).By = (ftype)(tByx + tByz);
	cube(i, j, k).Bz = (ftype)(tBzx + tBzy);
}
