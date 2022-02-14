#pragma once
#define _USE_MATH_DEFINES

#include "type data.h"

using namespace std;

template <typename type_data>
class data3d {
	type_data* data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	data3d(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		this->data = (type_data*)malloc((n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2)*sizeof(type_data));

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}
	inline type_data& operator()(int i, int j, int l) {
		return data[l * (m + 2 * delta_y + 2) * (n + 2 * delta_x + 2) + j * (n + 2 * delta_x + 2) + i];
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

template <class type>
type Kahan(type a, type b, type ca, type cb)
{
	return (((a - b) - ca) + cb);
}

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML_Kahan(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	data3d<ComponentSplit<ftypePML>>& compensatorPML, data3d<ComponentSplit<ftypePML>>& compensatorPML2,
	data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;
	ftypePML sum, y;
	//
	sum = cube_split(i, j, k).Exy * Coeff(i, j, k).Exy1;
	y = Kahan<ftypePML>((ftypePML)cube(i, j + 1, k).Bz, (ftypePML)cube(i, j, k).Bz, (ftypePML)compensator(i, j + 1, k).Bz, (ftypePML)compensator(i, j, k).Bz) * Coeff(i, j, k).Exy2 * (_1dy)
		-compensatorPML(i, j, k).Exy * Coeff(i, j, k).Exy1;
	tExy = sum + y;
	compensatorPML2(i, j, k).Exy = (tExy - sum) - y;
	//
	sum = cube_split(i, j, k).Exz * Coeff(i, j, k).Exz1;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k + 1).By, (ftypePML)cube(i, j, k).By, (ftypePML)compensator(i, j, k + 1).By, (ftypePML)compensator(i, j, k).By) * Coeff(i, j, k).Exz2 * (_1dz)
		-compensatorPML(i, j, k).Exz * Coeff(i, j, k).Exz1;
	tExz = sum + y;
	compensatorPML2(i, j, k).Exz = (tExz - sum) - y;
	//
	sum = cube_split(i, j, k).Eyx * Coeff(i, j, k).Eyx1;
	y = -Kahan<ftypePML>((ftypePML)cube(i + 1, j, k).Bz, (ftypePML)cube(i, j, k).Bz, (ftypePML)compensator(i + 1, j, k).Bz, (ftypePML)compensator(i, j, k).Bz) * Coeff(i, j, k).Eyx2 * (_1dx)
		-compensatorPML(i, j, k).Eyx * Coeff(i, j, k).Eyx1;
	tEyx = sum + y;
	compensatorPML2(i, j, k).Eyx = (tEyx - sum) - y;
	//
	sum = cube_split(i, j, k).Eyz * Coeff(i, j, k).Eyz1;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k + 1).Bx, (ftypePML)cube(i, j, k).Bx, (ftypePML)compensator(i, j, k + 1).Bx, (ftypePML)compensator(i, j, k).Bx) * Coeff(i, j, k).Eyz2 * (_1dz)
		-compensatorPML(i, j, k).Eyz * Coeff(i, j, k).Eyz1;
	tEyz = sum + y;
	compensatorPML2(i, j, k).Eyz = (tEyz - sum) - y;
	//
	sum = cube_split(i, j, k).Ezx * Coeff(i, j, k).Ezx1;
	y = Kahan<ftypePML>((ftypePML)cube(i + 1, j, k).By, (ftypePML)cube(i, j, k).By, (ftypePML)compensator(i + 1, j, k).By, (ftypePML)compensator(i, j, k).By) * Coeff(i, j, k).Ezx2 * (_1dx)
		-compensatorPML(i, j, k).Ezx * Coeff(i, j, k).Ezx1;
	tEzx = sum + y;
	compensatorPML2(i, j, k).Ezx = (tEzx - sum) - y;
	//
	sum = cube_split(i, j, k).Ezy * Coeff(i, j, k).Ezy1;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j + 1, k).Bx, (ftypePML)cube(i, j, k).Bx, (ftypePML)compensator(i, j + 1, k).Bx, (ftypePML)compensator(i, j, k).Bx) * Coeff(i, j, k).Ezy2 * (_1dy)
		-compensatorPML(i, j, k).Ezy * Coeff(i, j, k).Ezy1;
	tEzy = sum + y;
	compensatorPML2(i, j, k).Ezy = (tEzy - sum) - y;
	//
	cube_split(i, j, k).Exy = tExy;
	cube_split(i, j, k).Exz = tExz;
	cube_split(i, j, k).Eyx = tEyx;
	cube_split(i, j, k).Eyz = tEyz;
	cube_split(i, j, k).Ezx = tEzx;
	cube_split(i, j, k).Ezy = tEzy;

	cube(i, j, k).Ex = (ftype)(tExy + tExz);
	cube(i, j, k).Ey = (ftype)(tEyx + tEyz);
	cube(i, j, k).Ez = (ftype)(tEzx + tEzy);

	compensator2(i, j, k).Ex = (ftype)(((ftypePML)cube(i, j, k).Ex - tExy) - tExz);
	compensator2(i, j, k).Ey = (ftype)(((ftypePML)cube(i, j, k).Ey - tEyx) - tEyz);
	compensator2(i, j, k).Ez = (ftype)(((ftypePML)cube(i, j, k).Ez - tEzx) - tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML_Kahan(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	data3d<ComponentSplit<ftypePML>>& compensatorPML, data3d<ComponentSplit<ftypePML>>& compensatorPML2,
	data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;
	ftypePML sum, y;
	//
	sum = cube_split(i, j, k).Bxy * Coeff(i, j, k).Bxy1;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ez, (ftypePML)cube(i, j - 1, k).Ez, (ftypePML)compensator(i, j, k).Ez, (ftypePML)compensator(i, j - 1, k).Ez) * Coeff(i, j, k).Bxy2 * (_1dy)
		-compensatorPML(i, j, k).Bxy * Coeff(i, j, k).Bxy1;
	tBxy = sum + y;
	compensatorPML2(i, j, k).Bxy = (tBxy - sum) - y;
	//
	sum = cube_split(i, j, k).Bxz * Coeff(i, j, k).Bxz1;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ey, (ftypePML)cube(i, j, k - 1).Ey, (ftypePML)compensator(i, j, k).Ey, (ftypePML)compensator(i, j, k - 1).Ey) * Coeff(i, j, k).Bxz2 * (_1dz)
		-compensatorPML(i, j, k).Bxz * Coeff(i, j, k).Bxz1;
	tBxz = sum + y;
	compensatorPML2(i, j, k).Bxz = (tBxz - sum) - y;
	//
	sum = cube_split(i, j, k).Byx * Coeff(i, j, k).Byx1;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ez, (ftypePML)cube(i - 1, j, k).Ez, (ftypePML)compensator(i, j, k).Ez, (ftypePML)compensator(i - 1, j, k).Ez) * Coeff(i, j, k).Byx2 * (_1dx)
		-compensatorPML(i, j, k).Byx * Coeff(i, j, k).Byx1;
	tByx = sum + y;
	compensatorPML2(i, j, k).Byx = (tByx - sum) - y;
	//
	sum = cube_split(i, j, k).Byz * Coeff(i, j, k).Byz1;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ex, (ftypePML)cube(i, j, k - 1).Ex, (ftypePML)compensator(i, j, k).Ex, (ftypePML)compensator(i, j, k - 1).Ex) * Coeff(i, j, k).Byz2 * (_1dz)
		-compensatorPML(i, j, k).Byz * Coeff(i, j, k).Byz1;
	tByz = sum + y;
	compensatorPML2(i, j, k).Byz = (tByz - sum) - y;
	//
	sum = cube_split(i, j, k).Bzx * Coeff(i, j, k).Bzx1;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ey, (ftypePML)cube(i - 1, j, k).Ey, (ftypePML)compensator(i, j, k).Ey, (ftypePML)compensator(i - 1, j, k).Ey) * Coeff(i, j, k).Bzx2 * (_1dx)
		-compensatorPML(i, j, k).Bzx * Coeff(i, j, k).Bzx1;
	tBzx = sum + y;
	compensatorPML2(i, j, k).Bzx = (tBzx - sum) - y;
	//
	sum = cube_split(i, j, k).Bzy * Coeff(i, j, k).Bzy1;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ex, (ftypePML)cube(i, j - 1, k).Ex, (ftypePML)compensator(i, j, k).Ex, (ftypePML)compensator(i, j - 1, k).Ex) * Coeff(i, j, k).Bzy2 * (_1dy)
		-compensatorPML(i, j, k).Bzy * Coeff(i, j, k).Bzy1;
	tBzy = sum + y;
	compensatorPML2(i, j, k).Bzy = (tBzy - sum) - y;
	//
	cube_split(i, j, k).Bxy = tBxy;
	cube_split(i, j, k).Bxz = tBxz;
	cube_split(i, j, k).Byx = tByx;
	cube_split(i, j, k).Byz = tByz;
	cube_split(i, j, k).Bzx = tBzx;
	cube_split(i, j, k).Bzy = tBzy;

	cube(i, j, k).Bx = (ftype)(tBxy + tBxz);
	cube(i, j, k).By = (ftype)(tByx + tByz);
	cube(i, j, k).Bz = (ftype)(tBzx + tBzy);

	compensator2(i, j, k).Bx = (ftype)(((ftypePML)cube(i, j, k).Bx - tBxy) - tBxz);
	compensator2(i, j, k).By = (ftype)(((ftypePML)cube(i, j, k).By - tByx) - tByz);
	compensator2(i, j, k).Bz = (ftype)(((ftypePML)cube(i, j, k).Bz - tBzx) - tBzy);
}

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML_Kahan_2_0(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	data3d<ComponentSplit<ftypePML>>& compensatorPML, data3d<ComponentSplit<ftypePML>>& compensatorPML2,
	data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;
	ftypePML sum, y;
	//

	//c1*(A + b*c2)

	sum = cube_split(i, j, k).Exy;
	y = Kahan<ftypePML>((ftypePML)cube(i, j + 1, k).Bz, (ftypePML)cube(i, j, k).Bz, (ftypePML)compensator(i, j + 1, k).Bz, (ftypePML)compensator(i, j, k).Bz) * Coeff(i, j, k).Exy2 * (_1dy)
		-compensatorPML(i, j, k).Exy;
	tExy = sum + y;
	compensatorPML2(i, j, k).Exy = (tExy - sum) - y;
	//
	sum = cube_split(i, j, k).Exz;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k + 1).By, (ftypePML)cube(i, j, k).By, (ftypePML)compensator(i, j, k + 1).By, (ftypePML)compensator(i, j, k).By) * Coeff(i, j, k).Exz2 * (_1dz)
		-compensatorPML(i, j, k).Exz;
	tExz = sum + y;
	compensatorPML2(i, j, k).Exz = (tExz - sum) - y;
	//
	sum = cube_split(i, j, k).Eyx;
	y = -Kahan<ftypePML>((ftypePML)cube(i + 1, j, k).Bz, (ftypePML)cube(i, j, k).Bz, (ftypePML)compensator(i + 1, j, k).Bz, (ftypePML)compensator(i, j, k).Bz) * Coeff(i, j, k).Eyx2 * (_1dx)
		-compensatorPML(i, j, k).Eyx;
	tEyx = sum + y;
	compensatorPML2(i, j, k).Eyx = (tEyx - sum) - y;
	//
	sum = cube_split(i, j, k).Eyz;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k + 1).Bx, (ftypePML)cube(i, j, k).Bx, (ftypePML)compensator(i, j, k + 1).Bx, (ftypePML)compensator(i, j, k).Bx) * Coeff(i, j, k).Eyz2 * (_1dz)
		-compensatorPML(i, j, k).Eyz;
	tEyz = sum + y;
	compensatorPML2(i, j, k).Eyz = (tEyz - sum) - y;
	//
	sum = cube_split(i, j, k).Ezx;
	y = Kahan<ftypePML>((ftypePML)cube(i + 1, j, k).By, (ftypePML)cube(i, j, k).By, (ftypePML)compensator(i + 1, j, k).By, (ftypePML)compensator(i, j, k).By) * Coeff(i, j, k).Ezx2 * (_1dx)
		-compensatorPML(i, j, k).Ezx;
	tEzx = sum + y;
	compensatorPML2(i, j, k).Ezx = (tEzx - sum) - y;
	//
	sum = cube_split(i, j, k).Ezy;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j + 1, k).Bx, (ftypePML)cube(i, j, k).Bx, (ftypePML)compensator(i, j + 1, k).Bx, (ftypePML)compensator(i, j, k).Bx) * Coeff(i, j, k).Ezy2 * (_1dy)
		-compensatorPML(i, j, k).Ezy;
	tEzy = sum + y;
	compensatorPML2(i, j, k).Ezy = (tEzy - sum) - y;
	//
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

	compensator2(i, j, k).Ex = (ftype)(((ftypePML)cube(i, j, k).Ex - tExy) - tExz);
	compensator2(i, j, k).Ey = (ftype)(((ftypePML)cube(i, j, k).Ey - tEyx) - tEyz);
	compensator2(i, j, k).Ez = (ftype)(((ftypePML)cube(i, j, k).Ez - tEzx) - tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML_Kahan_2_0(data3d<Component<ftype>>& cube,
	data3d<ComponentSplit<ftypePML>>& cube_split, data3d<COEFF<ftypePML>>& Coeff,
	data3d<ComponentSplit<ftypePML>>& compensatorPML, data3d<ComponentSplit<ftypePML>>& compensatorPML2,
	data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;
	ftypePML sum, y;
	//
	sum = cube_split(i, j, k).Bxy;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ez, (ftypePML)cube(i, j - 1, k).Ez, (ftypePML)compensator(i, j, k).Ez, (ftypePML)compensator(i, j - 1, k).Ez) * Coeff(i, j, k).Bxy2 * (_1dy)
		-compensatorPML(i, j, k).Bxy;
	tBxy = sum + y;
	compensatorPML2(i, j, k).Bxy = (tBxy - sum) - y;
	//
	sum = cube_split(i, j, k).Bxz;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ey, (ftypePML)cube(i, j, k - 1).Ey, (ftypePML)compensator(i, j, k).Ey, (ftypePML)compensator(i, j, k - 1).Ey) * Coeff(i, j, k).Bxz2 * (_1dz)
		-compensatorPML(i, j, k).Bxz;
	tBxz = sum + y;
	compensatorPML2(i, j, k).Bxz = (tBxz - sum) - y;
	//
	sum = cube_split(i, j, k).Byx;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ez, (ftypePML)cube(i - 1, j, k).Ez, (ftypePML)compensator(i, j, k).Ez, (ftypePML)compensator(i - 1, j, k).Ez) * Coeff(i, j, k).Byx2 * (_1dx)
		-compensatorPML(i, j, k).Byx;
	tByx = sum + y;
	compensatorPML2(i, j, k).Byx = (tByx - sum) - y;
	//
	sum = cube_split(i, j, k).Byz;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ex, (ftypePML)cube(i, j, k - 1).Ex, (ftypePML)compensator(i, j, k).Ex, (ftypePML)compensator(i, j, k - 1).Ex) * Coeff(i, j, k).Byz2 * (_1dz)
		-compensatorPML(i, j, k).Byz;
	tByz = sum + y;
	compensatorPML2(i, j, k).Byz = (tByz - sum) - y;
	//
	sum = cube_split(i, j, k).Bzx;
	y = -Kahan<ftypePML>((ftypePML)cube(i, j, k).Ey, (ftypePML)cube(i - 1, j, k).Ey, (ftypePML)compensator(i, j, k).Ey, (ftypePML)compensator(i - 1, j, k).Ey) * Coeff(i, j, k).Bzx2 * (_1dx)
		-compensatorPML(i, j, k).Bzx;
	tBzx = sum + y;
	compensatorPML2(i, j, k).Bzx = (tBzx - sum) - y;
	//
	sum = cube_split(i, j, k).Bzy;
	y = Kahan<ftypePML>((ftypePML)cube(i, j, k).Ex, (ftypePML)cube(i, j - 1, k).Ex, (ftypePML)compensator(i, j, k).Ex, (ftypePML)compensator(i, j - 1, k).Ex) * Coeff(i, j, k).Bzy2 * (_1dy)
		-compensatorPML(i, j, k).Bzy;
	tBzy = sum + y;
	compensatorPML2(i, j, k).Bzy = (tBzy - sum) - y;
	//

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

	compensator2(i, j, k).Bx = (ftype)(((ftypePML)cube(i, j, k).Bx - tBxy) - tBxz);
	compensator2(i, j, k).By = (ftype)(((ftypePML)cube(i, j, k).By - tByx) - tByz);
	compensator2(i, j, k).Bz = (ftype)(((ftypePML)cube(i, j, k).Bz - tBzx) - tBzy);
}

template <class ftype>
void Update_electric_field_three_dimen_Kahan(data3d<Component<ftype>>& cube,
	data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;
	ftype sum, y;
	//исправленное значение
	sum = cube(i, j, k).Ex;
	y = Kahan<ftype>(cube(i, j + 1, k).Bz, cube(i, j, k).Bz, compensator(i, j + 1, k).Bz, compensator(i, j, k).Bz) * (dt_y)
		-Kahan<ftype>(cube(i, j + 1, k).By, cube(i, j, k).By, compensator(i, j + 1, k).By, compensator(i, j, k).By) * (dt_z)-compensator(i, j, k).Ex;
	tEx = sum + y;
	compensator2(i, j, k).Ex = (tEx - sum) - y;
	//
	sum = cube(i, j, k).Ey;
	y = Kahan<ftype>(cube(i, j, k + 1).Bx, cube(i, j, k).Bx, compensator(i, j, k + 1).Bx, compensator(i, j, k).Bx) * (dt_z)
		-Kahan<ftype>(cube(i + 1, j, k).Bz, cube(i, j, k).Bz, compensator(i + 1, j, k).Bz, compensator(i, j, k).Bz) * (dt_x)-compensator(i, j, k).Ey;
	tEy = sum + y;
	compensator2(i, j, k).Ey = (tEy - sum) - y;
	//
	sum = cube(i, j, k).Ez;
	y = Kahan<ftype>(cube(i + 1, j, k).By, cube(i, j, k).By, compensator(i + 1, j, k).By, compensator(i, j, k).By) * (dt_x)
		-Kahan<ftype>(cube(i, j + 1, k).Bx, cube(i, j, k).Bx, compensator(i, j + 1, k).Bx, compensator(i, j, k).Bx) * (dt_y)-compensator(i, j, k).Ez;
	tEz = sum + y;
	compensator2(i, j, k).Ez = (tEz - sum) - y;

	cube(i, j, k).Ex = tEx;
	cube(i, j, k).Ey = tEy;
	cube(i, j, k).Ez = tEz;
}

template <class ftype>
void Update_magnetic_field_three_dimen_Kahan(data3d<Component<ftype>>& cube, data3d<Component<ftype>>& compensator, data3d<Component<ftype>>& compensator2,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tBx, tBy, tBz;
	ftype sum, y;

	sum = cube(i, j, k).Bx;
	y = -Kahan<ftype>(cube(i, j, k).Ez, cube(i, j - 1, k).Ez, compensator(i, j, k).Ez, compensator(i, j - 1, k).Ez) * (dt_y)
		+Kahan<ftype>(cube(i, j, k).Ey, cube(i, j, k - 1).Ey, compensator(i, j, k).Ey, compensator(i, j, k - 1).Ey) * (dt_z)-compensator(i, j, k).Bx;
	tBx = sum + y;
	compensator2(i, j, k).Bx = (tBx - sum) - y;
	//
	sum = cube(i, j, k).By;
	y = -Kahan<ftype>(cube(i, j, k).Ex, cube(i, j, k - 1).Ex, compensator(i, j, k).Ex, compensator(i, j, k - 1).Ex) * (dt_z)
		+Kahan<ftype>(cube(i, j, k).Ez, cube(i - 1, j, k).Ez, compensator(i, j, k).Ez, compensator(i - 1, j, k).Ez) * (dt_x)-compensator(i, j, k).By;
	tBy = sum + y;
	compensator2(i, j, k).By = (tBy - sum) - y;
	//
	sum = cube(i, j, k).Bz;
	y = -Kahan<ftype>(cube(i, j, k).Ey, cube(i - 1, j, k).Ey, compensator(i, j, k).Ey, compensator(i - 1, j, k).Ey) * (dt_x)
		+Kahan<ftype>(cube(i, j, k).Ex, cube(i, j - 1, k).Ex, compensator(i, j, k).Ex, compensator(i, j - 1, k).Ex) * (dt_y)-compensator(i, j, k).Bz;
	tBz = sum + y;
	compensator2(i, j, k).Bz = (tBz - sum) - y;

	cube(i, j, k).Bx = tBx;
	cube(i, j, k).By = tBy;
	cube(i, j, k).Bz = tBz;
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
void Update_magnetic_field_three_dimen_PML(data3d<Component<ftype>>& cube,
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

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML_2_0(data3d<Component<ftype>>& cube,
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
void Update_magnetic_field_three_dimen_PML_2_0(data3d<Component<ftype>>& cube,
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
}
