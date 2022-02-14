#pragma once

template <class ftype>
class Component
{
public:
	ftype Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
};

template <class ftypePML>
class SIGMA
{
public:
	ftypePML sigmaE_x, sigmaE_y, sigmaE_z, sigmaH_x, sigmaH_y, sigmaH_z;
};
template <class ftypePML>
class COEFF
{
public:
	ftypePML Exy1, Exz1, Eyx1, Eyz1, Ezx1, Ezy1, Bxy1, Bxz1, Byx1, Byz1, Bzx1, Bzy1,
		Exy2, Exz2, Eyx2, Eyz2, Ezx2, Ezy2, Bxy2, Bxz2, Byx2, Byz2, Bzx2, Bzy2;
};
template <class ftypePML>
class ComponentSplit
{
public:
	ftypePML Exy, Exz, Eyx, Eyz, Ezx, Ezy,
		Bxy, Bxz, Byx, Byz, Bzx, Bzy;
};
enum Direction {
	x_comp_yz, x_comp_zy, y_comp_xz, y_comp_zx, z_comp_xy, z_comp_yx,
};