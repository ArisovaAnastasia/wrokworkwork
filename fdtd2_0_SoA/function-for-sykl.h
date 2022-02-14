#include <CL/sycl.hpp>
#include "type data.h"

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
template <class ftype, class ftypePML>
void Initializing_cube_split_3_dimen_PML_sycl(Fields<ftype>& cube,
	SplitFields<ftypePML>& cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	for (int i = 0; i < Nx + 2 * delta_x + 2; i++)
		for (int j = 0; j < Ny + 2 * delta_y + 2; j++)
			for (int k = 0; k < Nz + 2 * delta_z + 2; k++) {
				cube_split.Exy(i, j, k) = (ftypePML)0.0;
				cube_split.Exz(i, j, k) = (ftypePML)0.0;
				cube_split.Eyx(i, j, k) = (ftypePML)0.0;
				cube_split.Eyz(i, j, k) = (ftypePML)0.0;
				cube_split.Ezx(i, j, k) = (ftypePML)0.0;
				cube_split.Ezy(i, j, k) = (ftypePML)0.0;

				cube_split.Bxy(i, j, k) = (ftypePML)0.0;
				cube_split.Bxz(i, j, k) = (ftypePML)0.0;
				cube_split.Byx(i, j, k) = (ftypePML)0.0;
				cube_split.Byz(i, j, k) = (ftypePML)0.0;
				cube_split.Bzx(i, j, k) = (ftypePML)0.0;
				cube_split.Bzy(i, j, k) = (ftypePML)0.0;

				cube.Ex(i, j, k) = (ftype)0.0;
				cube.Ey(i, j, k) = (ftype)0.0;
				cube.Ez(i, j, k) = (ftype)0.0;
				cube.Bx(i, j, k) = (ftype)0.0;
				cube.By(i, j, k) = (ftype)0.0;
				cube.Bz(i, j, k) = (ftype)0.0;
			}
}

template <class ftypePML>
void Initializing_Sigma_3_dimen_PML_sycl(Sigma<ftypePML>& Sigma,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, double n,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++) {
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				Sigma.sigmaEx(i, j, k) = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta_x, i), n);
				Sigma.sigmaBx(i, j, k) = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta_x, i), n);

				Sigma.sigmaEy(i, j, k) = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta_y, j), n);
				Sigma.sigmaBy(i, j, k) = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta_y, j), n);

				Sigma.sigmaEz(i, j, k) = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta_z, k), n);
				Sigma.sigmaBz(i, j, k) = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
		}
	}
	//cout << sizeof(Sigma(10, 18, 11).sigmaH_x) << sizeof(double) << sizeof(float) << endl;
}

template <class ftypePML, class ftype = double>
void Initializing_Coeff_3_dimen_PML_correct_sycl(Coefficient<ftypePML>& Coeff, Sigma<ftypePML>& Sigma, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, ftype dt) {

	double Exy1, Exz1, Ezx1, Ezy1, Eyx1, Eyz1;
	double Bxy1, Bxz1, Bzx1, Bzy1, Byx1, Byz1;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
				{
				}
				else {

					Exy1 = exp(-dt * Sigma.sigmaEy(i, j, k));
					Exz1 = exp(-dt * Sigma.sigmaEz(i, j, k));
					Eyx1 = exp(-dt * Sigma.sigmaEx(i, j, k));
					Eyz1 = exp(-dt * Sigma.sigmaEz(i, j, k));
					Ezx1 = exp(-dt * Sigma.sigmaEx(i, j, k));
					Ezy1 = exp(-dt * Sigma.sigmaEy(i, j, k));

					Bxy1 = exp(-dt * Sigma.sigmaBy(i, j, k));
					Bxz1 = exp(-dt * Sigma.sigmaBz(i, j, k));
					Byx1 = exp(-dt * Sigma.sigmaBx(i, j, k));
					Byz1 = exp(-dt * Sigma.sigmaBz(i, j, k));
					Bzx1 = exp(-dt * Sigma.sigmaBx(i, j, k));
					Bzy1 = exp(-dt * Sigma.sigmaBy(i, j, k));

					Coeff.Exy1(i, j, k) = (ftypePML)Exy1;
					Coeff.Exz1(i, j, k) = (ftypePML)Exz1;
					Coeff.Eyx1(i, j, k) = (ftypePML)Eyx1;
					Coeff.Eyz1(i, j, k) = (ftypePML)Eyz1;
					Coeff.Ezx1(i, j, k) = (ftypePML)Ezx1;
					Coeff.Ezy1(i, j, k) = (ftypePML)Ezy1;

					Coeff.Bxy1(i, j, k) = (ftypePML)Bxy1;
					Coeff.Bxz1(i, j, k) = (ftypePML)Bxz1;
					Coeff.Byx1(i, j, k) = (ftypePML)Byx1;
					Coeff.Byz1(i, j, k) = (ftypePML)Byz1;
					Coeff.Bzx1(i, j, k) = (ftypePML)Bzx1;
					Coeff.Bzy1(i, j, k) = (ftypePML)Bzy1;

					if (Sigma.sigmaEx(i, j, k) != (ftypePML)0.0) {
						Coeff.Eyx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Eyx1 / Sigma.sigmaEx(i, j, k);
						Coeff.Ezx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Ezx1 / Sigma.sigmaEx(i, j, k);
					}
					else {
						Coeff.Eyx2(i, j, k) = dt;
						Coeff.Ezx2(i, j, k) = dt;
					}
					if (Sigma.sigmaEy(i, j, k) != (ftypePML)0.0) {
						Coeff.Exy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Exy1 / Sigma.sigmaEy(i, j, k);
						Coeff.Ezy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Ezy1 / Sigma.sigmaEy(i, j, k);
					}
					else {
						Coeff.Exy2(i, j, k) = dt;
						Coeff.Ezy2(i, j, k) = dt;
					}
					if (Sigma.sigmaEz(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Exz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Exz1 / Sigma.sigmaEz(i, j, k);
						Coeff.Eyz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Eyz1 / Sigma.sigmaEz(i, j, k);
					}
					else {
						Coeff.Exz2(i, j, k) = dt;
						Coeff.Eyz2(i, j, k) = dt;
					}
					if (Sigma.sigmaBx(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Byx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Byx1 / Sigma.sigmaBx(i, j, k);
						Coeff.Bzx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Bzx1 / Sigma.sigmaBx(i, j, k);
					}
					else {
						Coeff.Byx2(i, j, k) = dt;
						Coeff.Bzx2(i, j, k) = dt;
					}
					if (Sigma.sigmaBy(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Bxy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bxy1 / Sigma.sigmaBy(i, j, k);
						Coeff.Bzy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bzy1 / Sigma.sigmaBy(i, j, k);
					}
					else {
						Coeff.Bxy2(i, j, k) = dt;
						Coeff.Bzy2(i, j, k) = dt;
					}
					if (Sigma.sigmaBz(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Bxz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Bxz1 / Sigma.sigmaBz(i, j, k);
						Coeff.Byz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Byz1 / Sigma.sigmaBz(i, j, k);
					}
					else {
						Coeff.Bxz2(i, j, k) = dt;
						Coeff.Byz2(i, j, k) = dt;
					}
				}
			}

}

template <class ftype>
void Add_Currents_3_dimen_sycl(Fields<ftype>& cube,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, const ftype dt, int it,
	std::pair <ftype, ftype> ab, std::pair <ftype, ftype> cd, std::pair <ftype, ftype> fg)
{
	ftype x1, x2, y1, y2, z1, z2, t1, t2;

	ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.; // поперечное
	ftype ay_transverse = cd.second * (ftype)1. / (ftype)6.;
	ftype az_transverse = fg.second * (ftype)1. / (ftype)6.;

	//ftype ax_transverse = ab.second * (ftype)1. / (ftype)2.; // поперечное
	//ftype ay_transverse = cd.second * (ftype)1. / (ftype)2.;
	//ftype az_transverse = fg.second * (ftype)1. / (ftype)2.;

	ftype tp_x = ab.second / (ftype)6.0;
	ftype tp_y = cd.second / (ftype)6.0;
	ftype tp_z = fg.second / (ftype)6.0;

	ftype lambda_ = (ftype)2. * (ftype)M_PI / (ftype)4.; // wvelength               //?
	ftype k_ = (ftype)2. * (ftype)M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	ftype omega_ = (ftype)2. * (ftype)M_PI / lambda_; // 2 pi c /wavelength 
	ftype A = (ftype)1.; //amplitude

	t1 = dt * (ftype)it;
	t2 = t1 + (ftype)0.5 * dt;

	ftype x0 = (ab.second - ab.first) / (ftype)2.;
	ftype y0 = (cd.second - cd.first) / (ftype)2.;
	ftype z0 = (fg.second - fg.first) / (ftype)2.;
	ftype t0_x = (ftype)3. * tp_x;
	ftype t0_y = (ftype)3. * tp_y;
	ftype t0_z = (ftype)3. * tp_z;

	int index_start_x = delta_x + 1, index_start_y = delta_y + 1, index_start_z = delta_z + 1;

	int offset = 5;
	x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
	x2 = (ftype)(ab.first + dx * (ftype)offset);

	for (int j = 0; j < Ny; j++)
		for (int k = 0; k < Nz; k++) {

			y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
			y2 = (ftype)(cd.first + dy * (ftype)j);

			z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
			z2 = (ftype)(fg.first + dz * (ftype)k);

			// Ez and By has the same y coordinate
			cube.Ey(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			cube.Bz(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}

}

template <class ftype>
void Update_electric_field_three_dimen_sycl(Fields<ftype> cube,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;

	tEx = (cube.Bz(i, j + 1, k) - cube.Bz(i, j, k)) * (dt_y)
		-(cube.By(i, j, k + 1) - cube.By(i, j, k)) * (dt_z)+cube.Ex(i, j, k);

	tEy = (cube.Bx(i, j, k + 1) - cube.Bx(i, j, k)) * (dt_z)
		-(cube.Bz(i + 1, j, k) - cube.Bz(i, j, k)) * (dt_x)+cube.Ey(i, j, k);

	tEz = (cube.By(i + 1, j, k) - cube.By(i, j, k)) * (dt_x)
		-(cube.Bx(i, j + 1, k) - cube.Bx(i, j, k)) * (dt_y)+cube.Ez(i, j, k);

	cube.Ex(i, j, k) = tEx;
	cube.Ey(i, j, k) = tEy;
	cube.Ez(i, j, k) = tEz;
}

template <class ftype>
void Update_magnetic_field_three_dimen_sycl(Fields<ftype> cube,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tBx, tBy, tBz;

	tBx = -(cube.Ez(i, j, k) - cube.Ez(i, j - 1, k)) * (dt_y)
		+(cube.Ey(i, j, k) - cube.Ey(i, j, k - 1)) * (dt_z)+cube.Bx(i, j, k);

	tBy = -(cube.Ex(i, j, k) - cube.Ex(i, j, k - 1)) * (dt_z)
		+(cube.Ez(i, j, k) - cube.Ez(i - 1, j, k)) * (dt_x)+cube.By(i, j, k);

	tBz = -(cube.Ey(i, j, k) - cube.Ey(i - 1, j, k)) * (dt_x)
		+(cube.Ex(i, j, k) - cube.Ex(i, j - 1, k)) * (dt_y)+cube.Bz(i, j, k);

	cube.Bx(i, j, k) = tBx;
	cube.By(i, j, k) = tBy;
	cube.Bz(i, j, k) = tBz;
}

template <class ftype, class ftypePML>
inline void Update_electric_field_three_dimen_PML_sycl(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//c1*(A + b*c2)

	tExy = cube_split.Exy(i, j, k) * Coeff.Exy1(i, j, k) + ((ftypePML)cube.Bz(i, j + 1, k) - (ftypePML)cube.Bz(i, j, k)) * Coeff.Exy2(i, j, k) * (_1dy);

	tExz = cube_split.Exz(i, j, k) * Coeff.Exz1(i, j, k) - ((ftypePML)cube.By(i, j, k + 1) - (ftypePML)cube.By(i, j, k)) * Coeff.Exz2(i, j, k) * (_1dz);

	tEyx = cube_split.Eyx(i, j, k) * Coeff.Eyx1(i, j, k) - ((ftypePML)cube.Bz(i + 1, j, k) - (ftypePML)cube.Bz(i, j, k)) * Coeff.Eyx2(i, j, k) * (_1dx);

	tEyz = cube_split.Eyz(i, j, k) * Coeff.Eyz1(i, j, k) + ((ftypePML)cube.Bx(i, j, k + 1) - (ftypePML)cube.Bx(i, j, k)) * Coeff.Eyz2(i, j, k) * (_1dz);

	tEzx = cube_split.Ezx(i, j, k) * Coeff.Ezx1(i, j, k) + ((ftypePML)cube.By(i + 1, j, k) - (ftypePML)cube.By(i, j, k)) * Coeff.Ezx2(i, j, k) * (_1dx);

	tEzy = cube_split.Ezy(i, j, k) * Coeff.Ezy1(i, j, k) - ((ftypePML)cube.Bx(i, j + 1, k) - (ftypePML)cube.Bx(i, j, k)) * Coeff.Ezy2(i, j, k) * (_1dy);

	cube_split.Exy(i, j, k) = tExy;
	cube_split.Exz(i, j, k) = tExz;
	cube_split.Eyx(i, j, k) = tEyx;
	cube_split.Eyz(i, j, k) = tEyz;
	cube_split.Ezx(i, j, k) = tEzx;
	cube_split.Ezy(i, j, k) = tEzy;

	cube.Ex(i, j, k) = (ftype)(tExy + tExz);
	cube.Ey(i, j, k) = (ftype)(tEyx + tEyz);
	cube.Ez(i, j, k) = (ftype)(tEzx + tEzy);
}

template <class ftype, class ftypePML>
inline void Update_magnetic_field_three_dimen_PML_sycl(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - ((ftypePML)cube.Ez(i, j, k) - (ftypePML)cube.Ez(i, j - 1, k)) * Coeff.Bxy2(i, j, k) * (_1dy);

	tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + ((ftypePML)cube.Ey(i, j, k) - (ftypePML)cube.Ey(i, j, k - 1)) * Coeff.Bxz2(i, j, k) * (_1dz);

	tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + ((ftypePML)cube.Ez(i, j, k) - (ftypePML)cube.Ez(i - 1, j, k)) * Coeff.Byx2(i, j, k) * (_1dx);

	tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - ((ftypePML)cube.Ex(i, j, k) - (ftypePML)cube.Ex(i, j, k - 1)) * Coeff.Byz2(i, j, k) * (_1dz);

	tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - ((ftypePML)cube.Ey(i, j, k) - (ftypePML)cube.Ey(i - 1, j, k)) * Coeff.Bzx2(i, j, k) * (_1dx);

	tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + ((ftypePML)cube.Ex(i, j, k) - (ftypePML)cube.Ex(i, j - 1, k)) * Coeff.Bzy2(i, j, k) * (_1dy);

	cube_split.Bxy(i, j, k) = tBxy;
	cube_split.Bxz(i, j, k) = tBxz;
	cube_split.Byx(i, j, k) = tByx;
	cube_split.Byz(i, j, k) = tByz;
	cube_split.Bzx(i, j, k) = tBzx;
	cube_split.Bzy(i, j, k) = tBzy;

	cube.Bx(i, j, k) = (ftype)(tBxy + tBxz);
	cube.By(i, j, k) = (ftype)(tByx + tByz);
	cube.Bz(i, j, k) = (ftype)(tBzx + tBzy);
}


template <class ftype>
double CalculateEnergyPML_3dimen_sycl(Fields<ftype>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += (double)cube.Ex(i, j, k) * (double)cube.Ex(i, j, k)
					+ (double)cube.Ey(i, j, k) * (double)cube.Ey(i, j, k)
					+ (double)cube.Ez(i, j, k) * (double)cube.Ez(i, j, k);
				energy += (double)cube.Bx(i, j, k) * (double)cube.Bx(i, j, k)
					+ (double)cube.By(i, j, k) * (double)cube.By(i, j, k)
					+ (double)cube.Bz(i, j, k) * (double)cube.Bz(i, j, k);
			}
	return energy;
}
/*
template <class ftype>
void Graph_Solution_two_planes_3d_python_sycl(data3d_sycl<Component<ftype>>& cube, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, ftype dt, ftype T,
	int n, double R, string file_name1, string file_name2)
{
	ofstream numbEy(file_name1), numbBz(file_name2);

	numbEy << Nx + 2 * delta_x << ";" << Ny + 2 * delta_y << std::endl;
	numbBz << Nx + 2 * delta_x << ";" << Ny + 2 * delta_y << std::endl;

	for (int i = 1; i < Ny + 2 * delta_y + 1; i++)
	{
		ftype y = dy * (ftype)i;

		for (int j = 1; j < Nx + 2 * delta_x + 1; j++)
		{
			numbEy << cube(j, i, (Nz / 2 + delta_z) ? (Nz / 2 + delta_z) : 1).Ey << ";";
			numbBz << cube(j, i, (Nz / 2 + delta_z) ? (Nz / 2 + delta_z) : 1).Bz << ";";
		}
		numbEy << std::endl;
		numbBz << std::endl;

	}
	numbEy << std::endl;
	numbEy << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << std::endl;
	numbEy << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << std::endl;
	numbEy << "dx = " << dx << ";   dy = " << dy << ";   dz = " << dz << std::endl;
	numbEy << "T = " << T << ";   dt = " << dt << ";   n = " << n << std::endl;
	numbEy << "R = " << R;

	numbBz << std::endl;
	numbBz << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << std::endl;
	numbBz << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << std::endl;
	numbBz << "dx = " << dx << ";   dy = " << dy << ";   dz = " << dz << std::endl;
	numbBz << "T = " << T << ";   dt = " << dt << ";   n = " << n << std::endl;
	numbBz << "R = " << R;

	numbEy.close();
	numbBz.close();
}
*/