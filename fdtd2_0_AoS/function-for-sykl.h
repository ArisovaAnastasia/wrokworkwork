#include <CL/sycl.hpp>
#include "type data.h"

template <typename type_data>
class data3d_sycl {
	type_data* data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	data3d_sycl(sycl::queue q, int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		size_t array_size = (n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2);
		this->data = malloc_shared<type_data>(array_size, q);

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}
	inline type_data& operator()(int i, int j, int l) {
		return data[i * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2) + j * (k + 2 * delta_z + 2) + l];
	}
	type_data* ReturnData() {
		return data;
	}
};

template <class ftype, class ftypePML>
void Initializing_cube_split_3_dimen_PML_sycl(data3d_sycl<Component<ftype>>& cube, data3d_sycl<ComponentSplit<ftypePML>>& cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	for (int i = 0; i < Nx + 2 * delta_x + 2; i++)
		for (int j = 0; j < Ny + 2 * delta_y + 2; j++)
			for (int k = 0; k < Nz + 2 * delta_z + 2; k++) {
				cube_split(i, j, k).Exy = (ftypePML)0.0;
				cube_split(i, j, k).Exz = (ftypePML)0.0;
				cube_split(i, j, k).Eyx = (ftypePML)0.0;
				cube_split(i, j, k).Eyz = (ftypePML)0.0;
				cube_split(i, j, k).Ezx = (ftypePML)0.0;
				cube_split(i, j, k).Ezy = (ftypePML)0.0;

				cube_split(i, j, k).Bxy = (ftypePML)0.0;
				cube_split(i, j, k).Bxz = (ftypePML)0.0;
				cube_split(i, j, k).Byx = (ftypePML)0.0;
				cube_split(i, j, k).Byz = (ftypePML)0.0;
				cube_split(i, j, k).Bzx = (ftypePML)0.0;
				cube_split(i, j, k).Bzy = (ftypePML)0.0;

				cube(i, j, k).Ex = (ftype)0.0;
				cube(i, j, k).Ey = (ftype)0.0;
				cube(i, j, k).Ez = (ftype)0.0;
				cube(i, j, k).Bx = (ftype)0.0;
				cube(i, j, k).By = (ftype)0.0;
				cube(i, j, k).Bz = (ftype)0.0;		
			}
}

template <class ftypePML>
void Initializing_Sigma_3_dimen_PML_sycl(data3d_sycl<SIGMA<ftypePML>>& Sigma, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, double n,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;

	//cout << "delta_x= "<< delta_x << "  delta_y= " << delta_y << "   delta_z= " << delta_z << "    " << endl;
	//cout << "Nx= " << Nx << "  Ny= " << Ny << "   Nz= " << Nz << "    " << endl;
	//cout << "n = " << n << endl;
	cout << var_Sigma_max_x << "  " << var_Sigma_max_y << "   " << var_Sigma_max_z << "    " << std::endl;

	//cout << func_R(n, var_Sigma_max_x, delta, dx) << endl;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				Sigma(i, j, k).sigmaE_x = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta_x, i), n);
				Sigma(i, j, k).sigmaH_x = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta_x, i), n);

				Sigma(i, j, k).sigmaE_y = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta_y, j), n);
				Sigma(i, j, k).sigmaH_y = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta_y, j), n);

				Sigma(i, j, k).sigmaE_z = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta_z, k), n);
				Sigma(i, j, k).sigmaH_z = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
	//cout << sizeof(Sigma(10, 18, 11).sigmaH_x) << sizeof(double) << sizeof(float) << endl;
}

template <class ftypePML, class ftype = double>
void Initializing_Coeff_3_dimen_PML_correct_sycl(data3d_sycl<COEFF<ftypePML>>& Coeff, data3d_sycl<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz,
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

					Exy1 = exp(-dt * Sigma(i, j, k).sigmaE_y);
					Exz1 = exp(-dt * Sigma(i, j, k).sigmaE_z);
					Eyx1 = exp(-dt * Sigma(i, j, k).sigmaE_x);
					Eyz1 = exp(-dt * Sigma(i, j, k).sigmaE_z);
					Ezx1 = exp(-dt * Sigma(i, j, k).sigmaE_x);
					Ezy1 = exp(-dt * Sigma(i, j, k).sigmaE_y);

					Bxy1 = exp(-dt * Sigma(i, j, k).sigmaH_y);
					Bxz1 = exp(-dt * Sigma(i, j, k).sigmaH_z);
					Byx1 = exp(-dt * Sigma(i, j, k).sigmaH_x);
					Byz1 = exp(-dt * Sigma(i, j, k).sigmaH_z);
					Bzx1 = exp(-dt * Sigma(i, j, k).sigmaH_x);
					Bzy1 = exp(-dt * Sigma(i, j, k).sigmaH_y);

					Coeff(i, j, k).Exy1 = (ftypePML)Exy1;
					Coeff(i, j, k).Exz1 = (ftypePML)Exz1;
					Coeff(i, j, k).Eyx1 = (ftypePML)Eyx1;
					Coeff(i, j, k).Eyz1 = (ftypePML)Eyz1;
					Coeff(i, j, k).Ezx1 = (ftypePML)Ezx1;
					Coeff(i, j, k).Ezy1 = (ftypePML)Ezy1;

					Coeff(i, j, k).Bxy1 = (ftypePML)Bxy1;
					Coeff(i, j, k).Bxz1 = (ftypePML)Bxz1;
					Coeff(i, j, k).Byx1 = (ftypePML)Byx1;
					Coeff(i, j, k).Byz1 = (ftypePML)Byz1;
					Coeff(i, j, k).Bzx1 = (ftypePML)Bzx1;
					Coeff(i, j, k).Bzy1 = (ftypePML)Bzy1;

					if (Sigma(i, j, k).sigmaE_x != (ftypePML)0.0) {
						Coeff(i, j, k).Eyx2 = 1.0 / Sigma(i, j, k).sigmaE_x - Eyx1 / Sigma(i, j, k).sigmaE_x;
						Coeff(i, j, k).Ezx2 = 1.0 / Sigma(i, j, k).sigmaE_x - Ezx1 / Sigma(i, j, k).sigmaE_x;
					}
					else {
						Coeff(i, j, k).Eyx2 = dt;
						Coeff(i, j, k).Ezx2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_y != (ftypePML)0.0) {
						Coeff(i, j, k).Exy2 = 1.0 / Sigma(i, j, k).sigmaE_y - Exy1 / Sigma(i, j, k).sigmaE_y;
						Coeff(i, j, k).Ezy2 = 1.0 / Sigma(i, j, k).sigmaE_y - Ezy1 / Sigma(i, j, k).sigmaE_y;
					}
					else {
						Coeff(i, j, k).Exy2 = dt;
						Coeff(i, j, k).Ezy2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Exz2 = 1.0 / Sigma(i, j, k).sigmaE_z - Exz1 / Sigma(i, j, k).sigmaE_z;
						Coeff(i, j, k).Eyz2 = 1.0 / Sigma(i, j, k).sigmaE_z - Eyz1 / Sigma(i, j, k).sigmaE_z;
					}
					else {
						Coeff(i, j, k).Exz2 = dt;
						Coeff(i, j, k).Eyz2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_x != (ftypePML)0.0)
					{
						Coeff(i, j, k).Byx2 = 1.0 / Sigma(i, j, k).sigmaH_x - Byx1 / Sigma(i, j, k).sigmaH_x;
						Coeff(i, j, k).Bzx2 = 1.0 / Sigma(i, j, k).sigmaH_x - Bzx1 / Sigma(i, j, k).sigmaH_x;
					}
					else {
						Coeff(i, j, k).Byx2 = dt;
						Coeff(i, j, k).Bzx2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_y != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxy2 = 1.0 / Sigma(i, j, k).sigmaH_y - Bxy1 / Sigma(i, j, k).sigmaH_y;
						Coeff(i, j, k).Bzy2 = 1.0 / Sigma(i, j, k).sigmaH_y - Bzy1 / Sigma(i, j, k).sigmaH_y;
					}
					else {
						Coeff(i, j, k).Bxy2 = dt;
						Coeff(i, j, k).Bzy2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxz2 = 1.0 / Sigma(i, j, k).sigmaH_z - Bxz1 / Sigma(i, j, k).sigmaH_z;
						Coeff(i, j, k).Byz2 = 1.0 / Sigma(i, j, k).sigmaH_z - Byz1 / Sigma(i, j, k).sigmaH_z;
					}
					else {
						Coeff(i, j, k).Bxz2 = dt;
						Coeff(i, j, k).Byz2 = dt;
					}
				}
			}

}

template <class ftype>
void Add_Currents_3_dimen_sycl(data3d_sycl<Component<ftype>>& cube,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, const ftype dt, int it,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd, pair <ftype, ftype> fg)
{
	ftype x1, x2, y1, y2, z1, z2, t1, t2;

	ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.; // поперечное
	ftype ay_transverse = cd.second * (ftype)1. / (ftype)6.;
	ftype az_transverse = fg.second * (ftype)1. / (ftype)6.;

	//ftype ax_transverse = ab.second * (ftype)1. / (ftype)2.; // поперечное
	//ftype ay_transverse = cd.second * (ftype)1. / (ftype)2.;
	//ftype az_transverse = fg.second * (ftype)1. / (ftype)2.;

	ftype tp_x = ab.second / (ftype)6.;
	ftype tp_y = cd.second / (ftype)6.;
	ftype tp_z = fg.second / (ftype)6.;

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
			cube(offset + index_start_x, j + index_start_y, k + index_start_z).Ey += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			cube(offset + index_start_x, j + index_start_y, k + index_start_z).Bz += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}

}

template <class ftype>
inline void Update_electric_field_three_dimen_sycl(data3d_sycl<Component<ftype>> cube,
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
inline void Update_magnetic_field_three_dimen_sycl(data3d_sycl<Component<ftype>> cube,
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
inline void Update_electric_field_three_dimen_PML_sycl(data3d_sycl<Component<ftype>> cube,
	data3d_sycl<ComponentSplit<ftypePML>> cube_split, data3d_sycl<COEFF<ftypePML>> Coeff,
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
inline void Update_magnetic_field_three_dimen_PML_sycl(data3d_sycl<Component<ftype>> cube,
	data3d_sycl<ComponentSplit<ftypePML>> cube_split, data3d_sycl<COEFF<ftypePML>> Coeff,
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


template <class ftype>
double CalculateEnergyPML_3dimen_sycl(data3d_sycl<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += (double)cube(i, j, k).Ex * (double)cube(i, j, k).Ex
					+ (double)cube(i, j, k).Ey * (double)cube(i, j, k).Ey
					+ (double)cube(i, j, k).Ez * (double)cube(i, j, k).Ez;
				energy += (double)cube(i, j, k).Bx * (double)cube(i, j, k).Bx
					+ (double)cube(i, j, k).By * (double)cube(i, j, k).By
					+ (double)cube(i, j, k).Bz * (double)cube(i, j, k).Bz;
			}
	return energy;
}

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
