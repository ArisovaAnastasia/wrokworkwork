#pragma once
#include <string>
#include <iostream>
#include <chrono>
//#include <mpi.h>
#include "init_and_bound_conditions.h"
#include "update_formula_SoA.h"

template <class ftype>
void Check_Curant(ftype dx, ftype dy, ftype dz, ftype dt)
{
	double temp = dt * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy));
	std::cout << temp << " < 1  ";
	if (temp < 1)
		std::cout << "TRUE";
	else std::cout << "FALSE";
	std::cout << std::endl;
}

template <class ftype>
double CalculateEnergy_SoA(Fields<ftype>& field)
{
	int Nx = field.Ex.Get_Nx();
	int Ny = field.Ex.Get_Ny();
	int Nz = field.Ex.Get_Nz();
	int delta_x = field.Ex.Get_deltaX();
	int delta_y = field.Ex.Get_deltaY();
	int delta_z = field.Ex.Get_deltaZ();

	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += (double)field.Ex(i, j, k) * (double)field.Ex(i, j, k)
					+ (double)field.Ey(i, j, k) * (double)field.Ey(i, j, k)
					+ (double)field.Ez(i, j, k) * (double)field.Ez(i, j, k);
				energy += (double)field.Bx(i, j, k) * (double)field.Bx(i, j, k)
					+ (double)field.By(i, j, k) * (double)field.By(i, j, k)
					+ (double)field.Bz(i, j, k) * (double)field.Bz(i, j, k);
			}
	return energy;
}


template <class ftype, class ftypePML>
double run_fdtd(double n, int Nx, int Ny, int Nz, ftype T, ftype dt,
	int delta_x, int delta_y, int delta_z, double sigma_x, double sigma_y)
{
	//setlocale(LC_ALL, "Russian");

	int Nt = ceil(T / dt);
	std::cout << Nt << std::endl;

	Fields<ftype> field(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	SplitFields<ftypePML> split_field(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	Sigma<ftypePML> arr_sigma(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	Coefficient<ftypePML> arr_coef(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	std::pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny, dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;
	std::vector<double> vec_energy(Nt + 1);

	Init_Sigma_SoA<double>(arr_sigma, n, Nx, Ny, Nz, delta_x, delta_y, delta_z, sigma_x, sigma_y);
	Init_Coeff_SoA<ftypePML>(arr_coef, arr_sigma, dt);

	Check_Curant(dx, dy, dz, dt);

	auto start = std::chrono::steady_clock::now();

	for (int it = 0; it < Nt; ++it)
	{
		if (it % 100 == 0) {
			std::cout << it << std::endl;
		}

		//std::cout << it << std::endl;

		Update_boundary_condit<ftype>(field, dx, dy, dz, dt, it, ab, cd, fg);
		//vec_energy[it] = CalculateEnergy_SoA(field);
        //std::cout << vec_energy[it]<<std::endl;

#pragma omp parallel for collapse(2)
			//внутренняя область
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; ++j) {
				for (int k = delta_z + 1; k < Nz + delta_z + 1; ++k) {
					Update_electric_field_SoA<ftype>(field, dt_x, dt_y, dt_z, i, j, k);
				}
			}
		//грань от 0,0,0 по хy
#pragma omp parallel for collapse(2)
		for (int i = 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = 1; k < delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,0 по хz
#pragma omp parallel for collapse(2)
		for (int i = 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = 1; j < delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,0 по yz
#pragma omp parallel for collapse(2)
		for (int i = 1; i < delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от n,0,0 по yz
#pragma omp parallel for collapse(2)
		for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,m,0 по xz
#pragma omp parallel for collapse(2)
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,k по xy
#pragma omp parallel for collapse(2)
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; ++j)
				for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_electric_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}

#pragma omp parallel for collapse(2)
		//внутренняя область
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; ++k) {
					Update_magnetic_field_SoA<ftype>(field, dt_x, dt_y, dt_z, i, j, k);
				}
		//грань от 0,0,0 по хy
#pragma omp parallel for collapse(2)
		for (int i = 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = 1; k < delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,0 по хz
#pragma omp parallel for collapse(2)
		for (int i = 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = 1; j < delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,0 по yz
#pragma omp parallel for collapse(2)
		for (int i = 1; i < delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от n,0,0 по yz
#pragma omp parallel for collapse(2)
		for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,m,0 по xz
#pragma omp parallel for collapse(2)
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; ++j)
				for (int k = delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
		//грань от 0,0,k по xy
#pragma omp parallel for collapse(2)
		for (int i = delta_x + 1; i < Nx + delta_x + 1; ++i)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; ++j)
				for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; ++k) {
					Update_magnetic_field_PML_SoA<ftype, ftypePML>(field, split_field, arr_coef,
						_1dx, _1dy, _1dz, i, j, k);
				}
	}

	auto end = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	//vec_energy[Nt] = CalculateEnergy_SoA(field);
	//double result = vec_energy[Nt] / vec_energy[Nt/2];
	//std::cout << result << std::endl;

	std::cout << "elapsed time = " << elapsed_seconds.count() << "\n";


	return 0;
}

