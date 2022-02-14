//==============================================================
// Vector Add is the equivalent of a Hello, World! sample for data parallel
// programs. Building and running the sample verifies that your development
// environment is setup correctly and demonstrates the use of the core features
// of DPC++. This sample runs on both CPU and GPU (or FPGA). When run, it
// computes on both the CPU and offload device, then compares results. If the
// code executes on both CPU and offload device, the device name and a success
// message are displayed. And, your development environment is setup correctly!
//
// For comprehensive instructions regarding DPC++ Programming, go to
// https://software.intel.com/en-us/oneapi-programming-guide and search based on
// relevant terms noted in the comments.
//
// DPC++ material used in the code sample:
// •	A one dimensional array of data shared between CPU and offload device.
// •	A device queue and kernel.
//==============================================================
// Copyright © Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#define _USE_MATH_DEFINES
#include <CL/sycl.hpp>
#include <array>
#include <iostream>
#include <math.h>
#include <chrono>
#include "formula PML 3D 1_array.h"
#include "type data.h"
#include "function-for-sykl.h"

#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

using namespace sycl;
using namespace std;

// Parametres for this example.
#define ftype double
#define ftypePML double
double n = 4;
ftype T = 8.0 / 3.0 * M_PI, dt = 0.005;
int delta_x = 8, delta_y = 8, delta_z = 8;
double sigma_x = 46.5, sigma_y =46.5;
int Nx = 256, Ny = 128, Nz = 32;

pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.;
ftype ay_transverse = cd.second * (ftype)1. / (ftype)6.;
ftype az_transverse = fg.second * (ftype)1. / (ftype)6.;

ftype tp_x = ab.second / (ftype)6.;
ftype tp_y = cd.second / (ftype)6.;
ftype tp_z = fg.second / (ftype)6.;

ftype lambda_ = (ftype)2. * (ftype)M_PI / (ftype)4.;
ftype k_ = (ftype)2. * (ftype)M_PI / lambda_;
ftype omega_ = (ftype)2. * (ftype)M_PI / lambda_;
ftype A = (ftype)1.; //amplitude
ftype x00 = (ab.second - ab.first) / (ftype)2.;
ftype y00 = (cd.second - cd.first) / (ftype)2.;
ftype z00 = (fg.second - fg.first) / (ftype)2.;
ftype t0_x = (ftype)3. * tp_x;
ftype t0_y = (ftype)3. * tp_y;
ftype t0_z = (ftype)3. * tp_z;

int index_start_x = delta_x + 1;
int index_start_y = delta_y + 1;
int index_start_z = delta_z + 1;

int offset = 5;

template <typename type_data>
class data3d_sycl;

// Create an exception handler for asynchronous SYCL exceptions
static auto exception_handler = [](sycl::exception_list e_list) {
	for (std::exception_ptr const& e : e_list) {
		try {
			std::rethrow_exception(e);
		}
		catch (std::exception const& e) {
#if _DEBUG
			std::cout << "Failure" << std::endl;
#endif
			std::terminate();
		}
	}
};


//************************************
// Initialize the array from 0 to array_size - 1
//************************************
void InitializeArray(double* a, size_t size) {
	for (size_t i = 0; i < size; i++) a[i] = 0.0;
}

//************************************
// Demonstrate vector add both in sequential on CPU and in parallel on device.
//************************************
int main(int argc, char* argv[]) {
	// Change array_size if it was passed as argument
	//if (argc > 1) array_size = std::stoi(argv[1]);
	// Create device selector for the device of your interest.
#if FPGA_EMULATOR
  // DPC++ extension: FPGA emulator selector on systems without FPGA card.
	INTEL::fpga_emulator_selector d_selector;
#elif FPGA
  // DPC++ extension: FPGA selector on systems with FPGA card.
	INTEL::fpga_selector d_selector;
#else
  // The default device selector will select the most performant device.
	cpu_selector cpu_selector;
#endif
	int Nt = (T / dt);
	try {
		queue q(cpu_selector, exception_handler);

		// Print out the device information used for the kernel code.
		std::cout << "Running on device: "
			<< q.get_device().get_info<info::device::name>() << "\n";
		std::cout << "T =  " << T << "  dt = " << dt << "  Nt = " << Nt << "\n"
			<< "Nx = " << Nx << "  Ny = " << Ny << "  Nz = " << Nz << "\n"
			<< "delta_x = " << delta_x << "  delta_y = " << delta_y << "  delta_z = " << delta_z << "\n"
			<< "sigma_x = " << sigma_x << "  sigma_y = " << sigma_y << "  n = " << n << "\n";

		// Create grid with parametres to store input and output data. Allocate
		// unified shared memory so that both CPU and device can access them.
		data3d_sycl<Component<ftype>> cube(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		data3d_sycl<ComponentSplit<ftypePML>> cube_split(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		data3d_sycl<SIGMA<double>> Sigma(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		data3d_sycl<COEFF<ftypePML>> Coeff(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		double* vec_energy = malloc_shared<double>(size_t(Nt + 1), q);

		
		if ((cube.ReturnData() == nullptr) || (cube_split.ReturnData() == nullptr) || (Sigma.ReturnData() == nullptr) ||
			(Coeff.ReturnData() == nullptr) || (vec_energy == nullptr)) {
			if (cube.ReturnData() != nullptr) free(cube.ReturnData(), q);
			if (cube_split.ReturnData() != nullptr) free(cube_split.ReturnData(), q);
			if (Sigma.ReturnData() != nullptr) free(Sigma.ReturnData(), q);
			if (Coeff.ReturnData() != nullptr) free(Coeff.ReturnData(), q);
			//if (vec_energy != nullptr) free(vec_energy, q);

			std::cout << "Shared memory allocation failure.\n";
			return -1;
		}
		

		//Calculating the values of some parameters

		pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
		ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
			dz = (fg.second - fg.first) / (ftype)Nz;

		ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
		ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;

		// Initialize grid 
		Initializing_cube_split_3_dimen_PML_sycl<ftypePML>(cube, cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		Initializing_Sigma_3_dimen_PML_sycl<double>(Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, n, sigma_x, sigma_y);
		//InitializeArray(vec_energy, Nt + 1);

		

		Initializing_Coeff_3_dimen_PML_correct_sycl<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);

		// FDTD method in DPC++.
		std::cout << Nt << std::endl;
        auto start = std::chrono::steady_clock::now();
        
		for (int it = 0; it < Nt; ++it) {
			
			if (it % 100 == 0)
				cout << it << std::endl;
            
			//bound conditions
            ftype t1 = dt * (ftype)it;
            ftype t2 = t1 + (ftype)0.5 * dt;
            ftype x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
            ftype x2 = (ftype)(ab.first + dx * (ftype)offset);
            
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<2>(Nz, Ny), [=](auto index) {
					int j = index[1];
					int k = index[0];
                    auto kernel_cube = cube;
                    
                    ftype y1, y2, z1, z2;
                    y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
                    y2 = (ftype)(cd.first + dy * (ftype)j);

                    z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
                    z2 = (ftype)(fg.first + dz * (ftype)k);

                    // Ez and By has the same y coordinate
                    kernel_cube(offset + index_start_x, j + index_start_y, k + index_start_z).Ey += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
                        * exp(-(y2 - y00) * (y2 - y00) / (ay_transverse * ay_transverse)) *
                        sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
                    kernel_cube(offset + index_start_x, j + index_start_y, k + index_start_z).Bz += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
                        * exp(-(y2 - y00) * (y2 - y00) / (ay_transverse * ay_transverse)) *
                        sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
					});
				}).wait_and_throw();
            
			//vec_energy[it] = CalculateEnergyPML_3dimen_sycl(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//std::cout << vec_energy[it] << std::endl;
            
            //run calculate electric fields
            //internal area
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz, Ny, Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

                    Update_electric_field_three_dimen_sycl<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by XY
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(delta_z, Ny + 2 * delta_y, Nx + 2 * delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + 1;
					int k = index[0] + 1;

					Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by XZ
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz + delta_z, delta_y, Nx + 2 * delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + 1;
					int k = index[0] + delta_z + 1;

                    Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by YZ
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz + delta_z, Ny + delta_y, delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from n,0,0 by YZ
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz + delta_z, Ny + delta_y, delta_x), [=](auto index) {
					int i = index[2] + Nx + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,m,0 by XZ
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>( Nz + delta_z, delta_y, Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + Ny + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,k by XY
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(delta_z, Ny, Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + Nz + delta_z + 1;

					Update_electric_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
			//run calculate magnetic fields
            //internal area
            q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz, Ny, Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_magnetic_field_three_dimen_sycl<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by XY
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(delta_z, Ny + 2 * delta_y, Nx + 2 * delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + 1;
					int k = index[0] + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by XZ
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>( delta_z, delta_y, Nx + 2 * delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + 1;
					int k = index[0] + delta_z + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,0 by YZ
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(delta_z, delta_y, delta_x), [=](auto index) {
					int i = index[2] + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from n,0,0 by YZ
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(Nz + delta_z, Ny + delta_y,  delta_x), [=](auto index) {
					int i = index[2] + Nx + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,m,0 by XZ
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>( Nz + delta_z, delta_y, Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + Ny + delta_y + 1;
					int k = index[0] + delta_z + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
            //edge from 0,0,k by XY
			q.submit([&](sycl::handler& h) {
				h.parallel_for(range<3>(delta_z, Ny,  Nx), [=](auto index) {
					int i = index[2] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[0] + Nz + delta_z + 1;

					Update_magnetic_field_three_dimen_PML_sycl<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
					});
				}).wait_and_throw();
		}

		auto end = std::chrono::steady_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
/*
		vec_energy[Nt] = CalculateEnergyPML_3dimen_sycl(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		double result = vec_energy[Nt] / vec_energy[2512];

		// Print out the result,

		std::cout << "reflection coefficient = " << result <<" " << vec_energy[Nt] / vec_energy[Nt/2] << "\n";
        */
		std::cout << "elapsed time = " << elapsed_seconds.count() << "\n";


		free(cube.ReturnData(), q);
		free(cube_split.ReturnData(), q);
		free(Sigma.ReturnData(), q);
		free(Coeff.ReturnData(), q);
		free(vec_energy, q);
	}
	catch (sycl::exception const& e) {
		std::cout << "An exception is caught while adding two vectors.\n";
		std::terminate();
	}

	std::cout << "Vector add successfully completed on device.\n";
	return 0;
}
