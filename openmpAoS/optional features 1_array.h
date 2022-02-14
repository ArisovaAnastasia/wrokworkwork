#pragma once
#define _USE_MATH_DEFINES
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <clocale>
#include <vector>
#include <cmath>
#include "type data.h"

using namespace std;

template <class ftype>
void Check_Curant(ftype dx, ftype dy, ftype dz, ftype dt)
{
	//double temp = dt * sqrt(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)) ;
	double temp = dt * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy));
	cout << temp << " < 1  ";
	if (temp < 1)
		cout << "TRUE" << endl;
	else cout << "FALSE" << endl;
}

template <class ftypePML>
void Graph_for_Sigma_three_dimen(data3d<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftypePML dx, string file_sigma)
{
	//                         
	ofstream numb_sigma(file_sigma);

	numb_sigma << ";" << "value sigma_x" << endl;

	for (int s = 1; s < Nx + 2 * delta_x + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << Sigma(s, delta_y / 2, delta_z / 2 +1).sigmaH_x << endl;
		numb_sigma << dx * ((ftypePML)(s - 1) + 0.5) << ";" << Sigma(s, delta_y / 2, delta_z / 2+1).sigmaE_x << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
}

template <class ftypePML>
void Graph_for_Coeff_three_dimen(data3d<COEFF<ftypePML>>& Coeff, int Nx, int Ny, int Nz, int delta, ftypePML dx, string file_coeff)
{
	//                         
	ofstream numb_sigma(file_coeff);

	numb_sigma << ";" << "Coeff1" << ";;;" << "Coeff2" << ";;" << "Bxy" << endl;

	for (int s = 1; s < Ny + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << setprecision(16)<< Coeff(Nx + 1.5 * delta, s, 1).Bxy2 << ";;";
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << setprecision(16) << Coeff(Nx + 1.5 * delta, s, 1).Byx2 << ";;";
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16) << Coeff(Nx + 1.5 * delta, s, 1).Exy2 << ";;";
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16) << Coeff(Nx + 1.5*delta, s, 1).Eyx2 << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
}


template <class ftype>
double CalculateEnergyPML_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
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
double AverageValueMagnetic(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum += (double)cube(i, j, k).Bx + (double)cube(i, j, k).By + (double)cube(i, j, k).Bz;
			}
	return sum / (double)(3 * (Nx + 2 * delta_x) * (Ny + 2 * delta_y) * (Nz + 2 * delta_z));
}

template <class ftype>
double AverageValueElectric(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum += (double)cube(i, j, k).Ex + (double)cube(i, j, k).Ey + (double)cube(i, j, k).Ez;
			}
	return sum / (double)(3 * (Nx + 2 * delta_x) * (Ny + 2 * delta_y) * (Nz + 2 * delta_z));
}

template <class ftype>
double MaxValueElectric(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, abs((double)cube(i, j, k).Ex));
				sum = max(sum, abs((double)cube(i, j, k).Ey));
				sum = max(sum, abs((double)cube(i, j, k).Ez));
			}
	return sum;
}
template <class ftype>
double MaxValueCurrents(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, (double)cube(i, j, k).Jx);
				sum = max(sum, (double)cube(i, j, k).Jy);
				sum = max(sum, (double)cube(i, j, k).Jz);
			}
	return sum;
}

template <class ftype>
double MaxValueNe(data3d<ftype>& Ne, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, (double)Ne(i, j, k));
				sum = max(sum, (double)Ne(i, j, k));
				sum = max(sum, (double)Ne(i, j, k));
			}
	return sum;
}

template <class ftype>
double MinValueCurrents(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MAX;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = min(sum, (double)cube(i, j, k).Jx);
				sum = min(sum, (double)cube(i, j, k).Jy);
				sum = min(sum, (double)cube(i, j, k).Jz);
			}
	return sum;
}

template <class ftype>
double MaxValueElectricModule(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, abs((double)moduleE(cube, i, j, k)));
			}
	return sum;
}

template <class ftype>
double MaxValueMagnetic(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, abs((double)cube(i, j, k).Bx));
				sum = max(sum, abs((double)cube(i, j, k).By));
				sum = max(sum, abs((double)cube(i, j, k).Bz));
			}
	return sum;
}

template <class ftype>
double CalculateEnergyPML_3dimen_accurate(data3d<Component<ftype>>& cube, data3d<Component<ftype>>& compensator,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += ((double)cube(i, j, k).Ex - (double)compensator(i, j, k).Ex) * ((double)cube(i, j, k).Ex - (double)compensator(i, j, k).Ex)
					+ ((double)cube(i, j, k).Ey - (double)compensator(i, j, k).Ey) * ((double)cube(i, j, k).Ey - (double)compensator(i, j, k).Ey)
					+ ((double)cube(i, j, k).Ez - (double)compensator(i, j, k).Ez) * ((double)cube(i, j, k).Ez - (double)compensator(i, j, k).Ez);
				energy += ((double)cube(i, j, k).Bx - (double)compensator(i, j, k).Bx) * ((double)cube(i, j, k).Bx - (double)compensator(i, j, k).Bx)
					+ ((double)cube(i, j, k).By - (double)compensator(i, j, k).By) * ((double)cube(i, j, k).By - (double)compensator(i, j, k).By)
					+ ((double)cube(i, j, k).Bz - (double)compensator(i, j, k).Bz) * ((double)cube(i, j, k).Bz - (double)compensator(i, j, k).Bz);
			}
	return energy;
}

template <class ftype>
void Graph_Solution_in_two_planes_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction, string file_name)
{
	ofstream numbEx("graph_solution_E_x_float-float-sigma 49.6.csv"), numbBx("graph_solution_B_x_float-float-sigma 49.6.csv");
	ofstream numbEy("graph_solution_E_y_float-float-sigma 49.6.csv"), numbBy("graph_solution_B_y_float-float-sigma 49.6.csv");
	ofstream numbEz("graph_solution_E_z_float-float-sigma 49.6.csv"), numbBz("graph_solution_B_z_float-float-sigma 49.6.csv");

	numbEx << direction << ";" << "y" << endl;
	numbEx << "x" << ";" << ";";
	numbBx << direction << ";" << "y" << endl;
	numbBx << "x" << ";" << ";";
	numbEy << direction << ";" << "y" << endl;
	numbEy << "x" << ";" << ";";
	numbBy << direction << ";" << "y" << endl;
	numbBy << "x" << ";" << ";";
	numbEz << direction << ";" << "y" << endl;
	numbEz << "x" << ";" << ";";
	numbBz << direction << ";" << "y" << endl;
	numbBz << "x" << ";" << ";";
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbEx << dx * (ftype)i << ";";
		numbBx << dx * (ftype)i << ";";
		numbEy << dx * (ftype)i << ";";
		numbBy << dx * (ftype)i << ";";
		numbEz << dx * (ftype)i << ";";
		numbBz << dx * (ftype)i << ";";
	}
	numbEx << endl; 		numbBx << endl;	numbEy << endl; 		numbBy << endl;	numbEz << endl; 		numbBz << endl;

	for (int i = 1; i < Ny + 2 * delta + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numbEx << ";" << y << ";"; 			numbBx << ";" << y << ";";			numbEy << ";" << y << ";"; 			numbBy << ";" << y << ";";
		numbEz << ";" << y << ";"; 			numbBz << ";" << y << ";";

		for (int j = 1; j < Nx + 2 * delta + 1; j++)
		{
			numbEx << cube(j, i, Nz / 2 + delta).Ex << ";";
			numbBx << cube(j, i, Nz / 2 + delta).Bx << ";";
			numbEy << cube(j, i, Nz / 2 + delta).Ey << ";";
			numbBy << cube(j, i, Nz / 2 + delta).By << ";";
			numbEz << cube(j, i, Nz / 2 + delta).Ez << ";";
			numbBz << cube(j, i, Nz / 2 + delta).Bz << ";";
		}
		numbEx << endl; 			numbBx << endl;			numbEy << endl; 			numbBy << endl;
		numbEz << endl; 			numbBz << endl;

	}

	numbEx.close();
	numbBx.close();
	numbEy.close();
	numbBy.close();
	numbEz.close();
	numbBz.close();
}


//
template <class ftype>
void Graph_Solution_in_one_planes_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz)
{
	ofstream numbEx("E_x_1d_dd_4xgrid.csv"), numbBx("B_x_1d_dd_4xgrid.csv");
	ofstream numbEy("E_y_1d_dd_4xgrid.csv"), numbBy("B_y_1d_dd_4xgrid.csv");
	ofstream numbEz("E_z_1d_dd_4xgrid.csv"), numbBz("B_z_1d_dd_4xgrid.csv");

	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbEx << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).Ex << endl;
		numbEy << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).Ey << endl;
		numbEz << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).Ez << endl;

		numbBx << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).Bx << endl;
		numbBy << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).By << endl;
		numbBz << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz).Bz << endl;

	}

	numbEx.close();
	numbBx.close();
	numbEy.close();
	numbBy.close();
	numbEz.close();
	numbBz.close();
}

template <class ftype>
void Graph_Solution_two_planes_3d_python (data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz,
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
			numbEy << cube(j, i, 1).Ey << ";";
			numbBz << cube(j, i, 1).Bz << ";";
		}
		numbEy << endl;
		numbBz << endl;

	}
	numbEy << endl;
	numbBz << endl;

	numbEy << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << endl;
	numbEy << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << endl;
	numbEy << "dx = " << dx << ";   dy = " << dy << ";   dz = " << dz << endl;
	numbEy << "T = " << T << ";   dt = " << dt << ";   n = " << n << endl;
	numbEy << "Refl coeff = " <<R;

	numbBz << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << endl;
	numbBz << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << endl;
	numbBz << "dx = " << dx << ";   dy = " << dy << ";   dz = " << dz << endl;
	numbBz << "T = " << T << ";   dt = " << dt << ";   n = " << n << endl;
	numbBz << "Refl coeff = " << R;

	numbEy.close();
	numbBz.close();
}

