#include "chebyshevnetonsphere.h"
#include "chebyshevnet.h"

#include <iostream>
#include <cmath>

using namespace std;

ChebyshevNetOnSphere::ChebyshevNetOnSphere(double radius, int axis_size_arg)
	: ChebyshevNet(axis_size_arg), radius(radius), step_r(0.5 * radius * M_PI / axis_size)
{
}

bool ChebyshevNetOnSphere::compute(void)
{
	compute_axis_points();

	try
	{
		weave_net();
		cut_extra_points();
	}
	catch (exception &exc)
	{
		cerr << "Caught exception: " << exc.what() << endl;
		return false;
	}

	return true;
}

void ChebyshevNetOnSphere::compute_axis_points(void)
{
	for (int line = 0; line < axis_size; line++)
	{
		for (int col = 0; col < axis_size; col++)
		{
			if (line != median_index && col != median_index)
				continue;

			int offset_line_idx = line - median_index;
			int offset_col_idx = col - median_index;

			double x = radius * sin(0.5 * M_PI * offset_line_idx / axis_size);
			double y = radius * sin(0.5 * M_PI * offset_col_idx / axis_size);
			double z = radius * cos(0.5 * M_PI * offset_line_idx / axis_size);
			double U = offset_line_idx * step_r;
			double V = offset_col_idx * step_r;

			point_values_set(coordinates[line][col].first, x, y, z);
			point_values_set(coordinates[line][col].second, U, V);
		}
	}
}

void ChebyshevNetOnSphere::weave_net(void)
{
	int line, col, line_incr_sign, col_incr_sign;

	for (quarter = 1; quarter <= 4; quarter++)
	{
		tie(line, col) = get_start_indicies();
		tie(line_incr_sign, col_incr_sign) = get_incr_sign();

		for (int i = line; i < axis_size && i >= 0; i+= line_incr_sign)
			for (int j = col; j < axis_size && j >= 0; j += col_incr_sign)
			{
				Point3D &p00 = coordinates[i - line_incr_sign][j - col_incr_sign].first;
				Point3D &p10 = coordinates[i][j - col_incr_sign].first;
				Point3D &p01 = coordinates[i - line_incr_sign][j].first;
				Point3D &p_to_compute = coordinates[i][j].first;

				compute_next_point(p00, p10, p01, p_to_compute);

				double U = step_r * i * line_incr_sign;
				double V = step_r * j * col_incr_sign;

				point_values_set(coordinates[i][j].second, U, V);
			}
	}
}

//Based on paper-written formulas
void ChebyshevNetOnSphere::compute_next_point(const Point3D &p00, const Point3D &p10,
											  const Point3D &p01, Point3D &p_to_compute)
{
	double A1, B1, C1, A2, B2, C2;
	tie(A1, B1, C1) = p10;
	tie(A2, B2, C2) = p01;

	if (A1 == 0)
		return;

	double D1 = 0.5 * (pow2(radius) + pow2(A1) + pow2(B1) + pow2(C1) - pow2(step_r));
	double D2 = 0.5 * (pow2(radius) + pow2(A2) + pow2(B2) + pow2(C2) - pow2(step_r));

	double K = B2 - A2 * B1 / A1;
	double L = C2 - A2 * C1 / A1;
	double M = D2 - A2 * D1 / A1;

	if (K == 0)
		return;

	double P = M / K;
	double Q = L / K;

	double F = (D1 - B1 * P) / A1;
	double G = (B1 * Q - C1) / A1;

	double S = 1 + pow2(Q) + pow2(G);
	double T = 2 * (F * G - P * Q);
	double U = pow2(P) + pow2(F) - pow2(radius);

	if (S == 0)
		return;

	double x, y, z;
	z = (-1 * T + sqrt(pow2(T) - 4 * S * U)) / (2 * S);
	y = P - Q * z;
	x = (D1 - B1 * y - C1 * z) / A1;

	point_values_set(p_to_compute, x, y, z);

	int line_incr_sign, col_incr_sign;
	tie(line_incr_sign, col_incr_sign) = get_incr_sign();

	x *= line_incr_sign;
	y *=col_incr_sign;

	bool distance_check_failed = compute_distance_3D(p_to_compute, p00) < step_r;
	bool x_check_failed = (x < 0) || (x > radius);
	bool y_check_failed = (y < 0) || (y > radius);
	bool z_check_failed = (z < 0) || (z > radius);

	if (distance_check_failed || x_check_failed || y_check_failed || z_check_failed)
	{
		z = (-1 * T - sqrt(pow2(T) - 4 * S * U)) / (2 * S);
		y = P - Q * z;
		x = (D1 - B1 * y - C1 * z) / A1;

		point_values_set(p_to_compute, x, y, z);

		if (compute_distance_3D(p_to_compute, p00) < step_r)
			throw runtime_error("Bad computation");
	}
}

ChebyshevNetOnSphere::Cutter::Cutter(ChebyshevNetOnSphere *arg)
	: net(arg)
{}

void ChebyshevNetOnSphere::Cutter::operator ()(PointPair &point)
{
	double radius = this->net->radius;
	double x, y, z;
	tie(x, y, z) = point.first;

	if (abs(x) > radius || abs(y) > radius || z > radius || z < 0)
	{
		get<0>(point.first) = 0;
		get<1>(point.first) = 0;
		get<2>(point.first) = radius;
		get<0>(point.second) = 0;
		get<1>(point.second) = 0;
	}
}

void ChebyshevNetOnSphere::cut_extra_points(void)
{
	for_each_point(Cutter(this));
}

pair<int, int> ChebyshevNetOnSphere::get_incr_sign(void)
{
	int line_incr_sign = 0, col_incr_sign = 0;

	switch (quarter)
	{
		case 1:
			line_incr_sign = 1;
			col_incr_sign = 1;
			break;
		case 2:
			line_incr_sign = 1;
			col_incr_sign = -1;
			break;
		case 3:
			line_incr_sign = -1;
			col_incr_sign = -1;
			break;
		case 4:
			line_incr_sign = -1;
			col_incr_sign = 1;
			break;
		default:
			throw runtime_error("Bad quarter");
	}
	return make_pair(line_incr_sign, col_incr_sign);
}

pair<int, int> ChebyshevNetOnSphere::get_start_indicies(void)
{
	int line = 0, col = 0;

	switch (quarter)
	{
		case 1:
			line = median_index + 1;
			col = median_index + 1;
			break;
		case 2:
			line = median_index + 1;
			col = median_index - 1;
			break;
		case 3:
			line = median_index - 1;
			col = median_index - 1;
			break;
		case 4:
			line = median_index - 1;
			col = median_index + 1;
			break;
		default:
			throw runtime_error("Bad quarter");
	}

	return make_pair(line, col);
}
