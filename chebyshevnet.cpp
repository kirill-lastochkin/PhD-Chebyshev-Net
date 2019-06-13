#include "chebyshevnet.h"

#include <iomanip>
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

ChebyshevNet::ChebyshevNet(int axis_size_arg)
	: axis_size(axis_size_arg), quarter(0)
{
	assert(axis_size_arg > 3);

	if ((axis_size_arg & 0x1) == 0)
		axis_size++;     // odd value required

	coordinates.assign(axis_size, vector<PointPair>(0));

	for (auto &line : coordinates)
		line.assign(axis_size, PointPair(Point3D{0.0, 0.0, 0.0}, Point2D{0.0, 0.0}));

	median_index = axis_size / 2;
}

void ChebyshevNet::point_values_set(Point3D &point, double x, double y, double z)
{
	get<0>(point) = x;
	get<1>(point) = y;
	get<2>(point) = z;
}

void ChebyshevNet::point_values_set(Point2D &point, double U, double V)
{
	get<0>(point) = U;
	get<1>(point) = V;
}

double ChebyshevNet::compute_distance_3D(const Point3D &a, const Point3D &b)
{
	double x1, y1, z1, x2, y2, z2;
	tie(x1, y1, z1) = a;
	tie(x2, y2, z2) = b;
	return sqrt(pow2(x1 - x2) + pow2(y1 - y2) + pow2(z1 - z2));
}

void ChebyshevNet::print(void)
{
	cout.precision(2);
	cout << "Output 3D coordinates\n";
	for_each_point(Printer3D());

	cout << "Output 2D coordinates\n";
	for_each_point(Printer2D());
}

void ChebyshevNet::Printer2D::operator ()(const PointPair &point)
{
	double U, V;
	tie(U, V) = point.second;
	cout << "(" << setw(5) << U << " " << setw(5) << V << ")\n";
}

void ChebyshevNet::Printer3D::operator ()(const PointPair &point)
{
	double x, y, z;
	tie(x, y, z) = point.first;
	cout << "(" << setw(5) << x << " " << setw(5) << y << " " << setw(5) << z << ")\n";
}

void ChebyshevNet::save_to_file(const char *filename)
{
	for_each_point(FileSaver{filename});
}

ChebyshevNet::FileSaver::FileSaver(const char *filename)
	: file_out(filename)
{
}

void ChebyshevNet::FileSaver::operator ()(const PointPair &point)
{
	double x, y, z, U, V;
	tie(x, y, z) = point.first;
	tie(U, V) = point.second;

	file_out << x << " ";
	file_out << y << " ";
	file_out << z << " ";
	file_out << U << " ";
	file_out << V << "\n";
}

ChebyshevNet::FileSaver::~FileSaver()
{
	file_out.close();
}
