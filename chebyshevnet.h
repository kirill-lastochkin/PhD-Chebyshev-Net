#ifndef CHEBYSHEVNET_H
#define CHEBYSHEVNET_H

#include <tuple>
#include <vector>
#include <fstream>
#include <cmath>

class ChebyshevNet
{
	protected:
		typedef std::tuple<double, double, double> Point3D;
		typedef std::tuple<double, double> Point2D;
		typedef std::pair<Point3D, Point2D> PointPair;

		std::vector<std::vector<PointPair>> coordinates;
		int axis_size;
		int median_index;
		int quarter;

		double compute_distance_3D(const Point3D &a, const Point3D &b);
		void point_values_set(Point3D &point, double x, double y, double z);
		void point_values_set(Point2D &point, double U, double V);

		template<typename T>
		void for_each_point(T&& functor)
		{
			for (auto &line : coordinates)
				for (auto &point : line)
					functor(point);
		}

	public:
		ChebyshevNet() = delete;
		ChebyshevNet(int axis_size);
		virtual ~ChebyshevNet() = default;

		void print(void);
		void save_to_file(const char * filename);

		virtual bool compute(void) = 0;

	private:
		struct Printer2D
		{
			void operator ()(const PointPair &point);
		};

		struct Printer3D
		{
			void operator ()(const PointPair &point);
		};

		struct FileSaver
		{
			std::ofstream file_out;
			void operator ()(const PointPair &point);

			FileSaver() = delete;
			FileSaver(const char *filename);
			~FileSaver();
		};
};

template <typename T>
inline T pow2(T arg)
{
	return pow(arg, 2);
}

#endif // CHEBYSHEVNET_H
