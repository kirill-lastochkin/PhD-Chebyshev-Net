#ifndef CHEBYSHEVNETONSPHERE_H
#define CHEBYSHEVNETONSPHERE_H

#include "chebyshevnet.h"
#include "new_net_global.h"

class NEW_NETSHARED_EXPORT ChebyshevNetOnSphere : public ChebyshevNet
{
	friend struct Cutter;

	private:
		double radius;
		double step_r;

		void compute_axis_points(void);
		void weave_net(void);
		std::pair<int, int> get_start_indicies(void);
		std::pair<int, int> get_incr_sign(void);
		void compute_next_point(const Point3D &p00, const Point3D &p10,
								const Point3D &p01, Point3D &p_to_compute);
		void cut_extra_points(void);

		struct Cutter
		{
			ChebyshevNetOnSphere *net;

			void operator ()(PointPair &point);
			Cutter() = delete;
			Cutter(ChebyshevNetOnSphere *arg);
		};

	public:
		ChebyshevNetOnSphere() = delete;
		ChebyshevNetOnSphere(double radius, int axis_size);
		~ChebyshevNetOnSphere() = default;

		bool compute(void) override;
};

#endif // CHEBYSHEVNETONSPHERE_H
