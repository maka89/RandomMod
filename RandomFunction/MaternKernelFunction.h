#include <vector>
#include <cmath>

class MaternKernelFunction {
public:
	static vector<float> kernel_fn(const std::vector<float> dist, unsigned int smoothness, float scale = 1.0f);
	static inline float kernel_scalar(unsigned int smoothness, double d) {
		double v = smoothness + 0.5;
		float retval = 0.0f;
		switch (smoothness) {
		case 0:
			retval =(float)exp(-d);
			break;
		case 1:
			retval = (float)((1.0 + sqrt(3.0)*d)*exp(-sqrt(3.0)*d));
			break;
		case 2:
			retval = (float)((1.0 + sqrt(5.0)*d + (5.0 / 3.0)*d*d)*exp(-sqrt(5.0)*d));
			break;
		default:
			retval = (float)exp(-0.5*d*d);
		}
		return retval;
	}
};

// dist = x1-x2
// smoothness = matern v parameter
// scale = correlation length
std::vector<float> MaternKernelFunction::kernel_fn(const std::vector<float> dist, unsigned int smoothness, float scale)
{

	vector<float> out(dist.size());
	double s = (double)smoothness;
	double fac = pow(2.0, 1 - s) / tgamma(s);
	for (unsigned int i = 0; i < dist.size(); i++) {
		double d = (double)abs(dist[i]);
		out[i] = kernel_scalar(smoothness,d);
	}
	return out;

}