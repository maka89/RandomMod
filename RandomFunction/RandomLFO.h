#include "FFTConvolver.h"
#include <boost/circular_buffer.hpp>
#include <random>
using namespace fftconvolver;
using namespace std;
class RandomLFO {
public:
	RandomLFO();
	void setScale(float);
	void setSmoothness(float);
	void init(float sample_rate,float scale, float smoothness = 3.0, unsigned int n_samples_per_scale = 100, float n_scales = 7);
	float generateSamples(float *out, unsigned int N);
private:
	void generateFIRs();
	void calcFIR();
	void sampleNext();
	void updateConvolver();

	float scale;
	float smoothness;
	unsigned int n_samples_per_scale;
	float n_scales;
	float nextSample;
	float lastSample;
	float sample_rate;
	double current_time;
	double nextSampleTime;
	double ds;

	static const unsigned int n_samples_convolver = 128;
	static const unsigned int block_size_convolver = 2048;
	float conv_buf_in[n_samples_convolver];
	float conv_buf_out[n_samples_convolver];

	FFTConvolver convolver;
	vector<vector<float>> firs;
	vector<float> current_fir;
	default_random_engine generator;
	normal_distribution<float> distribution;
	boost::circular_buffer<float> rnd_buf;
	boost::circular_buffer<float> output_buf;


};

