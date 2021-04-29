#include "fftconvolver/FFTConvolver.h"
#include <boost/circular_buffer.hpp>

#include <random>
using namespace std;
using namespace fftconvolver;

class RandomLFO {
public:
	RandomLFO();

	// scale - Characteristic time-scale of the function.
	void setScale(float scale);

	// Smoothness - Smoothness of the function. Corresponds to v-parameter of the Matern Covariance kernel. For smoothness = 0,1,2 : v=smoothness+0.5. For higher smoothness, v=inf(gaussian covariance function).
	void setSmoothness(float smoothness);

	//
	// sample_rate - sample rate of audio. generateSamples will generate lfo-samples at the same rate.
	// scale - Characteristic time-scale of the Gaussian Function.
	// Smoothness - Smoothness of the function.
	// n_scale - Number of lengthscales(scale) to be included in the FIR. 
	// n_samples_per_scale - Number of samples per lengthscale in the FIR.
	void init(float sample_rate,float scale, float smoothness = 3.0, unsigned int n_samples_per_scale = 100, float n_scales = 7);

	//
	// out - buffer of size N to be filled
	// N - Number of samples to fill in "out".
	// output - Can set output = false if you want to progress the state of the RandomLFO without getting any samples.
	float generateSamples(float *out, unsigned int N,bool output=true);

	//
	// seed - Set seed of random number generator.
	void seed(unsigned int seed);


protected:
	float scale;
	float sample_rate;
	boost::circular_buffer<float> output_buf;
	double ds;
	static const unsigned int n_samples_convolver = 256;
private:
	void generateFIRs();
	void calcFIR();
	void sampleNext();
	void updateConvolver();

	
	float smoothness;
	unsigned int n_samples_per_scale;
	float n_scales;
	float nextSample;
	float lastSample;
	
	double current_time;
	double nextSampleTime;
	

	
	static const unsigned int block_size_convolver = 2048;
	float conv_buf_in[n_samples_convolver];
	float conv_buf_out[n_samples_convolver];

	FFTConvolver convolver;
	vector<vector<float>> firs;
	vector<float> current_fir;
	default_random_engine generator;
	normal_distribution<float> distribution;
	boost::circular_buffer<float> rnd_buf;
	

	bool is_init;


};

class RandomLFOSingle : public RandomLFO {
public:
	RandomLFOSingle();

	// Set output frequency of the RandomLFO.
	void setOutFreq(float freq);

	// Feed this function with samples of sample_frequency( See RandomLFO.init(...) ).
	// Will return True and overwrite "out" at the frequency specified in setOutFreq. Else returns False.
	bool processSingle(float *out);

private:
	double time;
	double out_dt;
};

