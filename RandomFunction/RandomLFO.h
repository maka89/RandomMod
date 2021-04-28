#include "fftconvolver/FFTConvolver.h"
#include <boost/circular_buffer.hpp>

#include <random>
using namespace std;
using namespace fftconvolver;

class RandomLFO {
public:
	RandomLFO();
	void setScale(float);
	void setSmoothness(float);
	void init(float sample_rate,float scale, float smoothness = 3.0, unsigned int n_samples_per_scale = 100, float n_scales = 7);
	float generateSamples(float *out, unsigned int N,bool output=true);
	void seed(unsigned int);


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
	RandomLFOSingle() :RandomLFO(), time(0.0), out_dt(1.0/1000.0) {}
	void setOutFreq(float freq) { this->out_dt = 1.0/freq; }
	bool processSingle(float *out) {


		bool retval = false;
		if (this->time >= out_dt) {
			this->generateSamples(out, 1, true);
			this->time = 0.0;
			retval = true;
		}
		else 
			this->generateSamples(NULL, 1, false);
		
		time += (1.0 / this->sample_rate) / this->scale ;
		return retval;

		
	}

private:
	double time;
	double out_dt;
};
