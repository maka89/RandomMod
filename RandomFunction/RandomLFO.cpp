#include "RandomLFO.h"
#include "MaternKernelFunction.h"

//Debug
#include <iostream>
#include <fstream>
#include <string>
#include <ctime> 
//Debug

void vec_out(vector<float> vec, string filename) {
	ofstream myfile;
	myfile.open(filename);
	myfile << "value\n";
	for (int i = 0; i < vec.size(); i++) {
		myfile << vec[i] << endl;
	}
	myfile.close();

}
RandomLFO::RandomLFO():distribution(0.0f, 1.0f) {
}

void RandomLFO::init(float sample_rate, float scale, float smoothness, unsigned int n_samples_per_scale, float n_scales) {
	// Set variables

	

	if (n_samples_per_scale < 1)
		throw "ERROR: n_samples_per_scale must be >=1";
	this->n_samples_per_scale = n_samples_per_scale;

	if (n_scales <= 0)
		throw "ERROR: n_scales must be > 0";
	this->n_scales = n_scales;

	this->sample_rate = sample_rate;
	this->generateFIRs();

	unsigned int N = (int)(this->n_scales*this->n_samples_per_scale);
	unsigned int num = 2 * N + 1;
	this->current_fir = vector<float>(num);
	
	this->ds = 2.0*n_scales / num;

	this->rnd_buf = boost::circular_buffer<float>(num);
	for (unsigned int i = 0; i < num; i++)
		this->rnd_buf.push_back(0.0);

	this->output_buf = boost::circular_buffer<float>(num);
	for (unsigned int i = 0; i < num; i++)
		this->output_buf.push_back(0.0);


	this->current_time = 0.0;
	this->nextSampleTime = 0.0;

	this->setScale(scale);
	this->setSmoothness(smoothness);
	this->generator.seed((unsigned int)time(NULL));
	
	vec_out(this->current_fir, "fir.csv");


}

float RandomLFO::generateSamples(float *out, unsigned int N) {

	double dt = (1.0 / this->sample_rate)/this->scale;

	size_t buflen = this->output_buf.size();

	unsigned int offset = 0;

	float last_sample;
	float this_sample;
	float ret = 0.0f;
	for (unsigned int i = 0; i < N; i++) {
		this->current_time += dt;
		if (this->current_time > this->nextSampleTime) {
			this->sampleNext();
			this->current_time -= this->nextSampleTime;
			this->nextSampleTime = this->ds*this->n_samples_convolver;

			offset = this->n_samples_convolver - (unsigned int)(this->current_time / this->ds);
			//cout << offset << endl;
			last_sample = this->output_buf[buflen - 2 - offset];
			this_sample = this->output_buf[buflen - 1 - offset];
		}

		offset = this->n_samples_convolver - (unsigned int)(this->current_time / this->ds);
		last_sample = this->output_buf[buflen - 2 - offset];
		this_sample = this->output_buf[buflen - 1 - offset];
		

		//Linear interpolation between sample_points.
		float intersample_time = (this->current_time / this->ds) - (unsigned int)(this->current_time / this->ds);
		out[i] = (1.0- intersample_time)*last_sample + this_sample*intersample_time;

	}
	return ret;

}


void RandomLFO::sampleNext() {
	
	for (int i = 0; i < this->n_samples_convolver; i++) {
		float eps = this->distribution(this->generator);
		this->rnd_buf.push_back(eps);
		this->conv_buf_in[i] = eps;
	}

	this->convolver.process(this->conv_buf_in, this->conv_buf_out, this->n_samples_convolver);

	for(int i = 0; i< this->n_samples_convolver;i++)
		this->output_buf.push_back(this->conv_buf_out[i]);
}

void RandomLFO::calcFIR() {

	if (this->smoothness >= this->firs.size()-1) {
		for (unsigned int i = 0; i < this->current_fir.size(); i++)
			this->current_fir[i] = this->firs[this->firs.size() - 1][i];
	}
	else if (this->smoothness <= 0) {
		for (unsigned int i = 0; i < this->current_fir.size(); i++)
			this->current_fir[i] = this->firs[0][i];
	}
	else {
		int low = (int)this->smoothness;
		int high = low+1;

		float w_low = (high - this->smoothness) / (high - low);
		float w_high = (this->smoothness - low) / (high - low);

		for (unsigned int i = 0; i < this->current_fir.size(); i++) 
			this->current_fir[i] = w_low * this->firs[low][i] + w_high * this->firs[high][i];
	}

	double sum1 = 0.0;
	for (unsigned int i = 0; i < this->current_fir.size(); i++)
		sum1 += this->current_fir[i]* this->current_fir[i];
	sum1 = sqrt(sum1);
	for (unsigned int i = 0; i < this->current_fir.size(); i++)
		this->current_fir[i] /= (float)sum1;

}


void RandomLFO::generateFIRs() {
	unsigned int N = (int)(n_scales*n_samples_per_scale);
	unsigned int num = 2 * N + 1;

	float start = -n_scales;
	float end = n_scales;
	
	float delta = (end - start) / (num - 1);
	vector<float> d(num);
	for (unsigned int i = 0; i < 2 * N + 1; i++)
		d[i] = start + i * delta;

	for (int i =0;i<4;i++)
		this->firs.push_back(MaternKernelFunction::kernel_fn(d, i));

	
}

void RandomLFO::updateConvolver() {
	this->convolver.reset();
	this->convolver.init(this->block_size_convolver, &this->current_fir[0], this->current_fir.size());

	float *outbuf = new float[this->sample_rate];
	
	unsigned int N = this->scale*this->n_scales * 2;

	for(unsigned int i = 0; i< N;i++)
		this->generateSamples(outbuf,this->sample_rate);
	delete outbuf;

}

void RandomLFO::setScale(float scale) {
	if (scale <= 0)
		throw "ERROR: scale must be > 0";
	this->scale = scale;
}
void RandomLFO::setSmoothness(float smoothness) {
	if (smoothness < 0.0)
		throw "ERROR: smoothness must be >= 0";
	this->smoothness = smoothness;
	this->calcFIR();
	this->updateConvolver();
}
