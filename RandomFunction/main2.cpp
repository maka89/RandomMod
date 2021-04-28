#include "RandomLFO.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>    
using namespace std::chrono;


using namespace std;
int main() {


	RandomLFOSingle rlfo;
	float sr = 44100.0f;
	float out_freq = 1000.0;

	rlfo.init(sr, 1.0);

	rlfo.setOutFreq(out_freq);
	int N = 64 * 10000;
	vector<float> samples = vector<float>(0);


	float scale = 1.0f;
	float smoothness = 2.0f;
	rlfo.seed(0);
	rlfo.setScale(scale);
	rlfo.setSmoothness(smoothness);


	ofstream myfile, plot;
	plot.open("plot_s2.csv");
	plot << "t;value\n" << endl;

	auto start = high_resolution_clock::now();

	float last_sample = 0.0f;
	for (int i = 0; i < N; i++) {
		float tmp;
		if (rlfo.processSingle(&tmp))
			samples.push_back(tmp);
	}
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);


	cout << (sr*((double)duration.count() / 1000000.0) / N) << endl;
	plot << std::setprecision(9);
	for (int i = 0; i < samples.size(); i++)
		plot << (float)i / out_freq << ";" << samples[i] << endl;
	plot.close();


	cout << "DONE" << endl;

}