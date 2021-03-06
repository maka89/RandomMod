#include "RandomLFO.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>    
using namespace std::chrono;


using namespace std;
void main_not() {


	RandomLFO rlfo;
	float sr = 44100.0f;
	rlfo.init(sr, 1.0);
	int N = 64*10000;
	float *samples = new float[N];
	for (int i = 0; i < N; i++)
		samples[i] = 0.0;

	
	float scale = 1.0f;
	float smoothness = 2.0f;
	rlfo.seed(0);
	rlfo.setScale(scale);
	rlfo.setSmoothness(smoothness);


	ofstream myfile,plot;
	plot.open("plot_s2.csv");
	plot << "t;value\n" << endl;

	auto start = high_resolution_clock::now();

	for (int i = 0; i < 10000; i++) 
		auto val = rlfo.generateSamples(&samples[i * 64], 64);
	
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);


	cout << (sr*((double)duration.count()/1000000.0)/N) << endl;
	plot << std::setprecision(9);
	for (int i = 0; i < 64 * 10000; i++)
		plot << (float)i/sr << ";" << samples[i] << endl;
	plot.close();
	
	delete samples;

	cout << "DONE" << endl;

}