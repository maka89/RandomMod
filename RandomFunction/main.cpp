#include "RandomLFO.h"
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std::chrono;


using namespace std;
int main() {


	RandomLFO rlfo;
	float sr = 44100.0f;
	rlfo.init(sr, 1.0);
	int N = 64*10000;
	float *samples = new float[N];
	for (int i = 0; i < N; i++)
		samples[i] = 0.0;

	
	float scale = 1.0f;
	float smoothness = 3.0f;
	rlfo.setScale(scale);
	rlfo.setSmoothness(smoothness);


	ofstream myfile,plot;
	plot.open("plot.csv");
	plot << "t;value\n" << endl;

	auto start = high_resolution_clock::now();

	for (int i = 0; i < 10000; i++) 
		auto val = rlfo.generateSamples(&samples[i * 64], 64);
	
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);


	cout << (sr*((double)duration.count()/1000000.0)/N) << endl;
	for (int i = 0; i < 64 * 10000; i++)
		plot << (float)i/sr << ";" << samples[i] << endl;
	plot.close();
	// To get the value of duration use the count()
	// member function on the duration object
	
	
	
	
	



	delete samples;

	cout << "DONE" << endl;

}