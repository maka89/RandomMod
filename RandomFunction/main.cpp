#include "RandomLFO.h"
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std::chrono;


using namespace std;
int main() {


	RandomLFO *rlfo = new RandomLFO();
	rlfo->init(44100.0f, 1.0);
	int N = 64*10000;
	float *samples = new float[N];
	for (int i = 0; i < N; i++)
		samples[i] = 0.0;

	
	float scale = 1.0f;// +1.0*((float)i) / 99.0;
	float smoothness = 3.0f;// 4.0 - 4.0*((float)i) / 99.0;
	rlfo->setScale(scale);
	rlfo->setSmoothness(smoothness);


	ofstream myfile,plot;
	plot.open("plot.csv");
	myfile.open("timing.csv");
	myfile << "value,value2\n";
	plot << "value\n" << endl;
	for (int i = 0; i < 10000; i++) {
		auto start = high_resolution_clock::now();
		auto val = rlfo->generateSamples(&samples[i * 64], 64);
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);
		myfile << duration.count() / 1000000.0 << ", " << val <<  endl;
		
	}
	
	for (int i = 0; i < 64 * 10000; i++)
		plot << samples[i] << endl;
	myfile.close();

	// To get the value of duration use the count()
	// member function on the duration object
	
	
	
	
	



	delete samples;
	delete rlfo;

	cout << "DONE" << endl;

}