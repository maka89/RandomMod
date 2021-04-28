#pragma once
#include "RandomLFO.h"
#include <dspfilters/so_lpf.h>

const float VeryAutoFilter::min_scale = 0.1f;
const float VeryAutoFilter::fc_max = 22000.0f;

class VeryAutoFilter {
public:
	VeryAutoFilter(float audio_sample_rate, float lfo_out_sample_rate):audio_sample_rate(audio_sample_rate),lfo_out_sample_rate(lfo_out_sample_rate),fc(20000.0),Q(0.707) {
		this->rlfo_fc.init(this->audio_sample_rate, 1.0);
		this->rlfo_fc.setSmoothness(3.0);
		rlfo_fc.setOutFreq(this->lfo_out_sample_rate);
	
	}
	void setNFilters(unsigned int n_filters ) {
		filters.clear();
		for (unsigned int i = 0; i < n_filters; i++)
			filters.push_back(SO_LPF());
		for (unsigned int i = 0; i < n_filters; i++)
			filters[i].calculate_coeffs(this->Q, this->fc, this->audio_sample_rate);
	}
	void setSmoothness(float smoothness) { this->rlfo_fc.setSmoothness(smoothness); }

	void setScale(float scale) { 
		if (scale >= min_scale)
			this->rlfo_fc.setScale(scale);
		else
			this->rlfo_fc.setScale(min_scale);
	}

	void setFcAmount(float amount) { this->fc_amount = amount; }
	void setFc(float fc) { this->fc_mean = fc; }
	void setQ(float Q) { this->Q_mean = Q; }


	void process(double *out, const double *in, unsigned int N) {


		float fc_val;
		for (unsigned int i = 0; i < N; i++) {
			if (rlfo_fc.processSingle(&fc_val))
				this->updateFilters(calcFc(this->fc_mean, fc_val*this->fc_amount),this->Q);

			out[i] = in[i];
			for (unsigned int j = 0; j < this->filters.size(); j++) 
				out[i] = this->filters[j].process(out[i]);

		}
	}



private:

	float calcFc(float fc_mean, float lfo_val) {return pow(2.0, log2(fc_mean) + lfo_val);}
	void updateFilters(float fc, float Q) {
		for (unsigned int i = 0; i < filters.size(); i++) {
			filters[i].calculate_coeffs(Q, fc, this->audio_sample_rate);
		}
	}

	RandomLFOSingle rlfo_res,rlfo_fc;
	vector<SO_LPF> filters;
	float fc, Q, audio_sample_rate, lfo_out_sample_rate;
	float fc_amount,res_amount;

	float fc_mean, Q_mean;

	static const float min_scale, fc_max;

};

