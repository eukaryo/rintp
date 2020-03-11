/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTP_REAL_LOGSUMEXP_H_
#define RINTP_REAL_LOGSUMEXP_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <functional>
#include <complex>
#include <random>
#include <chrono>
#include <sstream>

#include"interval_type.h"

namespace rintp {

template<typename RealScalar>class WideRealNumber {

private:

	static IntervalVar add_(const IntervalVar&, const IntervalVar&);
	static Floating add_(const Floating&, const Floating&);

	static IntervalVar sub_(const IntervalVar&, const IntervalVar&);
	static Floating sub_(const Floating&, const Floating&);

public:

	RealScalar log_scale;
	WideRealNumber() :log_scale(RealScalar(0.0)) {}
	WideRealNumber(const double number) {//�����̓���������Interval��Floating�ŕ�����͖̂ʓ|�������B
		assert(number >= 0.0);
		this->log_scale = log(RealScalar(number));
	}

	WideRealNumber& operator += (const WideRealNumber& obj) {
		this->log_scale = add_(this->log_scale, obj.log_scale);
		return *this;
	}
	WideRealNumber& operator -= (const WideRealNumber& obj) {
		
		assert(this->log_scale >= obj.log_scale);

		if (obj.log_scale == -std::numeric_limits<RealScalar>::infinity())return *this;

		const RealScalar big = this->log_scale;
		const RealScalar small = obj.log_scale;
		this->log_scale = big + log(RealScalar(1.0) - exp(small - big));
		return *this;
	}
	WideRealNumber& operator *= (const WideRealNumber& obj) {
		this->log_scale += obj.log_scale;
		return *this;
	}
	WideRealNumber& operator /= (const WideRealNumber& obj) {
		this->log_scale -= obj.log_scale;
		return *this;
	}
	WideRealNumber operator + (const WideRealNumber& obj)const { WideRealNumber re(*this); return re += obj; }
	WideRealNumber operator - (const WideRealNumber& obj)const { WideRealNumber re(*this); return re -= obj; }
	WideRealNumber operator * (const WideRealNumber& obj)const { WideRealNumber re(*this); return re *= obj; }
	WideRealNumber operator / (const WideRealNumber& obj)const { WideRealNumber re(*this); return re /= obj; }

	bool operator > (const WideRealNumber& obj)const {
		return this->log_scale > obj.log_scale;
	}
	bool operator >= (const WideRealNumber& obj)const {
		return this->log_scale >= obj.log_scale;
	}
	bool operator < (const WideRealNumber& obj)const {
		return this->log_scale < obj.log_scale;
	}
	bool operator <= (const WideRealNumber& obj)const {
		return this->log_scale <= obj.log_scale;
	}
	bool operator == (const WideRealNumber& obj)const {
		return this->log_scale == obj.log_scale;
	}
	RealScalar ToUsualReal()const {
		return exp(this->log_scale);
	}

	static WideRealNumber LogRealToWide(const RealScalar& x) {
		WideRealNumber ans;
		ans.log_scale = RealScalar(x);
		return ans;
	}
	static double Probability(const WideRealNumber& target, const WideRealNumber& total) {

		//return double(target/total)

		if (target.log_scale == -std::numeric_limits<RealScalar>::infinity())return 0.0;
		if (total.log_scale == -std::numeric_limits<RealScalar>::infinity())return 0.0;

		return exp(target.log_scale - total.log_scale);
	}
};

}

#endif//RINTP_REAL_LOGSUMEXP_H_
