/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"real_logsumexp.h"

namespace rintp {

//explicit instantiation

template IntervalVar WideRealNumber<IntervalVar>::add_(const IntervalVar&, const IntervalVar&);
template Floating WideRealNumber<Floating>::add_(const Floating&, const Floating&);

template IntervalVar WideRealNumber<IntervalVar>::sub_(const IntervalVar&, const IntervalVar&);
template Floating WideRealNumber<Floating>::sub_(const Floating&, const Floating&);

//function definition

template<typename RealScalar> IntervalVar WideRealNumber<RealScalar>::add_(const IntervalVar& x, const IntervalVar& y) {

	if (x.upper() == -std::numeric_limits<Floating>::infinity() || y.upper() == -std::numeric_limits<Floating>::infinity()) {
		return IntervalVar(-std::numeric_limits<Floating>::infinity(), -std::numeric_limits<Floating>::infinity());
	}
	if (x.lower() == std::numeric_limits<Floating>::infinity() || y.lower() == std::numeric_limits<Floating>::infinity()) {
		return IntervalVar(std::numeric_limits<Floating>::infinity(), std::numeric_limits<Floating>::infinity());
	}
	if (x.upper() == std::numeric_limits<Floating>::infinity() ||
		y.upper() == std::numeric_limits<Floating>::infinity() ||
		x.lower() == -std::numeric_limits<Floating>::infinity() ||
		y.lower() == -std::numeric_limits<Floating>::infinity()) {
		return x + y;
	}

	const IntervalVar p(x.upper() - mid(x) + y.upper() - mid(y));
	if (p.upper() == std::numeric_limits<Floating>::infinity()) {
		return IntervalVar(-std::numeric_limits<Floating>::infinity(), std::numeric_limits<Floating>::infinity());
	}
	return x + p + log(exp(-p) + exp(y - x - p));

	//   log(exp(x)+exp(y))
	// = log(exp(x)*(1+exp(y-x)))
	// = log(exp(x)*exp(p)*(exp(-p)+exp(y-x-p)))
	// = x+p+log(exp(-p)+exp(y-x-p))

	//こうすることで、yとxがほぼ同じ値でかつ区間幅が広いときにオーバーフローを防げる。

}
template<typename RealScalar> Floating WideRealNumber<RealScalar>::add_(const Floating& x, const Floating& y) {

	if (x == -std::numeric_limits<Floating>::infinity()) {
		return y;
	}
	if (y == -std::numeric_limits<Floating>::infinity()) {
		return x;
	}

	Floating big, small;
	if (x > y) {
		big = x;
		small = y;
	}
	else {
		big = y;
		small = x;
	}
	return big + log(Floating(1.0) + exp(small - big));
}

template<typename RealScalar> IntervalVar WideRealNumber<RealScalar>::sub_(const IntervalVar& x, const IntervalVar& y) {
	//精度保証演算＋区間型＋logsumexpで引き算をするとき、ナイーブにやると負の値にはみ出してエラーが出うる。
	//ここでは、「真の値どうしでは引き算が成立する」という仮定のもとでエラーが出ないようにする。

	assert(x.upper() >= y.lower());

	const IntervalVar p(x.upper() - mid(x) + y.upper() - mid(y));
	if (p.upper() == std::numeric_limits<Floating>::infinity()) {
		return IntervalVar(-std::numeric_limits<Floating>::infinity(), std::numeric_limits<Floating>::infinity());
	}
	IntervalVar t = exp(-p) - exp(y - x - p);
	t = /*kv::*/intersect(IntervalVar(0.0, std::numeric_limits<Floating>::infinity()), t);
	//t = /*kv::*/max(IntervalVar(0.0, 0.0), t);
	return x + p + log(t);

	//   log(exp(x)-exp(y))
	// = log(exp(x)*(1-exp(y-x)))
	// = log(exp(x)*exp(p)*(exp(-p)-exp(y-x-p)))
	// = x+p+log(exp(-p)-exp(y-x-p))

}
template<typename RealScalar> Floating WideRealNumber<RealScalar>::sub_(const Floating& x, const Floating& y) {
	assert(x >= y);

	if (y == -std::numeric_limits<RealScalar>::infinity())return x;

	const RealScalar big = x;
	const RealScalar small = y;
	return big + log(RealScalar(1.0) - exp(small - big));
}


}

