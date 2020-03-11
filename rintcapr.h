/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef RINTP_RINTCAPR_H_
#define RINTP_RINTCAPR_H_

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


#include"complex_number.h"
#include"interval_type.h"




namespace rintp {

template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintP1Dim(
	const std::string sequence,
	const std::vector<std::vector<int>> S,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop,
	const bool allow_fft);

std::pair<std::vector<Floating>, std::vector<std::vector<std::vector<Floating>>>>RegularizeRintP1Dim(
	const std::vector<WideComplexNumber<Floating>>& z,
	const std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>& profile);

std::pair<std::vector<IntervalVar>, std::vector<std::vector<std::vector<IntervalVar>>>>RegularizeRintP1Dim(
	const std::vector<WideComplexNumber<IntervalVar>>& z,
	const std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>& profile);

std::pair<std::vector<double>, std::vector<std::vector<std::vector<double>>>> BruteForceRintP1Dim(
	const std::string sequence,
	const std::string reference_structure,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop);


}


#endif//RINTP_RINTCAPR_H_