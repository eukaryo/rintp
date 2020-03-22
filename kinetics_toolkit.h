/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef RINTP_KINETICS_TOOLKIT_H_
#define RINTP_KINETICS_TOOLKIT_H_

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

std::vector<std::pair<std::string, double>>ComputeLocalMFE(
	const std::string& sequence,
	const std::string& initial_structure,
	const double temperature,
	const int max_span,
	const int max_loop);

std::vector<std::pair<std::string, double>> EnumerateStructureAndBoltzmannFactor(
	const std::string& sequence,
	const double temperature,
	const int max_span,
	const int max_loop);

std::vector<std::pair<std::string, double>> SimulateGillespie(
	const std::string& sequence,
	const std::string& initial_structure,
	const double temperature,
	const int max_span,
	const int max_loop,
	const int step,
	const uint64_t random_seed);
}

#endif//RINTP_KINETICS_TOOLKIT_H_

