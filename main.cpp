﻿/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

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
#include <thread>

#include <omp.h>

#include"complex_number.h"
#include"fourier_transform.h"
#include"misc.h"
#include"test.h"
#include"rintx.h"
#include"mccaskill_1990.h"
#include"sample_mccaskill.h"
#include"experiment.h"
#include"parameter.h"
#include"centroid_fold.h"
#include"interval_type.h"
#include"rintcapr.h"


namespace rintp {

void TestAll() {

	const int num = 100;

	TestMaxHamming(num);
	TestSampling(num);
	TestSimpleMcCaskillWide(num);
	TestRintD1Dim(num);
	TestRintD2Dim(num);
	TestRintW1Dim(num);
	TestRintW2Dim(num);
	TestCentroidFold(num);
	TestHagioNonFourier(num);
	TestMaxHammingPK(num);
	TestRintD1DimPK(num);

	TestRintP(num);
}



int verification(const std::string& structure, const std::string& sequence, const int max_span, const int max_loop) {

	for (const char c : sequence) {
		if (!(c == 'A' || c == 'U' || c == 'G' || c == 'C')) {
			std::cerr << "Error: The RNA sequence must be uppercase." << std::endl;
			return 1;
		}
	}
	for (const char c : structure) {
		if (!(c == '(' || c == '.' || c == ')')) {
			std::cerr << "Error: The RNA structure must consist of '(', '.' and ')' only." << std::endl;
			return 1;
		}
	}
	if (!(structure.size() == sequence.size())) {
		std::cerr << "Error: The RNA sequence and structure must be the same length." << std::endl;
		std::cerr << "Your sequence length  = " + std::to_string(sequence.size()) << std::endl;
		std::cerr << "Your structure length = " + std::to_string(structure.size()) << std::endl;
		return 1;
	}
	if (!(1 <= max_span)) {
		std::cerr << "Error: The value of max-span constraint must be a positive integer." << std::endl;
		return 1;
	}
	if (!(max_span <= int(structure.size()))) {
		std::cerr << "Error: The value of max-span constraint must be less than or equal to the length of the RNA sequence." << std::endl;
		std::cerr << "Your sequence length  = " + std::to_string(sequence.size()) << std::endl;
		std::cerr << "Your constraint value = " + std::to_string(max_span) << std::endl;
		return 1;
	}

	{
		const int n = int(sequence.size());
		const std::string bp = "AU UA GC CG GU UG";
		std::string query = "XX";
		std::vector<std::vector<int>>ans(n + 1, std::vector<int>(n + 1, 0));
		std::stack<int> bp_pos;
		for (int i = 1; i <= n; ++i) {
			switch (structure[i - 1]) {
			case '(':
				bp_pos.push(i);
				break;
			case ')':
				if (!(bp_pos.size() >= 1)) {
					std::cerr << "Error: The RNA structure is invalid." << std::endl;
					std::cerr << "')' of position " + std::to_string(i) + " cannot form a base pair." << std::endl;
					return 1;
				}
				if (!(TURN < (i - bp_pos.top()) && (i - bp_pos.top()) <= max_span)) {
					std::cerr << "Error: The RNA structure contains the base pair whose length is longer than the max-span constraint." << std::endl;
					return 1;
				}

				query[0] = sequence[bp_pos.top() - 1];
				query[1] = sequence[i - 1];
				if (!(bp.find(query) != std::string::npos)) {
					std::cerr << "Error: The RNA structure contains an illegal base pair." << std::endl;
					std::cerr << "Position " + std::to_string(bp_pos.top()) + " (the base is '" + query.substr(0, 1) + "') and" << std::endl;
					std::cerr << "position " + std::to_string(i) + " (the base is '" + query.substr(1, 1) + "') " << std::endl;
					return 1;
				}


				ans[bp_pos.top()][i] = 1;
				bp_pos.pop();
				break;
			case '.':
				break;
			default:
				assert(0);
				break;
			}
		}
		if (!(bp_pos.size() == 0)) {
			std::cerr << "Error: The RNA structure is invalid." << std::endl;
			std::cerr << "'(' is more than ')'." << std::endl;
			return 1;
		}
		if (!(ComputeMaxLoop(structure) <= max_loop)) {
			std::cerr << "Error: The RNA structure contains the base pair whose length is longer than the max-loop constraint." << std::endl;
			return 1;
		}
	}
	return 0;
}

int main_(int argc, char *argv[]) {

//#ifdef _WIN64 
//	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
//#endif

	if (argc == 2) {
		if (std::string(argv[1]) == std::string("test")) {
			TestAll();
			return 0;
		}
	}
	if (argc == 5) {

		const std::string sequence = std::string(argv[1]);
		const std::string structure = std::string(argv[2]);
		const int W = std::stoi(std::string(argv[3]));
		const std::string algo = std::string(argv[4]);
		const int n = sequence.length();
		const int max_loop = n < 30 ? n : 30;

		if (verification(structure, sequence, W, max_loop)) {
			return 1;
		}

		RintX1DOptions options;
		options.temperature = 37.0;
		options.param_file_name = std::string("Turner2004");
		options.max_span = W;
		options.max_loop = max_loop;
		options.sequence = sequence;
		options.reference_structure1 = structure;
		options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
		options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

		if (algo == std::string("RintPwithDFT")) {
			typedef WideComplexNumber<Floating> Comp;
			const auto ans1 = ComputeRintP1Dim<Comp>(sequence, options.S1, options.max_dim1, 37.0, options.max_span, options.max_loop, false);
			const auto ans2 = RegularizeRintP1Dim(ans1.first, ans1.second);
			for (int i = 0; i < ans2.first.size(); ++i) {
				std::cout << i << " " << ans2.first[i] << std::endl;
			}
			for (int i = 0; i < ans2.first.size(); ++i) {
				if (ans2.first[i] > 0.0) {
					for (int j = 1; j <= n; ++j)for (int k = 0; k < 6; ++k) {
						std::cout << i << " " << j - 1 << " " << k << " " << ans2.second[i][j][k] << std::endl;
					}
				}
			}
			return 0;
		}
		if (algo == std::string("RintPwithFFT")) {
			typedef WideComplexNumber<Floating> Comp;
			const auto ans1 = ComputeRintP1Dim<Comp>(sequence, options.S1, options.max_dim1, 37.0, options.max_span, options.max_loop, true);
			const auto ans2 = RegularizeRintP1Dim(ans1.first, ans1.second);
			for (int i = 0; i < ans2.first.size(); ++i) {
				std::cout << i << " " << ans2.first[i] << std::endl;
			}
			for (int i = 0; i < ans2.first.size(); ++i) {
				if (ans2.first[i] > 0.0) {
					for (int j = 1; j <= n; ++j)for (int k = 0; k < 6; ++k) {
						std::cout << i << " " << j - 1 << " " << k << " " << ans2.second[i][j][k] << std::endl;
					}
				}
			}
			return 0;
		}
		std::cout << "Error: invalid algo." << std::endl;
		return 1;
	}

	std::cout << "error" << std::endl;
	return 1;

}

}

int main(int argc, char *argv[]) {

	return rintp::main_(argc, argv);

}