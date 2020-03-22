/*
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

#include"kinetics_toolkit.h"


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

void OutputStructuralProfile(const std::pair<std::vector<IntervalVar>, std::vector<std::vector<std::vector<IntervalVar>>>>&answer, const std::string& sequence, const std::string& structure) {

	std::cout << ">" << sequence << std::endl;
	std::cout << ">" << structure << std::endl;
	std::cout << answer.first.size() << std::endl;
	for (int i = 0; i < answer.first.size(); ++i) {
		std::cout << i << " " << answer.first[i].lower() << " " << answer.first[i].upper() << std::endl;
	}

	std::vector<std::string>features{ "bulge","exterior","hairpin","internal","multi","stem" };

	for (int i = 0; i < answer.first.size(); ++i) {
		for (int j = 0; j < answer.second[i].size() - 1; ++j)for (int k = 0; k < 6; ++k) {
			if (answer.first[i] > 0.0) {
				std::cout << i << " " << j << " " << features[k] << " " << answer.second[i][j + 1][k].lower() << " " << answer.second[i][j + 1][k].upper() << std::endl;
			}
			else {
				std::cout << i << " " << j << " " << features[k] << " " << 0 << " " << 0 << std::endl;
			}
		}

	}
}
void OutputStructuralProfile(const std::pair<std::vector<Floating>, std::vector<std::vector<std::vector<Floating>>>>&answer, const std::string& sequence, const std::string& structure) {

	std::cout << ">" << sequence << std::endl;
	std::cout << ">" << structure << std::endl;
	std::cout << answer.first.size() << std::endl;
	for (int i = 0; i < answer.first.size(); ++i) {
		std::cout << i << " " << answer.first[i] << std::endl;
	}

	std::vector<std::string>features{ "bulge","exterior","hairpin","internal","multi","stem" };

	for (int i = 0; i < answer.first.size(); ++i) {
		for (int j = 0; j < answer.second[i].size() - 1; ++j)for (int k = 0; k < 6; ++k) {
			if (answer.first[i] > 0.0) {
				std::cout << i << " " << j << " " << features[k] << " " << answer.second[i][j + 1][k] << std::endl;
			}
			else {
				std::cout << i << " " << j << " " << features[k] << " " << 0 << std::endl;
			}
		}
	}
}



int main_(int argc, char *argv[]) {

	//#ifdef _WIN64 
	//	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	//#endif

	//{
	//	const std::string sequence = "AGGACCUACGCUGCCCUAGAGGUUUUGCUAGGGAGGAGACGUGUGUGGCUGUAGCCACCCGUCCCGGGUACAAGUCCCGGGUGGUGAGGACGGUGUCUGUGGUUGUCUUCCCAGACUCUGCUUUCUGCCGUCUUCGGUCAAGUACCAGCUGGUGGUCCGCAUGUUU";
	//	const std::string structure = ".(((((.((((((((((((........))))))....((((.(.(((((....))))))))))(((((.......))))).((.(((((((((((...((((..((((....)))).))))....)))))))))))..)).....))))..)))))))........";
	//	const int W = 166;
	//	const int n = sequence.length();
	//	const int max_loop = n < 30 ? n : 30;
	//	std::vector<std::pair<std::string, double>> answer = ComputeLocalMFE(sequence, structure, 37.0, W, max_loop);
	//	return 0;
	//}

	if (argc == 2) {
		if (std::string(argv[1]) == std::string("test")) {
			TestAll();
			return 0;
		}
	}
	if (argc == 4) {

		const std::string sequence = std::string(argv[1]);
		const int W = std::stoi(std::string(argv[2]));
		const std::string algo = std::string(argv[3]);
		const int n = sequence.length();
		const int max_loop = n < 30 ? n : 30;

		if (algo == std::string("CentroidFold")) {
			const auto bppm = SimpleMcCaskillWide(sequence, std::string("Turner2004"), 37.0, W, max_loop).first;
			const std::string answer = GetCentroidFoldMcCaskill(sequence, W, bppm, 1.0, max_loop);
			std::cout << answer << std::endl;
			return 0;
		}
		if (algo == std::string("EnumerateAll")) {
			std::vector<std::pair<std::string, double>>answer = EnumerateStructureAndBoltzmannFactor(sequence, 27.0, W, max_loop);

			std::cout << "sequence: " << sequence << std::endl;
			std::cout << "max_span: " << W << std::endl;
			std::cout << "max_loop: " << max_loop << std::endl;

			for (int i = 0; i < answer.size(); ++i) {
				std::cout << answer[i].first << " " << std::setprecision(20) << answer[i].second << std::endl;
			}

			return 0;
		}


		std::cout << "Error: invalid algo." << std::endl;
		return 1;
	}
	if (argc == 5) {

		const std::string algo = std::string(argv[4]);

		if (algo == std::string("BoltzmannSampling")) {
			const std::string sequence = std::string(argv[1]);
			const int num_sample = std::stoi(std::string(argv[2]));
			const int seed = std::stoi(std::string(argv[3]));
			const int n = sequence.length();
			const int max_loop = n < 30 ? n : 30;

			std::pair<std::vector<std::string>, WideFloating> result = SampleMcCaskillEnergyAware(
				sequence, std::string("Turner2004"), 37.0, num_sample, n, max_loop, false, seed);

			std::cout << ">" << sequence << std::endl;
			std::cout << num_sample << std::endl;
			for (int i = 0; i < num_sample; ++i) {
				std::cout << result.first[i] << std::endl;
			}
			return 0;
		}

		const std::string sequence = std::string(argv[1]);
		const std::string structure = std::string(argv[2]);
		const int W = std::stoi(std::string(argv[3]));
		const int n = sequence.length();
		const int max_loop = n < 30 ? n : 30;

		if (VerificateInput(sequence, structure, W, max_loop)) {
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
			OutputStructuralProfile(ans2, sequence, structure);
			return 0;
		}
		if (algo == std::string("RintPwithFFT")) {
			typedef WideComplexNumber<Floating> Comp;
			const auto ans1 = ComputeRintP1Dim<Comp>(sequence, options.S1, options.max_dim1, 37.0, options.max_span, options.max_loop, true);
			const auto ans2 = RegularizeRintP1Dim(ans1.first, ans1.second);
			OutputStructuralProfile(ans2, sequence, structure);
			return 0;
		}
		if (algo == std::string("RintPwithDFTInterval")) {
			typedef WideComplexNumber<IntervalVar> Comp;
			const auto ans1 = ComputeRintP1Dim<Comp>(sequence, options.S1, options.max_dim1, 37.0, options.max_span, options.max_loop, false);
			const auto ans2 = RegularizeRintP1Dim(ans1.first, ans1.second);
			OutputStructuralProfile(ans2, sequence, structure);
			return 0;
		}
		if (algo == std::string("RintPwithFFTInterval")) {
			typedef WideComplexNumber<IntervalVar> Comp;
			const auto ans1 = ComputeRintP1Dim<Comp>(sequence, options.S1, options.max_dim1, 37.0, options.max_span, options.max_loop, true);
			const auto ans2 = RegularizeRintP1Dim(ans1.first, ans1.second);
			OutputStructuralProfile(ans2, sequence, structure);
			return 0;
		}

		if (algo == std::string("RintDwithDFT")) {
			typedef WideComplexNumber<Floating> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintD1Dim<Comp>(options);
			const auto ans2 = RegularizeRintD1Dim(ans1);
			std::cout << ">" << sequence << std::endl;
			std::cout << ">" << structure << std::endl;
			std::cout << ans2.size() << std::endl;
			for (int i = 0; i < ans2.size(); ++i) {
				std::cout << i << " " << ans2[i] << std::endl;
			}

			return 0;
		}
		if (algo == std::string("RintDwithFFT")) {
			typedef WideComplexNumber<Floating> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintD1Dim<Comp>(options);
			const auto ans2 = RegularizeRintD1Dim(ans1);
			std::cout << ">" << sequence << std::endl;
			std::cout << ">" << structure << std::endl;
			std::cout << ans2.size() << std::endl;
			for (int i = 0; i < ans2.size(); ++i) {
				std::cout << i << " " << ans2[i] << std::endl;
			}

			return 0;
		}
		if (algo == std::string("LocalMFE")) {
			std::vector<std::pair<std::string, double>>answer = ComputeLocalMFE(sequence, structure, 37.0, W, max_loop);
			std::cout << ">" << sequence << std::endl;
			std::cout << ">" << structure << std::endl;
			std::cout << answer.size() << std::endl;
			for (int i = 0; i < answer.size(); ++i) {
				std::cout << answer[i].first << " " << answer[i].second << std::endl;
			}
			return 0;
		}


		std::cout << "Error: invalid algo." << std::endl;
		return 1;
	}

	if (argc == 7) {
		const std::string sequence = std::string(argv[1]);
		const std::string structure = std::string(argv[2]);
		const int W = std::stoi(std::string(argv[3]));
		const int step = std::stoi(std::string(argv[4]));
		const uint64_t seed = std::stoi(std::string(argv[5]));
		const std::string algo = std::string(argv[6]);
		const int n = sequence.length();
		const int max_loop = n < 30 ? n : 30;

		if (VerificateInput(sequence, structure, W, max_loop)) {
			return 1;
		}

		if (algo == std::string("Gillespie")) {

			std::vector<std::pair<std::string, double>> answer = SimulateGillespie(sequence, structure, 37.0, W, max_loop, step, seed);

			std::cout << "sequence: " << sequence << std::endl;
			std::cout << "max_span: " << W << std::endl;
			std::cout << "max_loop: " << max_loop << std::endl;
			std::cout << "step: " << step << std::endl;
			std::cout << "seed: " << seed << std::endl;

			for (int i = 0; i < answer.size(); ++i) {
				std::cout << answer[i].first << " " << std::setprecision(20) << answer[i].second << std::endl;
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
