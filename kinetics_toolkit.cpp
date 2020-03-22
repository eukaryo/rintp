/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include"kinetics_toolkit.h"

#include"parameter.h"
#include"misc.h"

namespace rintp {

std::vector<std::pair<std::string, double>> EnumerateStructureAndBoltzmannFactor(
	const std::string& sequence,
	const double temperature,
	const int max_span,
	const int max_loop) {

	parasor_param::InitializeParameter("Turner2004", temperature);

	const auto structures = EnumerateStructures(sequence, max_span, max_loop);

	std::vector<std::pair<std::string, double>> answer;

	for (const std::string s : structures) {
		const double F = EvalSpecificStructure(sequence, s);
		answer.push_back(std::make_pair(s, F));
	}

	return answer;
}

std::vector<std::string> EnumerateNeighbourStructure(
	const std::string& sequence,
	const std::string& structure,
	const int max_span,
	const int max_loop) {

	const int n = sequence.size();
	const std::string bp = "AU UA GC CG GU UG";
	std::string query = "XX";

	std::vector<std::string> answer;
	std::string new_structure = structure;
	std::vector<int> pv = DotNotationToPairVector(structure);

	const auto AddBpAtI = [&](const int i, const bool need_to_check_max_loop, const int forbid_j, const bool bigger_j_only) {
		query[0] = sequence[i];
		for (int j = (i + 1) % n; j != i; ++j) {

			const int span = std::abs(i - j);
			if (new_structure[j] == '.' && TURN < span && span <= max_span && j != forbid_j && !(bigger_j_only && j < i)) {
				query[1] = sequence[j];
				if (bp.find(query) != std::string::npos) {
					new_structure[std::min(i, j)] = '(';
					new_structure[std::max(i, j)] = ')';
					if ((!need_to_check_max_loop) || ComputeMaxLoop(new_structure) <= max_loop)answer.push_back(new_structure);
					new_structure[i] = '.';
					new_structure[j] = '.';
				}
			}

			if (pv[j] != -1) {
				j = pv[j];
			}
			if (j == n - 1) {
				j = -1;
			}
		}
	};

	//塩基を外す。
	//max-loop制約を破るような外し方はできない。

	int base_pair_count = 0;
	for (int i = 0; i < n; ++i) {
		if (i < pv[i]) {
			base_pair_count++;
			new_structure[i] = '.';
			new_structure[pv[i]] = '.';
			if (ComputeMaxLoop(new_structure) <= max_loop)answer.push_back(new_structure);
			new_structure[i] = '(';
			new_structure[pv[i]] = ')';
		}
	}

	//塩基を組む。
	//・exteriorに関して、max-span制約を破るような組み方はできない。
	//・pseudo-knotを作るような組み方はできない。つまり組むのは同じループ内の２塩基に限られる。
	//・TURN3以下のhairpinを作るような組み方はできない。

	if (base_pair_count == 0) {
		for (int i = 0; i < n; ++i) {
			query[0] = sequence[i];
			new_structure[i] = '(';
			for (int j = i + TURN + 1; j < n && j < i + max_span; ++j) {
				query[1] = sequence[j];
				if (bp.find(query) != std::string::npos) {
					new_structure[j] = ')';
					answer.push_back(new_structure);
					new_structure[j] = '.';
				}
			}
			new_structure[i] = '.';
		}
	}
	else {
		for (int i = 0; i < n; ++i) {
			if (structure[i] != '.')continue;
			AddBpAtI(i, false, -1, true);
		}
	}

	//MS2: 相異なるi,j,kについて、(i,j)が組まれている状態で、(i,j)が外れて(i,k)か(k,i)か(k,j)か(j,k)が組まれるような変化。
	//・塩基を外した瞬間にmax_loop制約が破られるのは構わないが、その場合は直後に組むことで制約が守られる必要がある。

	for (int i = 0; i < n; ++i) {
		if (i < pv[i]) {
			const int j = pv[i];
			pv[i] = -1;
			pv[j] = -1;
			new_structure[i] = '.';
			new_structure[j] = '.';

			const bool check_flag = max_loop < ComputeMaxLoop(new_structure);

			AddBpAtI(i, check_flag, j, false);
			AddBpAtI(j, check_flag, i, false);

			new_structure[i] = '(';
			new_structure[j] = ')';
			pv[i] = j;
			pv[j] = i;
		}
	}

	return answer;
}

std::vector<std::pair<std::string, double>>ComputeLocalMFE(
	const std::string& sequence,
	const std::string& initial_structure,
	const double temperature,
	const int max_span,
	const int max_loop) {

	//LocalMFEを探す。つまり、山登り法でボルツマン因子が高い構造を探す。

	parasor_param::InitializeParameter("Turner2004", temperature);

	std::vector<std::pair<std::string, double>>answer;
	answer.push_back(std::make_pair(initial_structure, EvalSpecificStructure(sequence, initial_structure)));

	while (true) {
		const std::vector<std::string> candidates = EnumerateNeighbourStructure(sequence, answer.back().first, max_span, max_loop);
		for (const auto c : candidates)VerificateInput(sequence, c, max_span, max_loop);
		double score = EvalSpecificStructure(sequence, candidates[0]);
		int champ = 0;
		for (int i = 1; i < candidates.size(); ++i) {
			
			const double s = EvalSpecificStructure(sequence, candidates[i]);
			if (score < s) {
				score = s;
				champ = i;
			}
		}
		if (answer.back().second < score)break;
		answer.push_back(std::make_pair(candidates[champ], score));
	}

	return answer;
}

std::vector<std::pair<std::string, double>> SimulateGillespie(
	const std::string& sequence,
	const std::string& initial_structure,
	const double temperature,
	const int max_span,
	const int max_loop,
	const int step,
	const uint64_t random_seed) {

	std::mt19937_64 rnd(random_seed);
	std::uniform_real_distribution<double>sampler01(0.0, 1.0);// [0.0,1.0)

	parasor_param::InitializeParameter("Turner2004", temperature);

	std::vector<std::pair<std::string, double>>answer;
	answer.push_back(std::make_pair(initial_structure, 0.0));

	for (int now_step = 0; now_step < step; ++now_step) {

		const std::string now_structure = answer.back().first;
		const double now_boltmann_factor = EvalSpecificStructure(sequence, now_structure);

		const std::vector<std::string> candidates = EnumerateNeighbourStructure(sequence, now_structure, max_span, max_loop);
		if (candidates.size() == 0)return answer;

		std::vector<double>candidate_kawasaki_flux(candidates.size());
		for (int i = 0; i < candidates.size(); ++i) {
			const double candidate_boltmann_factor = EvalSpecificStructure(sequence, candidates[i]);
			candidate_kawasaki_flux[i] = std::sqrt(candidate_boltmann_factor / now_boltmann_factor);
		}

		std::vector<double>candidate_kawasaki_flux_acc(candidates.size());
		candidate_kawasaki_flux_acc[0] = candidate_kawasaki_flux[0];
		for (int i = 1; i < candidates.size(); ++i) {
			candidate_kawasaki_flux_acc[i] = candidate_kawasaki_flux_acc[i - 1] + candidate_kawasaki_flux[i];
		}

		const double total_flux = candidate_kawasaki_flux_acc.back();
		const double r1 = 1.0 - sampler01(rnd);// (0.0,1.0]
		const double delta_t = -std::log(r1) / total_flux;
		const double next_time = answer.back().second + delta_t;

		const double r2 = sampler01(rnd) * total_flux;//[0.0, Φ)
		const int next_index = std::distance(candidate_kawasaki_flux_acc.begin(),
			std::lower_bound(
				candidate_kawasaki_flux_acc.begin(),
				candidate_kawasaki_flux_acc.end(), r2));

		answer.push_back(std::make_pair(candidates[next_index], next_time));
	}

	return answer;
}

}
