/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/


#include"rintcapr.h"

#include"parameter.h"
#include"misc.h"
#include"fourier_transform.h"
#include"complex_number.h"
#include"real_logsumexp.h"


namespace rintp {

//explicit instantiation

template std::pair<std::vector<WideComplexNumber<IntervalVar>>,
	std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>> ComputeRintP1Dim(
	const std::string sequence,
	const std::vector<std::vector<int>> S,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop,
	const bool allow_fft);

template std::pair<std::vector<WideComplexNumber<Floating>>,
	std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>> ComputeRintP1Dim(
	const std::string sequence,
	const std::vector<std::vector<int>> S,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop,
	const bool allow_fft);

//function definition

template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintP1Dim(
	const std::string sequence,
	const std::vector<std::vector<int>> S,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop,
	const bool allow_fft) {

	typedef decltype(Comp::real) RealScalar;

	const int fourier_dim = allow_fft ? Ceiling2Power(max_dim + 1) : (max_dim + 1);

	parasor_param::InitializeParameter("Turner2004", temperature);

	const int n = int(sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(S);

	const auto g_alpha_multiRclosing = [&](const int i, const int j) {return C[i][j] - C[i + 1][j]; };
	const auto g_alpha_multi2nLbranchRpaired = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto g_alpha_multi1Lbranch = [&](const int i, const int j) {return C[i][j] - C[i][j - 1]; };
	const auto g_alpha_stem = [&](const int i, const int j) {return C[i][j] - C[i + 1][j - 1] + 1 - 2 * S[i][j]; };
	const auto g_alpha_hairpin = [&](const int i, const int j) {return C[i][j] + 1 - 2 * S[i][j]; };
	const auto g_alpha_internal = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * S[i][j]; };

	const auto gz0 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C[i][j];
	};

	const auto g_alpha_all_1 = [&](const int i, const int j) {return C[i][j]; };
	const auto g_alpha_all_2 = [&](const int i, const int j, const int k) {return gz0(i, j) - gz0(i, k) - gz0(k + 1, j); };
	const auto g_alpha_all_3 = [&](const int i, const int j) {return gz0(i, j) - gz0(i + 1, j); };
	const auto g_alpha_all1_1 = g_alpha_multi1Lbranch;

	const auto g_beta_stem1 = [&](const int i, const int j) {return C[i][j] - C[i + 1][j - 1] + 1 - 2 * S[i][j]; };
	const auto g_beta_stem2 = [&](const int i, const int j) {return gz0(1, n) - gz0(i + 1, j - 1) - gz0(1, i - 1) - gz0(j + 1, n) + 1 - 2 * S[i][j]; };
	const auto g_beta_stem3 = [&](const int i, const int j) {return (1 - 2 * S[i][j] + fourier_dim) % fourier_dim; };
	const auto g_beta_stem4 = [&](const int i, const int j, const int h, const int l) {return gz0(h + 1, l - 1) - gz0(i + 1, j - 1) + 1 - 2 * S[i][j]; };
	const auto g_beta_multi1Lbranch1 = [&](const int i, const int j) {return C[i + 1][j] - C[i + 1][j - 1]; };
	const auto g_beta_multi1Lbranch2 = [&](const int i, const int j, const int k) {return (C[k + 1][j - 1] - C[k][i - 1] - C[i + 1][j - 1] + fourier_dim) % fourier_dim; };
	const auto g_beta_multi12nLbranch = [&](const int i, const int j, const int k) {return (C[i + 1][k - 1] - C[j + 1][k] - C[i + 1][j - 1] + fourier_dim) % fourier_dim; };
	const auto g_beta_multiRclosing1 = [&](const int i, const int j) {return C[i][j - 1] - C[i + 1][j - 1]; };
	const auto g_beta_multiRclosing2 = [&](const int i, const int j) {return C[i][j] - C[i + 1][j - 1]; };

	const auto g_pb1 = [&](const int p, const int j, const int k) {return C[j + 1][k - 1] - C[p][k - 1]; };
	const auto g_pb2 = [&](const int q, const int j, const int k) {return C[j + 1][k - 1] - C[j + 1][q]; };
	const auto g_pe1 = [&](const int i) {return gz0(1, n) - gz0(1, i - 1) - gz0(i + 1, n); };
	const auto g_ph1 = [&](const int j, const int k) {return C[j + 1][k - 1]; };
	const auto g_pi1 = [&](const int j, const int k, const int p, const int q) {return C[j + 1][k - 1] - C[p][q]; };
	const auto g_pm1 = [&](const int i, const int j) {return (C[i + 1][j - 1] - C[i + 1][j] + fourier_dim) % fourier_dim; };
	const auto g_pm2 = [&](const int y, const int i) {return (C[y + 1][i - 1] - C[y][i - 1] + fourier_dim) % fourier_dim; };

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	std::vector<Comp>zeta(fourier_dim, Comp(0.0, 0.0));
	std::vector<std::vector<std::vector<Comp>>>P(fourier_dim, std::vector<std::vector<Comp>>(n + 1, std::vector<Comp>(5, Comp(0.0, 0.0))));

#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < fourier_dim; ++x) {

		//1の(fourier_dim)乗根を全て求める。
		std::vector<Comp> root_of_unity(fourier_dim);
		for (int n1 = 0; n1 < fourier_dim; ++n1) {
			const int d = (x * n1) % fourier_dim;
			root_of_unity[n1] = Comp::GetPolar(
				RealScalar(1.0),
				RealScalar(2.0) * pi * RealScalar(double(d)) / RealScalar(double(fourier_dim)));
		}

		const auto GetRoot = [&](const int i) {
			assert(0 <= i && i < fourier_dim);
			return root_of_unity[i];
		};

		//メモ化再帰のためのデータ構造である。firstは初期化時はfalseで、値を計算したらtrueにして、secondに値を入れる。
		typedef std::pair<bool, Comp>MemComp;

		//[CapR論文]の11～12ページの式の左辺たちである。添字を1-originとする。
		std::vector<MemComp>AlphaMultiRclosing((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>AlphaMulti12nLbranchRanother((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>AlphaMulti2nLbranchRpaired((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>AlphaMulti1Lbranch((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>AlphaStem((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::vector<MemComp>AlphaAll((n + 1) * 2, std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>AlphaAll1((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::vector<MemComp>BetaStem((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>BetaMulti1Lbranch((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>BetaMulti2nLbranchRpaired((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>BetaMulti12nLbranchRanother((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>BetaMultiRclosing((n + 1) * (max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		//関数呼び出しの依存関係が循環しているので、先に全て宣言してから後で定義する。あとこのほうが論文の記述に沿っている

		//k<iなる(k,j+1)で閉じられるマルチループで、[k+1,i-1]がunpaired
		std::function<Comp(int, int)> GetAlphaMultiRclosing;

		//(i,j)より外側で閉じられるマルチループで、(i,j)内に2つ以上の分岐を含み、iがpairedでj+1がpaired
		std::function<Comp(int, int)> GetAlphaMulti2nLbranchRpaired;

		//(i,j)より外側で閉じられるマルチループで、
		//(i,j)内に1つ以上の分岐を含み、iがpairedでj+1がpairedで、j+1はclosingでない。かつ、i-1未満で最大のpaired baseがclosingである
		std::function<Comp(int, int)> GetAlphaMulti12nLbranchRanother;

		//(i,j)より外側で閉じられるマルチループで、(i,j)内にちょうど1つの分岐を含み、iがpaired
		std::function<Comp(int, int)> GetAlphaMulti1Lbranch;

		std::function<Comp(int, int)> GetAlphaStem;//(i,j)が組んでいる場合。
		std::function<Comp(int, int)> GetAlphaAll;//(i,j)内の分配関数、McCaskillのZ
		std::function<Comp(int, int)> GetAlphaAll1;//(i,j)内の分配関数のうちちょうど1つの分岐(i,k)を含む(ただしi<k<=j)。McCaskillのZ1

		std::function<Comp(int, int)> GetBetaStem;
		std::function<Comp(int, int)> GetBetaMulti1Lbranch;
		std::function<Comp(int, int)> GetBetaMulti12nLbranchRanother;
		std::function<Comp(int, int)> GetBetaMulti2nLbranchRpaired;
		std::function<Comp(int, int)> GetBetaMultiRclosing;

		std::function<Comp(int)> GetPB;
		std::function<Comp(int)> GetPE;
		std::function<Comp(int)> GetPH;
		std::function<Comp(int)> GetPI;
		std::function<Comp(int)> GetPM;

		GetAlphaMultiRclosing = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (AlphaMultiRclosing[at(i, j)].first)return AlphaMultiRclosing[at(i, j)].second;

			Comp ans(0.0, 0.0);
			ans += GetAlphaMulti2nLbranchRpaired(i, j);
			ans += GetAlphaMultiRclosing(i + 1, j)
				* GetRoot(g_alpha_multiRclosing(i, j));
			return (AlphaMultiRclosing[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetAlphaMulti2nLbranchRpaired = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (AlphaMulti2nLbranchRpaired[at(i, j)].first)return AlphaMulti2nLbranchRpaired[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = i + 1; k < j; ++k) {
				ans += GetAlphaMulti12nLbranchRanother(i, k)
					* GetAlphaMulti1Lbranch(k + 1, j)
					* GetRoot(g_alpha_multi2nLbranchRpaired(i, j, k));
			}
			return (AlphaMulti2nLbranchRpaired[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetAlphaMulti12nLbranchRanother = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (AlphaMulti12nLbranchRanother[at(i, j)].first)return AlphaMulti12nLbranchRanother[at(i, j)].second;

			Comp ans(0.0, 0.0);

			ans += GetAlphaMulti2nLbranchRpaired(i, j);
			ans += GetAlphaMulti1Lbranch(i, j);

			return (AlphaMulti12nLbranchRanother[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetAlphaMulti1Lbranch = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (AlphaMulti1Lbranch[at(i, j)].first)return AlphaMulti1Lbranch[at(i, j)].second;

			Comp ans(0.0, 0.0);

			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type != 0) {
				ans += GetAlphaStem(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
			ans += GetAlphaMulti1Lbranch(i, j - 1) * GetRoot(g_alpha_multi1Lbranch(i, j));
			return (AlphaMulti1Lbranch[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetAlphaStem = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			//iとjが塩基対を組み得ないならゼロを返す。
			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > max_span) {
				return Comp(0.0, 0.0);
			}

			if (AlphaStem[at(i, j)].first)return AlphaStem[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//hairpin
			ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence))) * GetRoot(g_alpha_hairpin(i, j));

			//internal, bulge
			for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
					if (type == 0)continue;
					const int unpaired_base2 = j - l - 1;
					if (unpaired_base1 + unpaired_base2 == 0)continue;
					ans += GetAlphaStem(k, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence)))
						* GetRoot(g_alpha_internal(i, j, k, l));
				}
			}

			//multi
			const int rtype = parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]);
			ans += GetAlphaMultiRclosing(i + 1, j - 1)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, sequence)))
				* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
				* GetRoot(g_alpha_stem(i, j));

			//stem
			if (i + 1 + TURN + 1 <= j - 1) {
				ans += GetAlphaStem(i + 1, j - 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, (i + 1) - 1, (j - 1) - 1, sequence)))
					* GetRoot(g_alpha_stem(i, j));
			}

			return (AlphaStem[at(i, j)] = std::make_pair(true, ans)).second;
		};

		GetAlphaAll = [&](const int i, const int j) {

			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(1.0, 0.0);
			if (j == i - 1)return Comp(1.0, 0.0);
			assert(i == 1 || j == n);
			const int index = (i == 1) ? j : (i + n + 1);
			assert(1 <= index && index < (n + 1) * 2);
			if (AlphaAll[index].first)return AlphaAll[index].second;

			Comp ans(0.0, 0.0);

			//左右端を含むケースだけを考慮すればよいことを利用して空間計算量を抑えた。
			//そのため、この関数の実装は[McCaskill, 1990]の記述とはかけ離れている。

			if (i == 1) {

				//一切塩基対を組まないケース
				ans += GetRoot(g_alpha_all_1(i, j));

				//最外側塩基対が1つだけで、かつそれが(i,*)であるケース
				ans += GetAlphaAll1(i, j);

				//最外側塩基対が1つ以上で、そのうち最も右側の塩基対が(k+1,*)であるケース
				for (int k = i; k <= j - 1; ++k) {
					ans += GetAlphaAll(i, k) * GetAlphaAll1(k + 1, j) * GetRoot(g_alpha_all_2(1, j, k));
				}
			}
			else {

				//iが最外側塩基対を組まないケース
				ans += GetAlphaAll(i + 1, j) * GetRoot(g_alpha_all_3(i, j));

				//(i,k)が最外側塩基対を組むケース
				for (int k = i + TURN + 1; k <= j && k - i <= max_span; ++k) {
					const int type = parasor_param::GetPairType(sequence[i - 1], sequence[k - 1]);
					if (type == 0)continue;
					ans += GetAlphaStem(i, k)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence)))
						* GetAlphaAll(k + 1, j)
						* GetRoot(g_alpha_all_2(i, j, k));
				}
			}

			return (AlphaAll[index] = std::make_pair(true, ans)).second;
		};
		GetAlphaAll1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j) {
				return Comp(0.0, 0.0);
			}

			if ((j - i) > max_span) {
				return GetAlphaAll1(i, i + max_span)*GetRoot(C[i][j] - C[i][i + max_span]);
			}

			if (AlphaAll1[at(i, j)].first)return AlphaAll1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type != 0) {
				ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)))
					* GetAlphaStem(i, j);
			}
			ans += GetAlphaAll1(i, j - 1) * GetRoot(g_alpha_all1_1(i, j));
			return (AlphaAll1[at(i, j)] = std::make_pair(true, ans)).second;
		};

		GetBetaStem = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			//iとjが塩基対を組み得ないならゼロを返す。
			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > max_span) {
				return Comp(0.0, 0.0);
			}

			if (BetaStem[at(i, j)].first)return BetaStem[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//(i,j)が最も外側の場合
			ans += GetAlphaAll(1, i - 1)
				* GetAlphaAll(j + 1, n)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)))
				* GetRoot(g_beta_stem2(i, j));

			//(i,j)がmultiloopの分岐の一つである場合
			ans += GetBetaMulti1Lbranch(i, j)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
				* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
				* GetRoot(g_beta_stem3(i, j));

			//(i,j)の外に塩基対(h,l)があってinternal loopかbulge loopを構成する場合(stem loopはダメ)
			for (int h = std::max(1, std::max(i - max_span + 1, i - max_loop - 1)); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= max_span && (i - h - 1) + (l - j - 1) <= max_loop; ++l) {
					if (h == i - 1 && l == j + 1)continue;
					ans += GetBetaStem(h, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence)))
						* GetRoot(g_beta_stem4(i, j, h, l));
				}
			}

			//stem
			if (2 <= i && j <= n - 1) {
				ans += GetBetaStem(i - 1, j + 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy((i - 1) - 1, (j + 1) - 1, i - 1, j - 1, sequence)))
					* GetRoot(g_beta_stem1(i, j));
			}

			return (BetaStem[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetBetaMulti1Lbranch = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || j == n/*定義よりjの右に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (BetaMulti1Lbranch[at(i, j)].first)return BetaMulti1Lbranch[at(i, j)].second;

			Comp ans(0.0, 0.0);
			ans += GetBetaMulti1Lbranch(i, j + 1)
				* GetRoot(g_beta_multi1Lbranch1(i, j));
			ans += GetBetaMulti12nLbranchRanother(i, j);
			for (int k = std::max(1, j - 1 - max_span); k < i; ++k) {
				ans += GetBetaMulti2nLbranchRpaired(k, j)
					* GetAlphaMulti12nLbranchRanother(k, i - 1)
					* GetRoot(g_beta_multi1Lbranch2(i, j, k));
			}
			return (BetaMulti1Lbranch[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetBetaMulti12nLbranchRanother = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (BetaMulti12nLbranchRanother[at(i, j)].first)return BetaMulti12nLbranchRanother[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = j + 1; k <= n - 1 && (k - i) <= max_span; k++) {
				ans += GetBetaMulti2nLbranchRpaired(i, k)
					* GetAlphaMulti1Lbranch(j + 1, k)
					* GetRoot(g_beta_multi12nLbranch(i, j, k));
			}
			return (BetaMulti12nLbranchRanother[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetBetaMulti2nLbranchRpaired = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (BetaMulti2nLbranchRpaired[at(i, j)].first)return BetaMulti2nLbranchRpaired[at(i, j)].second;

			Comp ans(0.0, 0.0);
			ans += GetBetaMulti12nLbranchRanother(i, j);
			ans += GetBetaMultiRclosing(i, j);
			return (BetaMulti2nLbranchRpaired[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetBetaMultiRclosing = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if (i + TURN >= j || (j - i) > max_span || i == 1 || j == n/*定義よりi,jの外に塩基が必要*/) {
				return Comp(0.0, 0.0);
			}

			if (BetaMultiRclosing[at(i, j)].first)return BetaMultiRclosing[at(i, j)].second;

			Comp ans(0.0, 0.0);
			ans += GetBetaMultiRclosing(i - 1, j)
				* GetRoot(g_beta_multiRclosing1(i, j));


			const int rtype = parasor_param::GetPairTypeReverse(sequence[(i - 1) - 1], sequence[(j + 1) - 1]);
			if (rtype != 0) {
				ans += GetBetaStem(i - 1, j + 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j + 1) - 1 - 1, (i - 1) - 1 + 1, false, sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* GetRoot(g_beta_multiRclosing2(i, j));
			}

			return (BetaMultiRclosing[at(i, j)] = std::make_pair(true, ans)).second;
		};

		GetPB = [&](const int i) {
			assert(1 <= i && i <= n);

			Comp ans(0.0, 0.0);

			for (int j = std::max(1, i - max_span); j <= i - 1; ++j)for (int k = i + 1; k <= std::min(n, j + max_span); ++k) {
				const int type = parasor_param::GetPairType(sequence[j - 1], sequence[k - 1]);
				if (type == 0)continue;
				for (int p = i + 1; p <= std::min(j + max_loop + 1, k - 1); ++p) {
					//j<i<p<kで、塩基対(j,k)と(p,k-1)がバルジを構成する場合
					if (!(1 <= j&&j < i&&i < p&&p + TURN < k - 1 && k <= n))continue;
					ans += GetBetaStem(j, k)
						* GetAlphaStem(p, k - 1)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(j - 1, k - 1, p - 1, (k - 1) - 1, sequence)))
						* GetRoot(g_pb1(p, j, k));
				}
				for (int q = std::max(j + 1, k - max_loop - 1); q < i; ++q) {
					//j<q<i<kで、塩基対(j,k)と(j+1,q)がバルジを構成する場合
					if (!(1 <= j&&j + 1 + TURN < q&&q < i&&i < k&&k <= n))continue;
					ans += GetBetaStem(j, k)
						* GetAlphaStem(j + 1, q)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(j - 1, k - 1, (j + 1) - 1, q - 1, sequence)))
						* GetRoot(g_pb2(q, j, k));
				}
			}
			return ans;
		};
		GetPE = [&](const int i) {
			assert(1 <= i && i <= n);

			Comp ans(0.0, 0.0);

			ans = GetAlphaAll(1, i - 1)
				* GetAlphaAll(i + 1, n)
				* GetRoot(g_pe1(i));

			return ans;
		};
		GetPH = [&](const int i) {
			assert(1 <= i && i <= n);

			Comp ans(0.0, 0.0);
			for (int j = std::max(1, i - max_span); j <= i - 1; ++j)for (int k = i + 1; k <= std::min(n, j + max_span); ++k) {
				const int type = parasor_param::GetPairType(sequence[j - 1], sequence[k - 1]);
				if (type == 0)continue;
				//j<i<kで、塩基対(j,k)がヘアピンを構成する場合
				if (j + TURN >= k)continue;
				ans += GetBetaStem(j, k)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(j - 1, k - 1, sequence)))
					* GetRoot(g_ph1(j, k));
			}
			return ans;
		};
		GetPI = [&](const int i) {
			assert(1 <= i && i <= n);

			Comp ans(0.0, 0.0);
			for (int j = std::max(1, i - max_span); j <= i - 1; ++j)for (int k = i + 1; k <= std::min(n, j + max_span); ++k) {
				const int type = parasor_param::GetPairType(sequence[j - 1], sequence[k - 1]);
				if (type == 0)continue;
				for (int p = i + 1; p <= k - 3 && (p - j - 1) <= max_loop; ++p) {
					for (int q = std::max(p + TURN + 1, (p - j - 1) + (k - max_loop - 1)); q <= k - 2; ++q) {
						//j<i<p<q<kで、塩基対(j,k)と(p,q)がinternal loopを構成する場合
						ans += GetBetaStem(j, k)
							* GetAlphaStem(p, q)
							* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(j - 1, k - 1, p - 1, q - 1, sequence)))
							* GetRoot(g_pi1(j, k, p, q));
					}
				}
				for (int p = j + 2; p <= i - 2 && (p - j - 1) <= max_loop; ++p) {
					for (int q = std::max(p + TURN + 1, (p - j - 1) + (k - max_loop - 1)); q <= i - 1; ++q) {
						//j<p<q<i<kで、塩基対(j,k)と(p,q)がinternal loopを構成する場合
						ans += GetBetaStem(j, k)
							* GetAlphaStem(p, q)
							* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(j - 1, k - 1, p - 1, q - 1, sequence)))
							* GetRoot(g_pi1(j, k, p, q));
					}
				}
			}
			return ans;
		};
		GetPM = [&](const int i) {
			assert(1 <= i && i <= n);

			Comp ans(0.0, 0.0);

			//iのすぐ下の塩基対がmultiloopを閉じる場合
			for (int j = i + 1; j <= std::min(n, i + max_span); ++j) {
				ans += GetBetaMultiRclosing(i, j)
					* GetAlphaMultiRclosing(i + 1, j)
					* GetRoot(g_pm1(i, j));
			}

			//iのすぐ下の塩基対がmultiloopを閉じるもの以外の分岐である場合
			for (int y = std::max(1, i - max_span); y <= i - 1; ++y) {
				ans += GetBetaMulti1Lbranch(y, i)
					* GetAlphaMulti1Lbranch(y, i - 1)
					* GetRoot(g_pm2(y, i));
			}


			return ans;
		};

		zeta[x] = GetAlphaAll(1, n);
		for (int i = 1; i <= n; ++i) {
			P[x][i][0] = GetPB(i);
			P[x][i][1] = GetPE(i);
			P[x][i][2] = GetPH(i);
			P[x][i][3] = GetPI(i);
			P[x][i][4] = GetPM(i);
		}
	}

	const std::vector<Comp>z_ = FourierTransform(zeta, allow_fft);

	std::vector<std::vector<std::vector<Comp>>>profile(max_dim + 1, std::vector<std::vector<Comp>>(n + 1, std::vector<Comp>(5, Comp(0.0, 0.0))));
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 1; i <= n; ++i) {
		for (int p = 0; p < 5; ++p) {
			std::vector<Comp>tmp_(max_dim + 1);
			for (int j = 0; j <= max_dim; ++j)tmp_[j] = P[j][i][p];
			const std::vector<Comp>tmp = FourierTransform(tmp_, allow_fft);
			for (int j = 0; j <= max_dim; ++j)profile[j][i][p] = tmp[j];
		}
	}

	if (allow_fft && (fourier_dim != max_dim)) {
		std::vector<Comp>zz_(max_dim + 1, Comp(0.0, 0.0));
		for (int d = 0; d <= max_dim; ++d)zz_[d] = z_[d];
		return std::make_pair(zz_, profile);
	}

	return std::make_pair(z_, profile);
}

std::pair<std::vector<Floating>, std::vector<std::vector<std::vector<Floating>>>>RegularizeRintP1Dim(
	const std::vector<WideComplexNumber<Floating>>& z,
	const std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>& profile) {

	//ComputeRintP1Dimの結果を実数にして、正規化して、logsumexpではない普通の数値型にして返す。

	typedef WideComplexNumber<Floating> Comp;
	const int max_dim = int(z.size());
	assert(profile.size() == max_dim);
	const int n = profile[0].size() - 1;

	//zについて
	std::vector<Floating> ans_z(max_dim, Floating(0.0));

	//内側分配関数の総和を計算する。
	Comp sum_comp(0.0, 0.0);
	for (int i = 0; i < max_dim; ++i) {
		sum_comp += z[i];
	}

	//内側分配関数の各値を総和で割ってans_zに入れる。
	for (int i = 0; i < max_dim; ++i) {
		ans_z[i] = z[i].IsZero() ? Floating(0.0) : exp(z[i].log_scale - sum_comp.log_scale);
		ans_z[i] = std::min<Floating>(Floating(1.0), ans_z[i]);
		ans_z[i] = std::max<Floating>(Floating(0.0), ans_z[i]);
	}

	//profileについて
	std::vector<std::vector<std::vector<Floating>>>ans_profile(
		max_dim, std::vector<std::vector<Floating>>(
			n + 1, std::vector<Floating>(6, Floating(0.0))));
	for (int i = 0; i < max_dim; ++i) {
		assert(profile[i].size() == n + 1);
		for (int j = 1; j <= n; ++j) {

			//この時点で、引数profile[i][j][k]は内側分配関数の値になっている。(k∈[0,5))
			//詳しく言うと、ハミング距離iの構造たちのうち塩基jのプロファイルがkであるような構造たちの内側分配関数である。
			//分母は引数z[i]である。

			//profile[i][j][k](k∈[0,5))の各値を総和で割ってans_profileに入れる。
			for (int k = 0; k < 5; ++k) {
				ans_profile[i][j][k] = (profile[i][j][k].IsZero() || z[i].IsZero()) ? Floating(0.0) : exp(profile[i][j][k].log_scale - z[i].log_scale);
				ans_profile[i][j][k] = std::min<Floating>(Floating(1.0), ans_profile[i][j][k]);
				ans_profile[i][j][k] = std::max<Floating>(Floating(0.0), ans_profile[i][j][k]);
			}

			//profile[i][j][5]はそれ以外からの引き算として処理する。
			ans_profile[i][j][5] = Floating(1.0);
			for (int k = 0; k < 5; ++k)ans_profile[i][j][5] -= ans_profile[i][j][k];
			ans_profile[i][j][5] = std::max(Floating(0.0), ans_profile[i][j][5]);
		}
	}

	return std::make_pair(ans_z, ans_profile);
}

std::pair<std::vector<IntervalVar>, std::vector<std::vector<std::vector<IntervalVar>>>>RegularizeRintP1Dim(
	const std::vector<WideComplexNumber<IntervalVar>>& z,
	const std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>& profile) {

	//ComputeRintP1Dimの結果を実数にして、正規化して、logsumexpではない普通の数値型にして返す。

	typedef WideComplexNumber<IntervalVar> Comp;
	const int max_dim = int(z.size());
	const Floating FINF = std::numeric_limits<Floating>::infinity();
	assert(profile.size() == max_dim);
	const int n = profile[0].size() - 1;

	//zについて
	std::vector<IntervalVar>ans_z(max_dim, IntervalVar(0.0, 0.0));

	//内側分配関数の総和を計算する。
	Comp sum_comp = z[0];
	for (int i = 1; i < max_dim; ++i)sum_comp += z[i];
	assert(sum_comp.real.lower() > Floating(0.0));
	Comp sum_comp_inverse = sum_comp;
	sum_comp_inverse.real = intersect(IntervalVar(0.0, std::numeric_limits<Floating>::infinity()), sum_comp_inverse.real);
	sum_comp_inverse.real = IntervalVar(1.0, 1.0) / sum_comp_inverse.real;
	assert(sum_comp_inverse.imag.lower() <= Floating(0.0) && Floating(0.0) <= sum_comp_inverse.imag.upper());
	sum_comp_inverse.log_scale *= IntervalVar(-1.0, -1.0);

	//内側分配関数の各値を総和で割ってans_zに入れる。
	for (int i = 0; i < max_dim; ++i) {
		ans_z[i] = (z[i] * sum_comp_inverse).ToUsualComp().real;
		ans_z[i] = /*kv::*/intersect(IntervalVar(0.0, 1.0), ans_z[i]);
		//ans_z[i] = std::min<IntervalVar>(IntervalVar(1.0, 1.0), ans_z[i]);
		//ans_z[i] = std::max<IntervalVar>(IntervalVar(0.0, 0.0), ans_z[i]);
	}

	//profileについて
	std::vector<std::vector<std::vector<IntervalVar>>>ans_profile(
		max_dim, std::vector<std::vector<IntervalVar>>(
			n + 1, std::vector<IntervalVar>(6, IntervalVar(0.0, 0.0))));
	{
		for (int i = 0; i < max_dim; ++i) {
			assert(profile[i].size() == n + 1);
			for (int j = 1; j <= n; ++j) {

				//profile[i][j][k](k∈[0,5))の各値を総和で割ってans_profileに入れる。
				for (int k = 0; k < 5; ++k) {
					Comp zi_comp_inverse = z[i];
					if (!z[i].IsZero()) {
						assert(zi_comp_inverse.real.upper() > Floating(0.0));
						bool flag = false;
						if (zi_comp_inverse.real.lower() < Floating(0.0)) {
							flag = true;
							zi_comp_inverse.real.upper() = std::max<Floating>(zi_comp_inverse.real.upper(), -zi_comp_inverse.real.lower());
							zi_comp_inverse.real.lower() = std::numeric_limits<Floating>::min();
						}
						zi_comp_inverse.real = IntervalVar(1.0, 1.0) / zi_comp_inverse.real;
						if(flag)zi_comp_inverse.real.upper() = std::numeric_limits<Floating>::infinity();
					}
					assert(zi_comp_inverse.imag.lower() <= Floating(0.0) && Floating(0.0) <= zi_comp_inverse.imag.upper());
					zi_comp_inverse.log_scale *= IntervalVar(-1.0, -1.0);
					ans_profile[i][j][k] = (profile[i][j][k] * zi_comp_inverse).ToUsualComp().real;					
					ans_profile[i][j][k] = /*kv::*/intersect(IntervalVar(0.0, 1.0), ans_profile[i][j][k]);
					//ans_profile[i][j][k] = /*kv::*/max(ans_profile[i][j][k], IntervalVar(0.0, 0.0));
					//ans_profile[i][j][k] = /*kv::*/min(ans_profile[i][j][k], IntervalVar(1.0, 1.0));
				}

				//profile[i][j][5]はそれ以外からの引き算として処理する。
				ans_profile[i][j][5] = IntervalVar(1.0, 1.0);
				for (int k = 0; k < 5; ++k)ans_profile[i][j][5] -= ans_profile[i][j][k];
				ans_profile[i][j][5] = /*kv::*/intersect(IntervalVar(0.0, 1.0), ans_profile[i][j][5]);
				//ans_profile[i][j][5] = /*kv::*/max(ans_profile[i][j][5], IntervalVar(0.0, 0.0));
				//ans_profile[i][j][5] = /*kv::*/min(ans_profile[i][j][5], IntervalVar(1.0, 1.0));

				//この時点で、引数profile[i][j][k]は内側分配関数の値になっている。(k∈[0,5))
				//詳しく言うと、ハミング距離iの構造たちのうち塩基jのプロファイルがkであるような構造たちの内側分配関数である。
				//分母は引数z[i]である。

				////まず内側分配関数は0以上の実数なので、そのように変換する。
				//std::vector<WideRealNumber<IntervalVar>>tmp_profile(5, 0.0);
				//for (int k = 0; k < 5; ++k) {
				//	tmp_profile[k] = WideRealNumber<IntervalVar>::LogRealToWide(profile[i][j][k].log_scale);
				//	WideRealNumber<IntervalVar> tmp;
				//	tmp.log_scale = log(/*kv::*/max(IntervalVar(0.0, 0.0), profile[i][j][k].real));
				//	tmp_profile[k] *= tmp;
				//}

				////profile[i][j][k](k∈[0,5))の各値を総和で割ってans_profileに入れる。
				//for (int k = 0; k < 5; ++k) {
				//	ans_profile[i][j][k] = (tmp_profile[k] / tmp_z[i]).ToUsualReal();
				//	ans_profile[i][j][k] = /*kv::*/max(ans_profile[i][j][k], IntervalVar(0.0, 0.0));
				//	ans_profile[i][j][k] = /*kv::*/min(ans_profile[i][j][k], IntervalVar(1.0, 1.0));
				//}

				////profile[i][j][5]はそれ以外からの引き算として処理する。
				//ans_profile[i][j][5] = IntervalVar(1.0, 1.0);
				//for (int k = 0; k < 5; ++k)ans_profile[i][j][5] -= ans_profile[i][j][k];
				//ans_profile[i][j][5] = /*kv::*/max(ans_profile[i][j][5], IntervalVar(0.0, 0.0));
				//ans_profile[i][j][5] = /*kv::*/min(ans_profile[i][j][5], IntervalVar(1.0, 1.0));
			}
		}
	}

	return std::make_pair(ans_z, ans_profile);
}


std::pair<std::vector<double>, std::vector<std::vector<std::vector<double>>>> BruteForceRintP1Dim(
	const std::string sequence,
	const std::string reference_structure,
	const int max_dim,
	const double temperature,
	const int max_span,
	const int max_loop) {

	//引数sequenceはAUCGから成るRNA配列とする。
	//引数reference_structureはリファレンス構造とする。
	//任意の二次構造について、リファレンス構造からのハミング距離はmax_dim以下と仮定する。

	//reference_structureからのハミング距離ごとにstructural profileを求めて返す。

	parasor_param::InitializeParameter("Turner2004", temperature);
	const int n = int(sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(sequence, max_span, max_loop);

	std::vector<double>z(max_dim + 1, 0);
	std::vector<std::vector<std::vector<double>>>ans(max_dim + 1, std::vector<std::vector<double>>(n + 1, std::vector<double>(6, 0.0)));
	for (const std::string s : structures) {
		const int d = ComputeHammingDistance(s, reference_structure);
		const double e = EvalSpecificStructure(sequence, s);
		z[d] += e;
		for (int i = 0; i < n; ++i) {
			const std::string context = ComputeStructuralContext(s, i);
			if (context == std::string("Bulge")) {
				ans[d][i + 1][0] += e;
			}
			else if (context == std::string("Exterior")) {
				ans[d][i + 1][1] += e;
			}
			else if (context == std::string("Hairpin")) {
				ans[d][i + 1][2] += e;
			}
			else if (context == std::string("Interior")) {
				ans[d][i + 1][3] += e;
			}
			else if (context == std::string("Multibranch")) {
				ans[d][i + 1][4] += e;
			}
			else if (context == std::string("Stem")) {
				ans[d][i + 1][5] += e;
			}
			else {
				assert(0);
			}
		}
	}

	for (int i = 0; i <= max_dim; ++i) {
		for (int j = 1; j <= n; ++j) {
			for (int k = 0; k < 6; ++k) {
				ans[i][j][k] /= z[i];
			}
		}
	}

	double z_sum = 0.0;
	for (int i = 0; i <= max_dim; ++i)z_sum += z[i];
	for (int i = 0; i <= max_dim; ++i)z[i] /= z_sum;

	return std::make_pair(z, ans);
}


}