/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"centroid_fold.h"

#include"interval_type.h"
#include"misc.h"
#include"parameter.h"

namespace rintp {

std::string GetCentroidFoldHamadaBook(
	const std::string& sequence,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop) {

	//����Ίm���s��̊e�v�f����͂Ƃ��āA��centroid�\�������߂ĕԂ��B
	//����bpp_matrix��(n+1)*(max_span+1)�s��Ƃ���B

	//�l�c�{��p140�ɂ��ƁAgamma < 1.0�̂Ƃ���2.(b)�̏����ł��悢���A(�������ʓ|�������̂�)�����ł͏��2.(a)���s���B

	//TODO:����̎��̗p�r�ɑ΂��ĕl�c�{�̕��@�͑S�R���߂��ƋC�t�����B
	//(1)max_loop�����max_span������l������Ȃ�
	//(2)�������\�\���́A���t�@�����X�\������̃n�~���O�������w��ł��Ȃ�
	//�̂ŁA���t�@�����X�\���Ƃ��Ă�Hagiotool�̉����Ƃ��Ă��g���Ȃ��B

	//�g���Ȃ��Ƃ������R�́AHagiotool�͏������l�̐��l�덷���������̂ŁA���A�ȃn�~���O�����̉���Ίm���s�񂪂߂��Ⴍ����Ȓl�ɂȂ邪
	//���̉���Ίm���s�񂩂��\�\�������߂悤�Ƃ���ƁAmax_loop�����max_span�������������Ȃ��\�������ۂ悭�o�邩��B

	//(2)�ɂ��Ă�RNAborMEA�̃A���S���Y���Ŏw��ł���B(1)�ɂ��Ă͂����炭McCaskill�^DP�Ńg���[�X�o�b�N����Ηǂ��B

	const int n = int(sequence.size());
	const int max_span = int(bpp_matrix[0].size()) - 1;

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(0.0 <= gamma && gamma <= 1000000000.0);
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(0.0 <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= 1.0);
	}

	const Floating G = (1.0 / (gamma + 1.0));

	//(�X�R�A, �g���[�X�o�b�N�p�̕ϐ�)
	typedef std::pair<Floating, int> Nus;

	//�������ċA�̂��߂̃f�[�^�\���ł���Bfirst�͏���������false�ŁA�l���v�Z������true�ɂ��āAsecond�ɒl������B
	typedef std::pair<bool, Nus> MemoNus;

	std::vector<std::vector<MemoNus>>MemoM(n + 1, std::vector<MemoNus>(n + 1, MemoNus(false, Nus(0.0, 0))));

	std::function<Nus(int, int)>Nussinov = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return std::make_pair(Floating(0.0), -9999);
		if (j == i - 1)return std::make_pair(Floating(0.0), -9999);
		if (MemoM[i][j].first)return MemoM[i][j].second;

		Floating max_score = -std::numeric_limits<Floating>::max();

		//�l�c�{p115�̃A���S���Y��3.3�̕ϐ�t�Ɠ����B
		//�������A�����l�͕l�c�{�ł�1�`3�̂Ƃ����-1�`-3�Ƃ��ɕς��Ă���B(���̂ق�������������)
		//�ł����������@�͑㐔�I�f�[�^�^���g�����̂����AC++�W���ɂ͖����B
		int trace = -9999;

		const auto update = [&](const Floating score, const int new_trace) {
			if (score > max_score) {
				max_score = score;
				trace = new_trace;
			}
		};

		update(Nussinov(i + 1, j).first, -1);//�l�c�{p115�̃A���S���Y��3.3��7�s��
		update(Nussinov(i, j - 1).first, -2);//�l�c�{p115�̃A���S���Y��3.3��9�s��

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0 && TURN + 1 <= (j - i)) {
			const Floating x = j - i <= max_span ? bpp_matrix[i][j - i] : Floating(0.0);
			update(Nussinov(i + 1, j - 1).first + x - G, -3);//�l�c�{p115�̃A���S���Y��3.3��12�s��
		}

		for (int k = i; k <= j - 1; ++k) {
			update(Nussinov(i, k).first + Nussinov(k + 1, j).first - G, k);//�l�c�{p115�̃A���S���Y��3.3��16�s��
		}

		return (MemoM[i][j] = std::make_pair(true, std::make_pair(max_score, trace))).second;
	};

	//��������͕l�c�{p116�̃A���S���Y��3.4�Ɠ������Ƃ�����B
	//�������A�l�c�{�ł̓X�^�b�N���g���Ă��邪�����ł͍ċA���g���B

	std::function<std::string(int, int)>Traceback = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);

		const int k = Nussinov(i, j).second;

		if (k == -9999) {
			std::string ans("");
			for (int x = i; x <= j; ++x)ans += std::string(".");
			return ans;
		}

		if (k == -1)return std::string(".") + Traceback(i + 1, j);
		if (k == -2)return Traceback(i, j - 1) + std::string(".");
		if (k == -3)return std::string("(") + Traceback(i + 1, j - 1) + std::string(")");
		return Traceback(i, k) + Traceback(k + 1, j);
	};

	return Traceback(1, n);
}

std::string GetCentroidFoldMcCaskill(
	const std::string& sequence,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop) {

	//max_span�����max_loop������l�����������ŁA
	//�S�Ă̓񎟍\���W���ɂ�����MEG�\�������߂ĕԂ��B

	//�l�c�{p135�̒藝3.6���g�����Ap136�̌n3.7�͎g��Ȃ��B
	//�Ȃ��Ȃ�A�n3.7�̏�����max_loop����ɖ������邩��ł���B

	const int n = int(sequence.size());

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(0.0 <= gamma && gamma <= 1000000000.0);
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(0.0 <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= 1.0);
	}

	const Floating G = (1.0 / (gamma + 1.0));
	const Floating inf = std::numeric_limits<Floating>::max() / Floating(1000000.0);

	//std::make_pair(�X�R�A, �g���[�X�o�b�N�p�̉ϒ��ϐ�)
	//�g���[�X�o�b�N�p�̕ϐ��Ƃ��đ㐔�I�f�[�^�^���g���������AC++�ɂ͖����B
	typedef std::pair<Floating, std::vector<int>> ScoreAndTracebackInfo;

	//�������ċA�̂��߂̃f�[�^�\���ł���Bfirst�͏���������false�ŁA�l���v�Z������true�ɂ��āAsecond�ɒl������B
	typedef std::pair<bool, ScoreAndTracebackInfo> MemoM;

	std::vector<MemoM>Z((n + 1), std::make_pair(false, std::make_pair(Floating(0.0), std::vector<int>{})));
	std::vector<MemoM>Z1((n + 1) * (max_span + 1), std::make_pair(false, std::make_pair(Floating(0.0), std::vector<int>{})));
	std::vector<MemoM>Zb((n + 1) * (max_span + 1), std::make_pair(false, std::make_pair(Floating(0.0), std::vector<int>{})));
	std::vector<MemoM>Zm((n + 1) * (max_span + 1), std::make_pair(false, std::make_pair(Floating(0.0), std::vector<int>{})));

	std::function<ScoreAndTracebackInfo(int)> ForwardZ;
	std::function<ScoreAndTracebackInfo(int, int)> ForwardZ1;
	std::function<ScoreAndTracebackInfo(int, int)> ForwardZb;
	std::function<ScoreAndTracebackInfo(int, int)> ForwardZm;

	std::function<std::string(int)> BackwardZ;
	std::function<std::string(int, int)> BackwardZ1;
	std::function<std::string(int, int)> BackwardZb;
	std::function<std::string(int, int)> BackwardZm;

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	ForwardZ = [&](const int j) {
		assert(0 <= j && j <= n);
		if (j <= 1)return std::make_pair(Floating(0.0), std::vector<int>{});
		if (Z[j].first)return Z[j].second;

		ScoreAndTracebackInfo ans = std::make_pair(Floating(0.0), std::vector<int>{});//0.0�͉���΂���ؑg�܂Ȃ��ꍇ�̃X�R�A

		for (int k = 0; k <= j - 1; ++k) {
			const Floating score = ForwardZ(k).first + ForwardZ1(k + 1, j).first;
			if (score > ans.first) ans = std::make_pair(score, std::vector<int>{k});
		}

		return (Z[j] = std::make_pair(true, ans)).second;
	};
	ForwardZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return std::make_pair(-inf, std::vector<int>{});

		if ((j - i) > max_span) {
			return ForwardZ1(i, i + max_span);
		}

		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		ScoreAndTracebackInfo ans = std::make_pair(-inf, std::vector<int>{});

		const Floating score1 = ForwardZb(i, j).first;
		if (score1 > ans.first) ans = std::make_pair(score1, std::vector<int>{1});

		const Floating score2 = ForwardZ1(i, j - 1).first;
		if (score2 > ans.first) ans = std::make_pair(score2, std::vector<int>{2});

		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	ForwardZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return std::make_pair(-inf, std::vector<int>{});

		//i��j������΂�g�ݓ��Ȃ��Ȃ�[����Ԃ��B
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > max_span) {
			return std::make_pair(-inf, std::vector<int>{});
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		ScoreAndTracebackInfo ans = std::make_pair(-inf, std::vector<int>{});

		//hairpin loop
		const Floating score1 = bpp_matrix[i][j - i] - G;
		if (score1 > ans.first) ans = std::make_pair(score1, std::vector<int>{1});

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				const Floating score2 = bpp_matrix[i][j - i] - G + ForwardZb(k, l).first;
				if (score2 > ans.first) ans = std::make_pair(score2, std::vector<int>{2, k, l});
			}
		}

		//multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			const Floating score3 = bpp_matrix[i][j - i] - G + ForwardZm(i + 1, k - 1).first + ForwardZ1(k, j - 1).first;
			if (score3 > ans.first) ans = std::make_pair(score3, std::vector<int>{3, k});
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	ForwardZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return std::make_pair(-inf, std::vector<int>{});
		if (j == i - 1)return std::make_pair(-inf, std::vector<int>{});
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		ScoreAndTracebackInfo ans = std::make_pair(-inf, std::vector<int>{});

		for (int k = i; k + TURN + 1 <= j; ++k) {
			const Floating score1 = ForwardZ1(k, j).first;
			if (score1 > ans.first) ans = std::make_pair(score1, std::vector<int>{1, k});

			const Floating score2 = ForwardZm(i, k - 1).first + ForwardZ1(k, j).first;
			if (score2 > ans.first) ans = std::make_pair(score2, std::vector<int>{2, k});
		}

		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};

	const auto EmptyStructure = [](const int i, const int j) {
		std::string ans("");
		for (int k = i; k <= j; ++k)ans += std::string(".");
		return ans;
	};

	BackwardZ = [&](const int j) {
		assert(0 <= j && j <= n);
		if (j == 1)return std::string(".");
		if (j == 0)return std::string("");

		const std::vector<int>trace = ForwardZ(j).second;
		if (trace.empty())return EmptyStructure(1, j);

		const int k = trace[0];
		return BackwardZ(k) + BackwardZ1(k + 1, j);
	};
	BackwardZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return std::string(".");

		if ((j - i) > max_span) {
			return BackwardZ1(i, i + max_span) + EmptyStructure(i + max_span + 1, j);
		}

		const std::vector<int>trace = ForwardZ1(i, j).second;
		if (trace.empty())return EmptyStructure(i, j);

		return trace[0] == 1 ? BackwardZb(i, j) : (BackwardZ1(i, j - 1) + std::string("."));
	};
	BackwardZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return std::string(".");

		//i��j������΂�g�ݓ��Ȃ��Ȃ�[����Ԃ��B
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > max_span) {
			return EmptyStructure(i, j);
		}

		const std::vector<int>trace = ForwardZb(i, j).second;
		if (trace.empty())return EmptyStructure(i, j);

		//hairpin loop
		if (trace[0] == 1) return std::string("(") + EmptyStructure(i + 1, j - 1) + std::string(")");

		//internal loop, stem, bulge
		if (trace[0] == 2) {
			const int k = trace[1];
			const int l = trace[2];
			return std::string("(") + EmptyStructure(i + 1, k - 1) + BackwardZb(k, l) + EmptyStructure(l + 1, j - 1) + std::string(")");
		}

		//multi loop
		if (trace[0] == 3) {
			const int k = trace[1];
			return std::string("(") + BackwardZm(i + 1, k - 1) + BackwardZ1(k, j - 1) + std::string(")");
		}

		assert(0);
		return EmptyStructure(i, j);
	};
	BackwardZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return std::string(".");
		if (j == i - 1)return std::string("");

		const std::vector<int>trace = ForwardZm(i, j).second;
		if (trace.empty())return EmptyStructure(i, j);

		const int k = trace[1];
		return (trace[0] == 1 ? EmptyStructure(i, k - 1) : BackwardZm(i, k - 1)) + BackwardZ1(k, j);
	};

	return BackwardZ(n);
}

std::vector<std::string> GetCentroidFoldForEachHammingDistance(
	const std::string& sequence,
	const std::vector<std::vector<int>>& S,
	const int max_dim,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop,
	const int distance) {

	//����Ίm���s��̊e�v�f����͂Ƃ��āA��centroid�\�������߂ĕԂ��B
	//����bpp_matrix��(n+1)*(max_span+1)�s��Ƃ���B

	//�ő�X�p������ƍő僋�[�v������l�����Čv�Z����B
	//�܂��A�S�Ă�k�ɂ��āA���t�@�����X�\������̃n�~���O������k�ł���悤�ȍ\���̒��ł�MEG�\�������߂�B
	//������distance!=-1�̂Ƃ������̓��t�@�����X�\������̃n�~���O������distance�ł���悤�ȍ\���̒��ł�MEG�\�����������߂�B

	const int n = int(sequence.size());

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(0.0 <= gamma && gamma <= 1000000000.0);
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(0.0 <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= 1.0);
	}

	//[Mori et al., 2014]��supp�̎�(S10)�̍��ӂ�C�ł���B
	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(S);

	//[Mori et al., 2014]��supp�̎�(S9)�A(S14)�`S(21)�Œ�`����Ă���g�ł���B
	const auto g1 = [&](const int i, const int j) {return C[i][j]; };
	const auto g2 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto g3 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k]; };
	const auto g4 = [&](const int i, const int j) {return C[i][j] + 1 - 2 * S[i][j]; };
	const auto g5 = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * S[i][j]; };
	const auto g6 = [&](const int i, const int j, const int k) {return C[i][j] - C[i + 1][k - 1] - C[k][j - 1] + 1 - 2 * S[i][j]; };
	const auto g7 = [&](const int i, const int j, const int k) {return C[i][j] - C[k][j]; };
	const auto g8 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k - 1] - C[k][j]; };

	const Floating G = (Floating(1.0) / (gamma + Floating(1.0)));
	const Floating inf = std::numeric_limits<Floating>::max() / Floating(1000000.0);

	//std::make_pair(�X�R�A, �g���[�X�o�b�N�p�̉ϒ��ϐ�)
	//�g���[�X�o�b�N�p�̕ϐ��Ƃ��đ㐔�I�f�[�^�^���g���������AC++�ɂ͖����B
	typedef struct infos {
		Floating score;
		int tag[3];
		infos() {
			score = Floating(0.0);
			for (int i = 0; i < 3; ++i)tag[i] = -1;
		}
		infos(Floating s) {
			score = s;
			for (int i = 0; i < 3; ++i)tag[i] = -1;
		}
		infos(Floating s, int i1) {
			score = s;
			tag[0] = i1;
			tag[1] = -1;
			tag[2] = -1;
		}
		infos(Floating s, int i1, int i2) {
			score = s;
			tag[0] = i1;
			tag[1] = i2;
			tag[2] = -1;
		}
		infos(Floating s, int i1, int i2, int i3) {
			score = s;
			tag[0] = i1;
			tag[1] = i2;
			tag[2] = i3;
		}
		void set(Floating s, int i1) {
			score = s;
			tag[0] = i1;
			tag[1] = -1;
			tag[2] = -1;
		}
		void set(Floating s, int i1, int i2) {
			score = s;
			tag[0] = i1;
			tag[1] = i2;
			tag[2] = -1;
		}
		void set(Floating s, int i1, int i2, int i3) {
			score = s;
			tag[0] = i1;
			tag[1] = i2;
			tag[2] = i3;
		}
		std::vector<int>get() {
			std::vector<int>ans;
			for (int i = 0; i < 3 && tag[i] != -1; ++i)ans.push_back(tag[i]);
			return ans;
		}
	}ScoreAndTracebackInfo;

	//�������ċA�̂��߂̃f�[�^�\���ł���Bfirst�͏���������false�ŁA�l���v�Z������true�ɂ��āAsecond�ɒl������B
	typedef std::pair<bool, ScoreAndTracebackInfo> MemoM;

	std::vector<std::vector<MemoM>>Z((distance == -1 ? max_dim : distance) + 1, std::vector<MemoM>((n + 1)));
	std::vector<std::vector<MemoM>>Z1((distance == -1 ? max_dim : distance) + 1, std::vector<MemoM>((n + 1) * (max_span + 1)));
	std::vector<std::vector<MemoM>>Zb((distance == -1 ? max_dim : distance) + 1, std::vector<MemoM>((n + 1) * (max_span + 1)));
	std::vector<std::vector<MemoM>>Zm((distance == -1 ? max_dim : distance) + 1, std::vector<MemoM>((n + 1) * (max_span + 1)));

	std::function<ScoreAndTracebackInfo(int, int)> ForwardZ;
	std::function<ScoreAndTracebackInfo(int, int, int)> ForwardZ1;
	std::function<ScoreAndTracebackInfo(int, int, int)> ForwardZb;
	std::function<ScoreAndTracebackInfo(int, int, int)> ForwardZm;

	std::function<std::string(int, int)> BackwardZ;
	std::function<std::string(int, int, int)> BackwardZ1;
	std::function<std::string(int, int, int)> BackwardZb;
	std::function<std::string(int, int, int)> BackwardZm;

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	ForwardZ = [&](const int j, const int x) {
		assert(0 <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (j <= 1)return ScoreAndTracebackInfo(x == 0 ? Floating(0.0) : -inf);
		if (Z[x][j].first)return Z[x][j].second;

		ScoreAndTracebackInfo ans = ScoreAndTracebackInfo(-inf);
		if (x == g1(1, j))ans.score = Floating(0.0);

		for (int k = 0; k <= j - 1; ++k) {
			const int y = x - g2(1, j, k);
			for (int z = 0; z <= y; ++z) {
				const Floating score = ForwardZ(k, z).score + ForwardZ1(k + 1, j, y - z).score;
				if (score > ans.score) ans.set(score, k, z);
			}
		}

		return (Z[x][j] = std::make_pair(true, ans)).second;
	};
	ForwardZ1 = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return ScoreAndTracebackInfo(-inf);

		if ((j - i) > max_span) {
			const int y = x - g3(i, j, i + max_span);
			return y >= 0 ? ForwardZ1(i, i + max_span, y) : ScoreAndTracebackInfo(-inf);
		}

		if (Z1[x][at(i, j)].first)return Z1[x][at(i, j)].second;

		ScoreAndTracebackInfo ans(-inf);

		const Floating score1 = ForwardZb(i, j, x).score;
		if (score1 > ans.score) ans.set(score1, 1);

		const int y = x - g3(i, j, j - 1);
		if (y >= 0) {
			const Floating score2 = ForwardZ1(i, j - 1, y).score;
			if (score2 > ans.score) ans.set(score2, 2);
		}

		return (Z1[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	ForwardZb = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return ScoreAndTracebackInfo(-inf);

		//i��j������΂�g�ݓ��Ȃ��Ȃ�[����Ԃ��B
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > max_span) {
			return ScoreAndTracebackInfo(-inf);
		}

		if (Zb[x][at(i, j)].first)return Zb[x][at(i, j)].second;

		ScoreAndTracebackInfo ans(-inf);

		//hairpin loop
		if (x == g4(i, j)) {
			const Floating score = bpp_matrix[i][j - i] - G;
			if (score > ans.score) ans.set(score, 1);
		}

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				const int y = x - g5(i, j, k, l);
				if (y >= 0) {
					const Floating score = bpp_matrix[i][j - i] - G + ForwardZb(k, l, y).score;
					if (score > ans.score) ans.set(score, 2, k, l);
				}
			}
		}

		//multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			const int y = x - g6(i, j, k);
			for (int z = 0; z <= y; ++z) {
				const Floating score = bpp_matrix[i][j - i] - G + ForwardZm(i + 1, k - 1, z).score + ForwardZ1(k, j - 1, y - z).score;
				if (score > ans.score) ans.set(score, 3, k, z);
			}
		}

		return (Zb[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	ForwardZm = [&](const int i, const int j, const int x) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return ScoreAndTracebackInfo(-inf);
		if (j == i - 1)return ScoreAndTracebackInfo(-inf);
		if (Zm[x][at(i, j)].first)return Zm[x][at(i, j)].second;

		ScoreAndTracebackInfo ans(-inf);

		for (int k = i; k + TURN + 1 <= j; ++k) {
			const int y1 = x - g7(i, j, k);
			if (y1 >= 0) {
				const Floating score = ForwardZ1(k, j, y1).score;
				if (score > ans.score) ans.set(score, 1, k);
			}
			const int y2 = x - g8(i, j, k);
			for (int z = 0; z <= y2; ++z) {
				const Floating score = ForwardZm(i, k - 1, z).score + ForwardZ1(k, j, y2 - z).score;
				if (score > ans.score) ans.set(score, 2, k, z);
			}
		}

		return (Zm[x][at(i, j)] = std::make_pair(true, ans)).second;
	};

	const auto EmptyStructure = [](const int i, const int j) {
		std::string ans("");
		for (int k = i; k <= j; ++k)ans += std::string(".");
		return ans;
	};

	BackwardZ = [&](const int j, const int x) {
		assert(0 <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (j == 1)return std::string(".");
		if (j == 0)return std::string("");

		const std::vector<int>trace = ForwardZ(j, x).get();
		if (trace.empty())return EmptyStructure(1, j);

		const int k = trace[0];
		const int z = trace[1];
		const int y = x - g2(1, j, k);
		return BackwardZ(k, z) + BackwardZ1(k + 1, j, y - z);
	};
	BackwardZ1 = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return std::string(".");

		if ((j - i) > max_span) {
			const int y = x - g3(i, j, i + max_span);
			return BackwardZ1(i, i + max_span, y) + EmptyStructure(i + max_span + 1, j);
		}

		const std::vector<int>trace = ForwardZ1(i, j, x).get();
		if (trace.empty())return EmptyStructure(i, j);

		assert(trace[0] == 1 || trace[0] == 2);
		return trace[0] == 1 ? BackwardZb(i, j, x) : BackwardZ1(i, j - 1, x - g3(i, j, j - 1)) + std::string(".");
	};
	BackwardZb = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return std::string(".");

		//i��j������΂�g�ݓ��Ȃ��Ȃ�[����Ԃ��B
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > max_span) {
			return EmptyStructure(i, j);
		}

		const std::vector<int>trace = ForwardZb(i, j, x).get();
		if (trace.empty())return EmptyStructure(i, j);

		//hairpin loop
		if (trace[0] == 1) return std::string("(") + EmptyStructure(i + 1, j - 1) + std::string(")");

		//internal loop, stem, bulge
		if (trace[0] == 2) {
			const int k = trace[1];
			const int l = trace[2];
			const int y = x - g5(i, j, k, l);
			return std::string("(") + EmptyStructure(i + 1, k - 1) + BackwardZb(k, l, y) + EmptyStructure(l + 1, j - 1) + std::string(")");
		}

		//multi loop
		if (trace[0] == 3) {
			const int k = trace[1];
			const int z = trace[2];
			const int y = x - g6(i, j, k);
			return std::string("(") + BackwardZm(i + 1, k - 1, z) + BackwardZ1(k, j - 1, y - z) + std::string(")");
		}

		assert(0);
		return EmptyStructure(i, j);
	};
	BackwardZm = [&](const int i, const int j, const int x) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		assert(0 <= x && x <= (distance == -1 ? max_dim : distance));
		if (i == j)return std::string(".");
		if (j == i - 1)return std::string("");

		const std::vector<int>trace = ForwardZm(i, j, x).get();
		if (trace.empty())return EmptyStructure(i, j);

		if (trace[0] == 1) {
			const int k = trace[1];
			const int y1 = x - g7(i, j, k);
			return EmptyStructure(i, k - 1) + BackwardZ1(k, j, y1);
		}

		if (trace[0] == 2) {
			const int k = trace[1];
			const int z = trace[2];
			const int y2 = x - g8(i, j, k);
			return BackwardZm(i, k - 1, z) + BackwardZ1(k, j, y2 - z);
		}

		assert(0);
		return EmptyStructure(i, j);
	};

	std::vector<std::string>ans;

	if (distance != -1) {
		assert(0 <= distance && distance <= max_dim);
		ans.push_back(BackwardZ(n, distance));
	}
	else {
		for (int x = 0; x <= max_dim; ++x) {
			ans.push_back(BackwardZ(n, x));
		}
	}

	return ans;
}

Floating ComputeGain(
	const std::string& structure,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop) {

	//����Ίm���s��bpp_matrix�ɑ΂��āA�񎟍\��structure���\�\���Ƃ����Ƃ���Gain�����߂ĕԂ��B
	//�l�c�{p139�̎�(3.28)���v�Z����B

	const int n = int(structure.size());
	const int max_span = int(bpp_matrix[0].size()) - 1;

	//std::cout << n << std::endl;
	//std::cout << max_loop << std::endl;
	//for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n; ++j)std::cout << i << " " << j << " " << bpp_matrix[i][j - i] << std::endl;
	//std::cout << gamma << std::endl;

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(Floating(0.0) <= gamma && gamma <= Floating(1000000000.0));
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(Floating(0.0) <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= Floating(1.0));
	}

	std::vector<std::vector<int>>pairing(n + 1, std::vector<int>(n + 1, 0));
	std::stack<int> bp_pos;
	for (int i = 1; i <= n; ++i) {
		switch (structure[i - 1]) {
		case '(':
			bp_pos.push(i);
			break;
		case ')':
			pairing[bp_pos.top()][i] = 1;
			bp_pos.pop();
			break;
		case '.':
			break;
		default:
			assert(0);
			break;
		}
	}

	Floating ans = Floating(0.0);

	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		if (pairing[i][j] == 0) {
			ans += 1.0 - bpp_matrix[i][j - i];
		}
		else {
			ans += gamma * bpp_matrix[i][j - i];
		}
	}

	return ans;
}

std::pair<std::string, Floating> BruteForceGetCentroidFoldMcCaskill(
	const std::string& sequence,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop) {
	//max_span�����max_loop������l�����������ŁA
	//�S�Ă̓񎟍\���W���ɂ�����MEG�\�������߂ĕԂ��B

	const int n = int(sequence.size());

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(0.0 <= gamma && gamma <= 1000000000.0);
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(0.0 <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= 1.0);
	}

	const std::vector<std::string> structures = EnumerateStructures(sequence, max_span, max_loop);

	std::pair<std::string, Floating> ans = std::make_pair(std::string(""), -std::numeric_limits<Floating>::max());
	for (const std::string s : structures) {
		const Floating g = ComputeGain(s, bpp_matrix, gamma, max_loop);
		if (g > ans.second)ans = std::make_pair(s, g);
	}
	return ans;
}

std::vector<std::pair<std::string, Floating>> BruteForceGetCentroidFoldForEachHammingDistance(
	const std::string& sequence,
	const std::string& reference_structure,
	const int max_dim,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop) {
	//����Ίm���s��̊e�v�f����͂Ƃ��āA��centroid�\�������߂ĕԂ��B
	//����bpp_matrix��(n+1)*(max_span+1)�s��Ƃ���B

	//�ő�X�p������ƍő僋�[�v������l�����Čv�Z����B
	//�܂��A�S�Ă�k�ɂ��āA���t�@�����X�\������̃n�~���O������k�ł���悤�ȍ\���̒��ł�MEG�\�������߂�B

	const int n = int(sequence.size());

	assert(int(bpp_matrix.size()) == n + 1);
	for (int i = 0; i <= n; ++i)assert(int(bpp_matrix[i].size()) == max_span + 1);
	assert(0.0 <= gamma && gamma <= 1000000000.0);
	for (int i = 1; i <= n; ++i)for (int j = i + TURN + 1; j <= n && j - i <= max_span; ++j) {
		assert(0.0 <= bpp_matrix[i][j - i] && bpp_matrix[i][j - i] <= 1.0);
	}

	const std::vector<std::string> structures = EnumerateStructures(sequence, max_span, max_loop);

	std::vector<std::pair<std::string, Floating>>ans(max_dim + 1, std::make_pair(std::string(""), -std::numeric_limits<Floating>::max()));

	for (const std::string s : structures) {
		const int d = ComputeHammingDistance(s, reference_structure);
		const Floating g = ComputeGain(s, bpp_matrix, gamma, max_loop);
		if (g > ans[d].second)ans[d] = std::make_pair(s, g);
	}
	return ans;
}

}
