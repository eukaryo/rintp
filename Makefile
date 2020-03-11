rintp: centroid_fold.cpp complex_number.cpp experiment.cpp fourier_transform.cpp main.cpp mccaskill_1990.cpp misc.cpp parameter.cpp real_logsumexp.cpp rintcapr.cpp rintx.cpp sample_mccaskill.cpp test.cpp
	g++ -I ./ -std=c++1y -fopenmp -static-libstdc++ -O2 centroid_fold.cpp complex_number.cpp experiment.cpp fourier_transform.cpp main.cpp mccaskill_1990.cpp misc.cpp parameter.cpp real_logsumexp.cpp rintcapr.cpp rintx.cpp sample_mccaskill.cpp test.cpp -o rintp
clean:
	rm -f *.o rintp
