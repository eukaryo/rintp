/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef RINTP_TEST_H_
#define RINTP_TEST_H_

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

namespace rintp {

void TestMaxHamming(const int num);
void TestSampling(const int num);
void TestSimpleMcCaskillWide(const int num);
void TestRintD1Dim(const int num);
void TestRintD2Dim(const int num);
void TestRintW1Dim(const int num);
void TestRintW2Dim(const int num);
void TestCentroidFold(const int num);
void TestHagioNonFourier(const int num);
void TestMaxHammingPK(const int num);
void TestRintD1DimPK(const int num);

void TestRintP(const int num);

}

#endif//RINTP_TEST_H_
