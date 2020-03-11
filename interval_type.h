/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef RINTDWR_INTERVAL_TYPE_H_
#define RINTDWR_INTERVAL_TYPE_H_

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

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

//#include <kv/dd.hpp>
//#include <kv/rdd.hpp>


namespace rintp {

//typedef kv::dd Floating;
typedef double Floating;
typedef kv::interval<Floating> IntervalVar;

}

#endif//RINTDWR_INTERVAL_TYPE_H_
