/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTP_INTERVAL_TYPE_H_
#define RINTP_INTERVAL_TYPE_H_

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

#endif//RINTP_INTERVAL_TYPE_H_
