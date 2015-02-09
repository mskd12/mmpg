#ifndef __ESTIMATION_HPP__
#define __ESTIMATION_HPP__

#include <assert.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <ctime>
#include <stack>

#include <limits.h>

#include "Value.hpp"

using namespace std;

namespace MeanPayoffGame {
	struct Estimation {
		map<long long int, Value> estimate;

		//Comstructor
		Estimation();
		//Estimation GetEstimate();

		bool operator==(Estimation);

		void AddEstimate(long long int id, Value val);
	};
}

#endif
