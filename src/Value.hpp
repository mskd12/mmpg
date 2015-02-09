#ifndef __VALUE_HPP__
#define __VALUE_HPP__

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

using namespace std;

namespace MeanPayoffGame {
	struct Value {
		long long int bestSum;

		/* Constructor */
		Value();
		Value(long long int n);

		bool operator==(Value v);
		Value operator+(const Value &v);
		Value& operator=(const Value &v);

	};
}


#endif /* __VALUE_HPP__ */
