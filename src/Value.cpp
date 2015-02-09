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
#include "CommonFile.hpp"

using namespace std;
using namespace MeanPayoffGame;


Value::Value() 
{
	bestSum = 0;
}

Value::Value(long long int n)
{
	bestSum = n;
}

bool Value::operator==(Value v) {
	return (bestSum == v.bestSum);
}

Value& Value::operator=(const Value &v) {
	bestSum = v.bestSum;

	return *this;
}

Value Value::operator+(const Value &v) {
	Value ret;
	if(this->bestSum == PINF || v.bestSum == PINF) {ret.bestSum = PINF;}
	else {ret.bestSum = this->bestSum + v.bestSum;}
	return ret;
}
