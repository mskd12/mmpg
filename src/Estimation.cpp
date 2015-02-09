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


#include "Estimation.hpp"
#include "Value.hpp"

using namespace std;
using namespace MeanPayoffGame;

Estimation::Estimation(){}

bool Estimation::operator==(Estimation e) {
  std::map<long long int, Value> lhs, rhs;
  lhs = estimate;
  rhs = e.estimate;
  
  if (lhs.size() != rhs.size()) return false; 
  
  map<long long int, Value>::const_iterator i, j;
  for(i = lhs.begin(), j = rhs.begin(); i != lhs.end(); ++i, ++j)
    {
      Value v1 = i-> second;
      Value v2 = j-> second; 

      if (i->first != j->first)  return false;
      //      if (i->second != j->second)  return false; 
      if (!(v1 == v2))  return false; 
    }

  return true;
}

void Estimation::AddEstimate(long long int id, Value val) {
  if(estimate.count(id)==0) {
    estimate[id]=val;
  }
  else {
    cout << "Estimation.cpp: Estimate Already Present" << endl;
  }
}
