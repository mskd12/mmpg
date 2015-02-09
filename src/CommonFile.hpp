
#ifndef __COMMON_FILE_HPP__
#define __COMMON_FILE_HPP__

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

#define PINF LLONG_MAX


using namespace std;

namespace MeanPayoffGame {
  
  enum PLAYER {even, odd};
  //enum STRAT_TYPE {s_even, s_odd, s_both, s_none};
  enum SymSIMExceptions {VertexNotFound, NoSuccessorFound, EndOfFile, FileNotFound};
  
  //  typedef vector<float> value;
  
  /*class SymSIMStat {
  public:
    std::set<int> even_winning_region;
    int iterations;
    double runningTimeInSecs;
  };*/
  

}

#endif // __COMMON_FILE_HPP__
