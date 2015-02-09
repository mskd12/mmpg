#ifndef __VERTEX_HPP__
#define __VERTEX_HPP__

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

#include "CommonFile.hpp"

using namespace std;

namespace MeanPayoffGame {
  
  class Vertex {

  public:
    long long int vid;
    long long int weight;
    vector<int> mmpgWeights;
    int mmpgPlayer;
    PLAYER player;

    std::vector<long long int> succ_list; // List of successors
    std::vector<long long int> pred_list; // List of predecessors
    
    //Constructor
    Vertex();
    Vertex(long long int id, long long int pr, PLAYER pl);  
    Vertex(long long int id, vector<int> w, int mmpgplayer);

    /* OVerloading assignment and comparison operators */
    Vertex& operator=(const Vertex &v);
    bool operator==(Vertex v);
    
    void add_succ(long long int succ_id); // Add a successor to this vertex
    void delete_succ(long long int succ_id);
    void add_pred(long long int pred_id); // Add a predecessor to this vertex
    void delete_pred(long long int succ_id);

    void show();
  };

}


#endif // __VERTEX_HPP__
