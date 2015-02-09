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
#include "Vertex.hpp"

using namespace std;
using namespace MeanPayoffGame;

Vertex::Vertex(){}

Vertex::Vertex(long long int id, long long int pr, PLAYER pl) 
{ 
  this->vid = id; 
  this->weight = pr;
  this->player = pl;
}

Vertex::Vertex(long long int id, vector<int> w, int mmpgplayer) 
{ 
  this->vid = id; 
  for(vector<int> :: iterator it=w.begin(); it!=w.end(); it++) {
    this->mmpgWeights.push_back(*it);
  }
  //cout<<weight0<<" "<<weight1<<" "<<weight2<<endl;
  //cout<<w1<<" "<<w2<<" "<<w3<<endl;
  this->mmpgPlayer = mmpgplayer;
}

Vertex& Vertex::operator=(const Vertex &v)
{
  this->vid = v.vid; 
  this->weight = v.weight;
  this->player = v.player;  
  
  return *this;
}

bool Vertex::operator==(Vertex v)
{
  return ((vid == v.vid) && (weight == v.weight) && (player == v.player));
}

void Vertex::add_succ(long long int succ_id) 
{
  //cout << "There" << endl;
  this->succ_list.push_back(succ_id); 

} 

void Vertex::delete_succ(long long int succ_id)
{
  /*for(vector<long long int>::iterator it = succ_list.begin(); it!=succ_list.end(); it++)
  {
    if(*it == succ_id) 
    {
      it=succ_list.erase(it);
      return;
    }
  }*/
  //cout << "NoSuccessorFound" << endl;
  //throw  NoSuccessorFound;

  for(int i=0;i<succ_list.size();i++) {
    if(succ_list[i] == succ_id) {
      succ_list.erase(succ_list.begin()+i);
      //cout << "SuccessorFound1" << endl;
      return;
    }
  }
  throw  NoSuccessorFound;  
}
     
void Vertex::add_pred(long long int pred_id) 
{
  this->pred_list.push_back(pred_id); 
} 

void Vertex::delete_pred(long long int pred_id)
{
  for(vector<long long int>::iterator it = pred_list.begin(); it!=pred_list.end(); it++)
  {
    if(*it == pred_id) 
    {
      pred_list.erase(it);
      return;
    }
  }
  cout << "NoPredecessorFound" << endl;
  //throw  NoSuccessorFound;
}

void Vertex::show() {
  cout << this->vid << " " << this->weight << " " << this->player << endl;
}