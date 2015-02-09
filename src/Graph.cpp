
#include "Graph.hpp" 

#include <bits/stdc++.h>
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
using namespace MeanPayoffGame;

Graph::Graph() {
  n_vertices = 0;
}

void Graph::add_vertex(Vertex v) 
{
  this->vertices.push_back(v);
  //this->IdToVertex.insert(pair<int, Vertex*>(v.vid, &v));
  if (v.player == even) {
      this->EvenVertices.push_back(n_vertices);
  }
  else {
    this->OddVertices.push_back(n_vertices);
  }
  
  n_vertices++; 
}

void Graph::add_vertex(long long int id, long long int priority, PLAYER player)
{
  Vertex v(id, priority, player);
  
  this->vertices.push_back(v);
  //this->IdToVertex.insert(pair<long long int, Vertex*>(v.vid, &v));
  if (v.player == even) {
    this->EvenVertices.push_back(n_vertices);
  }
  else {
    this->OddVertices.push_back(n_vertices);
  }
  
  n_vertices++; 
  //cout << "Out add_vertex" << endl;
}

void Graph::set_successor(long long int id, long long int succ)
{
  get_vertex(id).add_succ(succ);
  //cout << get_vertex(id).succ_list.size() << endl;
}

void Graph::fill_predecessors()
{
  vector<Vertex>::iterator i;

  // fill predecessors list
  for (i = vertices.begin(); i != vertices.end(); i++) {
    long long int vertex_id = (*i).vid;
    
    vector<long long int> succ_list = (*i).succ_list;
    vector<long long int>::iterator j = succ_list.begin();
    
    while (j != succ_list.end()) {
      get_vertex(*j).add_pred(vertex_id);
      j++;
    }
  }
}


Vertex& Graph::get_vertex(const long long int id)
{
  
  vector<Vertex>::iterator i;
  
  for (i = vertices.begin(); i != vertices.end(); i++) {
    if ((*i).vid == id) return (*i);
  }
  //if(IdToVertex.count(id)>0) return *IdToVertex[id];
  cout << "Error" << endl;
  throw VertexNotFound; 
}


Estimation Graph::UpdateEstimation(Estimation currentEstimate) {
  
  // Calculating Improvement Arena
  Graph ImprovementArena;
  //map<long long int, Vertex> IdToVertex;
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    Vertex temp(curr.vid, curr.weight, curr.player);
    ImprovementArena.add_vertex(temp);
    //IdToVertex.insert(pair<long long int,Vertex>(curr.vid, ImprovementArena.get_vertex(curr.vid)));   
  }
  
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    long long int CurrId = curr.vid;
    
    long long int count = 0;
    for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
      long long int SuccId = curr.succ_list[j];
      Vertex Succ = get_vertex(SuccId);
      if(currentEstimate.estimate[SuccId].bestSum + Succ.weight 
   >= currentEstimate.estimate[CurrId].bestSum || 
   currentEstimate.estimate[SuccId].bestSum==PINF) {
  count++;
  ImprovementArena.set_successor(CurrId, SuccId);
      }
    }
    //if(curr.player==odd) assert(count==(long long int)curr.succ_list.size());
  }
  ImprovementArena.fill_predecessors();
  
  //cout << "Improvement Arena" << endl;
  //ImprovementArena.show();
  
  //Update Calculation from Improvement Arena
  map<long long int, Value> update;
  map<long long int, bool> isEvaluated;
  Value tempValue;
  Value zero(0);

  update[0] = zero;
  isEvaluated[0] = true;
  //bool special = false;
  int CaseNo = 3;
  
  //vector<long long int> candidates = ImprovementArena.get_vertex(0).pred_list;
  while(true) { // Calculating Update
    //cout << "CaseNo" << CaseNo << endl;
    long long int size = vertices.size();
    bool ShouldExit = true;
    for(long long int i=0; i<size; i++) {
      Vertex curr = ImprovementArena.get_vertex(vertices[i].vid);
      long long int CurrId = curr.vid;
      //cout << "Vid " << CurrId << " " << isEvaluated[CurrId] << endl;
      if(!isEvaluated[CurrId]) ShouldExit = false;
    }
    
    if(ShouldExit) break;
    
    //vector<long long int> EvaluatedCandidates;
    if(CaseNo==1 || CaseNo==2) {
      //cout << "in" << CaseNo << endl;
      CaseNo = -1;
      for(long long int i=0; i<size; i++) {
  Vertex curr = ImprovementArena.get_vertex(vertices[i].vid);
  bool AllEvaluated = true; // Are all Successors Evaluated
  long long int CurrId = curr.vid;
  long long int min = -1;
  bool flag = true;
  if(isEvaluated[CurrId]) continue;
  if(curr.player == odd) {
    //cout << "CurrId : " << CurrId << endl;
    assert((long long int)curr.succ_list.size()!=0);
    for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
      long long int SuccId = curr.succ_list[j];
      //cout << "CurrId : " << CurrId <<  " SuccId : " << SuccId << endl;
      if(!isEvaluated[SuccId]) {
        AllEvaluated = false;
        //break;
      }
      else {
        long long int impr = currentEstimate.estimate[SuccId].bestSum + 
    get_vertex(SuccId).weight - currentEstimate.estimate[CurrId].bestSum;
        if(currentEstimate.estimate[SuccId].bestSum == PINF) {
    if(currentEstimate.estimate[CurrId].bestSum == PINF) {
      impr = 0;
    }
    else {
      impr = PINF;
    }
        }
        else if(currentEstimate.estimate[CurrId].bestSum == PINF) {
    assert(false);
        }
        
        long long int newVal;
        if(update[SuccId].bestSum==PINF || impr==PINF) {
    newVal = PINF;
        }
        else newVal = update[SuccId].bestSum + impr;
        
        /*cout << CurrId << " " << SuccId << endl;
    cout << currentEstimate.estimate[SuccId].bestSum << endl;
    cout << get_vertex(SuccId).weight << endl;
    cout << currentEstimate.estimate[CurrId].bestSum << endl;
        */
        //cout << "impr " << impr << "\nnewVal "  << newVal << endl;
        if(update[SuccId].bestSum == 0 && impr == 0) CaseNo=2;
        if(flag) {min=newVal; flag = false;}
        else if(newVal<min) {min=newVal;} 
      } 
    }
    //cout << "AllEvaluated " << AllEvaluated << endl;
    //cout << "CaseNo " << CaseNo << endl;
    if(AllEvaluated) {
      if(flag) {
        //cout << "Check" << endl;
        continue;
      }
      else {
        tempValue.bestSum = min;
        update[CurrId] = tempValue;
        isEvaluated[CurrId] = true;
        CaseNo = 1;
        //EvaluatedCandidates.push_back(CurrId);
        //cout << "Updated1 - " << CurrId << " " << min << endl;
        
      }
    }
    else if(CaseNo==2) {// Case 2
      update[CurrId] = zero;
      isEvaluated[CurrId] = true;
      //EvaluatedCandidates.push_back(CurrId);
      CaseNo = 1;
      //cout << "Updated2 - " << CurrId << " " << 0 << endl;
    }
  }
      }
      //if(CaseNo == 1) continue;
      CaseNo=3;
      continue;
    }
    else if(CaseNo==3) {
      //cout << "in" << CaseNo << endl;
      CaseNo = -1;
      for(long long int i=0; i<size; i++) {
  Vertex curr = ImprovementArena.get_vertex(vertices[i].vid);
  bool AllEvaluated = true; // Are all Successors Evaluated
  long long int CurrId = curr.vid;
  if(isEvaluated[CurrId]) continue;
  if(curr.player == even) {
    long long int max = -1;
    bool flag = true;
    for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
      long long int SuccId = curr.succ_list[j];
      //cout << "CurrId : " << CurrId << "SuccId : " << SuccId << endl;
      if(!isEvaluated[SuccId]) {
        AllEvaluated = false;
        break;
      }
      else {
        long long int impr = currentEstimate.estimate[SuccId].bestSum + 
    get_vertex(SuccId).weight - currentEstimate.estimate[CurrId].bestSum;
        
        if(currentEstimate.estimate[SuccId].bestSum == PINF) {
    if(currentEstimate.estimate[CurrId].bestSum == PINF) {
      impr = 0;
    }
    else {
      impr = PINF;
    }
        }
        else if(currentEstimate.estimate[CurrId].bestSum == PINF) {
    assert(false);
        }
        
        long long int newVal;
        if(update[SuccId].bestSum==PINF || impr==PINF) {
    newVal = PINF;
        }
        else newVal = update[SuccId].bestSum + impr;
        
        //cout << impr << " " << update[SuccId].bestSum << " " << newVal << endl;
        if(flag) {max = newVal; flag = false;}
        else {if(newVal>max) {max = newVal;}}
        //cout << "max " << max << endl;
      }
    }
    if(AllEvaluated) {
      tempValue.bestSum = max;
      update[CurrId] = tempValue;
      isEvaluated[CurrId] = true;
      CaseNo = 1;
      //cout << "Updated3 - " << CurrId << " " << tempValue.bestSum << endl;
    }
  }
      }
      
      if(CaseNo==1) continue;
      CaseNo = 4;
      continue;
    }
    else if(CaseNo==4) {
      //cout << "in" << CaseNo << endl;
      CaseNo = -1;
      //std::vector<long long int> ;
      map<long long int, bool> MinCandidates;
      vector<long long int> InfCandidates;
      map<long long int, long long int> InterImp;
      long long int MinImp = -1;
      bool flag1 = true;

      for(long long int i=0; i<size; i++) {
        Vertex curr = ImprovementArena.get_vertex(vertices[i].vid);
        bool AllEvaluated = true; // Are all Successors Evaluated
        long long int CurrId = curr.vid;
        long long int min = -1;
        bool flag = true;
        if(isEvaluated[CurrId]) continue;
        if(curr.player == odd) {
          for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
            long long int SuccId = curr.succ_list[j];
            //cout << "SuccId : " << SuccId << endl;
            if(!isEvaluated[SuccId]) {
              AllEvaluated = false;
              //break;
            }
            else {
              long long int impr = currentEstimate.estimate[SuccId].bestSum + 
          get_vertex(SuccId).weight - currentEstimate.estimate[CurrId].bestSum;
              if(currentEstimate.estimate[SuccId].bestSum == PINF) {
          if(currentEstimate.estimate[CurrId].bestSum == PINF) {
            impr = 0;
          }
          else {
            impr = PINF;
          }
              }
              else if(currentEstimate.estimate[CurrId].bestSum == PINF) {
          assert(false);
              }
              
              long long int newVal;
              if(update[SuccId].bestSum==PINF || impr==PINF) {
          newVal = PINF;
              }
              else newVal = update[SuccId].bestSum + impr;
              /*cout << CurrId << " " << SuccId << endl;
          cout << currentEstimate.estimate[SuccId].bestSum << endl;
          cout << get_vertex(SuccId).weight << endl;
          cout << currentEstimate.estimate[CurrId].bestSum << endl;
          cout << "impr " << impr << "\nnewVal "  << newVal << endl;
              */
              //if(update[SuccId].bestSum == 0 && impr == 0) CaseNo=2;
              if(flag) {min=newVal; flag = false;}
              else if(newVal<min) {min=newVal;} 
            }
          }
          if(AllEvaluated==true && curr.succ_list.size()!=0) {
            assert(false);
          }
          if(flag) {
            InfCandidates.push_back(CurrId);
              //tempValue.bestSum = PINF;
          }
          else {
            MinCandidates[CurrId] = true;
            InterImp[CurrId] = min;
            if(flag1) {
              MinImp = min;
              flag1 = false;
            }
            else if(min<MinImp) {
              MinImp = min;
            }
            //tempValue.bestSum = min;
          }
          //update[CurrId] = tempValue;
          //isEvaluated[CurrId] = true;
          //CaseNo = 1;
          //cout << "Updated4 - " << CurrId << " " << tempValue.bestSum << endl;
          //EvaluatedCandidates.push_back(CurrId);
        }
      }

      for(long long int i=0; i<size; i++) {
        Vertex curr = ImprovementArena.get_vertex(vertices[i].vid);
        long long int CurrId = curr.vid;
        if(MinCandidates[CurrId]) {
          if(InterImp[CurrId] == MinImp) {
            tempValue.bestSum = MinImp;
            update[CurrId] = tempValue;
            isEvaluated[CurrId] = true;
            CaseNo = 1;
            //cout << "Updated4 - " << CurrId << " " << tempValue.bestSum << endl;
          }
        }
      }
      
      if(CaseNo == -1) {
        for(int i=0; i<(int)InfCandidates.size(); i++) {
          long long int CurrId = InfCandidates[i];
          tempValue.bestSum = PINF;
          update[CurrId] = tempValue;
          isEvaluated[CurrId] = true;
          CaseNo = 1;
          //cout << "Updated4 - " << CurrId << " " << tempValue.bestSum << endl;          
        }
      }
      //assert(CaseNo!=-1);
      
    } 
  }
  Estimation ret;
  //cout << "Calculating New Estimation" << endl;
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    long long int CurrId = vertices[i].vid;
    if(isEvaluated[CurrId]) {
      //cout << update[CurrId].bestSum << " " << currentEstimate.estimate[CurrId].bestSum << endl;
      ret.estimate[CurrId] = update[CurrId] + currentEstimate.estimate[CurrId];
    }
    else {
      assert(false);
      ret.estimate[CurrId] = Value(PINF); 
      //cout << "What to do?" << endl;
    }
  }
  return ret;
}


vector<long long int> Graph::OptimalStrategyImprovement() {
  Estimation InitialEstimate, Prev, Next;
  Value zero(0);
  
  // Creating the first estimate
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex Curr = vertices[i];
    long long int CurrId = Curr.vid;
    if(Curr.player==even) {
      InitialEstimate.AddEstimate(CurrId, zero);
    }
    else {
      if(CurrId == 0) {  /* Why? */
  InitialEstimate.AddEstimate(0,zero);
  continue;
      }
      
      long long int min = -1;
      bool flag = true;
      for(long long int j=0; j<(long long int)Curr.succ_list.size(); j++) {
  long long int SuccId = Curr.succ_list[j];
  long long int temp = get_vertex(SuccId).weight;
  if(flag) { /* Why not use flags? This could be erroneous */
    min = temp;
    flag = false;
  }
  else {
    if(temp<min) min = temp;
  }
      }
      Value TempVal(min);
      InitialEstimate.AddEstimate(CurrId, TempVal);
    }
  }
  
  Prev = InitialEstimate;
  while(true) {
    //show(Prev);
    Next = UpdateEstimation(Prev);
    if(Next==Prev) break; // When the current estimation is equal to next, we stop
    Prev = Next;
  }
  //show(Next);
  
  vector<long long int> ret;
  
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    if(Next.estimate[i].bestSum == PINF) ret.push_back(i);
  }
  
  return ret;
}

void Graph::reduce(int &u, int &v) {
    int u1 = u;
    int v1 = v;
    while ( v1 != 0) {
        int r = u1 % v1;
        u1 = v1;
        v1 = r;
    }
    u = u/u1;
    v = v/u1;
}

bool Comaparator::operator ()(pair<int,int> a, pair<int,int> b){
  return(a.first*b.second < a.second*b.first);
}

map<int, pair<int, int> > Graph::GetValue(int a) {
    bool verbose = false;  
    int n = a;
    int num = n_vertices;
    map<int, pair<int, int> > VertexValue;
    
    vector<pair<int,int> > allFractions;
    for (int l = 1; l <=n; l++)
    { 
      for (int a = 0; a <= l ;a++)
      {
        pair<int,int> temp(a,l);
        reduce(temp.first, temp.second);
        allFractions.push_back(temp);
      }
    }
    Comaparator compare;
    sort(allFractions.begin(),allFractions.end(),compare);
    allFractions.erase(unique(allFractions.begin(),allFractions.end()), allFractions.end());

    if(verbose) {
      for(std::vector<pair<int,int> > :: iterator it=allFractions.begin();it!=allFractions.end();it++) {
        cout<<(*it).first<<"/"<<(*it).second<<endl;
      }
    }
    for(int i=0; i<num; i++) {
      Vertex curr = vertices[i];
      int CurrId = curr.vid;
      if(CurrId > n || CurrId == 0) continue;
      if(verbose) cout << "CurrId " << CurrId << endl;
      //assert(i==CurrId);
      int arg1 = 1;//a
      int arg2 = 2;//l
      int old_arg1 = 1;
      int old_arg2 = 2;
      bool t;
      while(true) {
        if(2*n < arg2) break;
        pair<int, int> p(arg1, arg2);
        t = isWinning(p,n,CurrId);
        old_arg1 = arg1;
        old_arg2 = arg2;
        if(t) {
          arg1 = 2*arg1 + 1;
        }
        else {
          arg1 = 2*arg1 - 1;
        }
        arg2 = 2*arg2;
        //if(2*arg2 > n) break;// TODO a or a+1 
        //int diff = //1/arg2;
      }
      pair<int, int> p1, p2;
      if(!t) {
        if(verbose) cout << "true" << endl;
        p1.first = (old_arg1-1)/2;
        p1.second = (old_arg2)/2;
        p2.first = old_arg1;
        p2.second = old_arg2;
      }
      else {
        if(verbose) cout << "false" << endl;
        p1.first = old_arg1;
        p1.second = old_arg2;
        p2.first = (old_arg1 + 1)/2;
        p2.second = old_arg2/2;
      }
      reduce(p1.first, p1.second);
      reduce(p2.first, p2.second);
      
      if(verbose) {
        cout << "Interval" << endl; 
        cout << p1.first << "/" << p1.second << endl;
        cout << p2.first << "/" << p2.second << endl;
      }
      int startIndex = -1;
      int endIndex = -1;
      for(int i=0;i<(int)allFractions.size();i++) {
        if(allFractions[i].first*p1.second >= allFractions[i].second*p1.first) {
          startIndex=i;
          break;
        }
      }
      endIndex = startIndex;
      for(int i=0;i<(int)allFractions.size();i++) {
        if(allFractions[i].first*p2.second > allFractions[i].second*p2.first) {
          endIndex=i-1;
          break;
        }
      }
      if(verbose) cout << "startIndex " << startIndex  << " endIndex " << endIndex<< endl;
      if(endIndex < startIndex || startIndex == -1 || endIndex == -1) {
        assert(false);
      }
      if(verbose) {
        cout<<"required : ";
        for(int i=startIndex;i<=endIndex;i++) {
          cout<<allFractions[i].first<<"/"<<allFractions[i].second<<" ";
        }
        cout<<endl;
      }

      // Calculating the value from the possible fractions in the interval
      int correctIndex = -1;
      int currIndex = (startIndex + endIndex)/2;
      pair<int, int> p;
      //bool t;
      bool diff = false;
      while(true) {
        if(endIndex == startIndex) {
          correctIndex = startIndex;
          break;
        }
        if(endIndex - startIndex == 1) { 
          diff = true;
          break;
        }
        p = allFractions[currIndex];
        t = isWinning(p,n,CurrId);
        if(t) {
          startIndex = currIndex;
          currIndex = (startIndex + endIndex)/2;
        }
        else {
          endIndex = currIndex;
          currIndex = (startIndex + endIndex)/2;
        } 
      }
      bool t1, t2;
      if(diff) {
        t1 = isWinning(allFractions[startIndex], n, CurrId);
        t2 = isWinning(allFractions[endIndex], n, CurrId);
        if(t2) {
          correctIndex = endIndex;
        }
        else {
          correctIndex = startIndex;
        }
        assert(t1 || t2);
      }
      assert (correctIndex!=-1);

      if(verbose) {
        cout << "correctIndex : " << correctIndex << " fraction is: ";
        cout << allFractions[correctIndex].first << "/" << allFractions[correctIndex].second << endl;
      }
      VertexValue[CurrId] = allFractions[correctIndex];
    }
    return VertexValue; 
}


bool Graph::isWinning(pair<int, int> w, int verticesInitial, int vid) {
  Graph g1;
  //map<long long int, Vertex> IdToVertex;
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    int weight;
    if(curr.weight == 0) {weight = -w.first;}
    else if(curr.weight == 1) {weight = w.second - w.first;}
    if(curr.vid == 0 || curr.vid > verticesInitial) weight = 0; 
    Vertex temp(curr.vid, weight, curr.player);
    g1.add_vertex(temp);
    //IdToVertex.insert(pair<long long int,Vertex>(curr.vid, ImprovementArena.get_vertex(curr.vid)));   
  }
  
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    long long int CurrId = curr.vid;
    
    for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
      long long int SuccId = curr.succ_list[j];
      Vertex Succ = get_vertex(SuccId);
      g1.set_successor(CurrId, SuccId);
    }
  }
    //if(curr.player==odd) assert(count==(long long int)curr.succ_list.size());
  g1.fill_predecessors();
  //g1.show();

  vector<long long int> Win = g1.OptimalStrategyImprovement();
  /*cout << w.first << "/" << w.second << endl;
  for(long long int i=0; i<(long long int)Win.size(); i++) {
    cout << Win[i] << " ";
  }
  cout << endl;
  */
  for(long long int i=0; i<(long long int)Win.size(); i++) {
    if(Win[i] == vid) return true;
  }
  return false;
}

vector<vector<int> > Graph::GenerateSCC() {
  stack<int> S;
  int size = vertices.size();
  int size1 = size;
  for(int i=0; i<size; i++) {
    Vertex curr = vertices[i];
    int CurrId = curr.vid;
    if(CurrId>size1) size1 = CurrId;
  }
  size1++;
  cout << "size " << size << endl;
  cout << "size1 " << size1 << endl;
  vector<bool> visited(size1, false);

  // Filling vertices in the stack as per their finishing order
  for(int i=0; i<size; i++) {
    Vertex curr = vertices[i];
    int CurrId = curr.vid;
    if(!visited[CurrId]) {
      DFSFinish(S,visited,curr);
    }
  }

  // Generating Transpose
  Graph Transpose;
  for(int i=0; i<(int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    Vertex temp(curr.vid, curr.weight, curr.player);
    Transpose.add_vertex(temp);
  }
  
  for(int i=0; i<(int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    int CurrId = curr.vid;
    
    for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
      long long int SuccId = curr.succ_list[j];
      Vertex Succ = get_vertex(SuccId);
      // Successor is inverted
      Transpose.set_successor(SuccId, CurrId);
    }
  }
  Transpose.fill_predecessors();
  //Transpose.showG();

  // reset the visited vector
  for(int i=0; i<size1; i++) {
    visited[i] = false;
  }

  std::vector<vector<int> > ret;
  int CurrId;
  while(S.empty() == false) {
    CurrId = S.top();
    S.pop();
    //cout << "CurrId : " << CurrId << endl;
    if(!visited[CurrId]) {
      //cout << "Not visited" << endl;  
      //cout << "CurrId : " << CurrId << endl;
      vector<int> ret1;
      Transpose.DFS(ret1, visited, Transpose.get_vertex(CurrId));
      ret.push_back(ret1);
      /*
      cout << "New SCC " << endl;
      for(int i=0; i<(int)ret1.size(); i++) {
        cout << ret1[i] << " ";
      }
      cout << endl;
      */
    } 
    //cout << CurrId << " ";
  }
  cout << endl;

  return ret;
}

void Graph::DFSFinish(stack<int> &S, std::vector<bool> &visited, Vertex curr) {
  int CurrId = curr.vid;
  visited[CurrId] = true;
  vector<long long int> succ = curr.succ_list;

  for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
    long long int SuccId = curr.succ_list[j];    
    if(!visited[SuccId]) DFSFinish(S,visited,get_vertex(SuccId));
  }
  // Pushing into stack after all children have returned
  S.push(CurrId);
}

void Graph::DFS(vector<int> &SCC,  std::vector<bool> &visited, Vertex curr) {
  //cout << "In DFS" << endl;
  int CurrId = curr.vid;
  visited[CurrId] = true;
  vector<long long int> succ = curr.succ_list;
  //cout << "CurrId " << CurrId << endl;
  SCC.push_back(CurrId);
  //cout << "Succ list size " << (long long int)curr.succ_list.size() << endl;
  for(long long int j=0; j<(long long int)curr.succ_list.size(); j++) {
    long long int SuccId = curr.succ_list[j];
    //cout << "Curr ID " << CurrId << " SuccId " << SuccId << " "<< visited[SuccId] << endl;
    if(!visited[SuccId]) DFS(SCC,visited,get_vertex(SuccId));
  }
}

Graph Graph::computeSubGame(vector<int> vert_set)
{
  //cout << "In computeSubGame" << endl;
  Graph result;
  
  for (vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); i++) {
    if (find(vert_set.begin(), vert_set.end(), (*i).vid) != vert_set.end()) {
      // If the vertex *i is in the vert_set 
      result.add_vertex((*i).vid, (*i).weight, (*i).player); 
    }
  }
  //cout << "Size: " << result.vertices.size() << endl;
  for (vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); i++) {
    if (find(vert_set.begin(), vert_set.end(), (*i).vid) != vert_set.end()) {
      vector<long long int> succ_list = (*i).succ_list;
      vector<long long int>::iterator j = succ_list.begin();
      while (j != succ_list.end()) {
        if(find(vert_set.begin(), vert_set.end(), (*j)) != vert_set.end()) {
          //cout << (*i).vid << " " << *j << endl;
          result.set_successor((*i).vid,(*j));
          //cout << "Te" << endl;
        }
        j++;
      }
    }
  }
  //cout << "CHeck"  << endl;
  result.fill_predecessors();
  return result;
}

void Graph::show() {
  cout << "Graph" << endl;
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    vertices[i].show();
    Vertex Curr = get_vertex(vertices[i].vid);
    cout << "Successors : ";
    for(long long int j=0; j<(long long int)Curr.succ_list.size(); j++) {
      cout << Curr.succ_list[j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void Graph::show(Estimation e) {
  cout << "id Est" << endl;
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    cout << vertices[i].vid << " " << e.estimate[vertices[i].vid].bestSum << endl;
  }
  cout << endl;
}


void Graph::output_dot(ostream& out)
{
  vector<Vertex>::iterator i;
  
  out << "digraph G { \n";
  for (i = vertices.begin(); i != vertices.end(); i++) {
    long long int temp = (*i).vid;
    //temp+=1;
    ostringstream ss;
    ss << temp;
    string v = ss.str();

    if ((*i).player == odd) {
      out << v <<"[label=\"\\N (" << (*i).weight << ")\", color=black, shape=polygon ]\n";
    }
    else {
      out << v <<"[label=\"\\N (" << (*i).weight << ")\" , color=black ]\n";
    }
    
    vector<long long int> succ_list = (*i).succ_list;
    vector<long long int>::iterator j = succ_list.begin();
    
    while (j != succ_list.end()) {
      out << v << "->"<< get_vertex(*j).vid <<"\n";
      j++;
    }
  }
  out << "}\n";
}


void Graph::showG()
{
  string dot_file = "all.dot";
  string pdf_file = "all.pdf";
  ofstream outs(dot_file);
    
  this->output_dot(outs);
  outs.close();

  string cmd1 = "dot -Tpdf " + dot_file + " -o" + pdf_file;
  string cmd2 = "xdg-open " + pdf_file + "&"; 
   
  system(cmd1.c_str());
  system(cmd2.c_str());
  
  getchar();
  
}
 


void Graph::output_dot_mmpg(ostream& out)
{
  vector<Vertex>::iterator i;
  
  out << "digraph G { \n";
  for (i = vertices.begin(); i != vertices.end(); i++) {
    long long int temp = (*i).vid;
    ostringstream ss;
    ss << temp;
    string v = ss.str();

    string weights;
    for(vector<int> :: iterator it=(*i).mmpgWeights.begin(); it!=(*i).mmpgWeights.end(); it++) {
      stringstream ss;
      int temp = *it;
      ss << temp;
      weights+=ss.str()+",";
    }
    stringstream ss1;
    ss1 << (*i).mmpgPlayer+3;
    out << v <<"[label=\"\\N (" << weights << ")\", color=black, shape=polygon, sides=" << ss1.str() <<"]\n";    

    vector<long long int> succ_list = (*i).succ_list;
    vector<long long int>::iterator j = succ_list.begin();
    
    while (j != succ_list.end()) {
      out << v << "->"<< get_vertex(*j).vid <<"\n";
      j++;
    }
  }
  out << "}\n";
}

void Graph::showGmmpg()
{
  string dot_file = "all.dot";
  string pdf_file = "all.pdf";
  ofstream outs(dot_file);
    
  this->output_dot_mmpg(outs);
  outs.close();

  string cmd1 = "dot -Tpdf " + dot_file + " -o" + pdf_file;
  string cmd2 = "xdg-open " + pdf_file + "&"; 
   
  system(cmd1.c_str());
  system(cmd2.c_str());
  
  getchar();
  
}

Graph Graph::createMpg(int p) {
  Graph g1;
  PLAYER pl = odd;
  Vertex v1(0,0,pl);
  g1.add_vertex(v1);
  vector<Vertex>::iterator i;
  for (i = vertices.begin(); i != vertices.end(); i++) {
    long long int temp = (*i).vid;
    PLAYER mpgplayer;
    if((*i).mmpgPlayer==p) mpgplayer = even;
    else mpgplayer = odd;
    long long int w = (*i).mmpgWeights[p];
    Vertex v(temp,w,mpgplayer);
    g1.add_vertex(v);
  }
  for(long long int i=0; i<(long long int)vertices.size(); i++) {
    Vertex Curr1 = g1.get_vertex(vertices[i].vid);
    if(Curr1.player==even) g1.set_successor(Curr1.vid,0);
    Vertex Curr = get_vertex(vertices[i].vid);
    for(long long int j=0; j<(long long int)Curr.succ_list.size(); j++) {
      g1.set_successor(Curr.vid,Curr.succ_list[j]);
    }
  }
  //g1.show();
  long long int next_vid=g1.vertices.size();
  for(long long int i=0; i<(long long int)g1.vertices.size(); i++) {
    Vertex Curr = g1.vertices[i];
    //cout<<"Currrent vertex id " << Curr.vid<<endl;
    vector<long long int> del;
    for(long long int j=0; j<(long long int)Curr.succ_list.size(); j++) {
      if(Curr.player==g1.get_vertex(Curr.succ_list[j]).player) {
        //cout<<"hello"<<endl;
        //cout<<Curr.vid<<" "<<Curr.succ_list[j]<<endl;
        del.push_back(j);
        PLAYER newplayer;
        if(Curr.player==even) newplayer = odd;
        else newplayer = even;
        Vertex v(next_vid,0,newplayer);
        g1.add_vertex(v);
        g1.set_successor(Curr.vid,next_vid);
        g1.set_successor(next_vid,Curr.succ_list[j]);
        g1.vertices[i].delete_succ(Curr.succ_list[j]);
        if(newplayer==even) g1.set_successor(next_vid,0);
        next_vid++;
      }
    }
  }
  g1.fill_predecessors();
  return g1;
}

void Graph::MakeNash() {
  int i;
  NumPlayers++;
  for (i = 0; i <(int)vertices.size(); i++) {
    Vertex curr = vertices[i];
    
    int NewPlayer = curr.mmpgPlayer+1;
    vector<int> NewWeight, mmpgWeight = curr.mmpgWeights;
    
    NewWeight.push_back(mmpgWeight[0]);
    for(int j=0; j<(int)mmpgWeight.size(); j++) {
      NewWeight.push_back(mmpgWeight[j]);
    }

    vertices[i].mmpgPlayer = NewPlayer;
    vertices[i].mmpgWeights = NewWeight;
  }
}
 
