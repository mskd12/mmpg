#ifndef __MEAN_PAYOFF_GAME_HPP
#define __MEAN_PAYOFF_GAME_HPP

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


#include "CommonFile.hpp"
#include "Vertex.hpp"
#include "Value.hpp"
#include "Estimation.hpp"

using namespace std;

namespace MeanPayoffGame {
	struct Graph {
		long long int n_vertices;
		vector<Vertex> vertices;
		int StartId;
		int NumPlayers;

		vector<long long int> EvenVertices;
		vector<long long int> OddVertices;
		map<long long int, Vertex*> IdToVertex;


		// Constructor
		Graph();
		
		/* Setter functions */
	    void add_vertex(Vertex v);
	    void add_vertex(long long int id, long long int weight, PLAYER player); 
	    void set_successor(long long int id, long long int succ);
	    void fill_predecessors();

	    /* Getter functions */
	    Vertex& get_vertex(const long long int id);

	    // Mpg functions
	    Estimation UpdateEstimation(Estimation);
	    vector<long long int> OptimalStrategyImprovement();

	    // Functions added for MMPG
	    bool isWinning(pair<int, int> w, int verticesInitial, int vid);
	    void reduce(int &u, int &v);
	    //bool compare(const pair<int,int> &a,const pair<int,int> &b);
	    map<int, pair<int, int> > GetValue(int a);
	    vector<vector<int> > GenerateSCC();
	    void DFSFinish(stack<int> &S, std::vector<bool> &visited, Vertex curr);
	    void DFS(vector<int> &SCC,  std::vector<bool> &visited, Vertex curr);
	    Graph computeSubGame(vector<int> vert_set);
		Graph createMpg(int p);
		void MakeNash();

	    //Functions for showing the graph
	    void show();
		void show(Estimation e);
		void output_dot(ostream& out);
		void output_dot_mmpg(ostream& out);
		void showG();
		void showGmmpg();

	};

	struct Comaparator
	{
	  bool operator ()(pair<int,int> a, pair<int,int> b);
	  /*bool compare(pair<int,int> a, pair<int,int> b) {
	    return a.released()>b.released();
		}
		*/
	};
}


#endif /* __PARITY_GAME_HPP__ */
