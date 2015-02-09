#include "Graph.hpp"
#include <bits/stdc++.h>
#include <math.h>
#include "../lp_dev/lp_lib.h"

using namespace std;
using namespace MeanPayoffGame;


int NumPlayers;
bool bribe;
bool nash;

void recurse(int current, vector<pair<int,int> > v, vector<vector<pair<int, int> > > &Iterate,
 vector<pair<int, int> > UniqueValue[]) 
{ 
  if (current==NumPlayers) {  
    //for (uint i=0;i<NumPlayers;i++) std::cout << levels[i] << " "; 
    //std::cout << std::endl; 
    Iterate.push_back(v);
  } 
  else {
  	vector<pair<int,int> > temp = UniqueValue[current];
    for (int i=0; i<(int)temp.size(); i++) {
      v.push_back(temp[i]); 
      recurse(current+1,v,Iterate,UniqueValue); 
      v.erase(v.end()-1);
    }
  } 
} 


int main(int argc, char* argv[]) {
	//PLAYER pl = odd;
	//Vertex v(0,0,pl);
	Graph g;
	bool GraphOutput = false;
	
	//command line arguments
	ifstream infile(argv[1]);
	
	string t(argv[2]);
	if(t != "-b") {
		cout << "Bribe should be given -b " << endl;
		return 0;
	}

	t.assign(argv[3]);
	if(t=="0") {bribe = false;}
	else if(t=="1") {bribe = true;}
	else {nash = true;}

	cout << "Bribe: " << bribe << endl;
	cout << "Nash: " << nash << endl;

	
	string t1(argv[4]);
	if(t1 != "-n") {
		cout << "Number of players should be given -n " << endl;
		return 0;
	}
	t1.assign(argv[5]);
	int n;
	stringstream(t1)>>n; 
	cout << "NumPlayers: " << n << endl;
	g.NumPlayers = n;

	t1.assign(argv[6]);
	if(t1 == "-out") {GraphOutput = true;}
	
	string line;
	getline(infile, line);
	istringstream iss(line);
	string word;
	iss>>word;
	string mode = word;
	//cout<<mode<<endl;
	long long int a;
	iss>>word;
	stringstream(word) >> a;
	a+=1;
	for(long long int i=0;i<a;i++) {
		getline(infile, line);
		//cout<<line<<endl;
		istringstream iss(line);
		long long int curr_vid,curr_pl,curr_pr1,curr_pr2,curr_pr3;
		iss>>word;
		stringstream(word) >> curr_vid;
		iss>>word;
		stringstream(word) >> curr_pr1;
		vector<int> pr_vec;
		pr_vec.push_back(curr_pr1);
		if(mode == "mmpg") {
			//cout<<"came1"<<endl;
			for(int k=0; k<n-1; k++) {
				iss>>word;
				stringstream(word) >> curr_pr1;
				pr_vec.push_back(curr_pr1);
			}
		}
		else {
			cout<<"give a valid file\n";
			return 0;
		}
		//cout<<"came2"<<endl;
		iss>>word;
		stringstream(word) >> curr_pl;
		iss>>word;
		curr_vid+=1;
		Vertex v(curr_vid,pr_vec,curr_pl);
		g.add_vertex(v);
		long long int prev=0;
		long long int pre=0;
		string suc;
		//if(curr_pl==0) g.set_successor(curr_vid, 0);
		
		while(true) {
			//cout<<"came3"<<endl;
			pre=prev+word.substr(prev).find(",");
			long long int temp;
			stringstream(word.substr(prev,pre-prev)) >> temp;
			g.set_successor(curr_vid,temp+1);
			if(prev>pre) break;
			prev=pre+1;
		}
		//cout<<curr_vid<<endl;
	}

	getline(infile, line);
	//cout<<line;

	//cout<<"came4"<<endl;
	string line1 = line;
	string word1;
	istringstream iss1(line1);
	iss1 >> word1;
	iss1 >> word1;
	stringstream(word1) >> g.StartId;
	
	g.StartId += 1;
	g.fill_predecessors();
	clock_t startTime = clock();	
	// Nash
	if(nash) {
		bribe = false;
		g.MakeNash();
		cout << "Nash Equilibrium " << endl;
		cout << "New Number of players : "<< g.NumPlayers << endl;
	}

	NumPlayers = g.NumPlayers - 1;
	if(GraphOutput) g.showGmmpg();
	//
	Graph G[NumPlayers]; 

	for(int pl=0; pl<NumPlayers; pl++) {
		G[pl] = g.createMpg(pl+1);
	}
	//G[1] = g.createMpg(2);
	//g1.showG();
	//g2.showG();
	map<int, pair<int, int> > VertexValue[NumPlayers];
	vector<vector<pair<int,int> > > Iterate;
	vector<pair<int, int> > UniqueValue[NumPlayers];
	std::vector<pair<int,int> > v;

	for(int pl=0; pl<NumPlayers; pl++) {
		VertexValue[pl] = G[pl].GetValue(a);
		cout << "VertexValue " << pl << endl;
		for(int i=0; i<(int)g.vertices.size(); i++) {
			if(g.vertices[i].vid > a || g.vertices[i].vid == 0) continue;
			cout << VertexValue[pl][g.vertices[i].vid].first << "/"; 
			cout << VertexValue[pl][g.vertices[i].vid].second << endl;
		}
		
		map<pair<int, int>, bool> temp;
		for(int i=0; i<(int)g.vertices.size(); i++) {
			if(g.vertices[i].vid > a || g.vertices[i].vid == 0) continue;
			pair<int, int> p = VertexValue[pl][g.vertices[i].vid];
			if(temp.count(p) == 0) {
				temp[p] = true;
				//count++;
				UniqueValue[pl].push_back(p);
			}
		}

		cout << "UniqueValue " << pl << endl;
		for(int i=0; i<(int)UniqueValue[pl].size(); i++) {
			cout << UniqueValue[pl][i].first << "/" << UniqueValue[pl][i].second << " " ;
		}
		cout << endl;
	}
	//g.showGmmpg();

	recurse(0,v,Iterate,UniqueValue);
	cout << "Iterate " << endl;
	for(int i=0; i<(int)Iterate.size(); i++) {
		vector<pair<int,int> > temp = Iterate[i];
		for(int j=0; j<(int)temp.size(); j++) {
			cout << temp[j].first << "/" << temp[j].second << " ";
		}
		cout << endl;
	}

	int StartId = g.StartId;
	cout << "StartId " << StartId << endl;
	Vertex Start = g.get_vertex(StartId);
	std::vector<bool> visited(a,false);
	vector<int> ReachableSet;
	map<int, bool> isReachable;
	for(int k=0; k<(int)g.vertices.size(); k++) {
		isReachable[g.vertices[k].vid] = false;
	}
	
	// Getting all the vertices reachable from the start vertex
	g.DFS(ReachableSet, visited, Start);

	for(int i=0; i<(int)ReachableSet.size(); i++) {
		isReachable[ReachableSet[i]] = true;
		cout << ReachableSet[i] << " " << isReachable[ReachableSet[i]] << endl;
	}
	cout << endl;
	//g.showG();
	float FinalAnswer;
	bool FinalFlag = false;

	for(int i1=0; i1<(int)Iterate.size(); i1++) {
		vector<pair<int,int> > li = Iterate[i1];
			for(int pl=0; pl<NumPlayers; pl++) {
				cout << "li" << pl << " " << li[pl].first << "/" << li[pl].second << endl;
			}
			vector<int> nodeSet, Q;
			for(int k=0; k<(int)g.vertices.size(); k++) {
				bool flag10 = false;
				for(int pl=0; pl<NumPlayers; pl++) {
					pair<int,int> CurrPair = VertexValue[pl][g.vertices[k].vid];
					if(CurrPair.first*li[pl].second>CurrPair.second*li[pl].first) {
						flag10 = true;
						break;
					}	
				}
				if(!flag10) {
					nodeSet.push_back(g.vertices[k].vid);
				}
			}

			cout << "nodeSet" << endl;
			for(int k=0; k<(int)nodeSet.size(); k++) {
				cout << nodeSet[k] << " ";
			}
			cout << endl;

			// check if start vertex is in nodeSet
			if(find(nodeSet.begin(), nodeSet.end(), StartId) == nodeSet.end()) continue;

			for(int k=0; k<(int)nodeSet.size(); k++) {
				if(isReachable[nodeSet[k]]) Q.push_back(nodeSet[k]);
			}

			cout << "Q" << endl;
			for(int k=0; k<(int)Q.size(); k++) {
				cout << Q[k] << " ";
			}
			cout << endl;
			
			pair<int,int> max[NumPlayers];//max value when player was maximiser among all v belongs to Q

			vector<bool> flag(NumPlayers, false);//1 = false, flag2 = false;
			for(int pl=0; pl<NumPlayers; pl++) {
				for(int k=0; k<(int)Q.size(); k++) {
					if(!flag[pl]) {
						max[pl]=VertexValue[pl][Q[k]];
						flag[pl] = true;
						continue;
					}
					else {
						pair<int,int> p = VertexValue[pl][Q[k]];
						if(max[pl].first*p.second < p.first*max[pl].second) max[pl] = p;
					}
				}
				cout << "max" << pl << " " << max[pl].first << "/" << max[pl].second << endl;
			}
			//g.showGmmpg();

			vector<vector<int> > SCCList;
			//g.show();
			Graph subGraph = g.computeSubGame(Q);
			//subGraph.showG();
			SCCList = subGraph.GenerateSCC();
			for(int k=0; k<(int)SCCList.size(); k++) {// iterate over all SCC's
				vector<int> SCC = SCCList[k];
				vector<int> reward[NumPlayers+1];
				Graph SCCGraph = subGraph.computeSubGame(SCC);
				cout << "SCC" << endl;
				map<pair<int,int> , int> EdgeIndex;
				map<int, pair<int,int> > EdgeIndexRev;				
				int numEdges = 0;
				for(int l=0; l<(int)SCC.size(); l++) {
					cout << SCC[l] << " " ;
					Vertex curr = SCCGraph.get_vertex(SCC[l]);
					for(long long int m=0; m<(long long int)curr.succ_list.size(); m++) {
						pair<int, int> p(curr.vid,curr.succ_list[m]);
						EdgeIndex[p] = numEdges;
						EdgeIndexRev[numEdges+SCC.size()+1] = p;
						numEdges++;
					}
					for(int pl=0; pl<=NumPlayers; pl++) {
						reward[pl].push_back(g.get_vertex(SCC[l]).mmpgWeights[pl]);
						cout << g.get_vertex(SCC[l]).mmpgWeights[pl] << " ";
					}
					cout << endl;
				}
				cout << endl;
				cout << "numEdges" << numEdges << endl;
				// Linear Programming
					lprec *lp;
					int Ncol, *colno = NULL, j, i, ret = 0;
					REAL *row = NULL;

					int n = SCC.size();//Should be number of vertices in S
					Ncol = n+NumPlayers+numEdges;
					if(true) cout << "Ncol " << Ncol << endl;
					if(true) cout << "n " << n << endl;
					lp = make_lp(0, Ncol);
					if(lp == NULL)
						ret = 1; /* couldn't construct a new model... */

					if(ret == 0) {
						int i;
						for(i=1; i<=n; i++) {
							string str = to_string(i);
							str = "f" + str; 
							char *cstr = new char[str.length() + 1];
							strcpy(cstr, str.c_str());
							set_col_name(lp, i, cstr);
						}

						for(i=i; i<=n+numEdges; i++) {
							pair<int,int> p = EdgeIndexRev[i];
							cout << i << " " << p.first << "," << p.second << endl; 
							stringstream ss, ss1;
							ss << p.first;
							string str = ss.str();
							ss1 << p.second;
							string str1 = ss1.str();
							str = "f " + str + "," + str1; 
							char *cstr = new char[str.length() + 1];
							strcpy(cstr, str.c_str());
							set_col_name(lp, i, cstr);
						}
						
						for(;i<=Ncol;i++) {
							string str = "b";
							stringstream ss;
							ss << i-(n+numEdges);
							str = str + ss.str();
							char* cstr = &str[0];
							set_col_name(lp, i, cstr);
						}
						
						colno = (int *) malloc(Ncol * sizeof(*colno));
						row = (REAL *) malloc(Ncol * sizeof(*row));
						if((colno == NULL) || (row == NULL))
							ret = 2;
					}

					set_add_rowmode(lp, TRUE); 


					for(j=1; j<=Ncol; j++) { // f(i) >= 0 for all i in S and b's >= 0
						i=0;
						while(true) {
							if(i>=Ncol) break;
						    colno[i] = i+1;
						    if(i!=j-1) {row[i++] = 0;}
						    else {row[i++] = 1;}
						}

						if(bribe) {
							if(!add_constraintex(lp, i, row, colno, GE, 0)) {
								cout << ret << endl;
								assert(false);
							};
						}
						else {
							if(j>Ncol-NumPlayers) {
								if(!add_constraintex(lp, i, row, colno, EQ, 0)) {
									cout << ret << endl;
									assert(false);
								}
							}
							else {
								if(!add_constraintex(lp, i, row, colno, GE, 0)) {
									cout << ret << endl;
									assert(false);
								}
							}
						}
					}

					// sigma f(i) = 1 for all i belongs to S
					i=0;
					while(true) {
						if(i>=Ncol) break;
					    colno[i] = i+1;
					    if(i>=n) {row[i++] = 0;}
					    else {row[i++] = 1;}
					}
					if(!add_constraintex(lp, i, row, colno, EQ, 1)) {
						cout << ret << endl;
						assert(false);
					};

					// sum of values of outgoing edges =  sum of values of incoming edges
					// sum of values of outgoing edges = value of vertex
					for(j=1; j<=n; j++) {
						int CurrId = SCC[j-1];
						//cout << "CurrId " << CurrId << endl;
						Vertex curr = SCCGraph.get_vertex(CurrId);
						// sum of successors f values is equal to f of vertex 
						vector<long long int> succ = curr.succ_list;
						vector<int> Edges;
						for(int z=0; z<(int)succ.size(); z++) {
							//cout << "CurrId " << CurrId << " SuccID " << succ[z] << endl;
							pair<int,int> p(CurrId,(int)succ[z]);
							Edges.push_back(EdgeIndex[p]+n);
						}
						i=0;
						while(true) {
							if(i>=Ncol) break;
						    colno[i] = i+1;
						    if(find(Edges.begin(),Edges.end(),i) != Edges.end()) {row[i++] = 1;}
						    else if(i==j-1) {row[i++] = -1;}
						    else {row[i++] = 0;}
						}

						if(!add_constraintex(lp, i, row, colno, EQ, 0)) {
							cout << ret << endl;
							assert(false);
						};

						// sum of predecessors f values is equal to f of vertex 
						vector<long long int> pred = curr.pred_list;
						vector<int> Edges1;
						for(int z=0; z<(int)pred.size(); z++) {
							//cout << "CurrId " << CurrId << " PredID " << pred[z] << endl;
							pair<int,int> p((int)pred[z],CurrId);
							Edges1.push_back(EdgeIndex[p]+n);
						}
						i=0;
						while(true) {
							if(i>=Ncol) break;
						    colno[i] = i+1;
						    if(find(Edges1.begin(),Edges1.end(),i) != Edges1.end()) {row[i++] = 1;}
						    else if(i==j-1) {row[i++] = -1;}
						    else {row[i++] = 0;}
						}

						if(!add_constraintex(lp, i, row, colno, EQ, 0)) {
							cout << ret << endl;
							assert(false);
						};
					}


					// Constraints saying that each player benefits with bribe
					for(int pl=0; pl<NumPlayers; pl++) {
						i=0;
						while(true) {
							if(i>=Ncol) break;
							colno[i] = i+1;
						    if(i==n+numEdges+pl) {row[i++] = max[pl].second;}
						    else if(i>=n) {row[i++] = 0;}
						    else {row[i++] = reward[pl+1][i-1]*(max[pl].second);}//weight of vertex w.r.t player 1
						}
						if(false) {
							for(i=0; i<Ncol; i++) {
						      cout << colno[i] << " " << row[i] << endl;
						    }
						}
						if(!add_constraintex(lp, i, row, colno, GE, max[pl].first)) {
							cout << ret << endl;
							assert(false);
						};						
					}

					set_add_rowmode(lp, FALSE);

					// Objective function - maximising player1's gains
					i=0;
					while(true) {
						if(i>=Ncol) break;
						colno[i] = i+1;
					    if(i>=n+numEdges) {row[i++] = -1;}
					    else if(i>=n) {row[i++] = 0;}
					    else {row[i++] = reward[0][i-1];}//weight of vertex 'i+1' w.r.t player 0
					}

					if(false) {
						for(i=0; i<Ncol; i++) {
					      cout << colno[i] << " " << row[i] << endl;
					    }
					}
					
					if(!set_obj_fnex(lp, i, row, colno)) {
						cout << "Error in Objective function part " << ret << endl;
						assert(false);
					}

					set_maxim(lp);

					write_LP(lp, stdout);

					set_verbose(lp, IMPORTANT);

					ret = solve(lp);
					if(ret == OPTIMAL)
						ret = 0;
					else {
						cout << "Error in solve part " << ret << endl;
					}

					printf("Objective value: %f\n", get_objective(lp));
					
					get_variables(lp, row);
					for(j = 0; j < Ncol; j++)
						printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);

					float Answer = (float)get_objective(lp);
					if(!FinalFlag) {
						FinalFlag = true;
						FinalAnswer = Answer;
					}
					else {
						if(FinalAnswer < Answer) FinalAnswer = Answer;
					}
					//cout << "Answer" << Answer << endl;
				 	if(row != NULL)
						free(row);
					if(colno != NULL)
						free(colno);

					if(lp != NULL) {
						/* clean up such that all used memory by lpsolve is freed */
						delete_lp(lp);
					}
				//lp
			}
		}
	cout << "FinalAnswer : " << FinalAnswer ;


	//g.show();
	//g.GenerateSCC();
	//
	
	/*vector<long long int> Win = g.OptimalStrategyImprovement();
	cout << "{";
	vector<long long int> WinMax;
	for(long long int i=0; i<(long long int)Win.size(); i++) {
		if(Win[i]<=a) WinMax.push_back(Win[i]-1);
	}
	for(long long int i=0; i<(long long int)WinMax.size(); i++) {
		if(i!=WinMax.size()-1) {cout << WinMax[i] << ",";}
		else {cout << WinMax[i];}
	}

	cout << "}";
	cout << endl;
*/
	if(bribe) cout << "Time: " <<  ((clock() - startTime)/(double) CLOCKS_PER_SEC) <<  endl; 
	else cout<<endl;
	
}

