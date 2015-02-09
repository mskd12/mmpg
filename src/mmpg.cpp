#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h> 
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
	string t(argv[1]);
	int n;
	stringstream(t)>>n;
	//cout<<"Players:"<<" "<<n<<endl;
	string line;
	int i=0;
	int max_id;
	srand(time(NULL));
	while(getline(cin,line))
	{
		istringstream iss(line);
    string word;
		int a;
		if(i==0) {
			iss>>word;
			cout<<"mmpg"<<" ";
			iss>>word;
			stringstream(word) >> a;
			max_id=a;
			cout<<a<<";"<<endl;
		}
		else {
			iss>>word;
			cout<<word<<" ";
			iss>>word;
			for(int l=0;l<n;l++) {
				cout<<rand()%2<<" ";
			}
			iss>>word;
			cout<<rand()%n;
			while(iss >> word) {
				cout<<" "<<word;
	    }
			cout<<endl;
		}
		i++;
	}
	cout<<"StartIndexis:"<<" "<<rand()%(max_id+1)<<endl;
}		
	
