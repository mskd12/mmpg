#include <bits/stdc++.h>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]) {
	string t1(argv[1]);
	int n;
	stringstream(t1)>>n; 
	string t2(argv[2]);
	int m;
	stringstream(t2)>>m; 
	cout << "mmpg " << m*n-1 << endl;
	for(int i=0; i<m*n; i+=m) {
		int l=1;
		cout<<i<<" "<<"1 ";
		while(l!=n) {
			if(l==(i/m)+1) cout<<"1 ";
			else cout<<"0 ";
			l++;
		}
		if(m!=1) cout<<((i/m)+1)%n<<" "<<(i+1)<<","<<(i+m)%(m*n)<<" "<<"\""<<i<<"\""<<endl;
		else cout<<((i/m)+1)%n<<" "<<(i+m)%(m*n)<<" "<<"\""<<i<<"\""<<endl;
		for(int j=1; j<m; j++) {
			cout<<i+j<<" ";
			l=n;
			while(l!=0) {
				cout<<"0 ";
				l--;
			}
			cout<<((i/m)+1)%n<<" ";
			if(j==m-1) cout<<i;
			else cout<< i+j+1;
			cout<<" "<<"\""<<i+j<<"\""<<endl;
		}
	}
}