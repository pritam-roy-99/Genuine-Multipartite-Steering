#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cstdint>
using namespace std; using U=uint32_t;
int rk(const vector<U>&g,U rows,U cols,int stop=100){U p[32]={};int r=0;while(rows){int i=__builtin_ctz(rows);rows&=rows-1;U v=g[i]&cols;while(v){int j=31-__builtin_clz(v);if(p[j])v^=p[j];else{p[j]=v;if(++r>=stop)return r;break;}}}return r;}
bool cert(const vector<U>&g,int a,int b){int n=g.size();U AB=(1u<<a)|(1u<<b),R=((1u<<n)-1)^AB;if(rk(g,AB,R)!=2)return false;for(U S=(R-1)&R;S;S=(S-1)&R){U T=R^S;if(!(S&(R&-R)))continue;int l=__builtin_popcount(S),r=__builtin_popcount(T);if(rk(g,S,T,3)>=3)continue;if(l+2<=r&&rk(g,AB|S,T)==l+2)continue;if(r+2<=l&&rk(g,AB|T,S)==r+2)continue;return false;}return true;}
int main(){mt19937_64 rng(20261127);ofstream f("random_summary.csv");f<<"n,samples,fixed_pair_pass,all_pair_samples,all_pair_pass\n";ofstream data("random_graphs.csv");data<<"n,index,adjacency_masks,fixed_pair_pass,all_pair_pass\n";for(int n=6;n<=20;n++){int N=100,Nall=n<=14?100:(n<=16?20:0),fixed=0,all=0;for(int it=0;it<N;it++){vector<U>g(n);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(rng()&1){g[i]|=1u<<j;g[j]|=1u<<i;}bool q=cert(g,0,1);fixed+=q;int z=-1;if(it<Nall){z=q;for(int a=0;a<n&&z;a++)for(int b=a+1;b<n&&z;b++)if(!cert(g,a,b))z=0;all+=z;}data<<n<<','<<it<<',';for(int i=0;i<n;i++){if(i)data<<';';data<<g[i];}data<<','<<q<<','<<z<<'\n';}f<<n<<','<<N<<','<<fixed<<','<<Nall<<','<<all<<'\n';f.flush();cerr<<n<<" "<<fixed<<" "<<all<<"\n";}}
