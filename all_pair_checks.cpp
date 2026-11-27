#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cstdint>
using namespace std; using U=uint32_t;
int rk(const vector<U>&g,U rows,U cols,int stop=100){U p[32]={};int r=0;while(rows){int i=__builtin_ctz(rows);rows&=rows-1;U v=g[i]&cols;while(v){int j=31-__builtin_clz(v);if(p[j])v^=p[j];else{p[j]=v;if(++r>=stop)return r;break;}}}return r;}
bool cert(const vector<U>&g,int a,int b){int n=g.size();U AB=(1u<<a)|(1u<<b),R=((1u<<n)-1)^AB;if(rk(g,AB,R)!=2)return false;for(U S=(R-1)&R;S;S=(S-1)&R){U T=R^S;if(!(S&(R&-R)))continue;int l=__builtin_popcount(S),r=__builtin_popcount(T);if(rk(g,S,T,3)>=3)continue;if(l+2<=r&&rk(g,AB|S,T)==l+2)continue;if(r+2<=l&&rk(g,AB|T,S)==r+2)continue;return false;}return true;}

bool cut(const vector<U>&g,U S,U T,int r){if(r>=5)return true;for(int k=0;k<2;k++,std::swap(S,T)){int s=__builtin_popcount(S),t=__builtin_popcount(T);if(s<3||r==s)continue;for(U a=S;a;a&=a-1){U bit=a&-a;for(U b=a^bit;b;b&=b-1){U AB=bit|(b&-b),R=S^AB;if(rk(g,R,T,3)>=3)continue;if(t+2<=s-2 && rk(g,T|AB,R)==t+2)continue;return false;}}}return true;}
bool allcert(const vector<U>&g){int n=g.size();U V=(1u<<n)-1;for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){U ab=(1u<<a)|(1u<<b);if(rk(g,ab,V^ab)!=2)return false;}
for(int s=1;s<=4;s++){U u=(1u<<s)-1;while(u<=V){if(!cut(g,u,V^u,rk(g,u,V^u,5)))return false;U c=u&-u,v=u+c;if(v>V)break;u=v|(((v^u)>>2)/c);}}
for(U u=1;u<(1u<<(n-1));u++){U v=V^u;if(__builtin_popcount(u)<=4||__builtin_popcount(v)<=4)continue;if(!cut(g,u,v,rk(g,u,v,5)))return false;}return true;}
int main(){// Check equivalence against original pair-first search on independent small graphs.
mt19937_64 checkrng(203);for(int n=6;n<=12;n++)for(int i=0;i<20;i++){vector<U>g(n);for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(checkrng()&1){g[a]|=1u<<b;g[b]|=1u<<a;}bool z=true;for(int a=0;a<n&&z;a++)for(int b=a+1;b<n&&z;b++)if(!cert(g,a,b))z=false;if(z!=allcert(g)){cerr<<"Self-test failed\n";return 1;}}cerr<<"140 pair-first/cut-first checks passed\n";mt19937_64 rng(20261128);ofstream f("all_pair_extension.csv");f<<"n,index,adjacency_masks,all_pair_pass\n";for(int n=17;n<=24;n++){int N=n<=20?20:10,pass=0;for(int i=0;i<N;i++){vector<U>g(n);for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(rng()&1){g[a]|=1u<<b;g[b]|=1u<<a;}bool c=allcert(g);pass+=c;f<<n<<','<<i<<',';for(int a=0;a<n;a++){if(a)f<<';';f<<g[a];}f<<','<<c<<'\n';f.flush();}cerr<<n<<" "<<pass<<"/"<<N<<'\n';}}
