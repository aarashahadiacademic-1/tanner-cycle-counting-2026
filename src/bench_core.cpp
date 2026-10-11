// Exact and KB baseline primitives; C++17 -O3.
// These are the primitives used by chunk_hybrid.cpp's timed runs.
#include <bits/stdc++.h>
using namespace std;
using ull=unsigned long long;
using Clock=chrono::steady_clock;
struct Graph {
 int n,m,M;vector<vector<int>> a,out;vector<int> to,pos;
 explicit Graph(const string&path){
  ifstream f(path);int E,u,w;if(!(f>>n>>m>>E))throw runtime_error("missing graph");
  M=n+m;a.resize(M);out.resize(M);
  while(f>>u>>w){w+=n;a[u].push_back(w);a[w].push_back(u);
   int j=to.size();pos.push_back(out[u].size());pos.push_back(out[w].size());
   to.push_back(w);to.push_back(u);out[u].push_back(j);out[w].push_back(j+1);}
  if((int)to.size()!=2*E)throw runtime_error("edge count mismatch");
 }
};
struct ExactRoots {
 const Graph&G;int g,T,a,R,epoch=0,root=0;
 unsigned long long total_credits=0,max_root_credits=0;
 vector<ull>C;vector<int>tag,dep,par,br,q;vector<unsigned char>used,blocked;
 explicit ExactRoots(const Graph&gg,int girth,int target):G(gg),g(girth),T(target),
 a(girth/2-1),R(target-girth+2),C(target+1),tag(gg.M),dep(gg.M),par(gg.M),
 br(gg.M),q(gg.M),used(gg.M),blocked(gg.M){}
 bool done()const{return root==G.n;}
 void next_root(){
  if(done())return;unsigned long long ops=1;int v=root++;++epoch;int h=0,tail=0;
  q[tail++]=v;tag[v]=epoch;dep[v]=0;par[v]=-1;
  while(h<tail){int z=q[h++];++ops;if(dep[z]==a)continue;
   for(int w:G.a[z])if((++ops,tag[w]!=epoch)){
    tag[w]=epoch;dep[w]=dep[z]+1;par[w]=z;br[w]=(z==v?w:br[z]);q[tail++]=w;}}
  for(int ii=0;ii<tail;ii++){int x=q[ii];++ops;if(dep[x]!=a)continue;
   for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=1;}used[x]=1;
   auto dfs=[&](auto&&self,int z,int l)->void{
    ++ops;
    if(l>=2&&!(l&1)&&tag[z]==epoch&&dep[z]==a&&x<z&&br[x]!=br[z]){
     bool ok=true;for(int w=par[z];w!=-1;w=par[w])if((++ops,used[w])){ok=false;break;}
     if(ok)++C[g-2+l];}
    if(l==R)return;
    for(int w:G.a[z])if((++ops,!used[w]&&!blocked[w])){
     used[w]=1;self(self,w,l+1);used[w]=0;}
   };dfs(dfs,x,0);used[x]=0;
   for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=0;}
  }
  total_credits+=ops;max_root_credits=max(max_root_credits,ops);
 }
 vector<ull> counts()const{
  if(!done())throw runtime_error("exact incomplete");auto c=C;
  for(int t=g;t<=T;t+=2){
   if((2*c[t])%t)throw runtime_error("nonintegral count");c[t]=2*c[t]/t;}
  return c;
 }
};
vector<ull> exact(const Graph&G,int g,int T){
 ExactRoots E(G,g,T);while(!E.done())E.next_root();return E.counts();
}
vector<ull> kb(const Graph&G,int g,int T){
 int D=0;for(int u=0;u<G.n;u++)D=max(D,(int)G.a[u].size());
 int E=G.to.size();vector<ull>C(T+1),cur(E*D),nxt(E*D),sum(G.M*D);
 vector<int>active,na,nodes;vector<char>touch(G.M);
 for(int root=0;root<G.n;root++){
  active.clear();
  for(int j=0;j<(int)G.out[root].size();j++){
   int e=G.out[root][j];cur[e*D+j]=1;active.push_back(e);}
  for(int l=1;l<T;l++){
   nodes.clear();
   for(int e:active){int w=G.to[e];if(!touch[w]){touch[w]=1;nodes.push_back(w);}
    for(int j=0;j<D;j++)sum[w*D+j]+=cur[e*D+j];}
   na.clear();
   for(int w:nodes){if(w==root)continue;for(int e:G.out[w]){
    int z=G.to[e];if(z<G.n&&z<root)continue;bool nonzero=false;
    for(int j=0;j<D;j++){ull val=sum[w*D+j]-cur[(e^1)*D+j];
     nxt[e*D+j]=val;nonzero|=val!=0;}
    if(nonzero)na.push_back(e);}}
   if((l+1)%2==0&&l+1>=g){
    for(int j=0;j<(int)G.out[root].size();j++){
     int e=G.out[root][j]^1;
     for(int k=0;k<D;k++)if(k!=j)C[l+1]+=nxt[e*D+k];}}
   for(int e:active)for(int j=0;j<D;j++)cur[e*D+j]=0;
   for(int w:nodes){touch[w]=0;for(int j=0;j<D;j++)sum[w*D+j]=0;}
   cur.swap(nxt);active.swap(na);
  }
  for(int e:active)for(int j=0;j<D;j++)cur[e*D+j]=0;
 }
 for(auto&x:C){if(x%2)throw runtime_error("odd KB count");x/=2;}return C;
}
