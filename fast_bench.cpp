#include <bits/stdc++.h>
using namespace std;
using ull=unsigned long long;
using Clock=chrono::steady_clock;
struct Graph {
 int n,m,M;vector<vector<int>> a,out; vector<int> to,pos;
 explicit Graph(const string&path){ifstream f(path);int E,u,w;if(!(f>>n>>m>>E))throw runtime_error("missing graph");M=n+m;a.resize(M);out.resize(M);while(f>>u>>w){w+=n;a[u].push_back(w);a[w].push_back(u);int j=to.size();pos.push_back(out[u].size());pos.push_back(out[w].size());to.push_back(w);to.push_back(u);out[u].push_back(j);out[w].push_back(j+1);}if((int)to.size()!=2*E)throw runtime_error("edge count mismatch");}
};
struct ExactRoots {
 const Graph&G; int g,T,a,R,epoch=0,root=0; unsigned long long total_credits=0,max_root_credits=0;
 vector<ull>C;vector<int>tag,dep,par,br,q;vector<unsigned char>used,blocked;
 explicit ExactRoots(const Graph&gg,int girth,int target):G(gg),g(girth),T(target),a(girth/2-1),R(target-girth+2),C(target+1),tag(gg.M),dep(gg.M),par(gg.M),br(gg.M),q(gg.M),used(gg.M),blocked(gg.M){}
 bool done()const{return root==G.n;}
 void next_root(){if(done())return;unsigned long long ops=1;int v=root++;++epoch;int h=0,tail=0;q[tail++]=v;tag[v]=epoch;dep[v]=0;par[v]=-1;
 while(h<tail){int z=q[h++];++ops;if(dep[z]==a)continue;for(int w:G.a[z])if((++ops,tag[w]!=epoch)){tag[w]=epoch;dep[w]=dep[z]+1;par[w]=z;br[w]=(z==v?w:br[z]);q[tail++]=w;}}
 for(int ii=0;ii<tail;ii++){int x=q[ii];++ops;if(dep[x]!=a)continue;
  for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=1;}used[x]=1;
  auto dfs=[&](auto&&self,int z,int l)->void{
   ++ops;
   if(l>=2&&!(l&1)&&tag[z]==epoch&&dep[z]==a&&x<z&&br[x]!=br[z]){
    bool ok=true;for(int w=par[z];w!=-1;w=par[w])if((++ops,used[w])){ok=false;break;}if(ok)++C[g-2+l];}
   if(l==R)return;
   for(int w:G.a[z])if((++ops,!used[w]&&!blocked[w])){used[w]=1;self(self,w,l+1);used[w]=0;}
  };dfs(dfs,x,0);used[x]=0;for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=0;}
 }
 total_credits+=ops;max_root_credits=max(max_root_credits,ops);
 }
 vector<ull> counts()const{if(!done())throw runtime_error("exact incomplete");auto c=C;for(int t=g;t<=T;t+=2){if((2*c[t])%t)throw runtime_error("nonintegral count");c[t]=2*c[t]/t;}return c;}
};
vector<ull> exact(const Graph&G,int g,int T){ExactRoots E(G,g,T);while(!E.done())E.next_root();return E.counts();}
vector<ull> kb(const Graph&G,int g,int T){
 int D=0;for(int u=0;u<G.n;u++)D=max(D,(int)G.a[u].size());int E=G.to.size();vector<ull>C(T+1),cur(E*D),nxt(E*D),sum(G.M*D);vector<int>active,na,nodes;vector<char>touch(G.M);
 for(int root=0;root<G.n;root++){active.clear();for(int j=0;j<(int)G.out[root].size();j++){int e=G.out[root][j];cur[e*D+j]=1;active.push_back(e);}for(int l=1;l<T;l++){nodes.clear();for(int e:active){int w=G.to[e];if(!touch[w]){touch[w]=1;nodes.push_back(w);}for(int j=0;j<D;j++)sum[w*D+j]+=cur[e*D+j];}na.clear();for(int w:nodes){if(w==root)continue;for(int e:G.out[w]){int z=G.to[e];if(z<G.n&&z<root)continue;bool nonzero=false;for(int j=0;j<D;j++){ull val=sum[w*D+j]-cur[(e^1)*D+j];nxt[e*D+j]=val;nonzero|=val!=0;}if(nonzero)na.push_back(e);}}if((l+1)%2==0&&l+1>=g){for(int j=0;j<(int)G.out[root].size();j++){int e=G.out[root][j]^1;for(int k=0;k<D;k++)if(k!=j)C[l+1]+=nxt[e*D+k];}}for(int e:active)for(int j=0;j<D;j++)cur[e*D+j]=0;for(int w:nodes){touch[w]=0;for(int j=0;j<D;j++)sum[w*D+j]=0;}cur.swap(nxt);active.swap(na);}for(int e:active)for(int j=0;j<D;j++)cur[e*D+j]=0;}for(auto&x:C){if(x%2)throw runtime_error("odd KB count");x/=2;}return C;}
struct Sampler{
 const Graph&G;int g,r,A=0;ull S=0;mt19937_64 rng;
 Sampler(const Graph&gr,int gg,int rr,uint64_t seed):G(gr),g(gg),r(rr),rng(seed){}
 void trial(){int origin=uniform_int_distribution<int>(0,G.n-1)(rng),z=origin,back=-1;for(int depth=0;depth<g;depth++){int d=G.out[z].size(),j=uniform_int_distribution<int>(0,d-1-(back!=-1))(rng);if(back!=-1&&j>=back)++j;int e=G.out[z][j];z=G.to[e];back=G.pos[e^1];}++S;A+=(z==origin);}
};
// Atomic roots avoid the heavy per-step resumable state machine.
// Costs are work credits. An atomic root may overshoot Q; the observed
// maximum is reported for every graph.
struct HybridOut{double result,time_ms;bool random;ull samples;int roots;ull exact_credits,random_credits,max_root;};
HybridOut hybrid(const Graph&G,int g,int r,uint64_t seed,ull Q,bool exact_only=false){
 ExactRoots E(G,g,g);Sampler R(G,g,r,seed);
 double B=G.a[0].size()*pow(G.a[0].size()-1,g/2-1)*pow(G.a[G.n].size()-1,g/2);
 ull rw=0;
 auto start=Clock::now();
 while(true){
  ull slice=0;
  do {ull before=E.total_credits;E.next_root();slice+=E.total_credits-before;
    if(E.done()){auto C=E.counts();return {(double)C[g],chrono::duration<double,milli>(Clock::now()-start).count(),false,R.S,E.root,E.total_credits,rw,E.max_root_credits};}
  }while(slice<Q);
  if(exact_only)continue;
  slice=0;
  do {R.trial();ull credits=1+8ULL*g;rw+=credits;slice+=credits;
    if(R.A>=r)return {((double)G.n*B/g)*(r-1)/(R.S-1),chrono::duration<double,milli>(Clock::now()-start).count(),true,R.S,E.root,E.total_credits,rw,E.max_root_credits};
  }while(slice<Q);
 }
}
int main(int argc,char**argv){
 if(argc<5){cerr<<"usage: fast_bench graph.edges g T exact_reps [hybrid_reps=0] [Q=4096] [seed=20261011]\n";return 2;}
 Graph G(argv[1]);int g=stoi(argv[2]),T=stoi(argv[3]),reps=stoi(argv[4]),hr=argc>5?stoi(argv[5]):0;
 ull Q=argc>6?stoull(argv[6]):4096;uint64_t seed=argc>7?stoull(argv[7]):20261011;
 auto ref=exact(G,g,T),kk=kb(G,g,T);if(ref!=kk){cerr<<"COUNT MISMATCH\n";return 3;}
 vector<double>et,kt;for(int k=0;k<reps;k++){auto st=Clock::now();auto c=exact(G,g,T);et.push_back(chrono::duration<double,milli>(Clock::now()-st).count());if(c!=ref)abort();st=Clock::now();c=kb(G,g,T);kt.push_back(chrono::duration<double,milli>(Clock::now()-st).count());if(c!=ref)abort();}
 sort(et.begin(),et.end());sort(kt.begin(),kt.end());
 cout<<setprecision(12)<<"{\"g\":"<<g<<",\"T\":"<<T<<",\"n\":"<<G.n<<",\"counts\":[";
 for(int t=g;t<=T;t+=2){if(t>g)cout<<",";cout<<ref[t];}
 cout<<"],\"exact_ms\":"<<et[reps/2]<<",\"kb_ms\":"<<kt[reps/2];
 if(hr>0){if(T!=g)throw runtime_error("hybrid only for T=g");
 auto alone=hybrid(G,g,101,seed,Q,true);
 if((ull)alone.result!=ref[g])throw runtime_error("resumable mismatch");
 double elapsed=0,err=0;int rwins=0;ull maxwork=0;double maxratio=0;
 for(int j=0;j<hr;j++){
  auto h=hybrid(G,g,101,seed+(uint64_t)j*0x9e3779b97f4a7c15ULL,Q);
  elapsed+=h.time_ms;err+=abs(h.result-(double)ref[g])/ref[g];rwins+=h.random;
  ull work=h.exact_credits+h.random_credits;maxwork=max(maxwork,work);
  maxratio=max(maxratio,double(work)/alone.exact_credits);
  if(work>2*alone.exact_credits+2*Q+2*alone.max_root)throw runtime_error("WORK RACE BOUND VIOLATED");
 }
 cout<<",\"hybrid_ms\":"<<elapsed/hr<<",\"MARE_percent\":"<<100*err/hr
 <<",\"random_wins\":"<<rwins<<",\"runs\":"<<hr<<",\"Q\":"<<Q
 <<",\"standalone_root_exact_ms\":"<<alone.time_ms<<",\"standalone_exact_credits\":"<<alone.exact_credits
 <<",\"max_root_credits\":"<<alone.max_root<<",\"max_total_credits\":"<<maxwork<<",\"max_credit_ratio\":"<<maxratio;
 }
 cout<<"}\n";
}