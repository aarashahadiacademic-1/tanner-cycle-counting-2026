// Bounded microbatch interleaving, IEEE ITW experiment, C++17 -O3.
// No per-operation step() calls. Work credits, not wall-clock quanta.
// The exact branch finishes whole BFS/x-DFS microtasks only after checking
// a worst-case bound that fits within the remaining Q credits.
// Compiles with: g++ -std=c++17 -O3 src/chunk_hybrid.cpp -o chunk_hybrid
#include "bench_core.cpp"
struct ChunkExact{
 const Graph&G;int g,a,root=0,epoch=0,tail=0,qi=0,delta;
 vector<int>tag,dep,par,br,que;vector<uint8_t>used,blocked;
 ull C=0,credits=0;bool needBFS=true,finished=false;
 ull bfsBound,xBound,maxObservedBfs=0,maxObservedX=0;
 ChunkExact(const Graph&gr,int gg):G(gr),g(gg),a(gg/2-1),tag(gr.M),dep(gr.M),par(gr.M),br(gr.M),que(gr.M),used(gr.M),blocked(gr.M){
  delta=0;for(auto &v:G.a)delta=max(delta,(int)v.size());
  ull power=1,pops=0,checks=0;
  for(int j=0;j<a;j++){pops+=power;power*=delta;checks+=power;}
  bfsBound=4+pops+checks;
  xBound=4+2*a+2*delta+(2+a)*(ull)delta*delta;
 }
 ull computeBFS(){
  int v=root;ull ops=1;++epoch;tail=1;int h=0;que[0]=v;
  tag[v]=epoch;dep[v]=0;par[v]=-1;
  while(h<tail){int z=que[h++];++ops;if(dep[z]==a)continue;
   for(int w:G.a[z]){++ops;if(tag[w]!=epoch){
    tag[w]=epoch;dep[w]=dep[z]+1;par[w]=z;
    br[w]=(z==v?w:br[z]);que[tail++]=w;}}
  }
  qi=0;maxObservedBfs=max(maxObservedBfs,ops);
  if(ops>bfsBound)throw runtime_error("BFS bound violated");
  return ops;
 }
 ull computeX(int x){
  ull ops=1;if(dep[x]!=a)return ops;
  for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=1;}
  used[x]=1;
  auto dfs=[&](auto&&self,int z,int l)->void{
   ++ops;
   if(l==2&&tag[z]==epoch&&dep[z]==a&&x<z&&br[x]!=br[z]){
    bool ok=true;
    for(int w=par[z];w!=-1;w=par[w])if((++ops,used[w])){ok=false;break;}
    if(ok)++C;
   }
   if(l==2)return;
   for(int w:G.a[z])if((++ops,!used[w]&&!blocked[w])){
    used[w]=1;self(self,w,l+1);used[w]=0;
   }
  };
  dfs(dfs,x,0);used[x]=0;
  for(int z=par[x];z!=-1;z=par[z]){++ops;blocked[z]=0;}
  maxObservedX=max(maxObservedX,ops);
  if(ops>xBound)throw runtime_error("X bound violated");
  return ops;
 }
 void advance(ull Q){
  ull remaining=Q;
  if(bfsBound>Q||xBound>Q)
   throw runtime_error("Q too small for guaranteed atomic microtask");
  while(remaining&&!finished){
   if(root==G.n){finished=true;break;}
   if(needBFS){
    if(remaining<bfsBound)break;
    ull c=computeBFS();remaining-=c;credits+=c;needBFS=false;continue;
   }
   if(qi==tail){++root;needBFS=true;--remaining;++credits;continue;}
   int x=que[qi];ull bound=dep[x]==a?xBound:1;
   if(remaining<bound)break;
   ull c=computeX(x);
   if(c>remaining)throw runtime_error("Q overshoot");
   remaining-=c;credits+=c;++qi;
  }
  if(root==G.n)finished=true;
 }
 ull result()const{
  if(!finished)throw runtime_error("incomplete exact");
  if(2*C%g)throw runtime_error("bad count");
  return 2*C/g;
 }
};
struct BatchSampler{
 const Graph&G;int g,r,A=0;ull S=0,credits=0;mt19937_64 rng;
 BatchSampler(const Graph&gr,int gg,int rr,uint64_t seed):G(gr),g(gg),r(rr),rng(seed){}
 void advance(ull Q){
  ull trialCost=1+8ULL*g,remain=Q;
  if(trialCost>Q)throw runtime_error("Q too small for atomic walk");
  while(remain>=trialCost&&A<r){
   int origin=uniform_int_distribution<int>(0,G.n-1)(rng),z=origin,back=-1;
   for(int i=0;i<g;i++){
    int d=G.out[z].size(),j=uniform_int_distribution<int>(0,d-1-(back!=-1))(rng);
    if(back!=-1&&j>=back)++j;
    int e=G.out[z][j];z=G.to[e];back=G.pos[e^1];
   }
   ++S;A+=(z==origin);remain-=trialCost;credits+=trialCost;
  }
 }
};
struct Result{double value,ms;bool random;ull exactCredits,samplerCredits,samples;};
Result runHybrid(const Graph&G,int g,int r,ull Q,uint64_t seed){
 ChunkExact E(G,g);BatchSampler S(G,g,r,seed);
 double B=G.a[0].size()*pow(G.a[0].size()-1,g/2-1)*pow(G.a[G.n].size()-1,g/2);
 auto start=Clock::now();
 for(;;){
  ull before=E.credits;E.advance(Q);
  if(E.credits-before>Q)throw runtime_error("EXACT BUDGET VIOLATION");
  if(E.finished)
   return {double(E.result()),chrono::duration<double,milli>(Clock::now()-start).count(),false,E.credits,S.credits,S.S};
  before=S.credits;S.advance(Q);
  if(S.credits-before>Q)throw runtime_error("SAMPLER BUDGET VIOLATION");
  if(S.A>=r)
   return {((double)G.n*B/g)*(r-1)/(S.S-1),chrono::duration<double,milli>(Clock::now()-start).count(),true,E.credits,S.credits,S.S};
 }
}
int main(int argc,char**argv){
 if(argc<4){
  cerr<<"usage: chunk_hybrid graph.edges girth hybrid_reps [Q=4096] [seed=20261011] [baseline_reps=7]\n";
  return 2;
 }
 Graph G(argv[1]);int g=stoi(argv[2]),runs=stoi(argv[3]);
 ull Q=argc>4?stoull(argv[4]):4096;
 uint64_t seed=argc>5?stoull(argv[5]):20261011;
 int reps=argc>6?stoi(argv[6]):7;
 auto ref=exact(G,g,g)[g];
 if(kb(G,g,g)[g]!=ref)throw runtime_error("KB count mismatch");
 ChunkExact standalone(G,g);
 while(!standalone.finished)standalone.advance(Q);
 if(standalone.result()!=ref)throw runtime_error("Chunk exact count mismatch");
 vector<double>et,kt;
 for(int j=0;j<reps;j++){
  auto start=Clock::now();auto c=exact(G,g,g);
  et.push_back(chrono::duration<double,milli>(Clock::now()-start).count());
  if(c[g]!=ref)throw runtime_error("Exact changed");
  start=Clock::now();auto k=kb(G,g,g);
  kt.push_back(chrono::duration<double,milli>(Clock::now()-start).count());
  if(k[g]!=ref)throw runtime_error("KB changed");
 }
 sort(et.begin(),et.end());sort(kt.begin(),kt.end());
 double total=0,err=0,maxRatio=0;int wins=0;
 ull maxExactCredit=0,maxSampleCredit=0;
 for(int j=0;j<runs;j++){
  auto v=runHybrid(G,g,101,Q,seed+(ull)j*0x9e3779b97f4a7c15ULL);
  if(!v.random&&v.value!=ref)throw runtime_error("Invalid exact hybrid output");
  total+=v.ms;err+=abs(v.value-double(ref))/ref;wins+=v.random;
  maxExactCredit=max(maxExactCredit,v.exactCredits);
  maxSampleCredit=max(maxSampleCredit,v.samplerCredits);
  maxRatio=max(maxRatio,double(v.exactCredits+v.samplerCredits)/standalone.credits);
 }
 cout<<setprecision(12)<<"{\"g\":"<<g<<",\"n\":"<<G.n<<",\"Ng\":"<<ref
 <<",\"exact_ms\":"<<et[reps/2]<<",\"kb_ms\":"<<kt[reps/2]
 <<",\"hybrid_ms\":"<<total/runs<<",\"MARE_percent\":"<<100*err/runs
 <<",\"random_wins\":"<<wins<<",\"runs\":"<<runs<<",\"Q\":"<<Q
 <<",\"bfs_bound\":"<<standalone.bfsBound<<",\"x_bound\":"<<standalone.xBound
 <<",\"observed_max_bfs\":"<<standalone.maxObservedBfs
 <<",\"observed_max_x\":"<<standalone.maxObservedX
 <<",\"exact_credits\":"<<standalone.credits
 <<",\"max_hybrid_exact_credits\":"<<maxExactCredit
 <<",\"max_hybrid_sample_credits\":"<<maxSampleCredit
 <<",\"max_credit_ratio\":"<<maxRatio<<"}\n";
}
