// Median-of-seven exact vs. KB, every even target g,...,T.
// g++ -std=c++17 -O3 src/exact_benchmark.cpp -o exact_benchmark
#include "bench_core.cpp"
int main(int argc,char**argv){
 if(argc<4){cerr<<"usage: exact_benchmark graph.edges g T [reps=7]\n";return 2;}
 Graph G(argv[1]);int g=stoi(argv[2]),T=stoi(argv[3]),reps=argc>4?stoi(argv[4]):7;
 if(reps<1||!(reps&1)||T<g||((T-g)&1)||T>=2*g)
  throw runtime_error("invalid target or repetitions");
 auto ref=exact(G,g,T),check=kb(G,g,T);
 if(ref!=check)throw runtime_error("Exact/KB counts disagree");
 vector<double>et,kt;
 for(int j=0;j<reps;j++){
  auto s=Clock::now();auto c=exact(G,g,T);
  et.push_back(chrono::duration<double,milli>(Clock::now()-s).count());
  if(c!=ref)throw runtime_error("Exact counts changed");
  s=Clock::now();c=kb(G,g,T);
  kt.push_back(chrono::duration<double,milli>(Clock::now()-s).count());
  if(c!=ref)throw runtime_error("KB counts changed");
 }
 sort(et.begin(),et.end());sort(kt.begin(),kt.end());
 cout<<setprecision(12)<<"{\"g\":"<<g<<",\"T\":"<<T<<",\"n\":"<<G.n<<",\"counts\":[";
 for(int t=g;t<=T;t+=2){if(t>g)cout<<",";cout<<ref[t];}
 cout<<"],\"exact_ms\":"<<et[reps/2]<<",\"kb_ms\":"<<kt[reps/2]<<"}\n";
}
