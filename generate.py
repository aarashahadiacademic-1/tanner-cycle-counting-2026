import random,json,itertools,pathlib
P=pathlib.Path(__file__).parent
(P/'graphs').mkdir(exist_ok=True)
manifest=[]
def save(name,n,m,edges,**meta):
 edges=sorted(edges)
 assert len(edges)==len(set(edges))
 with open(P/'graphs'/f'{name}.edges','w') as f:
  f.write(f'{n} {m} {len(edges)}\n')
  for u,w in edges:f.write(f'{u} {w}\n')
 manifest.append(dict(name=name,n=n,m=m,**meta))
def affine(q,dv,dc):
 return dc*q,dv*q,[(x*q+y,a*q+(y-a*x)%q) for x in range(dc) for y in range(q) for a in range(dv)]
def lift(n,m,edges,s,kind,seed):
 rng=random.Random(seed); out=[]
 for j,(u,w) in enumerate(edges):
  perm=list(range(s))
  if kind=='random':rng.shuffle(perm)
  elif j==0:perm=[(i+1)%s for i in range(s)]
  out.extend((u*s+i,w*s+perm[i]) for i in range(s))
 return n*s,m*s,out
for base,q,dv,dc in [('A',7,3,6),('B',5,3,4),('C',11,4,8)]:
 n,m,e=affine(q,dv,dc)
 for s in ([1,2,4,8,16,32,128,256] if base=='A' else [128]):
  for kind in (['cyclic','random'] if s>=128 else ['cyclic']):
   seed=20261008+s+q*100
   nn,mm,ee=lift(n,m,e,s,kind,seed)
   save(f'{base}{s}{kind[0]}',nn,mm,ee,dv=dv,dc=dc,base=base,lift=s,kind=kind,seed=seed,T=6)
# Symplectic generalized quadrangle W(5); split each degree-6 point into two degree-3 variables.
q=5
norm=lambda x: tuple((a*pow(next(t for t in x if t),-1,q))%q for a in x)
points=sorted(set(norm(x) for x in itertools.product(range(q),repeat=4) if any(x)))
idx={x:i for i,x in enumerate(points)}; lines=set()
for i,x in enumerate(points):
 for y in points[i+1:]:
  if (x[0]*y[2]+x[1]*y[3]-x[2]*y[0]-x[3]*y[1])%q:continue
  line=tuple(sorted({idx[norm(tuple((a*x[j]+b*y[j])%q for j in range(4)))] for a,b in itertools.product(range(q),repeat=2) if a or b}))
  assert len(line)==6
  lines.add(line)
lines=sorted(lines); inc=[[] for _ in points]
for j,line in enumerate(lines):
 for i in line:inc[i].append(j)
e=[]
for i,ls in enumerate(inc):
 assert len(ls)==6
 for k,j in enumerate(ls):e.append((2*i+k//3,j))
n,m=2*len(points),len(lines)
for s,kind in [(1,'cyclic'),(16,'random'),(16,'cyclic')]:
 seed=20261008+s+500
 nn,mm,ee=lift(n,m,e,s,kind,seed)
 save(f'D{s}{kind[0]}',nn,mm,ee,dv=3,dc=6,base='D',lift=s,kind=kind,seed=seed,T=8)
# Small overlapping-cycle correctness instance; targets exceed 2g.
n,m,e=affine(5,3,4);save('validation_B',n,m,e,dv=3,dc=4,T=14,validation=True)
# Irregular cactus unit tests, not coding-performance benchmarks.
for g in [10,12]:
 for k in [200,1000]:
  nv=0;ordinary=[];anchors=[]
  for block in range(k):
   last=None
   for length in range(g//2,g+2):
    vs=list(range(nv,nv+length));nv+=length
    ordinary.extend((vs[j],vs[(j+1)%length]) for j in range(length))
    if last is not None:ordinary.append((last,vs[0]))
    last=vs[0]
   anchors.append(last)
  ordinary.extend((anchors[i],anchors[i+1]) for i in range(k-1))
  e=[(v,j) for j,(u,w) in enumerate(ordinary) for v in (u,w)]
  save(f'cactus{g}_{k}',nv,len(ordinary),e,T=2*g+2,validation=True,expected=k)
json.dump(manifest,open(P/'manifest.json','w'),indent=2)
