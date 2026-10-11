import random,itertools,pathlib,hashlib
P=pathlib.Path(__file__).resolve().parent.parent/'graphs';P.mkdir(exist_ok=True)
def save(name,n,m,edges):
 edges=sorted(edges)
 with open(P/f'{name}.edges','w') as f:
  f.write(f'{n} {m} {len(edges)}\n')
  f.writelines(f'{u} {w}\n' for u,w in edges)
 print(name,n,m,hashlib.sha256((P/f'{name}.edges').read_bytes()).hexdigest(),flush=True)
def affine(q,dv,dc):return dc*q,dv*q,[(x*q+y,a*q+(y-a*x)%q) for x in range(dc) for y in range(q) for a in range(dv)]
def lift(n,m,e,s,kind,seed):
 rng=random.Random(seed);out=[]
 for j,(u,w) in enumerate(e):
  perm=list(range(s))
  if kind=='random':rng.shuffle(perm)
  elif j==0:perm=[(i+1)%s for i in range(s)]
  out.extend((u*s+i,w*s+perm[i]) for i in range(s))
 return n*s,m*s,out
for base,q,dv,dc,s,kind in [('A',7,3,6,16,'cyclic'),('A',7,3,6,128,'cyclic'),('A',7,3,6,128,'random'),('B',5,3,4,128,'cyclic'),('C',11,4,8,128,'cyclic'),('C',11,4,8,128,'random')]:
 n,m,e=affine(q,dv,dc);save(f'{base}{s}{kind[0]}',*lift(n,m,e,s,kind,20261008+s+q*100))
q=5
norm=lambda x:tuple((a*pow(next(t for t in x if t),-1,q))%q for a in x)
points=sorted(set(norm(x) for x in itertools.product(range(q),repeat=4) if any(x)))
idx={x:i for i,x in enumerate(points)};lines=set()
for i,x in enumerate(points):
 for y in points[i+1:]:
  if (x[0]*y[2]+x[1]*y[3]-x[2]*y[0]-x[3]*y[1])%q:continue
  line=tuple(sorted({idx[norm(tuple((a*x[j]+b*y[j])%q for j in range(4)))] for a,b in itertools.product(range(q),repeat=2) if a or b}))
  assert len(line)==6;lines.add(line)
lines=sorted(lines);inc=[[] for _ in points]
for j,line in enumerate(lines):
 for i in line:inc[i].append(j)
e=[]
for i,ls in enumerate(inc):
 assert len(ls)==6
 for k,j in enumerate(ls):e.append((2*i+k//3,j))
n,m=2*len(points),len(lines)
for s,kind in [(1,'cyclic'),(16,'cyclic')]:save(f'D{s}{kind[0]}',*lift(n,m,e,s,kind,20261008+s+500))
E=[[0]*6,[0,3,14,18,24,26],[0,19,62,107,170,224]]
for name,L,z,s,r,j in [('Q410',4,79,128,2,0),('Q412',4,142,64,0,3)]:
 Z=z*s;n=L*Z;m=3*Z;ex=[a[:L] for a in E];ex[r][j]+=z
 save(name,n,m,[(col*Z+i,row*Z+(i+ex[row][col])%Z) for col in range(L) for i in range(Z) for row in range(3)])
