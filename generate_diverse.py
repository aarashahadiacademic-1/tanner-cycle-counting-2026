import random,pathlib,json,collections,itertools
P=pathlib.Path(__file__).resolve().parent
(P/'graphs').mkdir(exist_ok=True)
def save(name,n,m,dv,dc,seed,family,U):
 edges=sorted((u,w) for u,ws in enumerate(U) for w in ws)
 assert len(edges)==n*dv and len(set(edges))==len(edges)
 W=[set() for _ in range(m)]
 for u,w in edges:W[w].add(u)
 assert all(len(x)==dv for x in U) and all(len(x)==dc for x in W)
 with (P/'graphs'/f'{name}.edges').open('w') as f:
  f.write(f'{n} {m} {len(edges)}\n');f.writelines(f'{u} {w}\n' for u,w in edges)
 return dict(name=name,n=n,m=m,dv=dv,dc=dc,seed=seed,family=family,T=6)
def bad(U):
 pairs={}
 for u,ws in enumerate(U):
  for w,z in itertools.combinations(sorted(ws),2):
   if (w,z) in pairs:return u,w
   pairs[w,z]=u
 return None
def random_regular(n,dv,dc,seed):
 m=n*dv//dc;r=random.Random(seed)
 while True:
  stubs=[w for w in range(m) for _ in range(dc)];r.shuffle(stubs)
  U=[set(stubs[u*dv:(u+1)*dv]) for u in range(n)]
  if all(len(ws)==dv for ws in U):break
 W=[set() for _ in range(m)]
 for u,ws in enumerate(U):
  for w in ws:W[w].add(u)
 def four(u,w):
  return any((U[a]-{w})&(U[u]-{w}) for a in W[w]-{u})
 switches=0
 while (e:=bad(U)):
  u,w=e
  for trial in range(100000):
   a=r.randrange(n);b=r.choice(sorted(U[a]))
   if a==u or b==w or b in U[u] or w in U[a]:continue
   U[u].remove(w);W[w].remove(u);U[a].remove(b);W[b].remove(a)
   U[u].add(b);W[b].add(u);U[a].add(w);W[w].add(a)
   if not four(u,b) and not four(a,w):switches+=1;break
   U[u].remove(b);W[b].remove(u);U[a].remove(w);W[w].remove(a)
   U[u].add(w);W[w].add(u);U[a].add(b);W[b].add(a)
  else:raise RuntimeError('repair failed')
 return U,switches
def peg(n,dv,dc,seed):
 m=n*dv//dc;r=random.Random(seed)
 for restart in range(100):
  U=[set() for _ in range(n)];W=[set() for _ in range(m)];ok=True
  for u in range(n):
   for k in range(dv):
    available=[w for w in range(m) if len(W[w])<dc and w not in U[u]]
    if k:
     dist={};seen_u={u};front={u};depth=0
     while front:
      nxt=set()
      for a in front:
       for w in U[a]:
        if w not in dist:
         dist[w]=depth+1
         for v in W[w]:
          if v not in seen_u:seen_u.add(v);nxt.add(v)
      front=nxt;depth+=2
     available=[w for w in available if dist.get(w,n+m)>3]
     if available:
      far=max(dist.get(w,n+m) for w in available);available=[w for w in available if dist.get(w,n+m)==far]
    if not available:ok=False;break
    least=min(len(W[w]) for w in available);available=[w for w in available if len(W[w])==least]
    w=r.choice(available);U[u].add(w);W[w].add(u)
   if not ok:break
  if ok:return U,restart
 raise RuntimeError('PEG saturated; all restarts failed')
def actual_girth(U,m):
 n=len(U);adj=[set(ws) for ws in U]+[set() for _ in range(m)]
 for u,ws in enumerate(U):
  adj[u]={n+w for w in ws}
  for w in ws:adj[n+w].add(u)
 best=n+m
 for root in range(n+m):
  dist={root:0};par={root:-1};q=collections.deque([root])
  while q:
   u=q.popleft()
   if 2*dist[u]>=best:continue
   for w in adj[u]:
    if w not in dist:dist[w]=dist[u]+1;par[w]=u;q.append(w)
    elif par[u]!=w and par[w]!=u:best=min(best,dist[u]+dist[w]+1)
  if best==6:break
 return best
manifest=[]
for name,n,dv,dc,seed,family in [('R1200',1200,4,8,202610091,'random-regular'),('R2400',2400,3,6,202610092,'random-regular'),('P1200',1200,3,6,202610093,'degree-constrained-PEG'),('P1600',1600,3,4,202610094,'degree-constrained-PEG')]:
 U,aux=random_regular(n,dv,dc,seed) if family=='random-regular' else peg(n,dv,dc,seed)
 row=save(name,n,n*dv//dc,dv,dc,seed,family,U);row['g']=actual_girth(U,n*dv//dc);row['T']=row['g'];row['repair_switches' if family=='random-regular' else 'restarts']=aux;manifest.append(row);print(row,flush=True)
(P/'manifest_diverse.json').write_text(json.dumps(manifest,indent=2))
