import random,pathlib,itertools,hashlib
P=pathlib.Path(__file__).resolve().parent.parent/'graphs';P.mkdir(exist_ok=True)
def save(name,U,dc):
 n=len(U);m=n*len(U[0])//dc;e=sorted((u,w) for u,ws in enumerate(U) for w in ws)
 with (P/f'{name}.edges').open('w') as f:
  f.write(f'{n} {m} {len(e)}\n');f.writelines(f'{u} {w}\n' for u,w in e)
 print(name,hashlib.sha256((P/f'{name}.edges').read_bytes()).hexdigest(),flush=True)
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
 def four(u,w):return any((U[a]-{w})&(U[u]-{w}) for a in W[w]-{u})
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
 raise RuntimeError('PEG saturated')
for name,n,dv,dc,seed,family in [('R1200',1200,4,8,202610091,'r'),('R2400',2400,3,6,202610092,'r'),('P1200',1200,3,6,202610093,'p'),('P1600',1600,3,4,202610094,'p')]:
 U,aux=random_regular(n,dv,dc,seed) if family=='r' else peg(n,dv,dc,seed)
 save(name,U,dc);print('aux',aux,flush=True)
