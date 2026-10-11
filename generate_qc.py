"""Reproduce column-weight-three QC benchmarks. No external dependencies.
The starting exponent matrix is Eq. (36) of Zhang & Wang, arXiv:1001.3916.
The two test instances below are our modified QC matrices,
not instances whose measured results are taken from that paper.
"""
import pathlib,json,hashlib
P=pathlib.Path(__file__).resolve().parent
E=[[0]*6,[0,3,14,18,24,26],[0,19,62,107,170,224]]
SPECS=[('Q410',4,10,79,128,2,0,40448),
       ('Q412',4,12,142,64,0,3,308992)]
def generate():
 (P/'graphs').mkdir(exist_ok=True);rows=[]
 for name,L,g,z,s,r,j,Ng in SPECS:
  Z=z*s;n=L*Z;m=3*Z;ex=[a[:L] for a in E];ex[r][j]+=z
  path=P/'graphs'/f'{name}.edges'
  with path.open('w') as f:
   f.write(f'{n} {m} {3*n}\n')
   for col in range(L):
    for i in range(Z):
     for row in range(3):f.write(f'{col*Z+i} {row*Z+(i+ex[row][col])%Z}\n')
  data=dict(name=name,n=n,m=m,dv=3,dc=L,g=g,T=g,Z=Z,base_Z=z,lift=s,
            changed_block=[r,j],exponents=ex,expected_Ng=Ng,
            lambda_expected=g*Ng/(3*2**(g//2-1)*(L-1)**(g//2)),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
  rows.append(data);print(name,n,g,Ng,flush=True)
 (P/'manifest_qc.json').write_text(json.dumps(rows,indent=2)+'\n')
 return rows
if __name__=='__main__':generate()
