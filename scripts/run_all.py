#!/usr/bin/env python3
"""Reproduce both tables from source, preserving independent fresh output.
Run from repository root after generating the 14 graph files and compiling C++.
Uses seven baseline runs, 100 hybrid runs, epsilon .08 (r=101), Q=4096.
"""
import csv, hashlib, json, pathlib, subprocess
ROOT=pathlib.Path(__file__).resolve().parent.parent
ROWS1=[('A128r',6,6),('A128c',6,8),('A128c',6,10),
 ('A128r',6,10),('C128r',6,6),('R1200',6,6),('R2400',6,6),
 ('P1200',8,8),('D1c',8,14),('P1600',10,10)]
ROWS2=[('A16c',6),('A128r',6),('A128c',6),('B128c',6),
 ('C128c',6),('D16c',8),('P1200',8),('R1200',6),
 ('P1600',10),('Q410',10),('Q412',12)]
def read_csv(name):
 with (ROOT/'results'/name).open(newline='') as f:return list(csv.DictReader(f))
source1=read_csv('table1_exact_7runs.csv')
source2=read_csv('table2_hybrid_100runs.csv')
expected_hash={r['graph']:r['graph_sha256'] for r in source1}
expected_hash.update({r['graph']:r['sha256'] for r in source2})
def check_graph(name):
 p=ROOT/'graphs'/(name+'.edges')
 if not p.is_file():raise FileNotFoundError(f'{p}; generate graphs first')
 sha=hashlib.sha256(p.read_bytes()).hexdigest()
 if sha!=expected_hash[name]:raise RuntimeError(f'{name} SHA256 mismatch: {sha}')
 return p
def run(args):
 return json.loads(subprocess.check_output([str(x) for x in args],text=True,timeout=600))
one=[];two=[]
for name,g,T in ROWS1:
 p=check_graph(name)
 print('Table I',name,g,T,flush=True)
 row=run([ROOT/'exact_benchmark',p,g,T,7])
 actual=';'.join(map(str,row['counts']))
 expected=next(x for x in source1 if x['graph']==name and int(x['T'])==T)
 if actual!=expected['counts']:raise RuntimeError(f'cycle tuple mismatch {name},T={T}')
 one.append(dict(graph=name,**row,sha256=expected_hash[name]))
for name,g in ROWS2:
 p=check_graph(name)
 print('Table II',name,g,flush=True)
 row=run([ROOT/'chunk_hybrid',p,g,100,4096,20261011,7])
 expected=next(x for x in source2 if x['graph']==name)
 if row['Ng']!=int(expected['Ng']):raise RuntimeError(f'girth count mismatch {name}')
 if row['observed_max_bfs']>row['bfs_bound'] or row['observed_max_x']>row['x_bound']:
  raise RuntimeError(f'microbatch credit-bound error {name}')
 two.append(dict(graph=name,**row,sha256=expected_hash[name]))
for name,items in [('recomputed_table1.json',one),('recomputed_table2.json',two)]:
 (ROOT/'results'/name).write_text(json.dumps(items,indent=2)+'\n')
print('PASS: all cycle counts and SHA256 values match; timings are fresh, not expected to match historical medians.')
