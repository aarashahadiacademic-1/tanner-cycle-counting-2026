#!/usr/bin/env python3
"""Rerun same-target exact/KB benchmarks and one-processor root-batched hybrid.
Requires g++ 14+, Python 3, and no third-party packages.
Recorded timings will depend on host and load.
"""
import csv, hashlib, json, pathlib, subprocess
P = pathlib.Path(__file__).resolve().parent
EXACT = [
 ("A128r",6,6),("A128c",6,8),("A128c",6,10),
 ("A128r",6,10),("C128r",6,6),("R1200",6,6),
 ("R2400",6,6),("P1200",8,8),("D1c",8,14),("P1600",10,10)]
HYBRID = [
 ("A16c",6),("A128r",6),("A128c",6),("B128c",6),
 ("C128c",6),("D16c",8),("R1200",6),("P1200",8),
 ("P1600",10),("Q410",10),("Q412",12)]
DEGREES={"A16c":(3,6),"A128r":(3,6),"A128c":(3,6),
 "B128c":(3,4),"C128c":(4,8),"C128r":(4,8),
 "D16c":(3,6),"D1c":(3,6),"R1200":(4,8),
 "R2400":(3,6),"P1200":(3,6),"P1600":(3,4),
 "Q410":(3,4),"Q412":(3,4)}
def run(name,g,T,hybrid):
    path=P/"graphs"/(name+".edges")
    if not path.exists():
        if name in ("Q410","Q412"):
            subprocess.run(["python3",str(P/"generate_qc.py")],check=True)
        if not path.exists(): raise FileNotFoundError(path)
    cmd=[str(P/"fast_bench"),str(path),str(g),str(T),"7",
         "100" if hybrid else "0","4096","20261011"]
    row=json.loads(subprocess.check_output(cmd,text=True))
    row.update(graph=name,dv=DEGREES[name][0],dc=DEGREES[name][1],
      graph_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return row
def write(name,cols,rows):
    with (P/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def main():
    subprocess.run(["g++","-std=c++17","-O3",str(P/"fast_bench.cpp"),"-o",str(P/"fast_bench")],check=True)
    table1=[]
    for name,g,T in EXACT:
        a=run(name,g,T,False)
        table1.append(dict(graph=name,dv=a["dv"],dc=a["dc"],g=g,n=a["n"],T=T,
            counts="("+", ".join(map(str,a["counts"]))+")",
            exact_ms=round(a["exact_ms"],3),KB_ms=round(a["kb_ms"],3),
            speedup=round(a["kb_ms"]/a["exact_ms"],2),graph_sha256=a["graph_sha256"]))
        print("Exact",name,g,T,flush=True)
    table2=[]
    for name,g in HYBRID:
        a=run(name,g,g,True)
        table2.append(dict(graph=name,dv=a["dv"],dc=a["dc"],g=g,n=a["n"],
            Ng=a["counts"][0],KB_ms=round(a["kb_ms"],3),
            Exact_ms=round(a["exact_ms"],3),
            Hybrid_ms=round(a["hybrid_ms"],3),
            MARE_percent=round(a["MARE_percent"],2),
            random_wins=a["random_wins"],runs=a["runs"],
            max_root_credits=a["max_root_credits"],
            max_credit_ratio=round(a["max_credit_ratio"],3),
            graph_sha256=a["graph_sha256"]))
        print("Hybrid",name,g,flush=True)
    write("table1_rerun.csv",list(table1[0]),table1)
    write("table2_rerun.csv",list(table2[0]),table2)
if __name__=="__main__":main()
