"""Analytic Pauli propagation for a five-CZ preparation, plus shot budgets."""
import csv,math
from scipy.optimize import brentq
edges=[(0,1),(2,3),(1,2),(3,4),(4,0)]
strings=['XZIIZ','ZXZII','IZXZI','IIZXZ','ZIIZX','XXYIY','XXXXX','XYZZY','XYIYX','YXYZZ','YXXYI','YYZIZ','YYIXI']
def info(s):
 x=[int(c in 'XY') for c in s];z=[int(c in 'ZY') for c in s];d=0
 for a,b in reversed(edges):
  d+=bool(x[a] or z[a] or x[b] or z[b]);z[a]^=x[b];z[b]^=x[a]
 assert not any(z),s
 return d,sum(c!='I' for c in s)
params=[info(s) for s in strings]
b=2*math.sqrt(2);lam=4-b
bounds=[3+math.sqrt(2),b,b];groups=[[2,2,1],[1,1,1,1+lam],[.5]*6+[1]]
weights=list(map(sum,groups))
def scores(q,r,v=1):
 p=[v*(1-q)**d*(1-2*r)**w for d,w in params];return [sum(p[:5]),sum(p[5::2])-lam*(1-p[3]),sum(p[5:])/2]
# p[5::2] indices 5 7 9 11

def shots(s,B,gs,delta=.01):
 return sum(math.ceil(8*math.log(1/delta)*sum(gs)*g/(s-B)**2) for g in gs) if s>B else None
if __name__=='__main__':
 for q,r,v in [(0,0,.9),(.01,.01,1),(.02,.02,1),(.03,.02,1)]:
  ss=scores(q,r,v);print(q,r,v,ss,[shots(s,B,L) for s,B,L in zip(ss,bounds,groups)])
 print('optimistic generator budgets', [shots(scores(q,r,v)[0],2+math.sqrt(5),groups[0]) for q,r,v in [(0,0,.9),(.01,.01,1)]])
 print('gate thresholds at r=.01',[brentq(lambda q:scores(q,.01)[i]-bounds[i],0,.2) for i in range(3)])
 with open('noise_budget.csv','w') as f:
  w=csv.writer(f);w.writerow(['gate_q','readout_r','generator','penalty','projected','generator_shots','penalty_shots','projected_shots']);
  for q in [i/1000 for i in range(101)]:
   ss=scores(q,.01);w.writerow([q,.01,*ss,*[shots(s,B,L) for s,B,L in zip(ss,bounds,groups)]])
