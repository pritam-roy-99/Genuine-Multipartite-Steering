"""Feasible product strategies for the five-generator benchmark; not an upper bound."""
import numpy as np,itertools,csv
I=np.eye(2); X=np.array([[0.,1],[1,0]]); Z=np.diag([1.,-1]); rng=np.random.default_rng(20261127)
strings=['XZIIZ','ZXZII','IZXZI','IIZXZ','ZIIZX']
def kron(ls):
 out=np.ones((1,1))
 for x in ls:out=np.kron(out,x)
 return out
def mats(s,theta):
 return np.array([kron([I if st[i]=='I' else Z if st[i]=='Z' else (np.cos(theta[i])*Z+np.sin(theta[i])*X if i<2 else X) for i in s]) for st in strings])
def ev(M,v):return np.einsum('i,kij,j->k',v,M,v)
def opt(M,w):return np.linalg.eigh(np.einsum('k,kij->ij',w,M))[1][:,-1]
rows=[]
for mask in range(1,16):
 S=[i for i in range(5) if mask>>i&1];T=[i for i in range(5) if i not in S];best=-1;record=None
 for run in range(12):
  theta=rng.uniform(0,2*np.pi,2); v=rng.normal(size=2**len(S));v/=np.linalg.norm(v);w=rng.normal(size=2**len(T));w/=np.linalg.norm(w);old=-100
  for step in range(80):
   A,B=mats(S,theta),mats(T,theta);v=opt(A,ev(B,w));w=opt(B,ev(A,v))
   for i in range(2):
    ts=theta.copy();ts[i]=0;z=ev(mats(S,ts),v)[i]*ev(mats(T,ts),w)[i];ts[i]=np.pi/2;x=ev(mats(S,ts),v)[i]*ev(mats(T,ts),w)[i];theta[i]=np.arctan2(x,z)
   score=ev(mats(S,theta),v)@ev(mats(T,theta),w)
   if abs(score-old)<1e-11:break
   old=score
  if score>best:best=score;record=(theta.copy(),v.copy(),w.copy())
 rows.append([mask,best,*[','.join(map(str,a)) for a in record]])
 print(mask,best,flush=True)
with open('generator_benchmark.csv','w') as f:
 writer=csv.writer(f);writer.writerow(['cut_mask','feasible_value','angles','state_left','state_right']);writer.writerows(rows)
