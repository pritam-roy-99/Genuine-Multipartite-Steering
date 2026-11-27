"""Search for smaller independent stabilizer sets for rank-certified pairs.
An even-mask span of rank m-1 proves case (i) on every trusted cut.
The output is an exact GF(2) certificate, not an optimality claim."""
import csv,random,itertools
rng=random.Random(20261127)
def basis(vs):
 p={}
 for v in vs:
  while v:
   i=v.bit_length()-1
   if i in p:v^=p[i]
   else:p[i]=v;break
 return list(p.values())
def zmap(s,g):
 z=0
 while s:
  t=s&-s;z^=g[t.bit_length()-1];s^=t
 return z
def kernel(g):
 # A=0,B=1; eliminate the 2-bit signatures of trusted generators
 p={};vs=[]
 for j in range(2,len(g)):
  s=1<<j; z=g[j]&3
  while z:
   i=z.bit_length()-1
   if i in p: z^=p[i][0];s^=p[i][1]
   else: p[i]=(z,s);break
  if not z:vs.append(s)
 return vs
def combo(coeff,K):
 x=0
 for i,v in enumerate(K):
  if coeff>>i&1:x^=v
 return x
def omrank(Q,g):
 Z=[zmap(s,g) for s in Q]
 return len(basis([(s&Z[j])^(Q[j]&Z[i]) for i,s in enumerate(Q) for j in range(i+1,len(Q))]))
seen=set(); out=[]
for row in csv.DictReader(open('random_graphs.csv')):
 n=int(row['n']);m=n-2
 if n in seen or row['fixed_pair_pass']!='1':continue
 seen.add(n);g=list(map(int,row['adjacency_masks'].split(';')));K=kernel(g)
 assert len(K)==n-4 and omrank(K,g)==m-1
 found=K
 for r in range(2,len(K)):
  if r*(r-1)//2<m-1:continue
  for j in range(3000):
   cand=basis([combo(rng.randrange(1,1<<len(K)),K) for _ in range(r)])
   if len(cand)==r and omrank(cand,g)==m-1:
    found=cand;break
  if len(found)==r:break
 out.append([n,row['index'],len(K),len(found),';'.join(map(str,found)),omrank(found,g)])
 print(out[-1])
with open('sparse_checks.csv','w') as f:
 w=csv.writer(f);w.writerow(['n','sample_index','kernel_rank','selected_checks','generator_labels','local_anticommutation_mask_rank']);w.writerows(out)
