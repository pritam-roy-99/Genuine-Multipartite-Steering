import networkx as nx,itertools,csv
from collections import defaultdict

def check(g,a,b):
 n=len(g); R=((1<<n)-1)^(1<<a)^(1<<b); adj=[sum(1<<j for j in g[i]) for i in range(n)]; z=[0]*(1<<n)
 for s in range(1,1<<n):
  l=s&-s; z[s]=z[s^l]^adj[l.bit_length()-1]
 def sig(s,i): return ((s>>i)&1)+2*((z[s]>>i)&1)
 reps={}; piv={}
 for s in range(1,1<<n):
  k=(sig(s,a),sig(s,b)); reps.setdefault(k,s)
  if k==(0,0):
   v=s
   while v:
    p=v.bit_length()-1
    if p in piv:v^=piv[p]
    else:piv[p]=v;break
 K=list(piv.values())
 def om(v,w):return (v&z[w])^(w&z[v])
 unresolved=[]; S=(R-1)&R
 while S:
  if S&(R&-R) and not any((om(v,w)&S).bit_count()%2 for v,w in itertools.combinations(K,2)):unresolved.append(S)
  S=(S-1)&R
 ok=[]; axes=[(1,2),(1,3),(2,3)]
 for ax,bx in itertools.product(axes,repeat=2):
  keys=list(itertools.product(ax,bx))
  if any(k not in reps for k in keys):continue
  ps=[reps[k] for k in keys]
  if all(any(sum((om(v,w)&S).bit_count()%2 for w in ps)>=2 for v in K) for S in unresolved):ok.append((ax,bx))
 return bool(ok),((1,3),(1,3)) in ok
if __name__=='__main__':
 stats=defaultdict(lambda:[0]*6)
 with open('atlas_results.csv','w') as f:
  w=csv.writer(f); w.writerow(['n','graph6','pairs','certified_some_frame','certified_XY'])
  for g in nx.graph_atlas_g():
   n=len(g)
   if n<3 or not nx.is_connected(g):continue
   r=[check(g,a,b) for a,b in itertools.combinations(g,2)]; v=[len(r),sum(x for x,y in r),sum(y for x,y in r)]; s=stats[n];s[0]+=1;s[1]+=any(x for x,y in r);s[2]+=all(x for x,y in r);s[3]+=v[0];s[4]+=v[1];s[5]+=v[2];w.writerow([n,nx.to_graph6_bytes(g,header=False).decode().strip(),*v])
 with open('atlas_summary.csv','w') as f:
  w=csv.writer(f);w.writerow(['n','graphs','at_least_one_pair','all_pairs','pairs','certified_some_frame','certified_XY']);w.writerows([n,*s] for n,s in sorted(stats.items()))
 print(dict(stats))
