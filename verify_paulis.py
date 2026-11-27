"""Check the finite Pauli identities used in the analytical examples.
This is not a graph-atlas enumeration or a proof of the general theorems.
"""
import numpy as np
from itertools import product
P={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
def kron(s):
 r=np.array([[1.]])
 for q in s: r=np.kron(r,P[q])
 return r
def graph_state(n,edges):
 return np.array([(-1)**sum(z[i]*z[j] for i,j in edges) for z in product([0,1],repeat=n)])/np.sqrt(2**n)
def op(s): return (1 if s[0]=='+' else -1)*kron(s[1:])
def check(n,edges,strings):
 v=graph_state(n,edges)
 for s in strings:
  assert np.linalg.norm(op(s)@v-v)<1e-12,s
 return v
c5=['+XXYIY','-XXXXX','-XYZZY','+XYIYX','-YXYZZ','+YXXYI','+YYZIZ','+YYIXI']
check(5,[(i,(i+1)%5) for i in range(5)],c5+['+IIZXZ'])
q=op('+IIZXZ')
for j in range(0,8,2): assert np.linalg.norm(op(c5[j])@q-op(c5[j+1]))<1e-12
assert sorted(sum(c!='I' for c in s[1:]) for s in c5)==[3,4,4,4,4,5,5,5]
prism=['-XXXZZZ','-XYIYIZ','-YXXYII','+YYIZZI','+IIXYYZ','+IIZZZX']
check(6,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)],prism)
def anti(s,t,subset): return sum(s[i]!=t[i] and s[i]!='I' and t[i]!='I' for i in subset)%2
for S in [(2,),(2,3),(2,4),(2,3,4),(2,5),(2,3,5),(2,4,5)]:
 if anti(prism[4][1:],prism[5][1:],S): continue
 assert any(sum(anti(q[1:],s[1:],S) for s in prism[:4])>=2 for q in prism[4:]),S
for S in [(2,),(3,),(4,)]: assert sum(anti('IIZXZ',c5[j][1:],S) for j in range(0,8,2))==2
print('C5 eight strings, projector products, weights, anticommutation patterns; prism stabilizers and cut certificates: PASS')
from scipy.optimize import brentq
p=brentq(lambda p: ((1-p)**3+4*(1-p)**4+3*(1-p)**5)/2-2*np.sqrt(2),0,.2)
print(f'Local depolarization threshold p={p:.12f}')

# Independent check of the explicit saturating strategy in the Q=J=+1 space.
Q=kron('ZXZ'); J=kron('XIX'); PP=(np.eye(8)+Q)@(np.eye(8)+J)/4
e,V=np.linalg.eigh(PP); V=V[:,e>.5]; assert V.shape==(8,2)
T=[V.conj().T@op(s)@V for s in ['+YIY','-ZZY','-YZZ','+ZIZ']]
assert np.allclose(T[2],T[1]) and np.allclose(T[3],-T[0])
assert np.allclose(T[0]@T[1]+T[1]@T[0],0)
L0=(T[0]+T[2])/np.sqrt(2); L1=(T[1]+T[3])/np.sqrt(2)
for L in [L0,L1]:assert np.allclose(L@L,np.eye(2))
phi=np.array([1,0,0,1])/np.sqrt(2)
score=phi.conj()@(np.kron(L0.T,T[0]+T[2])+np.kron(L1.T,T[1]+T[3]))@phi
assert np.allclose(score,2*np.sqrt(2));print('Explicit biseparable saturator:',score.real)
# Exact rational check of the fixed-frame counterexample.
from fractions import Fraction as F
h=F(27,32);u=h*h;d=F(7,10);e=F(79,80)
assert 2*u*u>1 and d+e==2*h and d+e<=1+u and abs(d-e)<=1-u and d*d+u*u<=1
ss=list(product([1,-1],repeat=2))
def tau(b,y,s):return (1+b*s[0]*u,0,d*s[0]+e*b) if y==0 else (1,u*s[1]*b,d*s[0]) # 32 tau in I,X,Z coordinates
for a,b,x,y in product([1,-1],[1,-1],range(2),range(2)):
 arr=[tau(b,y,s) for s in ss if s[x]==a]+[tau(a,x,s) for s in ss if s[y]==b]
 lhs=tuple(sum(q[i] for q in arr) for i in range(3))
 rhs=(1+a*b*u,0,(a+b)*h) if x==y==0 else (1,0,a*h) if x==0 else (1,0,b*h) if y==0 else (1,a*b*u,0)
 assert lhs==tuple(4*v for v in rhs)
print('Fixed-frame biseparable certificate: exact rational checks PASS')
