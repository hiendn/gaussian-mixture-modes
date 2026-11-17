"""Bounded checks for the polynomial semidefinite mixture certificate.
New work only. Does not read earlier reports or modify source notes.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from mpmath import iv

OUT=Path(__file__).resolve().parent
iv.dps=70
D=iv.mpf('2.8'); r=iv.mpf('.68'); beta=iv.mpf('.99'); tau=iv.mpf('.01'); q=iv.mpf(2)
A=iv.exp(-D*D/2+r*D); a=iv.sqrt(3)*D/2; b=D/2; c=iv.mpf('.51')
u=iv.mpf('2.72'); v=iv.mpf('2.7'); z=iv.mpf('2.7')
H00=beta+u-(a-c)**2; H11=beta+v-b*b; H01=z-(a-c)*b
ival={
 'separation_slack':D*D-(q+1+tau+r*r+q/r**2),
 'invariance_slack':r-2*A*(D-r),
 'lambda_determinant':u*v-z*z,
 'main_slack_x':beta-c*c-2*A*u,
 'main_slack_y':beta-2*A*v,
 'second_slack_x':H00,
 'second_slack_y':H11,
 'second_determinant':H00*H11-H01**2,
}
assert all(x.a>0 for x in ival.values())
rr=D*D-q-1-tau
rmin=iv.sqrt(2*q/(rr+iv.sqrt(rr*rr-4*q)))
Amin=iv.exp(-D*D/2+D*rmin)
trace_min=2*Amin*D*D/(1+2*Amin)-3*Amin**2*D*D/(1+2*Amin)**2
assert trace_min.a>1
interval_results={k:str(value) for k,value in ival.items()}
interval_results.update({'A':str(A),'rmin':str(rmin),'trace_at_rmin':str(trace_min)})

# Actual convex-inner EM/MM step within a symmetric subfamily. The general
# derivation permits all triangular B, all v and all logits. Here B=b I,
# equal weights and an equilateral configuration make every component's
# semidefinite inequalities identical and only 2 by 2.
rng=np.random.default_rng(20261117)
angles=np.array([0,2*np.pi/3,4*np.pi/3])
V=np.c_[np.cos(angles),np.sin(angles)]/np.sqrt(3) # unit side length
base_offsets=rng.normal(size=(16,2))*.65
X=[]
for j,angle in enumerate(angles):
 R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
 X.extend(3.2*V[j]+base_offsets@R.T)
X=np.array(X)
D0=2.8;r0=.68;beta0=.99;tau0=.01
A0=np.exp(-D0**2/2+r0*D0)
z0=np.array([1.,D0,r0,.51,A0,2.72,2.7,2.7,D0-r0])
# Order: b,d,r,c,a,u,v,z,t, with [[u,z],[z,v]]=Lambda_+.
logkern=-.5*np.sum((X[:,None,:]-D0*V[None,:,:])**2,axis=2)
resp=np.exp(logkern-logsumexp(logkern,axis=1)[:,None])
gamma=.2

def Q(w):
 b0,d=w[:2]
 err=b0*X[:,None,:]-d*V[None,:,:]
 return 2*np.log(b0)-.5*np.mean(np.sum(resp*np.sum(err*err,axis=2),axis=1))-np.log(3)

def LL(w):
 b0,d=w[:2]
 err=b0*X[:,None,:]-d*V[None,:,:]
 return np.mean(logsumexp(-.5*np.sum(err*err,axis=2),axis=1)-np.log(3)+2*np.log(b0)-np.log(2*np.pi))

def Uscalar(a,b,at,bt):
 return a*b+.25*((a-at)-(b-bt))**2

def constraints(w):
 b0,d,rad,cc,aa,uu,vv,zz,tt=w
 _,dt,rt,ct,at,ut,vt,zt,tt0=z0
 s=dt-rt
 low=s*(d-rad)-.5*s*s-.5*rad*rad
 major11=aa*uu+.25*((aa-at-(uu-ut))**2+(zz-zt)**2)
 major22=aa*vv+.25*((aa-at-(vv-vt))**2+(zz-zt)**2)
 lam=np.array([[uu,zz],[zz,vv]])
 delta_c=np.array([np.sqrt(3)*d/2-cc,d/2])
 return np.r_[
  d-np.sqrt(3+tau0+rad*rad+2/(rad*rad)),
  aa-np.exp(-low),
  tt-(d-rad),
  rad-2*Uscalar(aa,tt,at,tt0),
  beta0-cc*cc-2*major11,
  beta0-2*major22,
  np.linalg.eigvalsh(lam),
  np.linalg.eigvalsh(beta0*np.eye(2)+lam-np.outer(delta_c,delta_c))]

def obj(w):
 return -Q(w)+gamma*.5*np.sum((w-z0)**2)

sol=minimize(obj,z0,method='SLSQP',bounds=[(.3,3),(2.5,7),(.4,1),(-1,1),(1e-6,1),(0,20),(0,20),(-20,20),(0,10)],
 constraints={'type':'ineq','fun':constraints},options={'ftol':2e-12,'maxiter':2000,'disp':False})
# SLSQP can finish several units in its stopping tolerance outside a cone.
# Retreat toward the strictly feasible certificate to get a checked feasible
# point. The current odds and distance epigraphs start on their boundaries;
# those are evaluated explicitly as well.
step=.999
w=z0+step*(sol.x-z0)
while constraints(w).min() < -2e-12 and step>1e-8:
 step*=.99
 w=z0+step*(sol.x-z0)
incQ=float(Q(w)-Q(z0));prox=float(gamma*.5*np.sum((w-z0)**2));incLL=float(LL(w)-LL(z0))
assert constraints(w).min()>-2e-12
assert incQ>prox and incLL+1e-10>=incQ
b1,d1,r1,c1,a1,u1,v1,z1,t1=w
Lam=np.array([[u1,z1],[z1,v1]])
Acap=np.exp(-d1*d1/2+r1*d1)
actual_main_slack=np.array([beta0-c1*c1-2*Acap*u1,beta0-2*Acap*v1])
actual_inv_slack=r1-2*Acap*(d1-r1)
second_vec=np.array([np.sqrt(3)*d1/2-c1,d1/2])
actual_second_eig=np.linalg.eigvalsh(beta0*np.eye(2)+Lam-np.outer(second_vec,second_vec))
assert actual_main_slack.min()>0 and actual_inv_slack>0 and actual_second_eig.min()>0
# Bound checks for U(a,Lambda)-a Lambda=(Y-Yt)^2/4 in random dimensions.
worst_residual=0.
for p in [1,2,4,7]:
 for _ in range(100):
  W=rng.normal(size=(p,p));L=W@W.T
  Wt=rng.normal(size=(p,p));Lt=Wt@Wt.T
  aa=float(np.exp(rng.normal()));at=float(np.exp(rng.normal()))
  Xmat=aa*np.eye(p)+L;Y=aa*np.eye(p)-L;Yt=at*np.eye(p)-Lt
  Up=.25*(Xmat@Xmat-(Yt@Y+Y@Yt-Yt@Yt))
  residual=Up-aa*L-.25*(Y-Yt)@(Y-Yt)
  worst_residual=max(worst_residual,float(np.max(np.abs(residual))))

result={
 'interval_checks':interval_results,
 'convex_inner_demo':{
  'n':len(X),'solver':'scipy SLSQP on the displayed convex inner constraints',
  'solver_success':bool(sol.success),'solver_message':str(sol.message),
  'iterations':int(sol.nit),'retained_step_fraction':step,
  'start':z0.tolist(),'next':w.tolist(),
  'min_inner_slack':float(constraints(w).min()),
  'Q_increase_per_observation':incQ,'proximal_penalty':prox,
  'log_likelihood_increase_per_observation':incLL,
  'original_main_slack':actual_main_slack.tolist(),'original_invariance_slack':float(actual_inv_slack),
  'original_second_slack_eigenvalues':actual_second_eig.tolist(),
  'end_shared_covariance_diagonal':float(b1**-2),
  'density_parameter_scope':'B=b I, equilateral means, equal weights; full algorithm derived but not numerically implemented here'
 },
 'random_matrix_identity':{'cases':400,'dimensions':[1,2,4,7],'largest_absolute_identity_residual':worst_residual}
}
(OUT/'operator_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
