"""One bounded joint-M / variable-nu conic step; no convergence simulation.

The curvature formula is justified in initial_report.md / followup_report.md.
Floating-point checks here are implementation evidence, not interval proofs.
All output remains in this task's algorithm working directory.
"""
from pathlib import Path
import json
import time
import numpy as np
import cvxpy as cp
from scipy.special import logsumexp

OUT=Path(__file__).resolve().parent
P=2; K=3; Q=2; BETA=.95; TAU=.01
GAMMA=.05; KAPPA=25.; RHO=.12; RHO_S=.1; CURVATURE_FLOOR=1e-8
PAIRS=[(j,h) for j in range(K) for h in range(K) if h!=j]
I=np.eye(P)
v0=np.array([[0.,0.],[2.3,0.],[0.,20.]])
state={'B':I.copy(),'v':v0.copy(),'alpha':np.zeros(K),
       'r':np.array([.6,.6,.06]),'nu':np.full(K,1.001),
       'M':np.array([np.diag([1.,0.]),np.diag([1.,0.]),np.diag([0.,1.])]),
       'c':np.zeros((K,P)),'a':np.zeros(len(PAIRS)),
       'lam':np.zeros((len(PAIRS),P,P))}
INDEX={pair:k for k,pair in enumerate(PAIRS)}
for j in range(K):
    ids=[INDEX[j,h] for h in range(K) if h!=j]
    deltas=np.array([v0[h]-v0[j] for h in range(K) if h!=j])
    ds=np.linalg.norm(deltas,axis=1)
    caps=np.exp(-.5*ds**2+state['r'][j]*ds)
    centre=(caps[:,None]*deltas).sum(axis=0)/(1+caps.sum())
    state['c'][j]=centre
    for k,delta,cap in zip(ids,deltas,caps):
        state['a'][k]=cap
        scalar=np.sum((delta-centre)**2)-BETA+1e-4
        state['lam'][k]=scalar*I

offsets=np.array([[-.20,.10],[.15,-.20],[.10,.25],[-.05,-.10]])
centres=v0+np.array([[.02,0.],[.04,.05],[.12,-.02]])
X=np.vstack([centre+offsets for centre in centres]); N=len(X)

def ll(z):
    err=X@z['B'].T
    err=err[:,None,:]-z['v'][None,:,:]
    logk=z['alpha'][None,:]-.5*np.sum(err**2,axis=2)
    return float(np.mean(logsumexp(logk,axis=1))-logsumexp(z['alpha'])
                 +np.log(np.diag(z['B'])).sum()-P*.5*np.log(2*np.pi))

def responsibilities(z):
    err=X@z['B'].T
    err=err[:,None,:]-z['v'][None,:,:]
    logk=z['alpha'][None,:]-.5*np.sum(err**2,axis=2)
    return np.exp(logk-logsumexp(logk,axis=1)[:,None])

def auxiliary(z,R):
    err=X@z['B'].T
    err=err[:,None,:]-z['v'][None,:,:]
    return float(np.log(np.diag(z['B'])).sum()-.5*np.sum(R*np.sum(err**2,axis=2))/N
                 +(R.sum(axis=0)/N)@z['alpha']-logsumexp(z['alpha']))

def curvature_bound(d,s):
    # Both displacement and log-odds trust constraints are in the subproblem.
    lower=d-RHO; upper=d+RHO
    peak=float(np.clip(np.sqrt(Q),lower,upper))
    Amax=np.exp(s+RHO_S-.5*peak**2+np.sqrt(Q)*peak)
    G=np.sqrt(1+(upper+np.sqrt(Q))**2)
    H=1+np.sqrt(Q)/lower
    B1=np.sqrt(4*upper**2+upper**4)
    B2=2+2*upper
    raw=Amax*(upper**2*(G**2+H)+2*G*B1+B2)
    return max(CURVATURE_FLOOR,float(raw))

def exact_checks(z):
    vals={name:[] for name in ['separation','odds','first_matrix','second_matrix',
                               'lambda_psd','M_psd','I_minus_M_psd','budget','nu_lower']}
    for k,(j,h) in enumerate(PAIRS):
        delta=z['v'][h]-z['v'][j]; d=np.linalg.norm(delta); s=z['alpha'][h]-z['alpha'][j]
        A=np.exp(s-.5*d*d+z['r'][j]*d)
        vals['separation'].append(d*d-(1+TAU+z['nu'][j]+z['r'][j]**2+z['nu'][j]/z['r'][j]**2))
        vals['odds'].append(z['a'][k]-A)
        dc=delta-z['c'][j]
        vals['second_matrix'].append(np.linalg.eigvalsh(BETA*I+z['lam'][k]-np.outer(dc,dc)).min())
        vals['lambda_psd'].append(np.linalg.eigvalsh(z['lam'][k]).min())
    for j in range(K):
        ids=[INDEX[j,h] for h in range(K) if h!=j]
        vals['first_matrix'].append(np.linalg.eigvalsh(BETA*I-np.outer(z['c'][j],z['c'][j])-
                                        sum(z['a'][k]*z['lam'][k] for k in ids)).min())
        vals['M_psd'].append(np.linalg.eigvalsh(z['M'][j]).min())
        vals['I_minus_M_psd'].append(np.linalg.eigvalsh(I-z['M'][j]).min())
        residual=0.
        for h in range(K):
            if h==j: continue
            delta=z['v'][h]-z['v'][j]; d=np.linalg.norm(delta); s=z['alpha'][h]-z['alpha'][j]
            residual+=np.exp(s-.5*d*d+np.sqrt(Q)*d)*(delta@(I-z['M'][j])@delta)
        vals['budget'].append(z['nu'][j]-np.trace(z['M'][j])-residual/BETA)
        vals['nu_lower'].append(z['nu'][j]-1.)
    return {name:float(np.min(values)) for name,values in vals.items()}

B=cp.Variable((P,P)); v=cp.Variable((K,P)); al=cp.Variable(K)
r=cp.Variable(K,pos=True); nu=cp.Variable(K); c=cp.Variable((K,P)); a=cp.Variable(len(PAIRS))
lam=[cp.Variable((P,P),symmetric=True) for _ in PAIRS]
Ms=[cp.Variable((P,P),symmetric=True) for _ in range(K)]
constraints=[B[1,0]==0,cp.diag(B)>=.85,cp.diag(B)<=1.15,B[0,1]>=-.15,B[0,1]<=.15,
             v>=-3.,v<=23.,al>=-.25,al<=.25,al[-1]==0,r>=.02,r<=1.,nu>=1.,nu<=1e4,a>=0]
budget_terms=[[] for _ in range(K)]
curvatures=[]
for j in range(K):
    constraints.extend([Ms[j]>>0,I-Ms[j]>>0])
for k,(j,h) in enumerate(PAIRS):
    dt=state['v'][h]-state['v'][j]; d0=np.linalg.norm(dt); st=state['alpha'][h]-state['alpha'][j]
    delta=v[h]-v[j]; s=al[h]-al[j]
    constraints.extend([cp.norm(delta-dt)<=RHO,cp.abs(s-st)<=RHO_S,lam[k]>>0])
    rt=state['r'][j]; nt=state['nu'][j]; tangent=d0-rt
    low=tangent*((dt/d0)@delta-r[j])-.5*tangent**2-.5*cp.square(r[j])
    constraints.append(cp.exp(s-low)<=a[k])
    separation=1+TAU+nu[j]+cp.square(r[j])+cp.square(nu[j])/(2*nt*rt**2)+.5*nt*rt**2*cp.power(r[j],-4)
    constraints.append(separation<=2*dt@delta-d0*d0)
    dc=cp.reshape(delta-c[j],(P,1),order='F')
    constraints.append(cp.bmat([[BETA*I+lam[k],dc],[dc.T,np.ones((1,1))]])>>0)
    Pm=I-state['M'][j]; b0=dt@Pm@dt; A0=np.exp(st-.5*d0*d0+np.sqrt(Q)*d0)
    g_delta=A0*(2*Pm@dt+b0*(-dt+np.sqrt(Q)*dt/d0))
    g_s=A0*b0; g_M=-A0*np.outer(dt,dt)
    L=curvature_bound(d0,st); curvatures.append(L)
    delta_M=Ms[j]-state['M'][j]
    term=A0*b0+g_delta@(delta-dt)+g_s*(s-st)+cp.sum(cp.multiply(g_M,delta_M))
    term+=.5*L*(cp.sum_squares(delta-dt)+cp.square(s-st)+cp.sum_squares(delta_M))
    budget_terms[j].append(term)
for j in range(K):
    constraints.append(cp.trace(Ms[j])+sum(budget_terms[j])/BETA<=nu[j])
    linear=cp.Constant(np.zeros((P,P))); Ws=[]
    for h in range(K):
        if h==j: continue
        k=INDEX[j,h]; at=state['a'][k]; Lt=state['lam'][k]
        linear+=at*lam[k]+a[k]*Lt-at*Lt
        Ws.append(.5*(np.sqrt(KAPPA)*(a[k]-at)*I+(lam[k]-Lt)/np.sqrt(KAPPA)))
    cc=cp.reshape(c[j],(P,1),order='F')
    blocks=[[BETA*I-linear,cc]+Ws,[cc.T,np.ones((1,1))]+[np.zeros((1,P)) for _ in Ws]]
    for i,W in enumerate(Ws):
        blocks.append([W,np.zeros((P,1))]+[I if i==h else np.zeros((P,P)) for h in range(len(Ws))])
    constraints.append(cp.bmat(blocks)>>0)
R=responsibilities(state); BX=X@B.T
weighted_error=sum(cp.sum_squares(cp.multiply(np.sqrt(R[:,j])[:,None],BX-v[j])) for j in range(K))
Qt=cp.sum(cp.log(cp.diag(B)))-.5*weighted_error/N+(R.sum(axis=0)/N)@al-cp.log_sum_exp(al)
variables={'B':B,'v':v,'alpha':al,'r':r,'nu':nu,'c':c,'a':a}
prox=sum(cp.sum_squares(variables[key]-state[key]) for key in variables)
prox+=sum(cp.sum_squares(lam[k]-state['lam'][k]) for k in range(len(PAIRS)))
prox+=sum(cp.sum_squares(Ms[j]-state['M'][j]) for j in range(K))
problem=cp.Problem(cp.Maximize(Qt-.5*GAMMA*prox),constraints)
assert problem.is_dcp()
start=time.perf_counter()
problem.solve(solver=cp.CLARABEL,tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=250)
elapsed=time.perf_counter()-start
if problem.status not in ['optimal','optimal_inaccurate']:
    raise RuntimeError(problem.status)
raw={key:np.array(variable.value) for key,variable in variables.items()}
raw['lam']=np.array([variable.value for variable in lam]); raw['M']=np.array([variable.value for variable in Ms])
fraction=.999
new={key:state[key]+fraction*(raw[key]-state[key]) for key in state}
for key,variable in variables.items(): variable.value=new[key]
for k,variable in enumerate(lam): variable.value=new['lam'][k]
for j,variable in enumerate(Ms): variable.value=new['M'][j]
inner_violation=float(max(np.max(con.violation()) for con in constraints))
checks=exact_checks(new)
qgain=auxiliary(new,R)-auxiliary(state,R)
penalty=.5*GAMMA*sum(np.sum((new[key]-state[key])**2) for key in state)
lgain=ll(new)-ll(state)
result={'scope':'One bounded conic step; M and nu are joint free variables; no convergence claim from this computation.',
        'status':problem.status,'is_dcp':problem.is_dcp(),'observations':N,'elapsed_seconds':elapsed,
        'proximal_coefficient':GAMMA,'displacement_trust_radius':RHO,'log_odds_trust_radius':RHO_S,
        'pair_curvature_bounds':curvatures,'retained_fraction':fraction,
        'initial_exact_checks':exact_checks(state),'accepted_exact_checks':checks,
        'max_inner_constraint_violation':inner_violation,'likelihood_increase_per_observation':lgain,
        'auxiliary_increase_per_observation':qgain,'proximal_penalty':penalty,
        'change_norms':{key:float(np.linalg.norm(new[key]-state[key])) for key in state},
        'M':new['M'].tolist(),'nu':new['nu'].tolist(),'B':new['B'].tolist(),
        'shared_covariance':(np.linalg.inv(new['B'])@np.linalg.inv(new['B']).T).tolist(),
        'checks_are_floating_point':True}
assert min(checks.values())>=-1e-7,checks
assert inner_violation<=1e-7,inner_violation
assert qgain>=penalty-1e-7 and lgain>=qgain-1e-7,(qgain,penalty,lgain)
(OUT/'joint_budget_step_check.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
