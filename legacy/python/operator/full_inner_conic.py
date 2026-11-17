"""Bounded nonsymmetric implementation check, three convex-inner EM/MM steps.

Run with the task-local CVXPY/Clarabel dependency directory on PYTHONPATH.
This checks feasibility, the exact matrix certificate and likelihood ascent.
It does not diagnose convergence to a global or local optimum.
"""
from pathlib import Path
import json
import time
import numpy as np
import cvxpy as cp
import clarabel
from scipy.special import logsumexp

OUT=Path(__file__).resolve().parent
P=2; K=3; Q=2; NU=2.; BETA=.99; TAU=.01
PROX=.15; SCALE=25.; STEPS=3
PAIRS=[(j,h) for j in range(K) for h in range(K) if h!=j]
INDEX={pair:i for i,pair in enumerate(PAIRS)}
I=np.eye(P)
SEED=2026111701
TOL={'tol_gap_abs':1e-9,'tol_gap_rel':1e-9,'tol_feas':1e-9,'max_iter':250}

# Equilateral start with the independently interval-checked rational witness.
angles=np.arange(K)*2*np.pi/K
v0=2.8/np.sqrt(3)*np.c_[np.cos(angles),np.sin(angles)]
state={'B':I.copy(),'v':v0.copy(),'alpha':np.zeros(K),'r':np.full(K,.68),
       'c':np.zeros((K,P)),'a':np.full(len(PAIRS),np.exp(-2.016)),
       'lam':np.zeros((len(PAIRS),P,P))}
for j in range(K):
 e1=-v0[j]/np.linalg.norm(v0[j]);e2=np.array([-e1[1],e1[0]])
 frame=np.column_stack([e1,e2]);state['c'][j]=.51*e1
 for h in range(K):
  if h==j:continue
  sign=np.sign((v0[h]-v0[j])@e2)
  state['lam'][INDEX[j,h]]=frame@np.array([[2.72,sign*2.70],[sign*2.70,2.70]])@frame.T

# Nonsymmetric observations: different component counts and perturbed means,
# shared nonspherical covariance. Labels are used only to generate the data.
rng=np.random.default_rng(SEED)
mu_data=1.08*v0+np.array([[.18,-.10],[-.13,.21],[.06,-.16]])
C_data=np.array([[.65,.19],[.19,.88]])
chol=np.linalg.cholesky(C_data)
counts=[42,33,21]
X=np.concatenate([mu_data[j]+rng.normal(size=(counts[j],P))@chol.T for j in range(K)])
rng.shuffle(X)
np.savetxt(OUT/'full_inner_data.csv',X,delimiter=',',header='x1,x2',comments='')
N=len(X)

def rho(D,nu=NU):
 z=D*D-nu-1
 return np.sqrt(2*nu/(z+np.sqrt(z*z-4*nu)))

def logkern(z):
 err=X@z['B'].T
 err=err[:,None,:]-z['v'][None,:,:]
 return z['alpha'][None,:]-.5*np.sum(err*err,axis=2)

def ll(z):
 return float(np.mean(logsumexp(logkern(z),axis=1)-logsumexp(z['alpha']))
              +np.linalg.slogdet(z['B'])[1]-P*.5*np.log(2*np.pi))

def resp(z):
 L=logkern(z)
 return np.exp(L-logsumexp(L,axis=1)[:,None])

def aux(z,R):
 err=X@z['B'].T
 err=err[:,None,:]-z['v'][None,:,:]
 return float(np.linalg.slogdet(z['B'])[1]-.5*np.sum(R*np.sum(err*err,axis=2))/N
              +(R.sum(axis=0)/N)@z['alpha']-logsumexp(z['alpha']))

def norm2(z,zt):
 return float(sum(np.sum((z[key]-zt[key])**2) for key in z))

def exact_checks(z):
 vals={'separation':[],'odds_epigraph':[],'lambda_psd':[],
       'first_lifted_matrix':[],'first_actual_matrix':[],'second_matrix':[],
       'automatic_invariance':[],'dominance_bound':[]}
 Ds=[];Ss=[];Pi=np.exp(z['alpha']-logsumexp(z['alpha']))
 A=np.zeros(len(PAIRS))
 for k,(j,h) in enumerate(PAIRS):
  delta=z['v'][h]-z['v'][j];d=np.linalg.norm(delta)
  A[k]=np.exp(z['alpha'][h]-z['alpha'][j]-.5*d*d+z['r'][j]*d)
  vals['separation'].append(d*d-(NU+1+TAU+z['r'][j]**2+NU/z['r'][j]**2))
  vals['odds_epigraph'].append(z['a'][k]-A[k])
  vals['lambda_psd'].append(np.linalg.eigvalsh(z['lam'][k]).min())
  dc=delta-z['c'][j]
  vals['second_matrix'].append(np.linalg.eigvalsh(BETA*I+z['lam'][k]-np.outer(dc,dc)).min())
 for j in range(K):
  ids=[INDEX[j,h] for h in range(K) if h!=j]
  cc=np.outer(z['c'][j],z['c'][j])
  vals['first_lifted_matrix'].append(np.linalg.eigvalsh(BETA*I-cc-sum(z['a'][k]*z['lam'][k] for k in ids)).min())
  vals['first_actual_matrix'].append(np.linalg.eigvalsh(BETA*I-cc-sum(A[k]*z['lam'][k] for k in ids)).min())
  D=min(np.linalg.norm(z['v'][h]-z['v'][j]) for h in range(K) if h!=j)
  sroot=rho(D/np.sqrt(BETA))**2
  vals['automatic_invariance'].append(z['r'][j]-np.sqrt(BETA*sroot))
  S=float(sum(A[k] for k in ids));Ss.append(S);Ds.append(D)
  vals['dominance_bound'].append(sroot-S)
 out={name:float(np.min(v)) for name,v in vals.items()}
 out.update({'min_weight':float(Pi.min()),'weights':Pi.tolist(),'radii':z['r'].tolist(),
             'minimum_local_distances':Ds,'S':Ss,'minimum_B_diagonal':float(np.diag(z['B']).min())})
 return out

def inner_checks(z,zt):
 vals={'separation':[],'odds_epigraph':[],'first_matrix':[]}
 for k,(j,h) in enumerate(PAIRS):
  dt=zt['v'][h]-zt['v'][j];d0=np.linalg.norm(dt);u0=dt/d0;s0=d0-zt['r'][j]
  delta=z['v'][h]-z['v'][j];rad=z['r'][j]
  low=s0*(u0@delta-rad)-.5*s0*s0-.5*rad*rad
  vals['odds_epigraph'].append(z['a'][k]-np.exp(z['alpha'][h]-z['alpha'][j]-low))
  vals['separation'].append(u0@delta-np.sqrt(NU+1+TAU+rad*rad+NU/(rad*rad)))
 for j in range(K):
  U=np.zeros((P,P))
  for h in range(K):
   if h==j:continue
   k=INDEX[j,h];a=z['a'][k];at=zt['a'][k];L=z['lam'][k];Lt=zt['lam'][k]
   W=np.sqrt(SCALE)*(a-at)*I+(L-Lt)/np.sqrt(SCALE)
   U+=at*L+a*Lt-at*Lt+.25*W@W
  vals['first_matrix'].append(np.linalg.eigvalsh(BETA*I-np.outer(z['c'][j],z['c'][j])-U).min())
 return {name:float(np.min(v)) for name,v in vals.items()}

def solve_step(zt):
 R=resp(zt);nj=R.sum(axis=0)
 B=cp.Variable((P,P),name='B');v=cp.Variable((K,P),name='v')
 al=cp.Variable(K,name='alpha');r=cp.Variable(K,nonneg=True,name='r')
 c=cp.Variable((K,P),name='c');a=cp.Variable(len(PAIRS),nonneg=True,name='a')
 lam=[cp.Variable((P,P),symmetric=True,name=f'Lambda_{j}_{h}') for j,h in PAIRS]
 constraints=[B[1,0]==0,cp.diag(B)>=.3,cp.diag(B)<=3.,B[0,1]>=-2.,B[0,1]<=2.,
              v>=-8.,v<=8.,al>=-3.,al<=3.,al[-1]==0,r>=.01,r<=1.]
 for k,(j,h) in enumerate(PAIRS):
  delta=v[h,:]-v[j,:]
  dt=zt['v'][h]-zt['v'][j];d0=np.linalg.norm(dt);u0=dt/d0;s0=d0-zt['r'][j]
  H=cp.norm(cp.hstack([np.sqrt(NU+1+TAU),r[j],np.sqrt(NU)*cp.inv_pos(r[j])]),2)
  low=s0*(u0@delta-r[j])-.5*s0*s0-.5*cp.square(r[j])
  constraints.extend([H<=u0@delta,cp.exp(al[h]-al[j]-low)<=a[k],lam[k]>>0])
  dc=cp.reshape(delta-c[j,:],(P,1),order='F')
  constraints.append(cp.bmat([[BETA*I+lam[k],dc],[dc.T,np.ones((1,1))]])>>0)
 for j in range(K):
  linear=cp.Constant(np.zeros((P,P)));Ws=[]
  for h in range(K):
   if h==j:continue
   k=INDEX[j,h];at=zt['a'][k];Lt=zt['lam'][k]
   linear=linear+at*lam[k]+a[k]*Lt-at*Lt
   Ws.append(.5*(np.sqrt(SCALE)*(a[k]-at)*I+(lam[k]-Lt)/np.sqrt(SCALE)))
  cc=cp.reshape(c[j,:],(P,1),order='F')
  blocks=[[BETA*I-linear,cc]+Ws,
          [cc.T,np.ones((1,1))]+[np.zeros((1,P)) for _ in Ws]]
  for i,W in enumerate(Ws):
   blocks.append([W,np.zeros((P,1))]+[I if i==h else np.zeros((P,P)) for h in range(len(Ws))])
  constraints.append(cp.bmat(blocks)>>0)
 BX=X@B.T
 weighted_error=sum(cp.sum_squares(cp.multiply(np.sqrt(R[:,j])[:,None],BX-v[j,:])) for j in range(K))
 Qt=cp.sum(cp.log(cp.diag(B)))-.5*weighted_error/N+(nj/N)@al-cp.log_sum_exp(al)
 terms=[cp.sum_squares(B-zt['B']),cp.sum_squares(v-zt['v']),cp.sum_squares(al-zt['alpha']),
        cp.sum_squares(r-zt['r']),cp.sum_squares(c-zt['c']),cp.sum_squares(a-zt['a'])]
 terms.extend(cp.sum_squares(lam[k]-zt['lam'][k]) for k in range(len(PAIRS)))
 prob=cp.Problem(cp.Maximize(Qt-.5*PROX*sum(terms)),constraints)
 assert prob.is_dcp()
 begin=time.perf_counter()
 data,_,_=prob.get_problem_data(cp.CLARABEL)
 canonicalisation_seconds=time.perf_counter()-begin
 dims=data['dims']
 size={'independent_model_coordinates':44,'cvxpy_scalar_variables':int(prob.size_metrics.num_scalar_variables),
       'cvxpy_scalar_equalities':int(prob.size_metrics.num_scalar_eq_constr),
       'cvxpy_scalar_inequalities':int(prob.size_metrics.num_scalar_leq_constr),
       'canonical_variables':int(data['A'].shape[1]),'canonical_rows':int(data['A'].shape[0]),
       'canonical_A_nonzeros':int(data['A'].nnz),'zero_cone_dim':int(dims.zero),'nonnegative_cone_dim':int(dims.nonneg),
       'exp_cones':int(dims.exp),'soc_orders':[int(x) for x in dims.soc],
       'psd_orders':[int(x) for x in dims.psd]}
 begin=time.perf_counter()
 value=prob.solve(solver=cp.CLARABEL,verbose=False,**TOL)
 solve_wall_seconds=time.perf_counter()-begin
 if prob.status not in ['optimal','optimal_inaccurate']:
  raise RuntimeError(f'Unexpected solver status: {prob.status}')
 raw={'B':np.array(B.value),'v':np.array(v.value),'alpha':np.array(al.value),'r':np.array(r.value),
      'c':np.array(c.value),'a':np.array(a.value),'lam':np.array([L.value for L in lam])}
 # A fixed small retreat protects the checked point from cone-solver boundary
 # tolerances. It is a convex combination inside the same inner problem.
 frac=.999
 new={key:zt[key]+frac*(raw[key]-zt[key]) for key in zt}
 orig=exact_checks(new);inner=inner_checks(new,zt)
 for _ in range(12):
  core=[orig[key] for key in ['separation','odds_epigraph','lambda_psd','first_lifted_matrix','first_actual_matrix','second_matrix']]+list(inner.values())
  if min(core)>=-1e-10:break
  frac*=.8
  new={key:zt[key]+frac*(raw[key]-zt[key]) for key in zt}
  orig=exact_checks(new);inner=inner_checks(new,zt)
 incQ=aux(new,R)-aux(zt,R);penalty=.5*PROX*norm2(new,zt);incLL=ll(new)-ll(zt)
 assert min(core)>=-1e-10,(orig,inner)
 assert incQ>=penalty-1e-9 and incLL>=incQ-1e-9 and incLL>0
 result={'solver_status':prob.status,'solver_iterations':int(prob.solver_stats.num_iters),
         'solver_reported_seconds':float(prob.solver_stats.solve_time),
         'canonicalisation_seconds':canonicalisation_seconds,'solve_wall_seconds':solve_wall_seconds,
         'dcp_verified':True,'size':size,'retained_step_fraction':frac,
         'log_likelihood_per_observation':ll(new),'log_likelihood_increase_per_observation':incLL,
         'EM_auxiliary_increase_per_observation':incQ,'proximal_penalty':penalty,
         'exact_checks':orig,'inner_checks':inner,
         'B':new['B'].tolist(),'whitened_means':new['v'].tolist(),'logits':new['alpha'].tolist()}
 return new,result

initial_checks=exact_checks(state)
initial_ll=ll(state)
results=[]
for step in range(STEPS):
 state,res=solve_step(state);res['step']=step+1;results.append(res)
 print(json.dumps({'step':step+1,'status':res['solver_status'],
  'likelihood_increase':res['log_likelihood_increase_per_observation'],
  'Q_increase':res['EM_auxiliary_increase_per_observation'],'penalty':res['proximal_penalty'],
  'seconds':res['solve_wall_seconds'],'min_original_matrix_slack':min(res['exact_checks'][k] for k in ['lambda_psd','first_lifted_matrix','first_actual_matrix','second_matrix']),
  'min_separation_slack':res['exact_checks']['separation']}),flush=True)

# Fixed-density SDP check of the analytic optimum, not an additional fitting step.
D=2.8;rfix=.68;Afix=np.exp(-D*D/2+rfix*D)
deltas=np.array([[np.sqrt(3)*D/2,D/2],[np.sqrt(3)*D/2,-D/2]])
cc=cp.Variable(P);lev=cp.Variable();Ls=[cp.Variable((P,P),symmetric=True) for _ in range(2)]
col=cp.reshape(cc,(P,1),order='F')
cons=[lev>=0,cp.bmat([[lev*I-Afix*sum(Ls),col],[col.T,np.ones((1,1))]])>>0]
for h in range(2):
 dc=cp.reshape(deltas[h]-cc,(P,1),order='F')
 cons.extend([Ls[h]>>0,cp.bmat([[lev*I+Ls[h],dc],[dc.T,np.ones((1,1))]])>>0])
fp=cp.Problem(cp.Minimize(lev),cons)
fp.solve(solver=cp.CLARABEL,verbose=False,tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=250)
analytic=lambda rr:(2*np.exp(-D*D/2+D*rr))*(np.sqrt(3)*D/2)**2/(1+2*np.exp(-D*D/2+D*rr))**2
rzero=rho(D)
z=D*D-Q-1-TAU
rmin=np.sqrt(2*Q/(z+np.sqrt(z*z-4*Q)))
fixed={'solver_status':fp.status,'sdp_value':float(fp.value),'analytic_value':float(analytic(.68)),
       'absolute_difference':float(abs(fp.value-analytic(.68))),
       'radius_zero_tau':float(rzero),'value_zero_tau':float(analytic(rzero)),
       'radius_min_tau_01':float(rmin),'value_min_tau_01':float(analytic(rmin)),
       'radius_used_for_rational_witness':.68,'value_at_068':float(analytic(.68)),
       'theoretical_status':'The analytic value is the proved SDP optimum in the stated symmetric regime; solver agreement is a numerical check.'}
assert fixed['absolute_difference']<2e-8
Binv=np.linalg.inv(state['B']);Cfinal=Binv@Binv.T
payload={'scope':'Three nonsymmetric full-triangular-whitener convex-inner steps; no stationarity or global-argmax diagnosis.',
 'versions':{'cvxpy':cp.__version__,'clarabel':clarabel.__version__,'numpy':np.__version__},
 'settings':{'n':N,'p':P,'K':K,'q':Q,'nu':NU,'M':'identity','beta':BETA,'tau':TAU,'proximal_gamma':PROX,
             'matrix_majorant_scale':SCALE,'fitting_steps':STEPS,'seed':SEED,'solver_tolerances':TOL,
             'scalar_invariance_constraint_imposed':False},
 'initial_log_likelihood_per_observation':initial_ll,'initial_checks':initial_checks,'steps':results,
 'total_log_likelihood_increase_per_observation':ll(state)-initial_ll,
 'final_shared_covariance':Cfinal.tolist(),'final_original_coordinate_means':np.linalg.solve(state['B'],state['v'].T).T.tolist(),
 'fixed_density_sdp':fixed,'final_lift':{key:value.tolist() for key,value in state.items()}}
(OUT/'full_inner_results.json').write_text(json.dumps(payload,indent=2))
print(json.dumps({'fixed_density_sdp':fixed,'final_covariance':Cfinal.tolist(),
                  'total_log_likelihood_increase_per_observation':payload['total_log_likelihood_increase_per_observation'],
                  'solver_size':results[0]['size']},indent=2))
