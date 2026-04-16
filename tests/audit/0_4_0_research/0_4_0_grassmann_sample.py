"""
Quick Grassmannian optimisation on 5 matched molecules, 20 workers.
Tests: canonical, adapted, Grassmann(vf), Grassmann(v_min), Grassmann(combined)
at m = 2, 3, 5.

Usage: python 0_4_0_grassmann_sample.py
"""
import sys, time, numpy as np
from numpy.linalg import inv, eigh, svd
from scipy.linalg import null_space, expm
from scipy.optimize import minimize
from pathlib import Path
from multiprocessing import Pool

sys.path.insert(0, str(Path("C:/observer_geometry_workspace_v0.3.2/nomogeo/src")))
DATA_ROOT = Path("C:/observer_geometry_workspace_v0.3.2/datasets/hessian_qm9_DatasetDict/hessian_qm9_DatasetDict")
ATOMIC_MASS = {1:1.00782503223, 6:12.0, 7:14.00307400443, 8:15.99491461957, 9:18.99840316273}

def sym(M): return 0.5*(M+M.T)

def build_hvib(row):
    an = np.asarray(row["atomic_numbers"],dtype=int); pos = np.asarray(row["positions"],dtype=float)
    na = len(an); masses = np.array([ATOMIC_MASS[int(z)] for z in an])
    h_raw = np.asarray(row["hessian"],dtype=float)
    if h_raw.shape==(na,3,na,3): hc=sym(h_raw.reshape(3*na,3*na))
    else: hc=sym(h_raw.transpose(0,2,1,3).reshape(3*na,3*na))
    dm=np.repeat(masses,3); ism=1.0/np.sqrt(dm); hmw=sym((ism[:,None]*hc)*ism[None,:])
    w=np.sqrt(masses); ctr=np.average(pos,axis=0,weights=masses); ctrd=pos-ctr
    cols=[]
    for ax in range(3): v=np.zeros((na,3)); v[:,ax]=w; cols.append(v.reshape(-1))
    for ax in np.eye(3): v=np.cross(ax[None,:],ctrd)*w[:,None]; cols.append(v.reshape(-1))
    raw=np.column_stack(cols); u,sv,_=np.linalg.svd(raw,full_matrices=False)
    rk=int(np.sum(sv>1e-10*max(1,float(np.max(sv))))); vb=null_space(u[:,:rk].T)
    return sym(vb.T @ hmw @ vb)

def compute_diag(H, C, Hdot):
    n=H.shape[0]; m=C.shape[0]; Hinv=inv(H); Phi=inv(C@Hinv@C.T)
    L=Hinv@C.T@Phi; _,_,Vt=svd(C); Z=Vt[m:].T; R=Z.T@H@Z; Rinv=inv(R)
    V=L.T@Hdot@L; U_h=Z.T@Hdot@Z
    dHinv=-Hinv@Hdot@Hinv; dPhi=sym(-Phi@(C@dHinv@C.T)@Phi)
    vis=float(np.trace(inv(Phi)@dPhi)); hid=float(np.trace(Rinv@U_h))
    amb=float(np.trace(Hinv@Hdot)); vm=float(np.min(eigh(V)[0]))
    P=L@C; leak=float(np.linalg.norm(Hdot@P-P@Hdot,'fro'))
    vf=vis/amb if abs(amb)>1e-10 else float('nan')
    return {'vf':vf,'v_min':vm,'leak':leak,'vis':vis,'hid':hid,'amb':amb}

def grassmann_opt(H,Hdot,m,objective,n_restarts=15,max_iter=150):
    n=H.shape[0]; Hinv=inv(H); amb=float(np.trace(Hinv@Hdot))
    if abs(amb)<1e-10: return None,None
    np_=n*(n-1)//2; rng=np.random.default_rng(42)
    def neg_obj(Af,Q0):
        A=np.zeros((n,n)); idx=0
        for i in range(n):
            for j in range(i+1,n): A[i,j]=Af[idx]; A[j,i]=-Af[idx]; idx+=1
        Q=expm(A)@Q0; C=Q[:m,:]
        try:
            d=compute_diag(H,C,Hdot)
            if objective=='vf': return -d['vf'] if not np.isnan(d['vf']) else 1e10
            elif objective=='v_min': return -d['v_min']
            elif objective=='combined': return -(d['vf']+0.3*max(d['v_min'],0)-0.1*d['leak'])
        except: return 1e10
    best=1e10; best_C=None
    for _ in range(n_restarts):
        Q0=np.linalg.qr(rng.standard_normal((n,n)))[0]; A0=rng.standard_normal(np_)*0.05
        try:
            res=minimize(neg_obj,A0,args=(Q0,),method='L-BFGS-B',options={'maxiter':max_iter,'ftol':1e-12})
            if res.fun<best:
                best=res.fun; A=np.zeros((n,n)); idx=0
                for i in range(n):
                    for j in range(i+1,n): A[i,j]=res.x[idx]; A[j,i]=-res.x[idx]; idx+=1
                best_C=(expm(A)@Q0)[:m,:].copy()
        except: pass
    if best_C is None: return None,None
    return best_C,compute_diag(H,best_C,Hdot)

def process_one(args):
    label, hv, hw, n = args
    from nomogeo import closure_adapted_observer
    Hdot=hw-hv; H=hv; Hinv=inv(H); amb=float(np.trace(Hinv@Hdot))
    if abs(amb)<0.01: return None
    rows = []
    for m in [2,3,5]:
        if m>=n: continue
        C_c=np.zeros((m,n))
        for i in range(m): C_c[i,i]=1.0
        d_c=compute_diag(H,C_c,Hdot)
        try:
            res=closure_adapted_observer(H,[Hdot],m); d_a=compute_diag(H,res.C,Hdot)
        except: d_a=None
        _,d_gvf=grassmann_opt(H,Hdot,m,'vf')
        _,d_gvm=grassmann_opt(H,Hdot,m,'v_min')
        _,d_gcomb=grassmann_opt(H,Hdot,m,'combined')
        rows.append({'label':label,'n':n,'m':m,'canon':d_c,'adapted':d_a,
                     'gr_vf':d_gvf,'gr_vm':d_gvm,'gr_comb':d_gcomb})
    return rows

def main():
    import pyarrow.ipc as ipc
    t0=time.time()
    vac={}; wat={}
    for shard in range(5):
        p=DATA_ROOT/'vacuum'/f'data-{shard:05d}-of-00005.arrow'
        if not p.exists(): continue
        seen=0
        with ipc.open_stream(p) as r:
            for b in r:
                for row in b.to_pylist():
                    try:
                        hv=build_hvib(row)
                        if np.min(eigh(hv)[0])>1e-6: vac[row['label']]=(hv,hv.shape[0])
                    except: pass
                    seen+=1
                    if seen>=200: break
                if seen>=200: break
    for shard in range(5):
        p=DATA_ROOT/'water'/f'data-{shard:05d}-of-00005.arrow'
        if not p.exists(): continue
        seen=0
        with ipc.open_stream(p) as r:
            for b in r:
                for row in b.to_pylist():
                    try:
                        hw=build_hvib(row)
                        if np.min(eigh(hw)[0])>1e-6: wat[row['label']]=(hw,hw.shape[0])
                    except: pass
                    seen+=1
                    if seen>=200: break
                if seen>=200: break
    matched=set(vac.keys())&set(wat.keys())
    tasks=[]
    for l in sorted(matched):
        hv,nv=vac[l]; hw,nw=wat[l]
        if nv==nw and nv>=6: tasks.append((l,hv,hw,nv))
        if len(tasks)>=5: break
    print(f'Processing {len(tasks)} molecules on {20} workers...')
    with Pool(20) as pool:
        all_rows=[r for r in pool.imap_unordered(process_one,tasks) if r]
    for mol_rows in all_rows:
        if not mol_rows: continue
        label=mol_rows[0]['label']; n=mol_rows[0]['n']
        print(f'\n=== {label} (n_vib={n}) ===')
        for row in mol_rows:
            m=row['m']
            print(f'  m={m}:')
            print(f'    {"Observer":>20s}  {"vf":>8s}  {"v_min":>8s}  {"leak":>8s}  {"V>0":>4s}')
            def pr(nm,d):
                if d is None: print(f'    {nm:>20s}  {"failed":>8s}'); return
                vp='Y' if d['v_min']>1e-6 else 'N'
                print(f'    {nm:>20s}  {d["vf"]:8.4f}  {d["v_min"]:8.4f}  {d["leak"]:8.4f}  {vp:>4s}')
            pr('Canonical',row['canon']); pr('Adapted',row['adapted'])
            pr('Grassmann(vf)',row['gr_vf']); pr('Grassmann(v_min)',row['gr_vm'])
            pr('Grassmann(combined)',row['gr_comb'])
    print(f'\nElapsed: {time.time()-t0:.1f}s')

if __name__=='__main__': main()
