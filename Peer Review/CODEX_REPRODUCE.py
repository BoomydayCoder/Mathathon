from pathlib import Path
import concurrent.futures, subprocess, time, json, re
root=Path(__file__).resolve().parent
jobs=[]
for readme in sorted(root.glob('Codex-*/README.md')):
    for f in re.findall(r'python3 (computations/\S+\.py)',readme.read_text()): jobs.append((readme.parent,f))
def run(job):
    d,f=job; start=time.monotonic()
    try:
        r=subprocess.run(['python3',f],cwd=d,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
        result={'script':f,'exit_code':r.returncode,'seconds':round(time.monotonic()-start,2),'output':r.stdout}
    except subprocess.TimeoutExpired as e:
        result={'script':f,'exit_code':None,'seconds':180,'output':'TIMEOUT; this is not a passed check.\n'+str(e.stdout or '')}
    out=d/'verification'; out.mkdir(exist_ok=True)
    (out/(Path(f).stem+'.txt')).write_text(result['output'])
    return d.name,result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    results=list(pool.map(run,jobs))
for folder in sorted({d for d,r in results}):
    subset=[{k:v for k,v in r.items() if k!='output'} for d,r in results if d==folder]
    (root/folder/'verification'/'results.json').write_text(json.dumps({'checked_at_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'checks':subset},indent=2)+'\n')
for d,r in results: print(d,Path(r['script']).name,r['exit_code'],r['seconds'])
assert all(r['exit_code']==0 for d,r in results),'Some checks did not pass'
