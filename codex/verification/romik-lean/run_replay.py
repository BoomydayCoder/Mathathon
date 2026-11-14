"""Durably record a bounded Lean kernel replay, with line-buffered PTY output."""
from pathlib import Path
import argparse,datetime,json,os,pty,resource,select,signal,subprocess,sys,time
parser=argparse.ArgumentParser()
parser.add_argument('--max-seconds',type=int,default=2700)
args=parser.parse_args()
base=Path(__file__).resolve().parent
repo=base.parents[1]/'.cache/lean-romik/repo'
assert (repo/'lean-toolchain').read_text().strip()=='leanprover/lean4:v4.34.0'
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
assert commit=='fc2fefa77076e22e7e80b4c101ded5a98d53e14d'
status={'status':'starting','command':['lake','env','leanchecker','--verbose','Zeta5'],'commit':commit,
        'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'threads':os.getenv('LEAN_NUM_THREADS'),
        'max_seconds':args.max_seconds,'scope':'project modules replayed against imported environments, not --fresh replay of all Mathlib constants'}
status_path=base/'kernel-replay-result.json'
def save():
    tmp=status_path.with_suffix('.json.tmp');tmp.write_text(json.dumps(status,indent=2)+'\n');tmp.replace(status_path)
stop=[]
def interrupted(signum,frame):stop.append(signum)
signal.signal(signal.SIGTERM,interrupted)
signal.signal(signal.SIGINT,interrupted)
start=time.monotonic();save()
master,slave=pty.openpty()
with (base/'kernel-replay.log').open('wb',buffering=0) as log:
    proc=subprocess.Popen(status['command'],cwd=repo,stdout=slave,stderr=slave,start_new_session=True)
    os.close(slave)
    status.update(status='running',pid=proc.pid);save()
    banner=('Kernel replay started '+status['started_utc']+'; PID '+str(proc.pid)+'\n').encode()
    log.write(banner);sys.stdout.buffer.write(banner);sys.stdout.buffer.flush()
    last_save=start;ended=False
    while not ended:
        now=time.monotonic()
        if stop or now-start>args.max_seconds:
            status['termination_reason']='signal '+str(stop[0]) if stop else 'bounded timeout'
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=15)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
            break
        ready,_,_=select.select([master],[],[],2)
        if ready:
            try:data=os.read(master,65536)
            except OSError:data=b''
            if not data:ended=True
            else:
                log.write(data);sys.stdout.buffer.write(data);sys.stdout.buffer.flush()
        if now-last_save>15:
            status['elapsed_seconds']=now-start;status['last_checkpoint_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save();last_save=now
    proc.wait();os.close(master)
    status.update(status='passed' if proc.returncode==0 and 'termination_reason' not in status else 'failed_or_interrupted',
                  exit_code=proc.returncode,elapsed_seconds=time.monotonic()-start,
                  finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  max_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    save()
    result=(json.dumps(status,indent=2)+'\n').encode();log.write(result);sys.stdout.buffer.write(result);sys.stdout.buffer.flush()
sys.exit(0 if status['status']=='passed' else 1)
