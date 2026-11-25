// Independent verifier enumerator (written for R8 verification).
// Counts/stores all nonzero integer x (one per +-pair) with sum_i r[i]*(x_i + sum_{j>i} mu[j][i] x_j)^2 <= B.
#include <math.h>
#include <string.h>
typedef struct {int n; const double*mu; const double*r; double B; long cap; long cnt; int*out; long nodes; long maxnodes; int st; int x[64];} E;
static void rec(E*e,int i,double part,int allzero){
  int n=e->n;
  double c=0.0; for(int j=i+1;j<n;j++) c+= e->mu[j*n+i]*(double)e->x[j];
  double rem=e->B-part; if(rem<0) return;
  double rad=sqrt(rem/e->r[i]);
  long lo=(long)floor(-c-rad), hi=(long)ceil(-c+rad);
  if(allzero && lo<0) lo=0;
  for(long v=lo; v<=hi; v++){
    double y=(double)v+c; double p=part+e->r[i]*y*y;
    if(p>e->B) continue;
    if(++e->nodes>e->maxnodes){e->st=1;return;}
    e->x[i]=(int)v;
    int az=allzero&&(v==0);
    if(i==0){ if(az) continue; if(e->cnt<e->cap) memcpy(e->out+e->cnt*n,e->x,n*sizeof(int)); e->cnt++; }
    else { rec(e,i-1,p,az); if(e->st) return; }
  }
  e->x[i]=0;
}
long my_enum(int n,const double*mu,const double*r,double B,long cap,int*out,long maxnodes,long*nodes,int*status){
  E e; e.n=n; e.mu=mu; e.r=r; e.B=B; e.cap=cap; e.cnt=0; e.out=out; e.nodes=0; e.maxnodes=maxnodes; e.st=0; memset(e.x,0,sizeof(e.x));
  rec(&e,n-1,0.0,1);
  *nodes=e.nodes; *status=e.st; return e.cnt;
}
