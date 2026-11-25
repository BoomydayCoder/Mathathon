// LLL + BKZ (floating point, basis rows b in R^m, integer transformation U tracked exactly)
#include <math.h>
#include <stdlib.h>
#include <string.h>
static double *BS=0; static size_t BSsz=0;
static void gso_row(int n,int m,const double*b,double*mu,double*r,int k){
    // modified Gram-Schmidt with stored b* (numerically stable for nearly dependent bases)
    if(BSsz<(size_t)n*m){ BS=(double*)realloc(BS,sizeof(double)*(size_t)n*m); BSsz=(size_t)n*m; }
    const double*bk=b+(size_t)k*m; double*sk=BS+(size_t)k*m;
    for(int t=0;t<m;t++) sk[t]=bk[t];
    for(int j=0;j<k;j++){
        const double*sj=BS+(size_t)j*m; double s=0;
        for(int t=0;t<m;t++) s+=sk[t]*sj[t];
        double c=s/r[j]; mu[k*n+j]=c;
        for(int t=0;t<m;t++) sk[t]-=c*sj[t];
    }
    double s=0; for(int t=0;t<m;t++) s+=sk[t]*sk[t];
    r[k]=s;
}
int lll_c(int n,int m,double*b,long*U,double delta,double*mu,double*r){
    int k=1; long iters=0;
    gso_row(n,m,b,mu,r,0);
    while(k<n){
        if(++iters>200000000L) return -1;
        gso_row(n,m,b,mu,r,k);
        for(int rep=0;rep<6;rep++){
            int ch=0;
            for(int j=k-1;j>=0;j--){
                double q=rint(mu[k*n+j]);
                if(q!=0.0){
                    ch=1;
                    double*bk=b+(size_t)k*m,*bj=b+(size_t)j*m;
                    for(int t=0;t<m;t++) bk[t]-=q*bj[t];
                    long ql=(long)q; long*Uk=U+(size_t)k*n,*Uj=U+(size_t)j*n;
                    for(int t=0;t<n;t++) Uk[t]-=ql*Uj[t];
                    for(int i=0;i<j;i++) mu[k*n+i]-=q*mu[j*n+i];
                    mu[k*n+j]-=q;
                }
            }
            if(!ch) break;
            gso_row(n,m,b,mu,r,k);
        }
        double mm=mu[k*n+k-1];
        if(r[k]>=(delta-mm*mm)*r[k-1]) k++;
        else{
            double*bk=b+(size_t)k*m,*bj=b+(size_t)(k-1)*m;
            for(int t=0;t<m;t++){double tmp=bk[t];bk[t]=bj[t];bj[t]=tmp;}
            long*Uk=U+(size_t)k*n,*Uj=U+(size_t)(k-1)*n;
            for(int t=0;t<n;t++){long tmp=Uk[t];Uk[t]=Uj[t];Uj[t]=tmp;}
            if(k>1) k--; else gso_row(n,m,b,mu,r,0);
        }
    }
    return 0;
}
typedef struct{int n,jj,bs;const double*mu,*r;double R2;long x[200],best[200];long nodes,maxnodes;int found;}ectx;
static void erec(ectx*c,int t,double partial,int az){
    if(c->nodes>c->maxnodes) return;
    int n=c->n,jj=c->jj;
    double cen=0;
    for(int i=t+1;i<c->bs;i++) cen-=c->x[i]*c->mu[(jj+i)*n+(jj+t)];
    double rt=c->r[jj+t];
    double rem=c->R2-partial; if(rem<=0) return;
    double rad=sqrt(rem/rt);
    long lo=(long)ceil(cen-rad),hi=(long)floor(cen+rad);
    if(az&&lo<0) lo=0;
    long x0=(long)rint(cen); if(x0<lo) x0=lo; if(x0>hi) x0=hi;
    long a=x0,b=x0+1;
    while(a>=lo||b<=hi){
        long v;
        if(a>=lo&&(b>hi||fabs(a-cen)<=fabs(b-cen))) v=a--; else v=b++;
        double d=v-cen; double p=partial+d*d*rt;
        if(p>=c->R2) continue;
        c->nodes++;
        c->x[t]=v;
        int az2=az&&(v==0);
        if(t==0){ if(!az2){ c->R2=p*0.999999; memcpy(c->best,c->x,sizeof(long)*c->bs); c->found=1; } }
        else erec(c,t-1,p,az2);
        if(c->nodes>c->maxnodes) return;
    }
    c->x[t]=0;
}
static void xgcd(long a,long b,long*g,long*s,long*t){
    long s0=1,s1=0,t0=0,t1=1;
    while(b){long q=a/b;long tmp=a-q*b;a=b;b=tmp;tmp=s0-q*s1;s0=s1;s1=tmp;tmp=t0-q*t1;t0=t1;t1=tmp;}
    if(a<0){a=-a;s0=-s0;t0=-t0;}
    *g=a;*s=s0;*t=t0;
}
static void recompute(int n,int m,const double*b0,double*b,const long*U){
    for(int i=0;i<n;i++) for(int t=0;t<m;t++){
        double s=0; for(int k=0;k<n;k++) s+=(double)U[(size_t)i*n+k]*b0[(size_t)k*m+t];
        b[(size_t)i*m+t]=s;
    }
}
// BKZ: b0 = original basis (n x m), U in/out (n x n), returns number of tours done
int bkz_c(int n,int m,const double*b0,long*U,int beta,int maxtours,double delta,long maxnodes,double*r_out){
    double*mu=calloc((size_t)n*n,8),*r=calloc(n,8),*b=calloc((size_t)n*m,8);
    recompute(n,m,b0,b,U);
    lll_c(n,m,b,U,delta,mu,r);
    int tour;
    for(tour=0;tour<maxtours;tour++){
        int changed=0;
        for(int j=0;j<n-1;j++){
            int bs=beta<n-j?beta:n-j;
            for(int k=0;k<n;k++) gso_row(n,m,b,mu,r,k);
            ectx c; c.n=n;c.jj=j;c.bs=bs;c.mu=mu;c.r=r;c.R2=0.99*r[j];c.nodes=0;c.maxnodes=maxnodes;c.found=0;
            memset(c.x,0,sizeof c.x);
            erec(&c,bs-1,0.0,1);
            if(c.found){
                long*x=c.best;
                for(int i=bs-1;i>=1;i--){
                    if(x[i]==0) continue;
                    long a=x[i-1],bb=x[i],g,s,t; xgcd(a,bb,&g,&s,&t);
                    long ca=a/g,cb=bb/g;
                    double*bp=b+(size_t)(j+i-1)*m,*bq=b+(size_t)(j+i)*m;
                    for(int t2=0;t2<m;t2++){double u=bp[t2],v=bq[t2];bp[t2]=ca*u+cb*v;bq[t2]=-t*u+s*v;}
                    long*Up=U+(size_t)(j+i-1)*n,*Uq=U+(size_t)(j+i)*n;
                    for(int t2=0;t2<n;t2++){long u=Up[t2],v=Uq[t2];Up[t2]=ca*u+cb*v;Uq[t2]=-t*u+s*v;}
                    x[i-1]=g;x[i]=0;
                }
                lll_c(n,m,b,U,delta,mu,r);
                changed=1;
            }
        }
        recompute(n,m,b0,b,U);
        if(!changed) break;
    }
    for(int k=0;k<n;k++) gso_row(n,m,b,mu,r,k);
    if(r_out) memcpy(r_out,r,sizeof(double)*n);
    free(mu);free(r);free(b);
    return tour;
}
