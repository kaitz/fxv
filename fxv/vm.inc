// VM for fxv

// ilog2
// returns floor(log2(x)), e.g. 30->4  31->4  32->5,  33->5
#ifdef _MSC_VER
#include <intrin.h>
inline U32 ilog2(U32 x) {
  DWORD tmp=0;
  if(x!=0)_BitScanReverse(&tmp,x);
  return tmp;
}
#elif __GNUC__
inline U32 ilog2(U32 x) {
  if(x!=0)x=31-__builtin_clz(x);
  return x;
}
#else
inline U32 ilog2(U32 x) {
  //copy the leading "1" bit to its left (0x03000000 -> 0x03ffffff)
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >>16);
  //how many trailing bits do we have (except the first)? 
  return BitCount(x >> 1);
}
#endif
#include "vector.inc"


#include "vmh.inc"
#include "mmap.inc"
#include "distance.inc"
#include "column.inc"
#include "table.inc"
#include "record.inc"
#include "dictionary.inc"

void VM::initComp(){
    //printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",ds,mx,smc,apm1,rcm,scm,cm,st,av,mm,dhs,sm,sk,apm2,em,uas,lmx,statec,byt);
    
    doGetComponents=false;
    if (totalc) alloc(comp,totalc);
    if (ds) alloc(dsA,ds);
    if (mx) alloc(mxA,mx);
    if (smc) alloc(smA,smc);
    if (apm1) alloc(apm1A,apm1);
    if (rcm) alloc(rcmA,rcm);
    if (scm) alloc(scmA,scm);
    if (cm) alloc(cmC,cm);
    if (st) alloc(stA,st);
    if (av) alloc(avA,av);
    if (mm) alloc(mmA,mm);
    if (dhs) alloc(dhsA,dhs);
    if (sm) alloc(smcA,sm);
    if (sk) alloc(skA,sk);
    if (apm2) alloc(apm2A,apm2);
    if (em) alloc(emA,em);
    if (uas) alloc(uasA,uas);
    if (lmx) alloc(lmxA,lmx);
    statec++; // for first statetable
    if (statec) alloc(vmstate,statec);
    if (byt) alloc(bmA,byt);
    vmstate[0].Init(22,56,16,5,1,18,30); // first default, needs better parameters 13 states at the end not reachable
}


VM::~VM()  {
    if (smc>0) {
       for (int i=0;i<smc;i++)  smA[i].Free(); 
       free(smA);
    }
    if (apm1>0) {
        for (int i=0;i<apm1;i++) apm1A[i].Free();
        free(apm1A);
    }
    if (apm2>0) {
        for (int i=0;i<apm2;i++) apm2A[i].Free();
        free(apm2A);
    }
    if (scm>0) {
        for (int i=0;i<scm;i++)  scmA[i].Free();
        free(scmA);
    }
    if (cm>0) {
        for (int i=0;i<cm;i++)  cmC[i].Free(); 
        free(cmC);
    }
   if (mx>0) {
        for (int i=0;i<mx;i++) mxA[i].Free();
        free(mxA);
    }
    if (ds>0) {
        for (int i=0;i<ds;i++)  dsA[i].Free();
        free(dsA);
    }
    if (dhs>0) {
        for (int i=0;i<dhs;i++)  dhsA[i].Free();
        free(dhsA);
    }
    if (sm>0) {
        for (int i=0;i<sm;i++)  smcA[i].Free();
        free(smcA);
    }
    if (rcm>0) {
        for (int i=0;i<rcm;i++)  rcmA[i].Free();
        free(rcmA);
    }
    if (uas>0) {
        for (int i=0;i<uas;i++)  uasA[i].Free();
        free(uasA);
    }
    if (lmx>0) {
        free(lmxA);
    }
    if (em>0) {
        free(emA);
    }
    if (byt>0) {
        free(bmA);
    }
    if (st>0) {
        free(stA);
    }
    if (av>0) {
        free(avA);
    }
    if (sk>0) {
        free(skA);
    }
    if (mm>0) {
        free(mmA);
    }
    if (statec>0) {
        free(vmstate);
    }
    if (totalc) free(comp);
    
    // free memory allocated by vmmalloc
    if (mindex){
        for (int i=0;i<mindex;i++){
            //printf("Free %x\n",mem[i]);
            free(mem[i]);
        }
    }
    smc=apm1=apm2=rcm=scm=cm=mx=mm=st=av=ds=dhs=sm=sk=lmx=currentc=totalc=mindex=compIdx=statec=byt=uas=em=0;
    if (sym) free(sym),sym=0;

    if (text) free(text),text=0;
    if (sp0) free(sp0),sp0=sp=0;
    if (vmMode==VMCOMPRESS && doJIT) {
        UnmapViewOfFile(jitmem);
        UnmapViewOfFile(jitmempr);
        UnmapViewOfFile(jitmemup);
    }
}

void printcomponent(int component) {
    printf("%s",cNames[component-1]);
}
void setComponents(VM* v) {
    v->totalc=v->smc+v->apm1+v->ds+v->av+v->scm+v->rcm+v->cm+v->mx+v->st+v->mm+v->dhs+v->sm+v->sk+v->apm2+v->em+v->uas+v->lmx+v->statec+v->byt;
    v->mcomp.resize(v->totalc); 
    if (v->totalc==0 && v->mx>0) quit("No inputs for mixer defined VM\n");
}
void setParameter(Parameter *p, int id, int idx, int val, int max, int min){
    p->ID=id;
    p->IDX=idx;
    p->param=val;
    p->max=max;
    p->min=min;
}
#ifdef PRINTTUNECOMPONENTS
void setTuneComponents(COMPONENT *p, int id, int idx){
    p->id=id;
    p->idx=idx;
    p->flags=0;
}
#endif
// vmi - init components
void initcomponent(VM* v, int component, int componentIndex, int f, int d, int indexOfInputs, int iParm6) {
    // Collect info about components (stage 1)
    if (v->doGetComponents){
        if (component > vmBYT)  {kprintf("VM vmi error: unknown component id %d\n ",component);quit();}
        if (componentIndex > 255)  {kprintf("VM vmi error: %s(%d) > 255\n",cNames[component],componentIndex);quit();}
        // does component has inputs for mixers
        bool isInputs = (component==vmAPM1 || component==vmAPM2 || component==vmUAS||component==vmDS||component==vmERR||component==vmBYT||component==vmSTA|| component==vmDHS && indexOfInputs==-1 || component==vmAVG ||
         (component==vmST && indexOfInputs==-1)||(component==vmSK && indexOfInputs==-1)||(component==vmSMC && indexOfInputs==-1));
        if (indexOfInputs>=0 &&  v->x.cInputs < indexOfInputs && isInputs==false){// input sets for mixers
           v->x.cInputs++;
        }
        //printf("%s(%d) indexOfInputs=%d\n ",cNames[component-1],componentIndex,indexOfInputs);
        switch (component) {
            case vmSMC: {
                v->smc++;
                if (indexOfInputs==-1)v->totalPR++; //uses direct
                else if (indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount++; // use mixer
                v->posParamCount=v->posParamCount+1;
                break; 
            }
            case vmAPM1:{
                v->totalPR++; 
                v->apm1++;
                v->posParamCount=v->posParamCount+1;
                break; 
            }
            case vmAPM2: {
                v->totalPR++; 
                v->apm2++;
                v->posParamCount=v->posParamCount+2;
                break; 
            }
            case vmUAS: {
                v->totalPR++;
                v->uas++;
                v->posParamCount=v->posParamCount+1;
                break;
            }   
            case vmDS: {
                v->ds=v->ds+indexOfInputs; // we can have multiple predictions
                v->totalPR+=indexOfInputs; 
                if (f<1) {kprintf("VM vmi error:ds(%d) memory bits must be larger then 0.",v->ds);quit();}
                v->posParamCount=v->posParamCount+1;
                break;      
            }
            case vmDHS: {
                if (indexOfInputs==-1) v->dhs=v->dhs+iParm6; // we can have multiple predictions
                else {
                    v->x.mxInputs[indexOfInputs].ncount=v->x.mxInputs[indexOfInputs].ncount+iParm6;
                    v->dhs=v->dhs+1;
                }
                if (indexOfInputs==-1) v->totalPR+=iParm6;
                if (f<1) {kprintf("VM vmi error:dhs(%d) memory bits must be larger then 0.",v->dhs);quit();}
                break;      
            }
            case vmRCM: {
                v->rcm++;
                if (indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount++; // use mixer
                v->posParamCount=v->posParamCount+1;
                break;  
            }
            case vmSCM: {
                v->scm++;
                if (indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount+=2; // use mixer
                break;  
            }
            case vmAVG: {
                v->av++;
                v->totalPR++;
                v->posParamCount=v->posParamCount+2;
                break;  
            }
            case vmLMX: {
                v->lmx++;
                v->totalPR++;
                v->posParamCount=v->posParamCount+1;
                break;  
            }
            case vmCM: {
                v->cm++;
                if (indexOfInputs>=0) v->x.mxInputs[(indexOfInputs)].ncount+=6 *(d&255); // 6* num of contexts
                v->posParamCount=v->posParamCount+6;
                break;  
            }
            case vmMX: {
                v->mx++;
                v->totalPR++;
                v->posParamCount=v->posParamCount+2;//3;
                break;  
            }
            case vmST: {
                v->st++;
                if (indexOfInputs==-1) v->totalPR++; //uses direct
                else if (indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount++; // use mixer
                v->posParamCount=v->posParamCount+1;
                break;  
            }
            case vmMM: {
                v->mm++;
                v->x.mxInputs[indexOfInputs].ncount++;
                break;  
            }
            case vmSM: {
                v->sm++;
                if ( indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount+=2;
                v->posParamCount=v->posParamCount+1;
                break;  
            }
            case vmSK: {
                v->sk++;
                if (indexOfInputs==-1) v->totalPR++; 
                if (indexOfInputs>=0) v->x.mxInputs[indexOfInputs].ncount+=1;
                break;  
            }
            case vmERR: {
                v->em++;
                v->posParamCount=v->posParamCount+2;
                break;  
            } 
            case vmSTA: {
                if (v->statec<componentIndex) quit("VM stage 1. bad STA id\n");
                v->statec++;
                v->posParamCount=v->posParamCount+7;
                break;  
            }     
            case vmBYT: {
                v->byt++;
                v->posParamCount=v->posParamCount+1;
                break;  
            } 
            default: 
                quit("VM stage 1. Bad component ID\n");
        }
        // realloc component memory
        setComponents(v);
        
        int prindex=0;
        // Get next free PR index if any
        if (component==vmAPM1 || component==vmAPM2 || component==vmUAS || component==vmDS || component==vmDHS && indexOfInputs==-1 || component==vmAVG || component==vmLMX ||
            (component==vmST && indexOfInputs==-1)||(component==vmSK && indexOfInputs==-1)||(component==vmSMC && indexOfInputs==-1) || component==vmMX) prindex=v->totalPR;
        
        // Get mixer input index
        int m=indexOfInputs;
        if (component==vmAVG || component==vmLMX || component==vmAPM1 || component==vmAPM2 ||  component==vmERR || component==vmBYT || component==vmSTA || component==vmUAS) m=0;
        if (indexOfInputs==-1) m=0;
        // Look for same component and complain to user if values are wrong
        const int reqMixer=(component==vmDS)?0:m;
        for (int j=v->currentc;j>=0 ;j--) {
            if (v->mcomp[j].component == component){
                if ((v->mcomp[j].compIndex+1)!=componentIndex ) {
                    printf("Error: Component %s(%d) redefined or not in order. Index must be %d, defined %d.\n", cNames[component-1],v->mcomp[j].compIndex,v->mcomp[j].compIndex+1,componentIndex);
                    exit(1);
                } else {
                    break;
                }
            }
        }
        // DS and DHS output only to PR index
        if (component==vmDS ||  component==vmDHS && indexOfInputs==-1 ) {
            int totalcount=0;
            if (component==vmDS) totalcount=indexOfInputs;
            else if (component==vmDHS) totalcount=iParm6;
            for (int j=0;j< totalcount;j++) {
                v->mcomp[v->currentc].mixerIndex = reqMixer;
                v->mcomp[v->currentc].component = component;
                v->mcomp[v->currentc].compIndex = componentIndex;
                v->mcomp[v->currentc].prIndex = (prindex-totalcount+j+1);
                v->currentc++;
            }
        }else {
            // 0xppiiccmm pr index,component index, component, mixer index
            v->mcomp[v->currentc].mixerIndex = reqMixer;
            v->mcomp[v->currentc].component = component;
            v->mcomp[v->currentc].compIndex = componentIndex;
            v->mcomp[v->currentc].prIndex = prindex;
            v->currentc++;
        }
        if (doDebugInfo==true) {
            int pri=v->mcomp[v->currentc-1].prIndex;
            printf("0x%02x%02x%02x%02x ", v->mcomp[v->currentc-1].prIndex,v->mcomp[v->currentc-1].compIndex,v->mcomp[v->currentc-1].component,v->mcomp[v->currentc-1].mixerIndex);
            if (pri==0) printf("       ("); else printf(" pr[%d](",pri-1);
            printcomponent(v->mcomp[v->currentc-1].component);
            printf("[%3d]) input[%d] ",v->mcomp[v->currentc-1].compIndex,v->mcomp[v->currentc-1].mixerIndex);
            //if (v->mcomp[v->currentc-1]&0xff=vmMM) printf("pr[%d] ",v->mcomp[v->currentc-1]>>16&0xff);
            printf("\n");
        }
        return;
    }
#ifdef FXTUNE    
    // Alloc parameter memory
    if (v->parm){
        if (v->parm->isactive){
            if (v->parm->parameters==0){
                alloc(v->parm->parameters,sizeof(Parameter)*v->posParamCount);
                v->parm->total=v->posParamCount;
            }
#ifdef PRINTTUNECOMPONENTS
            if (v->parm->components==0){
                alloc(v->parm->components,sizeof(COMPONENT)*v->totalc);
                v->parm->totalComponents=v->totalc;
            }
#endif
        }
    }
#endif
    // Set components (stage 2)
    // If Autotune then ignore model parameters, first run is allways with model parameters.
    switch (component) {
    case vmSMC:{
        int smc_l=d;
        if (smc_l==0) smc_l=1023;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmSMC,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmSMC,componentIndex,smc_l,1023,1);
                } else {
                    smc_l=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        U8 *n=&v->vmstate[0].ns[0];
        v->smA[componentIndex].Init(&v->x,f, smc_l,n);
        break;  
    } 
    case vmAPM1: {   
        int apm_l=d;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmAPM1,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmAPM1,componentIndex,apm_l,16,1);
                } else {
                    apm_l=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        v->apm1A[componentIndex].Init(f,apm_l,indexOfInputs);
        break;
    }
    case vmAPM2: {   
        int apm_l=d>>8;
        int apm_s=d&255;
        if (apm_l==0) apm_l=0;
        if (apm_s==0) apm_s=24;
#ifdef FXTUNE        
        if (v->parm){
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmAPM2,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmAPM2,componentIndex,apm_l,254,0);
                    setParameter(&parm->parameters[parm->current++],vmAPM2,componentIndex,apm_s,31,2);
                } else {
                    apm_l=parm->parameters[parm->current++].param;
                    apm_s=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        v->apm2A[componentIndex].Init(f,(apm_l<<8)+apm_s,indexOfInputs);
        break;
    }
    case vmSTA: {   
        int sta_l0=f&0xffff;
        int sta_l1=f>>16;
        int sta_l2=d&0xffff;
        int sta_l3=(d>>16);
        int sta_l4=indexOfInputs&0xffff;
        int sta_l5=(indexOfInputs>>16)&255;
        int sta_l6=(indexOfInputs>>24)&255;
        int isSta=sta_l0|sta_l1|sta_l2|sta_l3|sta_l4|sta_l5|sta_l6;
        if (isSta==0)sta_l0=22,sta_l1=56,sta_l2=16,sta_l3=5,sta_l4=1,sta_l5=18,sta_l6=30; // use STA 0
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmSTA,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l0,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l1,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l2,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l3,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l4,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l5,63,0);
                    setParameter(&parm->parameters[parm->current++],vmSTA,componentIndex,sta_l6,63,0);
                } else {
                    sta_l0=parm->parameters[parm->current++].param;
                    sta_l1=parm->parameters[parm->current++].param;
                    sta_l2=parm->parameters[parm->current++].param;
                    sta_l3=parm->parameters[parm->current++].param;
                    sta_l4=parm->parameters[parm->current++].param;
                    sta_l5=parm->parameters[parm->current++].param;
                    sta_l6=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        assert((componentIndex+1) < v->statec);
        v->vmstate[componentIndex+1].Init(sta_l0,sta_l1,sta_l2,sta_l3,sta_l4,sta_l5,sta_l6);
        break;
    }
    case vmUAS: {   
        int bits=f;
        int mask=(d==0)?((1<<bits)-1):d;
        //bool domask=false;
        int rate=indexOfInputs;
        if (rate==0) rate=5;
#ifdef FXTUNE
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmUAS,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmUAS,componentIndex,bits,24,1);
                } else {
                    bits=parm->parameters[parm->current++].param;
                    mask=((1<<bits)-1);
                }
            }
            /*if (parm->vm_uasm) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmUAS,componentIndex,mask,(1<<bits)-1,1);
                } else {
                    mask=parm->parameters[parm->current++].param;
                }
                //domask=true;
                //printf("%x\n",mask);
            }*/
            // needs work
            /*if (v->parm->vm_uasm[componentIndex]){
                if (v->parm->isactive==true) {
                    v->parm->vm_uas_mask[componentIndex]=(1<<bits)-1;
                }
                mask=v->parm->vm_uas_mask[componentIndex];
                domask=true;
            }
            if (v->parm->vm_uasr[componentIndex]){
                if (v->parm->isactive==true) {
                    v->parm->vm_uas_rate[componentIndex]=rate;
                }
                rate=v->parm->vm_uas_rate[componentIndex];
            }*/
        }
#endif
        //vm_uas_mask_max[componentIndex]=(1<<bits)-1; //set max mask
        v->uasA[componentIndex].Init(bits,mask,rate);
        break;
    }
    case vmDS: {
        int ds_l=U32(d)&0xffff;
        int stable=f>>16;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmDS,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmDS,componentIndex,ds_l,1023,1);
                } else {
                    ds_l=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        U8 *n=&v->vmstate[stable].ns[0];
        v->dsA[componentIndex].Init(&v->x,f&255,ds_l,indexOfInputs,n);
        break; 
    }
    case vmDHS: {
        int stable=(f>>16)&255;
        int slimit=4;
        if (iParm6) slimit=iParm6;
        U8 *n=&v->vmstate[stable].ns[0];
        v->dhsA[componentIndex].Init(&v->x,f&255,d,iParm6,slimit,n);
        break; 
    }
    case vmRCM: {    
        int rcm_ml=d&255;
        if (rcm_ml==0) rcm_ml=8;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmRCM,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmRCM,componentIndex,rcm_ml,30,1);
                } else {
                    rcm_ml=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        v->rcmA[componentIndex].Init(&v->x,f<=0?4096:f*4096,rcm_ml);
        break;
    }
    case vmSCM: {
#ifdef FXTUNE    
        /*if (v->parm){
            VMParam *parm=v->parm;
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmSCM,componentIndex);
        }*/
#endif
        v->scmA[componentIndex].Init(f); 
        break;
    }
    case vmAVG: {
        int avg_l0=f&255;
        int avg_l1=(f>>8)&255;  
        if (avg_l0==0) avg_l0=1;
        if (avg_l1==0) avg_l1=1;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmAVG,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (avg_l0==avg_l1 && avg_l0>1) avg_l0++;
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmAVG,componentIndex,avg_l0,63,1);
                    setParameter(&parm->parameters[parm->current++],vmAVG,componentIndex,avg_l1,63,1);
                } else {
                    avg_l0=parm->parameters[parm->current++].param;
                    avg_l1=parm->parameters[parm->current++].param;
                }
                
                U32 f0;
                if (ispowerof2((avg_l0+avg_l1))) {
                }
                else {
                f0= nextPOTwo(avg_l0+avg_l1);
                avg_l1=  f0  -avg_l0;
                }
            }
            /*
            if (v->parm->vm_avg[componentIndex]){
                if (v->parm->isactive==true) v->parm->vm_avg_limit0[componentIndex]=avg_l0;
                if (v->parm->isactive==true) v->parm->vm_avg_limit1[componentIndex]=avg_l1;
                avg_l0=v->parm->vm_avg_limit0[componentIndex];
                avg_l1=v->parm->vm_avg_limit1[componentIndex];
                if (avg_l0==avg_l1 && avg_l0>1) avg_l0++,v->parm->vm_avg_limit0[componentIndex]=avg_l0;
                U32 f0;
                if (ispowerof2((avg_l0+avg_l1))) {
                }
                else {
                f0= nextPOTwo(avg_l0+avg_l1) ;
                avg_l1=  f0  -avg_l0;
                }
            }*/
        }
#endif
        v->avA[componentIndex].Init(indexOfInputs&0xff,(indexOfInputs>>8)&255,f,d);
        break;
    }
    case vmLMX: {
        int lmx_l0=f&255;
        int lmx_l1=(f>>8)&255;  
        int w=d;
        if (w==0) w=2048;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmLMX,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmLMX,componentIndex,w,4095,1);
                } else {
                    w=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        v->lmxA[componentIndex].Init(lmx_l0,lmx_l1,w);
        break;
    }
    case vmERR: {
        int e_h=(U32(f)>>16);
        if (e_h==0) e_h=4095;
        int e_l=U32(f)&0xffff;
        if (e_l==0) e_l=2047;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmERR,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmERR,componentIndex,e_l,2047,1);
                    setParameter(&parm->parameters[parm->current++],vmERR,componentIndex,e_h,4094,1);
                } else {
                    e_l=parm->parameters[parm->current++].param;
                    e_h=parm->parameters[parm->current++].param;
                }
            }
            //printf("%d %d\n",e_l,e_h);
        }
#endif
        v->emA[componentIndex].Init(&v->x,e_l,e_h);
        break;
    } 
    case vmBYT: {
        int b_v=U32(f);
        int b_l=U32(d);
#ifdef FXTUNE        
        if (v->parm ){
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmBYT,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmBYT,componentIndex,b_v,b_l,0);
                } else {
                    b_v=parm->parameters[parm->current++].param;
                    //printf("%d\n",b_v);
                }
            }
        }
#endif        
        v->bmA[componentIndex].Init(b_v);
        break;
    }       
    case vmCM: {
        int cms3=(iParm6>>8)&255;
        int cms4=(iParm6>>16)&255;
        int cm_l=(U32(d)>>8)&255;
        if (cm_l==0) cm_l=4;
        int cms_l=(U32(d)>>16)&255;
        if (cms_l==0) cms_l=32;
        int cms2_l=(U32(d)>>24)&255;
        if (cms2_l==0) cms2_l=12;
        if (cms3==0) cms3=32;
        if (cms4==0) cms4=12;
        int mem=U32(f)&0xffffff;
        int stindex=U32(f)>>24;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmCM,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmCM,componentIndex,cm_l,30,0);
                    setParameter(&parm->parameters[parm->current++],vmCM,componentIndex,cms_l,60,1);
                    setParameter(&parm->parameters[parm->current++],vmCM,componentIndex,cms2_l,32,1);
                    setParameter(&parm->parameters[parm->current++],vmCM,componentIndex,cms3,60,1);
                    setParameter(&parm->parameters[parm->current++],vmCM,componentIndex,cms4,60,1);
                } else {
                    cm_l=parm->parameters[parm->current++].param;
                    cms_l=parm->parameters[parm->current++].param;
                    cms2_l=parm->parameters[parm->current++].param;
                    cms3=parm->parameters[parm->current++].param;
                    cms4=parm->parameters[parm->current++].param;
                }
            }
        }
#endif
        U8 *n=&v->vmstate[stindex].ns[0];U8 *n1=&v->vmstate[stindex].nn01[0];
        v->cmC[componentIndex].Init(mem<=0?4096:mem*4096,d&255,v->x,cm_l,cms_l,cms2_l,cms3,n,n1,cms4);
        break;
    }
    case vmMX: {
        int mx_err=(f>>8)&0xffff;
        int mx_sh=f&255;
        if (mx_sh==0)mx_sh=64;
        int mx_ue=f>>24;
        if (mx_ue==0)mx_ue=28;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmMX,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    //setParameter(&parm->parameters[parm->current++],vmMX,componentIndex,mx_err,32768,0);
                    setParameter(&parm->parameters[parm->current++],vmMX,componentIndex,mx_sh,255,1);
                    setParameter(&parm->parameters[parm->current++],vmMX,componentIndex,mx_ue,112,1);
                }else{
                   // mx_err=parm->parameters[parm->current++].param;
                    mx_sh=parm->parameters[parm->current++].param;
                    mx_ue=parm->parameters[parm->current++].param;
                }
            }
        }
#endif        
        v->mxA[componentIndex].Init(d,mx_sh,mx_err,mx_ue);
        break;
    }
    case vmST: {
        int range=f;
        if (f==0) range=128;
#ifdef FXTUNE        
        if (v->parm) {
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmST,componentIndex);
#endif
            if (parm->vm_comp[component][componentIndex]) {
                if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmST,componentIndex,range,255,0);
                } else {
                    range=parm->parameters[parm->current++].param;
                }
            }
        }
#endif        
        v->stA[componentIndex].Init(range,indexOfInputs);
        break;
    }
    case vmMM: {
#ifdef FXTUNE    
        /*if (v->parm){
            VMParam *parm=v->parm;
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmMM,componentIndex);
        }*/
#endif        
        v->mmA[componentIndex].Init(d,f);
        break;
    }
    case vmSM: {
        int sm_l=(d>>8)&255;
        int sm_b=d&255;  
        if (sm_l==0)sm_l=8<<1;
#ifdef FXTUNE        
        if (v->parm){
            VMParam *parm=v->parm;
#ifdef PRINTTUNECOMPONENTS
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmSM,componentIndex);
#endif
            if (v->parm->vm_comp[component][componentIndex]) {
                   if (v->parm->isactive==true) {
                    setParameter(&parm->parameters[parm->current++],vmSM,componentIndex,sm_l,32,1);
                } else {
                    sm_l=parm->parameters[parm->current++].param;
                }
            }
        }
#endif        
        v->smcA[componentIndex].Init(&v->x,f,sm_b,sm_l); 
        break;
    }
    case vmSK: {
#ifdef FXTUNE    
        /*if (v->parm){
            VMParam *parm=v->parm;
            if (parm->cmpCurrent < v->totalc) setTuneComponents(&parm->components[parm->cmpCurrent++],vmSK,componentIndex);
        }*/
#endif
        v->skA[componentIndex].Init(); 
        break;
    }
    default:
        quit("VM vmi error\n");
        break;
    }
}

void initcomponent1(const int iParm6, const int indexOfInputs, const int d, const int f,const int componentIndex, const int component, VM* v) {
    initcomponent(v, component, componentIndex, f, d, indexOfInputs, iParm6); //reverse order
}

//set context to component
int setcomponent(VM* v, int c, int i, U32 f) {
    int a=0;
    switch (c&255) {
        case vmSMC: {
            v->smA[i].set(f);
            break;
        }
        case vmAPM1: {
            v->apm1A[i].set(f);
            break;
        }
        case vmAPM2: {
            v->apm2A[i].set(f);
            break;
        }
        case vmDS: {
            v->dsA[i].set(f);
            break;
        }
        case vmDHS: {
            if (((c>>8)&255)==1)
                  a= v->dhsA[i].setbit(f);
                else
            a=v->dhsA[i].set(f);
            break;
        }
        case vmRCM: {
            v->rcmA[i].set(f);
            break;
        }
        case vmSCM: {
            v->scmA[i].set(f);
            break;
        }
        case vmCM: {
            
            if (v->isDualMode==true){
                if (v->isBitUpdate==true) {
                    a=v->cmC[i].get();
                }else{
                    if (v->x.bpos==0) {
                        if (((c>>8)&255)==1)
                            v->cmC[i].sets();
                        else
                            v->cmC[i].set(f);
                    }
                }
            }else{
                if (((c>>8)&255)==1)
                   v->cmC[i].sets();
                else
                    v->cmC[i].set(f);
            }
            break;
        }
        case vmUAS: {
            v->uasA[i].set(f);
            break;
        }
        case vmMX: {
            v->mxA[i].set(f);
            break;
        }
        case vmST: {
            v->stA[i].set(f);
            break;
        }
        case vmSM: {
            v->smcA[i].set(f);
            break;}
        case vmSK: {
            v->skA[i].set(f);
            break;
        }
        case vmERR: {
            a=v->emA[i].set(v->prSize[v->totalPR]);
            break;
        }
        case vmBYT: {
            a=v->bmA[i].set(f);
            break;
        }
        case vmAVG:{
            v->avA[i].set(f);
            break;
        }
        default:{
             quit("VM vmx error\n");
             break;}
    }
    return a;
}

int setcomponent1( U32 f, int i, int c, VM* v) {
    return setcomponent(v, c, i, f);
}

//  i    size
//  pos     -2  - Seek to pos
//  0       -2  - Seek to end
//  0       -1  - Seek to start
int readfile(VM* v, U8 *i, int size) {
    assert(size>-3); 
    assert(v->inFile!=NULL);
    if (size>0) return fread (i, 1, size,v->inFile);//,printf("Read\n");
    if (size==-2) fseek (v->inFile , 0 , SEEK_END),printf("Read Seek end\n");
    else if (size==-1) fseek(v->inFile, v->inpos, SEEK_SET),printf("Read Seek set %d\n",v->inpos); // set to block start pos not file start pos
    else return -1;
    return 0;
}

int readfile1(int size, U8 *i, VM* v) {
    return readfile(v, i, size);
}
int writefile(VM* v,U8 *i,int size){
    assert(size>-3);
    assert(v->outFile!=NULL);
    if (size==-2){ 
        fseek (v->outFile, 0, SEEK_END),printf("Write Seek end\n");
        return ftell(v->outFile);
    }
    else if (size==-1) fseek(v->outFile, 0, SEEK_SET),printf("Write Seek set\n");
    else return fwrite (i, 1, size, v->outFile);//,printf("Write\n");
    return -1;
}

int writefile1(int size, U8 *i,VM* v) {
    return writefile(v, i, size);
}
VM::VM(object &obj, BlockData& bd, int mode, VMParam *p):
    data1(2024*1024),mod(obj),x(bd),vmMode(mode),mem(0),memSize(0),membound(0),prSize(0),mcomp(0),parm(p) {    
    // Init blockdata
    x.Init(); 
    data=&data1[0];
    smc=apm1=apm2=uas=rcm=scm=cm=mx=st=av=mm=ds=dhs=sm=sk=statec=byt=lmx=currentc=totalc=mindex=totalPR=em=compIdx=posParamCount=0;
    debug=0; 
    doGetComponents=false;
    isDualMode=false;
    isBitUpdate=false;
    // Load cfg file, if error then exit
    if (initvm()==-1) {
        exit(1);
    }

    if (totalc!=currentc) printf("WARNING: Component counts differ for init and runtime: total %d current %d\n",totalc,currentc);
    if (vmMode==VMCOMPRESS && totalPR<1) printf("Error: No predictors provided.\n"),quit("");
    totalc=currentc; //update total count to current count 
    prSize.resize(totalPR--);

    // if mixer is used parse all input arrays
    if(x.cInputs>=0 && x.cInputs<256&& (idupdate->Val || isDualMode==true)) {
        // init resize inputs[x] arrays to correct size
        for (int j=0;j<x.cInputs+1;j++)  {
            x.mxInputs[j].ncount=(x.mxInputs[j].ncount+15)&-16;
            x.mxInputs[j].Init(x.mxInputs[j].ncount);
        }     
        // provide inputs array info to mixers
        for (int i=0;i<totalc;i++){
            // individual components
            if (mcomp[i].prIndex>0 && mcomp[i].component==vmMX){
                int index=mcomp[i].compIndex;
                int input=mcomp[i].mixerIndex;
                // set mixer input size and let it alloc weight memory 
                mxA[index].setTxWx(x.mxInputs[input].size,&x.mxInputs[input].n[0]);
            }
        }
    }

    if (vmMode==VMCOMPRESS && doJIT) {
        if (isDualMode==true)
            getBytePredictionJIT();
        getPredictionJIT();
    }
}


void VM::next() {
    char *pp;
    int nu;
    while (tk = *p) {
        ++p;
        if (tk == '\n') {
            ++line;
        }
        else if (tk == '#') {
            while (*p != 0 && *p != '\n') ++p;
        } else if ((tk >= 'a' && tk <= 'z') || (tk >= 'A' && tk <= 'Z') || tk == '_') {
            pp = p - 1;
            while ((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z') || (*p >= '0' && *p <= '9') || *p == '_')
                tk = tk * 147 + *p++;
            tk = (tk << 6) + (p - pp);
            id = sym;
            while (id->Tk) {
                if (tk == id->Hash && !memcmp(id->Name, pp, p - pp)) { tk = id->Tk; id->ePTR=int(e-text);return; }
                
                id = id + 1;
                idCount++;
            }
            id->Name = pp;
            id->Hash = tk;
            // id[IDLen] = p - pp;
            tk = id->Tk = Id;
            return;
        } else if (tk == '0' && *(p)== 'x') { //Hexadecimal numbers
            p++;
            for (ival = 0; '\0' != (nu = *p); p++) {
                if ( nu >= 'a' && nu <= 'f') {
                    nu = nu - 'a' + 10;
                } else if (nu >= 'A' && nu <= 'F') {
                    nu = nu - 'A' + 10;
                } else if (nu >= '0' && nu <= '9') {
                    nu = nu - '0';
                } else {
                    tk = Num;
                    return;
                }
                ival = ival<<4;
                ival = ival + nu;
            }
        } else if (tk >= '0' && tk <= '9') { //numbers
            ival = tk - '0';
            while (*p >= '0' && *p <= '9') ival = ival * 10 + *p++ - '0';
            tk = Num;
            return;
        } else if (tk == '/') { //comment
            if (*p == '/') {
                ++p;
                while (*p != 0 && *p != '\n') ++p;
            } else {
                tk = Div;
                return;
            }
        } else if (tk == '\'' || tk == '"') {
            pp = data;
            while (*p != 0 && *p != tk) {
                if ((ival = *p++) == '\\') {
                    if ((ival = *p++) == 'n') ival = '\n';
                }
                if (tk == '"') *data++ = ival;
            }
            ++p;
            if (tk == '"') ival = (int)pp; else tk = Num;
            return;
        } else if (tk == '=') { if (*p == '=') { ++p; tk = Eq; } else tk = Assign; return; }
        else if (tk == '+') { if (*p == '+') { ++p; tk = Inc; } else tk = Add; return; }
        else if (tk == '-') { if (*p == '-') { ++p; tk = Dec; } else tk = Sub; return; }
        else if (tk == '!') { if (*p == '=') { ++p; tk = Ne; } return; }
        else if (tk == '<') { if (*p == '=') { ++p; tk = Le; } else if (*p == '<') { ++p; tk = Shl; } else tk = Lt; return; }
        else if (tk == '>') { if (*p == '=') { ++p; tk = Ge; } else if (*p == '>') { ++p; tk = Shr; } else tk = Gt; return; }
        else if (tk == '|') { if (*p == '|') { ++p; tk = Lor; } else tk = Or; return; }
        else if (tk == '&') { if (*p == '&') { ++p; tk = Lan; } else tk = And; return; }
        else if (tk == '^') { tk = Xor; return; }
        else if (tk == '%') { tk = Mod; return; }
        else if (tk == '*') { tk = Mul; return; }
        else if (tk == '[') { tk = Brak; return; }
        else if (tk == '?') { tk = Cond; return; }
        else if (tk == ',') { tk = Comma; return;}
        else if (tk == '~' || tk == ';' || tk == '{' || tk == '}' || tk == '(' || tk == ')' || tk == ']' || tk == ':') return;
    }
}

void VM::expr(int lev){
    int t, *b, fc;
    struct ident_s *d;
    if (!tk) { kprintf("%d: unexpected eof in expression\n", line); exit(-1); }
    else if (tk == Num) { *++e = IMM; *++e = ival; next(); ty = iINT; }
    else if (tk == '"') {
        *++e = IMM; *++e = ival;  next(); 
        while (tk == '"') next();
        data = (char *)(((int)data + sizeof(int)) & -sizeof(int)); ty = PTR;
    }
    else if (tk == Sizeof) {
        next(); if (tk == '(') next(); else { kprintf("%d: open paren expected in sizeof\n", line); exit(-1); }
        ty = iINT; if (tk == Int) next(); else if (tk == Char) { next(); ty = rCHAR; } else if (tk == Short) { next(); ty = sSHORT; }
        while (tk == Mul) { next(); ty = ty + PTR; }
        if (tk == ')') next(); else { kprintf("%d: close paren expected in sizeof\n", line); exit(-1); }
        *++e = IMM; *++e = (ty == rCHAR) ? sizeof(char) :(ty == sSHORT) ? sizeof(short) : sizeof(int);
        ty = iINT;
    }
    else if (tk == Id) {
        d = id; 
        /* nam=(char*)id[Name];
    printf("\n//");
        for (int y=0;y<id[IDLen];y++)printf("%c",nam[y]);
        if (id[UBound]>0)
        printf("// id-val 0x%08x, max size %d",id[Val], id[UBound]);
        printf("\n");*/
        next();
        if (tk == '(') {      
            if ((d->Val>=VMI &&  d->Val<READ) || d->Val==READ|| d->Val==WRTE){//for special functions in vm
                *++e = VTHIS;
                next();
                t = 1; //adjust stack
            }
            else{
                next();
                t=0;
            }
            fc=0;
            if ( (d->Val >= VMI &&  d->Val<READ) && (vmMode==VMDECODE || vmMode==VMENCODE)) {
                printf("VMI, VMX allowed only in (de)comprassion stage.");
                exit(-1);
            }
            if (d->Val ==VMI  ) {fc=7 ;}
            else if (d->Val == VMX  ) {fc=4 ;}
            else if (d->Val == READ || d->Val == WRTE  ) {
                fc=3;
                if (!(vmMode==VMDECODE || vmMode==VMENCODE)) {
                    printf("read/write allowed only in decode or encode stage.");
                    exit(-1);
                }
            }
            int val1=-1,val2=-1;
            int *curE;
            curE=e;
            while (tk != ')') { 
                expr(Assign);
                if (doJIT){
                    if (d->Val == VMX && *(e-1)==IMM && t==1 && *(e)<64) val1=*e; // parameter 1 // fixme
                    if (d->Val == VMX && *(e-1)==IMM && t==2 && *(e)<64 && ((int)(e-curE))==5 ) {
                    val2=*e;
                    }
                }
                *++e = PSH; ++t; 
                if (doJIT) {
                    // allow APM1, CM, MX, DHS
                    if (d->Val == VMX && val1>=0 && val2>=0 && (val1==1 ||val1==2 ||val1==3 ||val1==7 && idbitupdate->Val==0 ||val1==8||val1==11||val1==12||val1==14) && t==3 && idbitupdate->Val==0) {
                        e=e-7,*++e = NOP, *++e = NOP, *++e = NOP, *++e = NOP, *++e = NOP, *++e = NOP, *++e = NOP;
                        //printf("%s %d\n",cNames1[val1-1],val2);
                    }
                }
                if (tk == Comma) next(); 
            }
            
            next();
            if (d->Class == Sys) {
                if (doJIT && d->Val == VMX && val1>=0 && val2>=0 &&  (val1==1 ||val1==2 ||val1==3 ||val1==7&&idbitupdate->Val==0 ||val1==8||val1==11||val1==12||val1==14)) {
                    *++e =VMXD,*++e =val1+val2*256,t=0;
                }
                else {
                *++e =d->Val;
                
                if (t!=fc && fc!=0){ kprintf("%d: wrong number of arguments, (%d) expected %d\n", line,t,fc); exit(-1);}
                }
            }
            else if (d->Class == Fun) { *++e = JSR; *++e = d->Val; }
            else { 
                kprintf("%d: bad function call\n", line);
                exit(-1);
            }
            if (t) { *++e = ADJ; *++e = t; }
            ty = d->Type; 
            
        }
        else if (d->Class == Num) { *++e = IMM; *++e = d->Val; ty = iINT; }
        else {
            if (d->Class == Loc) { *++e = LEA; *++e = loc - d->Val; }
            else if (d->Class == Glo) { *++e = IMM; *++e = d->Val; }
            else { kprintf("%d: undefined variable\n", line); exit(-1); }
            *++e = ((ty = d->Type) == rCHAR) ? LC : ((ty = d->Type) == sSHORT) ? LS : LI;
            // add dDST
        }
    }
    else if (tk == '(') {
        next();
        if (tk == Int || tk == Char || tk == Short) {
            t = (tk == Int) ? iINT : (tk == Short) ? sSHORT : rCHAR; next();
            while (tk == Mul) { next(); t = t + PTR; }
            if (tk == ')') next(); else { kprintf("%d: bad cast\n", line); exit(-1); }
            expr(Inc);
            ty = t;
        }
        else {
            expr(Assign);
            if (tk == ')') next(); else { kprintf("%d: close paren expected\n", line); exit(-1); }
        }
    }
    /*else if (tk == Mul) { // *name
    next(); expr(Inc);
    if (ty > iINT) ty = ty - PTR; else { kprintf("%d: bad dereference\n", line); exit(-1); }
    emit((ty == rCHAR) ? LC : (ty == sSHORT) ? LS : LI);// *++e = (ty == rCHAR) ? LC : (ty == sSHORT) ? LS : LI;
}*/
    /*else if (tk == And) { // &name
    next(); expr(Inc);
    if (*e == LC || *e == LI || *e == LS) --e; else { kprintf("%d: bad address-of\n", line); exit(-1); }
    ty = ty + PTR;
}*/
    else if (tk == '!') { 
        next(); *++e = IMM;
        if (tk == Num) {*++e = !ival;next(); }else { *++e = 0; *++e = PSH; expr(Inc); *++e = EQ; }
        ty = iINT;
    }
    else if (tk == '~') {
        next(); *++e = IMM;
        if (tk == Num) {*++e = ~ival;next(); }else { *++e = -1; *++e = PSH; expr(Inc); *++e = XOR; }
        ty = iINT;
    }
    else if (tk == Add) { next(); expr(Inc); ty = iINT; }
    else if (tk == Sub) {
        next(); *++e = IMM;
        if (tk == Num) {  *++e = -ival; next(); } else {  *++e = -1; *++e = PSH; expr(Inc); *++e = MUL; }
        ty = iINT;
    }
    else if (tk == Inc || tk == Dec) {
        t = tk; next(); expr(Inc);
        if (*e == LC) { *e = PSH; *++e = LC; }
        else if (*e == LI) { *e = PSH; *++e = LI; }
        else if (*e == LS) { *e = PSH; *++e = LS; }
        else { kprintf("%d: bad lvalue in pre-increment\n", line); exit(-1); }
        *++e = PSH;
        *++e = IMM; *++e = (ty > PTR) ? sizeof(int) : (ty > iINT) ?  sizeof(short) : sizeof(char);
        *++e = (t == Inc) ? ADD : SUB;
        *++e = (ty == rCHAR) ? SC : (ty == sSHORT) ? SS : SI;
    }

    else { kprintf("%d: bad expression\n", line); exit(-1); }

    while (tk >= lev) { // "precedence climbing" or "Top Down Operator Precedence" method
        t = ty;
        if (tk == Comma) { 
            next(); expr(Assign);
        }
        else if (tk == Assign) {
            next();
            if (ty > dDST) {  printf("%d: cant assign to pointer\n", line); exit(-1); }
            else if (ty==dDST){
                int op=0;
                if (d->HClass==Distance) op=DISTANCE;
                else if (d->HClass==Column) op=COLUMN;
                else if (d->HClass==Table) op=TABLE;
                else if (d->HClass==Record) op=RECORD;
                else if (d->HClass==Dictionary) op=DICTIONARY;
                if (tk==Num && ival==0 && op!=DICTIONARY) {
                    next();
                    *e = TYPR;
                    *++e=op; 
                } else if ( op==DICTIONARY) {
                    if (tk==Id && id->Type==dDST && id->Hash==d->Hash && id->Class==Glo) printf("Bad dict assign\n");
                    //next();
                    tk == Assign;
                    int addr=*(e-1);
                    if (*e == LC || *e == LI || *e == LS) *e = PSH; else { kprintf("%d: bad lvalue in assignment\n", line); exit(-1); }
                    expr(Assign); //parse assign
                    *++e = SI;
                    *++e = TYPR;
                    *++e=op;
                    *++e = addr; // our dict address
                } else {
                    // make sure that dest = dest + val;
                    if (tk==Id && id->Type==dDST && id->Hash==d->Hash && id->Class==Glo) {
                        next(); // skip add
                        if (tk !=Add) {  printf("%d: must be add\n", line); exit(-1); }
                        tk == Assign;
                        if (*e == LC || *e == LI || *e == LS) *e = PSH; else { kprintf("%d: bad lvalue in assignment\n", line); exit(-1); }
                        expr(Assign); //parse assign
                        // SI; // TYPA will do SI itself
                        *++e=TYPA;
                        *++e=op;
                    }else {
                        printf("%d: distance bad assign\n", line); exit(-1); 
                    }
                }
                
            }else{
                
                if (*e == LC || *e == LI || *e == LS) *e = PSH; else { kprintf("%d: bad lvalue in assignment\n", line); exit(-1); }
                expr(Assign); *++e = ((ty = t) == rCHAR) ? SC : ((ty = t) == sSHORT) ? SS : SI;
            }
        }
        else if (tk == Cond) {
            next();
            *++e=BZ; b = ++e;
            expr(Assign);
            if (tk == ':') next(); else { kprintf("%d: conditional missing colon\n", line); exit(-1); }
            *b = (int)(e + 3); *++e = JMP; b = ++e;
            expr(Cond);
            *b = (int)(e + 1);
        }
        else if (tk == Lor) { next(); *++e=(LBNZ) ;  b = ++e; expr(Lan); *b = (int)(e + 1); ty = iINT; }
        else if (tk == Lan) { next(); *++e=(LBZ) ;  b = ++e; expr(Or);  *b = (int)(e + 1); ty = iINT; }
        else if (tk == Or)  { next();*++e=(PSH); expr(Xor); *++e=( OR);  ty = iINT; }
        else if (tk == Xor) { next(); *++e=(PSH); expr(And);*++e=( XOR); ty = iINT; }
        else if (tk == And) { next(); *++e=(PSH); expr(Eq);  *++e=( AND); ty = iINT; }
        else if (tk == Eq)  { next(); *++e=(PSH); expr(Lt);  *++e=( EQ);  ty = iINT; }
        else if (tk == Ne)  { next(); *++e=(PSH); expr(Lt);  *++e=( NE);  ty = iINT; }
        else if (tk == Lt)  { next(); *++e=(PSH); expr(Shl); *++e=( LT);  ty = iINT; }
        else if (tk == Gt)  { next(); *++e=(PSH); expr(Shl); *++e=( GT);  ty = iINT; }
        else if (tk == Le)  { next(); *++e=(PSH); expr(Shl); *++e=( LE);  ty = iINT; }
        else if (tk == Ge)  { next(); *++e=(PSH); expr(Shl); *++e=( GE);  ty = iINT; }
        else if (tk == Shl) { next();*++e=(PSH); expr(Add); *++e=( SHL); ty = iINT; }
        else if (tk == Shr) { next(); *++e=(PSH); expr(Add); *++e=( SHR); ty = iINT; }
        else if (tk == Add) {
            next(); *++e=(PSH); expr(Mul);
            if ((ty = t) > PTR) {*++e=(PSH); *++e=(IMM);*++e=sizeof(int); *++e=(MUL);  }//pointer
            *++e=(ADD);
        }
        else if (tk == Sub) {
            next(); *++e = PSH; expr(Mul);
            if (t > PTR && t == ty) { *++e = SUB; *++e = PSH; *++e=(IMM);*++e=(sizeof(int));*++e=(DIV); ty = iINT; }
            else if ((ty = t) > PTR) { *++e = PSH; *++e = IMM; *++e = sizeof(int); *++e=(MUL); *++e=(SUB);   }
            else *++e = SUB;
        }
        else if (tk == Mul) { next(); *++e = PSH; expr(Inc); *++e=(MUL);   ty = iINT; }
        else if (tk == Div) { next(); *++e = PSH; expr(Inc);*++e=(DIV);   ty = iINT; }
        else if (tk == Mod) { next(); *++e = PSH; expr(Inc); *++e=(MOD);   ty = iINT; }
        else if (tk == Inc || tk == Dec) {
            if (ty > PTR) {  printf("%d: cant assign to pointer\n", line); exit(-1); }
            if (*e == LC) { *e = PSH; *++e = LC; }
            else if (*e == LI) { *e = PSH; *++e = LI; }
            else if (*e == LS) { *e = PSH; *++e = LS; }
            else { kprintf("%d: bad lvalue in post-increment\n", line); exit(-1); }
            *++e = PSH; *++e = IMM; *++e = (ty > PTR) ?  sizeof(int) : (ty > iINT) ?  sizeof(short) : sizeof(char);;
            *++e = (tk == Inc) ? ADD : SUB;
            *++e = (ty == rCHAR) ? SC : (ty == sSHORT) ? SS : SI;;
            *++e = PSH; *++e = IMM; *++e = (ty > PTR) ?  sizeof(int) : (ty > iINT) ?  sizeof(short) : sizeof(char);;
            *++e = (tk == Inc) ? SUB : ADD;
            next();
        }
        else if (tk == Brak) {
            // Bounds
            int directarray=(int)e;
            unsigned int upperbound=id->UBound-1;
            int *boundptr;

            next(); 
            *++e=(PSH);
            expr(Assign);     
            if (doBounds==true){
                if (((int)e-directarray)==(4*3) && *(e-2)==PSH && *(e-1)==IMM && (unsigned int)*e>upperbound) printf("Array out of bounds: defined %d used %d, line %d\n",upperbound,*e, line), exit(-1);
            }
            if (tk == ']') next(); else { kprintf("%d: close bracket expected\n", line); exit(-1); }
            if (t > PTR) { 
                *++e=(PSH); 
                if (doBoundsRun==true){     
                    // Runtime bounds check
                    *++e=(PSH);                                              // push index value again
                    *++e=(IMM);*++e=(upperbound);*++e=(GT);  *++e=(BZ);boundptr=++e;  // compare index>upperbound
                    *++e=(IMM);*++e=(line); *++e=(PSH); *++e=(BOUND); *++e=(ADJ);*++e=(1); // fail if larger
                    *boundptr = (int)(e + 1);
                }
                *++e=(IMM);*++e=(((ty = t - PTR) == rCHAR) ? 1 : ((ty = t - PTR) == sSHORT) ? 2 : 4);
                *++e=(MUL);
            } //fixed to int !!!
            else if (t < PTR) { kprintf("%d: pointer type expected\n", line); exit(-1); }
            *++e=(ADD);
            char aa=((ty = t - PTR) == rCHAR) ? LC : ((ty = t - PTR) == sSHORT) ? LS : LI;//9 i 10 s 11 c
            *++e=(aa);
        }
        else { kprintf("%d: compiler error tk=%d\n", line, tk); exit(-1); }
    }
}

void VM::stmt() {
    int *a, *b, *c,*d;
    switch (tk) {
    case If:  
        next();
        if (tk == '(') next(); else { kprintf("%d: open paren expected\n", line); exit(-1); }
        expr(Assign);
        if (tk == ')') next(); else { kprintf("%d: close paren expected\n", line); exit(-1); }
        *++e=(BZ) ; b = ++e;
        stmt();
        if (tk == Else) {
            *b = (int)(e + 3); *++e = JMP; b = ++e;
            next();
            stmt();
        }
        *b = (int)(e + 1);
        break;
    case While:
        next();
        a = e + 1;
        if (tk == '(') next(); else { kprintf("%d: open paren expected\n", line); exit(-1); }
        expr(Comma);
        if (tk == ')') next(); else { kprintf("%d: close paren expected\n", line); exit(-1); }
        *++e=(BZ) ; b = ++e;
        stmt();
        *++e = JMP; *++e = (int)a;
        *b = (int)(e + 1);
        break;
    case For:
        next();
        if (tk == '(') next(); else { kprintf("%d: open paren expected\n", line); exit(-1); }
        if (tk != ';') expr(Comma);  //1
        next();
        a = e + 1;
        if (tk != ';') expr(Comma); //2
        *++e=(BZ) ; b = ++e;
        next();
        if (tk != ')') {
            *++e = JMP; c=e+1;*++e = (int)0; //j1
            d=e+1;
            expr(Comma); //3
            *++e = JMP;  ;*++e = (int)a;//j2
            *c = (int)(e + 1); //patch j1
            a=d; //replace jmp to //3
        }
        if (tk == ')') next(); else { kprintf("%d: close paren expected\n", line); exit(-1); }
        stmt();
        *++e = JMP; *++e = (int)a;
        *b = (int)(e + 1);
        break; 
    case Return:
        next();
        if (tk != ';') expr(Comma);
        *++e = LEV;
        if (tk == ';') next(); else { kprintf("%d: semicolon expected\n", line); exit(-1); }
        break; 
    case '{':
        next();
        while (tk != '}') stmt();
        next();
        break; 
    case ';': 
        next();
        break; 
        default :
        expr(Comma);
        if (tk == ';') next(); else { kprintf("%d: semicolon expected\n", line); exit(-1); }
    }  
}
 
int vmbound(int line){
    if (line!=0)printf("Bounds error line: %d\n",line);
    exit(-1);
}

int VM::dovm(int *ttt) {
    int i;
    if (!(pc = ttt)) { kprintf("main() not defined\n"); return -1; }
    pc0=pc;
 
    cycle = 0;
    a=0;
    while (1) {
    i = *pc++; ++cycle;
 /*if (debug) {
      kprintf("%d>%x  %.4s", cycle,pc-pc0,
        &"LEA ,IMM ,JMP ,JSR ,BZ  ,BNZ ,ENT ,ADJ ,LEV ,LI  ,LS  ,LC  ,SI  ,SS  ,SC  ,PSH ,"
         "OR  ,XOR ,AND ,EQ  ,NE  ,LT  ,GT  ,LE  ,GE  ,SHL ,SHR ,ADD ,SUB ,MUL ,DIV ,MOD ,BOUN,THIS,"
         "PRTF,VMI ,VMX ,READ,WRTE,EXIT,TYPA,TYPR,VMXD"[i * 5]);
    if (i < JMP) kprintf(" %d\n",*pc); //? +1MALC,MSET,MCMP,MCPY,
     else if (i <= ADJ) kprintf(" %x\n",(int *)*pc-pc0+1); else kprintf("\n");
    }*/
    if      (i == LEA) a = (unsigned int)(bp + *pc++);                    // load local address
    else if (i == IMM) a = (unsigned int)*pc++;                           // load global address or immediate
    else if (i == JMP) pc = (int *)*pc;                                   // jump
    else if (i == JSR) { *--sp = (int)(pc + 1); pc = (int *)*pc; }        // jump to subroutine
    else if (i == BZ||i == LBZ)  pc = a ? pc + 1 : (int *)*pc;            // branch if zero
    else if (i == BNZ ||i == LBNZ) pc = a ? (int *)*pc : pc + 1;          // branch if not zero
    else if (i == ENT) { *--sp = (int)bp; bp = sp; sp = sp - *pc++; }     // enter subroutine
    else if (i == ADJ) sp = sp + *pc++;                                   // stack adjust
    else if (i == LEV) { sp = bp; bp = (int *)*sp++; pc = (int *)*sp++; } // leave subroutine
    else if (i == LI)  a = *(unsigned int *)a;                            // load int
    else if (i == LS)  a = *(unsigned short *)a;                          // load short
    else if (i == LC)  a = *(unsigned char *)a;                           // load char
    else if (i == SI) *(unsigned int *)*sp++ = (unsigned int)a;           // store int
    else if (i == SC)   a = *(unsigned char *)*sp++ = (unsigned char)a;   // store char
    else if (i == SS)   a = *(unsigned short *)*sp++ =(unsigned short) a; // store short
    else if (i == PSH) *--sp = a;                                         // push

    else if (i == OR)  a = (unsigned int)*sp++ |  (unsigned int)a;
    else if (i == XOR)  a = (unsigned int)*sp++ ^  (unsigned int)a;
    else if (i == AND)   a = (unsigned int)*sp++ & (unsigned int)a;
    else if (i == EQ)  a = *sp++ == a;
    else if (i == NE)   a = *sp++ !=a;
    else if (i == LT)   a = *sp++ <  a;
    else if (i == GT)   a = *sp++ >  a;
    else if (i == LE)   a = *sp++ <= a;
    else if (i == GE)   a = *sp++ >= a;
    else if (i == SHL) a = (unsigned int)*sp++ << (unsigned int)a;
    else if (i == SHR) a = (unsigned int)*sp++ >> (unsigned int)a;
    else if (i == ADD)  a =*sp++ +  a;
    else if (i == SUB) a = *sp++ -  a;
    else if (i == MUL) a = *sp++ * a;
    else if (i == DIV) a = *sp++ /   a;
    else if (i == MOD) a = *sp++ % a;

    else if (i == PRTF) { t = sp + pc[1]; a = printf((char *)t[-1], t[-2], t[-3], t[-4], t[-5], t[-6]); }
    else if (i == VMI) a=0, initcomponent(this, sp[5], sp[4], sp[3],sp[2], sp[1],*sp);
    else if (i == VMX) a=setcomponent(this, sp[2], sp[1],*sp);
    else if (i == VTHIS)  --sp;  //ignore
    else if (i == BOUND) { return vmbound(*sp); }
    else if (i == TYPA) { 
        if (*pc==DISTANCE) {
            pc++;
            // our BracketContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)*sp+4;
            BracketContext<U8> *brcxt1;
            brcxt1=reinterpret_cast<BracketContext<U8>*>(brptr);
            brcxt1->Update(a);
            a=brcxt1->getcxt();
            *(unsigned int *)*sp++ = (unsigned int)a; 
        } else if (*pc==COLUMN) {
            pc++;
            // our ColumnContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)*sp+4;
            ColumnContext *brcxt1;
            brcxt1=reinterpret_cast<ColumnContext*>(brptr);
            brcxt1->Update(a);
            a=brcxt1->context;
            *(unsigned int *)*sp++ = (unsigned int)a; 
        } else if (*pc==TABLE) {
            pc++;
            // our TableContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)*sp+4;
            TableContext *brcxt1;
            brcxt1=reinterpret_cast<TableContext*>(brptr);
            brcxt1->Update(a);
            a=brcxt1->context;
            *(unsigned int *)*sp++ = (unsigned int)a; 
        } else if (*pc==RECORD) {
            pc++;
            // our RecordContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)*sp+4;
            RecordContext *brcxt1;
            brcxt1=reinterpret_cast<RecordContext*>(brptr);
            brcxt1->Update(a);
            a=brcxt1->context;
            *(unsigned int *)*sp++ = (unsigned int)a; 
        } else if (*pc==DICTIONARY) {
            pc++;
            // our RecordContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)*sp+4;
            DictionaryContext *brcxt1;
            brcxt1=reinterpret_cast<DictionaryContext*>(brptr);
            //brcxt1->Update(a);
            a=brcxt1->context;
            *(unsigned int *)*sp++ = (unsigned int)a; 
        } else { 
            kprintf("unknown TYPA \n"); return -1; 
        }
    } else if (i == TYPR) {
        int context=*pc&255;
        int operation=*pc>>8;
        if (context==DISTANCE){
            pc++;
            // our BracketContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)a+4;
            BracketContext<U8> *brcxt1;
            brcxt1=reinterpret_cast<BracketContext<U8>*>(brptr);
            brcxt1->Reset();
        } else if (context==TABLE){
            pc++;
            // our TableContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)a+4;
            TableContext *brcxt1;
            brcxt1=reinterpret_cast<TableContext*>(brptr);
            brcxt1->Reset();
        } else if (context==RECORD){
            pc++;
            // our BracketContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)a+4;
            RecordContext *brcxt1;
            brcxt1=reinterpret_cast<RecordContext*>(brptr);
            brcxt1->Reset();
        } else if (context==DICTIONARY){
            pc++;
            int *vala=(int *)*pc; 
            int val=*vala; // our value in dict.
             //printf("%d\n",val);
            // our BracketContext is at sp+4, we store our result at sp
            U8 *brptr=(unsigned char *)vala+4;
            DictionaryContext *brcxt1;
            brcxt1=reinterpret_cast<DictionaryContext*>(brptr);
           
            a=brcxt1->get(val); 
            //printf("d %d\n",a);
            *vala=a;
            pc++;
            
        } else if (context==COLUMN){
            pc++;
        } else { kprintf("unknown TYPR %d\n",context); return -1; }
    }
    else if (i == EXIT) { /*printf("exit(%d) cycle = %d\n", *sp, cycle);*/ return *sp; }
    else if (i == READ) a = (int)readfile(this,(U8 *)sp[1], *sp); //pointer,lenght
    else if (i == WRTE) a = (int)writefile(this,(U8 *)sp[1], *sp); //pointer,lenght
    else { kprintf("unknown instruction = %d! cycle = %d\n", i, cycle); return -1; }
    // if (debug) printf("a=%d ",a);
  }

}

void printOC(int a){
     printf("   %.4s",  
        &"LEA ,IMM ,JMP ,JSR ,BZ  ,BNZ ,LBZ ,LBNZ,ENT ,ADJ ,LEV ,LI  ,LS  ,LC  ,SI  ,SS  ,SC  ,PSH ,"
         "OR  ,XOR ,AND ,EQ  ,NE  ,LT  ,GT  ,LE  ,GE  ,SHL ,SHR ,ADD ,SUB ,MUL ,DIV ,MOD ,BOUN,THIS,"
         "PRTF,VMI ,VMX ,READ,WRTE,EXIT,TYPA,TYPR,NOP ,VMXD"[a * 5]);
}

//jitadr
#define REGEAX 0
#define REGECX 1
#define jmov(a,c) ( *jitadr++ = 0xb8+c,  *(int*)jitadr =a, jitadr = jitadr+4) //mov eax,x or mov ecx,x c=0=eax c=1=ecx
#define calleax() ( *jitadr++ = 0xff, *jitadr++ = 0xd0) //call eax
#define pusheax() ( *jitadr++ = 0x50)// push eax
#define getatreg(a) (*jitadr++ = 0x8B, *jitadr++ = (0x00+a)) //mov eax,[eax/ecx] c=0=eax c=1=ecx
#define movatecx() (*jitadr++ = 0x89, *jitadr++ = 0x01) // mov [ecx], eax
#define enters() (*(int *)jitadr++ = 0x51,*(int *)jitadr++ = 0x52) // push ecx edx
#define leaves() (*(int *)jitadr++ = 0x5a,*(int *)jitadr++ = 0x59, *jitadr++ = 0xc3) //pop edx ecx, ret 

#define xoreax() (*jitadr++ = 0x31, *jitadr++ = 0xc0) //xor eax,eax   
#define push(a) ( *jitadr++ = 0x50+a)// push eax/ecx
#define pop(a) ( *jitadr++ = 0x58+a)// pop eax/ecx
#define jitadr je

int VM::dojit() {
    int u,i;
    // setup jit memory
    jitmem = (char*)mmap(0, poolsz, PROT_EXEC | PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (!jitmem) { dprintf("could not mmap(%d) jit executable memory\n", poolsz); return -1; }

    // first pass: emit native code
    pc = text + 1; 
    je = jitmem; 
    line = 0;
    //int lastcode=0;
    while (pc <= e) {
        //lastcode=i;
        i = *pc;
        //if (debug) {
        //dprintf("// %x: ",pc);
        /* printf("%x: ",int(pc-text));
    printOC(i);
    printf("\n");
    */
        dprintf("   %.4s",  
        &"LEA ,IMM ,JMP ,JSR ,BZ  ,BNZ ,LBZ ,LBNZ,ENT ,ADJ ,LEV ,LI  ,LS  ,LC  ,SI  ,SS  ,SC  ,PSH ,"
        "OR  ,XOR ,AND ,EQ  ,NE  ,LT  ,GT  ,LE  ,GE  ,SHL ,SHR ,ADD ,SUB ,MUL ,DIV ,MOD ,BOUN,THIS,"
        "PRTF,VMI ,VMX ,READ,WRTE,EXIT,TYPA,TYPR,NOP ,VMXD"[i * 5]);
        if (i < JMP) dprintf(" 0x%x\n",*(pc+1)); //? +1
        else if (i <= ADJ) dprintf(" 0x%x\n",(int *)*pc); else dprintf("\n");
        // }
        *pc++ = ((int)je << 8) | i; // for later relocation of JMP/JSR/BZ/BNZ

        if (i == LEA) {
            if (*(pc+1)==LI /*&& i>=0*/){
                i = 4 * *pc++; if (i < -128 || i > 127) { kprintf("jit: LEA out of bounds\n"); return -1; }
                dprintf("   LI\n");//dprintf("// %x:    LI\n",pc);
                *je++ =  0x8B, *je++ = 0x45,*je++ = i;dprintf("\tmov eax,DWORD PTR [ebp%s%d]\n",i>=0?"+":"",i);
                i = *pc;*pc++ = ((int)je << 8) | i; 
            } else {
                i = 4 * *pc++; if (i < -128 || i > 127) { kprintf("jit: LEA out of bounds\n"); return -1; }
                *(int*)je = 0x458d; je = je + 2; *je++ = i;
                dprintf("\tlea eax,[ebp%s%d]\n",i>=0?"+":"",i);
            }
        } else if (i == ENT) {
            i = 4 * *pc++; if (i < -128 || i > 127) { kprintf("jit: ENT out of bounds\n"); return -1; }
            if (*(pc)==LEV && i==0){
                *je++ = 0xc3; //ret
                dprintf("\tBlank proc\n\tret\n"); 
                pc++;
            }else{
                *(int *)je = 0xe58955; je = je + 3;
                dprintf("\tpush ebp\n\tmov ebp, esp\n",i);
                if (i > 0) { 
                    *(int *)je = 0xec83; je = je + 2; 
                    *(int*)je++ = i; 
                    dprintf("\tsub esp,BYTE %x\n",i); // fix me, zero mem
                }
                dprintf("\tpush esi push ecx push edx  push ebx push   edi\n"); 
                *(int *)je++ = 0x56;*(int *)je++ = 0x51;*(int *)je++ = 0x52;//*(int *)je++ =0x53;//*(int *)je++ =0x57;
            }
        } else if (i == IMM) { 
            if (*(pc+1)==LI){
                dprintf("   LI\n");
                *je++ = 0xa1; *(int *)je = i=*pc++; je = je + 4; dprintf("\tmov eax,DWORD [0x%x]\n",i);
                i = *pc;*pc++ = ((int)je << 8) | i; 
            } else if (*(pc+1)==LS) {
                dprintf("   LS\n");
                *je++ =0x0F,  *je++ =0xB7,  *je++ =0x05;
                *(int *)je = i=*pc++; je = je + 4; dprintf("\tmovzx eax,WORD [0x%x]\n",i);
                i = *pc;*pc++ = ((int)je << 8) | i; 
            } else if (*(pc+1)==LC){
                dprintf("   LC\n");
                *je++ =0x0F,  *je++ =0xB6,  *je++ =0x05;
                *(int *)je = i=*pc++; je = je + 4; dprintf("\tmovzx eax,BYTE [0x%x]\n",i);
                i = *pc;*pc++ = ((int)je << 8) | i; 
            } else {
                *je++ = 0xb8; *(int *)je = i=*pc++; je = je + 4; dprintf("\tmov eax,DWORD 0x%x\n",i);
            }
        } else if (i == NOP) { 
        } else if (i == ADJ) { 
            i = 4 * *pc++; 
            if (i) {
                *(int *)je = 0xc483; je = je + 2; *(int *)je = i; je++;  // add esp,BYTE (n * 4)
            }
        } else if (i == PSH && *(pc)==IMM && (*(pc+1)==4|| *(pc+1)==2 )&& *(pc+2)==MUL && *(pc+3)==ADD){    //array index*int
            *(int*)je = 0x59; je = je + 1; // pop ecx 
            if (*(pc+1)==2)
            *(int *)je =0x41048d;            // lea    eax,[ecx+eax*2] 8d 04 41
            else if (*(pc+1)==4) {
                i=*(pc+4);
                if (*(pc+4)==LI) {
                    *(int *)je = 0x81048b; //0x8B, 0x04, 0x81  mov    eax,DWORD PTR [ecx+eax*4] 
                } else {
                    *(int *)je = 0x81048d;           // lea    eax,[ecx+eax*4] 0x8D, 0x04, 0x81
                }
            }
            dprintf("array index* optimized\n");
            dprintf("   IMM 0x4\n");
            dprintf("   MUL\n");
            dprintf("   ADD\n");
            dprintf("\tpop ecx\n");
            dprintf("\tlea    eax,[ecx+eax*4]\n");
            i = *pc;//IMM
            *pc++ = ((int)je << 8) | i; pc++;
            i = *pc;//mul
            *pc++ = ((int)je << 8) | i;
            i = *pc;//add
            *pc++ = ((int)je << 8) | i;
            je = je + 3; 
            if (*(pc)==LI) {
                i = *pc;//LI
                *pc++ = ((int)je << 8) | i;
            }
        } else if (i == PSH && *pc==IMM && *(pc+2)==AND) {
            pc++;
            dprintf("   AND\n");
            *je++ =0x25;
            *(int *)je = i=*pc++; je = je + 4; dprintf("\tand eax, 0x%x\n",i);
            i = *pc;*pc++ = ((int)je << 8) | i;
        } else if (i == PSH && *pc==IMM && *(pc+2)==SHR) {
            pc++;
            dprintf("   SHR\n");
            *je++ =0xc1;*je++ =0xe8,*je++ =*pc++;//, 0x11, 0x11, 0xFF, 0xFF
            i = *pc;*pc++ = ((int)je << 8) | i;
        }
        else if (i == PSH && *pc==IMM && *(pc+2)==SHL) {
            pc++;
            dprintf("   SHL\n");
            *je++ =0xc1;*je++ =0xe0,*je++ =*pc++;
            i = *pc;*pc++ = ((int)je << 8) | i;
        }  else if (i == PSH && *pc==IMM && *(pc+2)==SUB) {
            pc++;
            dprintf("   SUB\n");
            *je++ =0x2d;
            *(int *)je = i=*pc++; je = je + 4; dprintf("\tsub eax, 0x%x\n",i);
            i = *pc;*pc++ = ((int)je << 8) | i;
        } else if (i == PSH && *pc==IMM && *(pc+2)==ADD) {
            pc++;
            dprintf("   ADD\n");
            *je++ =0x05;
            *(int *)je = i=*pc++; je = je + 4; dprintf("\tadd eax, 0x%x\n",i);
            i = *pc;*pc++ = ((int)je << 8) | i;
        } else if (i == PSH && *pc==IMM && *(pc+2)==MUL) {
            pc++;
            dprintf("   MUL\n");
            *je++ =0xB9;
            *(int *)je = i=*pc++; je = je + 4; dprintf("\tmov ecx, 0x%x\n",i);
            *je++ =0x0F, *je++ =0xAF, *je++ =0xC1;// imul eax,ecx
            i = *pc;*pc++ = ((int)je << 8) | i;
        } else if (i == PSH && *pc==IMM && *(pc+2)==DIV){
            pc++;
            dprintf("DIV optimized\n");
            i=*(pc); pc++; // imm value
            *je++ = 0xb9; *(int *)je =i; je = je + 4; dprintf("\tmov ecx,DWORD 0x%x\n",i);
            *je++=0x99;         dprintf("\tcdq\n");    //cdq
            *je++=0xF7; *je++=0xF9; dprintf("\tidiv eax, ecx\n");//idiv eax, ecx
            i = *pc;//Div
            *pc++ = ((int)je << 8) | i;
        }
        else if (i == PSH) {  *(int *)je++ = 0x50;dprintf("\tpush eax\n"); }
        else if (i == LEV ) {
            /*if (*(pc+1)==LEV) {
            i = *pc;//LEV
            *pc++ = ((int)je << 8) | i; 
            }else*/{
                //*(int *)je++ = 0x5f;
                //*(int *)je++ =0x5b;
                *(int *)je++ = 0x5a;*(int *)je++ = 0x59; *(int *)je++ = 0x5e; 
                dprintf("\tpop edi pop ebx pop edx pop ecx pop esi\n");
                *(int *)je = 0xc35dec89; je = je + 4; 
                dprintf("\tmov esp, ebp\n\tpop ebp\n\tret\n");
            }
        }
        else if ((i == LI ||i == LC ||i == LS )&& *(pc)==PSH && *(pc+1)==IMM && 
                (*(pc+3)==SHR || *(pc+3)==SHL || *(pc+3)==SUB|| *(pc+3)==ADD|| 
                    (*(pc+3)==MUL && *(pc+4)!=ADD) ||
                    *(pc+3)==OR || *(pc+3)==XOR || *(pc+3)==AND    )){
            
            if (i == LI )*(int *)je = 0x008b,     je = je + 2; //LI
            else if (i == LC )*(int *)je = 0x00b60f,     je = je + 3; //LC
            else if (i == LS )*(int *)je = 0x00B70F,     je = je + 3; //LS
            i = *pc;//PSH
            *pc++ = ((int)je << 8) | i;    
            i = *pc;//IMM
            *pc++ = ((int)je << 8) | i; 
            *je++ = 0xb9; *(int *)je = i=*pc++; je = je + 4; 
            dprintf("    PSH\n");
            dprintf("    IMM 0x%x\n",i);
            dprintf("    %s\n",*pc==SHR?"SHR":*pc==SHL?"SHL":*pc==SUB?"SUB":*pc==ADD?"ADD":*pc==MUL?"MUL":*pc==OR?"OR":*pc==XOR?"XOR":"AND");
            dprintf("\tmov eax,DWORD PTR [eax]\n");
            dprintf("\tmov ecx,DWORD 0x%x\n",i);
            dprintf("\t%s eax,ecx\n",*pc==SHR?"SHR":*pc==SHL?"SHL":*pc==SUB?"SUB":*pc==ADD?"ADD":*pc==MUL?"MUL":*pc==OR?"OR":*pc==XOR?"XOR":"AND");
            i = *pc;//XXX
            *pc++ = ((int)je << 8) | i; 
            
            *(int*)je = i==SHR?0xe8d3:i==SHL?0xe0d3:i==SUB?0xc829:i==ADD?0xc801:i==MUL?0xc1af0f:i==OR?0xc809:i==XOR?0xc831:0xc821; //AND
            je = je + 2;//dprintf("\tXXX  ;optimized\n");  
            if (i==MUL) je=je+1;
        }
        else if (i == LI)  { *(int *)je = 0x008b;     je = je + 2; dprintf("\tmov eax,DWORD PTR [eax]\n");} 
        else if (i == LC)  { *(int *)je = 0x00b60f;   je = je + 3; dprintf("\tmovzx eax,BYTE PTR [eax]\n"); } 
        else if (i == LS)  { *(int *)je = 0x00B70F;   je = je + 3; dprintf("\tmovzx eax,WORD PTR [eax]\n"); } 
        else if (i == SI)  { *(int *)je = 0x018959;   je = je + 3; dprintf("\tpop ecx\n\tmov DWORD PTR [ecx],eax    \n");}    
        else if (i == SC)  { *(int *)je = 0x018859;   je = je + 3; dprintf("\tpop ecx\n\tmov BYTE PTR [ecx],al    \n"); }
        else if (i == SS)  { *(int *)je = 0x01896659; je = je + 4; dprintf("\tpop ecx\n\tmov WORD PTR [ecx],ax    \n"); }
        else if (i == OR)  { *(int *)je = 0xc80959;   je = je + 3; dprintf("\tpop ecx\n\tor eax, ecx\n"); }
        else if (i == XOR) { *(int *)je = 0xc83159;   je = je + 3; dprintf("\tpop ecx\n\txor eax, ecx\n"); }
        else if (i == AND) { *(int *)je = 0xc82159;   je = je + 3; dprintf("\tpop ecx\n\tand eax, ecx\n"); }
        else if (EQ <= i && i <= GE) {
            if (i == LT && *(pc)==BZ){// logical brancj problem
                *je++ =  0x59, *je++ =  0x39,  *je++ = 0xC1;dprintf("\tpop ecx\n\tcmp ecx, eax");
                i = *pc;
                *pc++ = ((int)je << 8) | i;   ++pc; 
                *je++ =0x90;*je++ =0x90; 
                *je++ = 0x0f , *je++ =0x8d;//jge
                je = je + 4;
            }else if (i == GT && *(pc)==BZ){// logical brancj problem
                *je++ =  0x59, *je++ =  0x39,  *je++ = 0xC1;dprintf("\tpop ecx\n\tcmp ecx, eax");
                i = *pc;
                *pc++ = ((int)je << 8) | i;   ++pc; 
                *je++ =0x90;*je++ =0x90;
                *je++ = 0x0f , *je++ =0x8e;// jle
                je = je + 4;
            }else if (i == NE && *(pc)==BZ){// logical brancj problem
                *je++ =  0x59, *je++ =  0x39,  *je++ = 0xC1;dprintf("\tpop ecx\n\tcmp ecx, eax");
                i = *pc;
                *pc++ = ((int)je << 8) | i;   ++pc; 
                *je++ =0x90;*je++ =0x90;  *je++ = 0x0f , *je++ =0x84;// je
                je = je + 4;
            }else{
                
                *(int*)je=0x0fc13959; je+=4; 
                dprintf("\tpop ecx\n\tcmp ecx, eax");
                if      (i == EQ)  { *je++ = 0x94; dprintf("\n\tsete al");} // setne al
                else if (i == NE)  { *je++ = 0x95; dprintf("\n\tsetne al");} // setne al
                else if (i == LT)  { *je++ = 0x9c; dprintf("\n\tsetb al");} // setb al
                else if (i == GT)  { *je++ = 0x9f; dprintf("\n\tseta al");} // seta al
                else if (i == LE)  { *je++ = 0x9e; dprintf("\n\tsetbe al");} // setbe al
                else if (i == GE)  { *je++ = 0x9d; dprintf("\n\tsetae al");} // setae al
                *je++=0xc0;
                dprintf("\n\tmovsx  eax,al\n");
                *je++= 0x0F, *je++=0xBE,*je++= 0xC0;
            }
        }
        else if (i == SHL) { *(int*)je = 0xe0d39159; je = je + 4;dprintf("\tpop ecx\n\txchg ecx, eax\n\tshl eax, cl\n"); } // pop ecx; xchg ecx, eax; shl eax, cl
        else if (i == SHR) { *(int*)je = 0xe8d39159; je = je + 4;dprintf("\tpop ecx\n\txchg ecx, eax\n\tshr eax, cl\n");  } // pop ecx; xchg ecx, eax; shr eax, cl
        else if (i == ADD) { *(int*)je = 0xc80159;   je = je + 3;dprintf("\tpop ecx\n\tadd eax, ecx\n"); } // pop ecx; add eax, ecx
        else if (i == SUB) { *(int*)je = 0xc8299159; je = je + 4;dprintf("\tpop ecx\n\txchg ecx, eax\n\tsub eax, ecx\n");  } // pop ecx; xchg ecx,eax; sub eax,ecx
        else if (i == MUL) { *(int*)je = 0xc1af0f59; je = je + 4;dprintf("\tpop ecx\n\txchg ecx, eax\n\t imul eax, ecx\n");  } // pop ecx; imul eax,ecx
        else if (i == DIV) { 
            *je++=0x59; *(int*)je = 0xF9F79991; je = je + 4;dprintf("\txor edx,edx\tpop ecx\n\txchg ecx, eax\n\tdiv eax, ecx\n"); 
        } // pop ecx; xchg ecx,eax; idiv eax,ecx
        else if (i == MOD) { *(int*)je = 0x999159; je += 3; *(int *)je = 0x92f9f7; je += 3; dprintf("\txor edx,edx\n\tpop ecx\n\txchg ecx,eax\n\tdiv ecx\n\txchg   edx,eax" ); }
        else if (i == JMP) { ++pc; *je       = 0xe9;     je = je + 5; dprintf("\tjmp  %x\n", *(pc-1) ); } // jmp <off32>
        else if (i == JSR) { ++pc; *je       = 0xe8;     je = je + 5; dprintf("\tcall  %x\n", *(pc-1) ); } // call <off32>
        else if (i == BZ||i == LBZ)  { ++pc; *(int*)je = 0x840fc085; je = je + 8;dprintf("\ttest eax, eax\n\tjz  %x\n", *(pc-1) ); } // test %eax, %eax; jz <off32>
        else if (i == BNZ ||i == LBNZ) { ++pc; *(int*)je = 0x850fc085; je = je + 8;dprintf("\ttest eax, eax\n\tjnz  %x\n", *(pc-1)  );  } // test %eax, %eax; jnz <off32>
        else if (i == VTHIS) { 
            *je++ = 0x68; //push this 
            *(int*)je =i=(unsigned int)(size_t(this));
            je += 4; dprintf("\tpush %x    ;this\n",i); 
        } else if (i==VMXD) {
            i = *pc; //PSH
            *pc++ = ((int)je << 8) | VMXD;
            if ((i&255)==7){ //CM
                int cmIDX=i>>8;
                int cmCN=(int)&cmC[cmIDX].cn;
                int cmCXT=(int)&cmC[cmIDX].cxt[0];
                int cmCXMASK=(int)&cmC[cmIDX].cxtMask;
                // use this or one below
                /* jmov(cmC[cmIDX].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &cmC[cmIDX].set),REGEAX);
                calleax();
                */
                jmov(cmCN,REGECX);  // mov ecx,&cn
                getatreg(REGECX);   // mov eax,[ecx]
                *je++=0x8D, *je++=0x50, *je++=0x01;    // lea edx, [eax+1]
                *je++=0x89, *je++=0x11;                // mov [ecx], edx
                *je++=0x5A;                            // pop edx -  our value 
                jmov(cmCXT,REGECX); //// mov ecx,cxt
                *je++=0x89, *je++=0x14, *je++=0x81; //mov  [ecx+eax*4],edx
                //  get mask
                jmov(cmCXMASK,REGECX);  // mov  ecx,cmCXMASK 
                getatreg(REGECX);       // mov  edx,[ecx]
                *je++ = 0x8B, *je++ =0x51, *je++ =0x04; //mov    edx,DWORD PTR [ecx+0x4]
                // shld    edx, eax, 1
                // add     eax, eax
                *je++ =0x0F,  *je++ =0xA4,  *je++ =0xC2,  *je++ =0x01;
                *je++ =0x01,  *je++ =0xC0;
                // store cxtMask (64 bit)
                movatecx();                            // mov    [ecx],eax
                *je++ =0x89,*je++ = 0x41,*je++ = 0x04; // mov [ecx+0x4],eax
            } else if ((i&255)==8){ // MX
                int mxIDX=i>>8;
                int mxcxt=(int)&mxA[mxIDX].cxt;
                je--;// remove push eax
                jmov(mxcxt,REGECX);
                movatecx();
            } else if ((i&255)==2){ // APM1
                int mxIDX=i>>8;
                int mxcxt=(int)&apm1A[mxIDX].cxt;
                je--;// remove push eax
                jmov(mxcxt,REGECX);
                movatecx();
            } else if ((i&255)==14){ // APM2
                int mxIDX=i>>8;
                int mxcxt=(int)&apm2A[mxIDX].cx;
                je--;// remove push eax
                jmov(mxcxt,REGECX);
                movatecx();
            } else if ((i&255)==11){ // DHS
                int mxIDX=i>>8;
                jmov(dhsA[mxIDX].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &dhsA[mxIDX].set),REGEAX);
                calleax();
            } else if ((i&255)==3){ // DS
                int mxIDX=i>>8;
                jmov(dsA[mxIDX].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &dsA[mxIDX].set),REGEAX);
                calleax();
            } else if ((i&255)==12){ // SM
                int mxIDX=i>>8;
                jmov(smcA[mxIDX].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &smcA[mxIDX].set),REGEAX);
                calleax();
            } else if ((i&255)==1){ // SMC
                int mxIDX=i>>8;
                jmov(smA[mxIDX].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &smA[mxIDX].set),REGEAX);
                calleax();
            }     
        } else if (i==VMX || i==VMI || i==READ || i==WRTE || (i == BOUND)) {
            //do vmx(x,y) or vmi(x,y,z,w) or ...
            //e8 x
            //83 c4, xx <- esp
            if (i==VMX) {
                tmp=(int)setcomponent1; 
                dprintf("\tcall vmx\n"); 
            } else if (i==VMI) {
                tmp=(int)initcomponent1; 
                dprintf("\tcall vmi\n"); 
            } else if (i==READ) {
                tmp=(int)readfile1;
                dprintf("\tcall read\n");
            } else if (i==WRTE) {
                tmp=(int)writefile1;  
                dprintf("\tcall write\n");
            } else if (i == BOUND) { 
                tmp = (int)vmbound;
                dprintf("\tcall bounds\n");
            } else quit("Bad function.\n");
            if (*pc++ == ADJ) { i = *pc++; } else { kprintf("no ADJ after native proc!\n"); exit(2); }
            *je++ = 0xe8;
            *(int*)je = tmp - (int)(je + 4); je = je + 4; // <*tmp offset>;
            *je++ = 0x83;  *je++ = 0xc4; *je++ =  i << 2;  dprintf("\tsub  esp, 0x%x\n", i << 2 );  // movl $(4 * n), %ecx;
        } else if (i==TYPA) {
            if (*pc==DISTANCE) { // not identical between jit and none
                // last value in stack is our BracketContext address
                pop(REGECX);  // address
                push(REGECX);// save
                // ecx= address of brcxt+4
                // eax == new value
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                push(REGECX);//  this
                BracketContext<U8> *brcxt1;
                char *brptr; // dummy ptr
                brcxt1=reinterpret_cast<BracketContext<U8>*>(brptr); // any valid address
                push(REGEAX);// our new byte value
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Update),REGEAX);
                calleax(); 
                pop(REGECX);
                jmov((int)&brcxt1->context,REGEAX);// mov eax,brcxt1->context
                *je++ =0x01,*je++ = 0xC1; //add    ecx,eax
                getatreg(REGECX);
                pop(REGECX);// our dest variable
                movatecx(); // store
            } else if (*pc==COLUMN) {
                // last value in stack is our BracketContext address
                pop(REGECX);  // address
                push(REGECX);// save
                // ecx= address of brcxt+4
                // eax == new value
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                push(REGECX);//  this
                ColumnContext *colcxt1;
                char *colptr; // dummy ptr
                colcxt1=reinterpret_cast<ColumnContext*>(colptr); // any valid address
                push(REGEAX);// our new byte value
                jmov((int)reinterpret_cast<int>((void *) &colcxt1->Update),REGEAX);
                calleax(); 
                pop(REGECX);
                jmov((int)&colcxt1->context,REGEAX);// mov eax,brcxt1->context
                *je++ =0x01,*je++ = 0xC1; //add    ecx,eax
                getatreg(REGECX);
                pop(REGECX);// our dest variable
                movatecx(); // store
            } else if (*pc==TABLE) {
                // last value in stack is our BracketContext address
                pop(REGECX);  // address
                push(REGECX);// save
                // ecx= address of brcxt+4
                // eax == new value
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                push(REGECX);//  this
                TableContext *brcxt1;
                char *brptr; // dummy ptr
                brcxt1=reinterpret_cast<TableContext*>(brptr); // any valid address
                push(REGEAX);// our new byte value
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Update),REGEAX);
                calleax(); 
                pop(REGECX);
                jmov((int)&brcxt1->context,REGEAX);// mov eax,brcxt1->context
                *je++ =0x01,*je++ = 0xC1; //add    ecx,eax
                getatreg(REGECX);
                pop(REGECX);// our dest variable
                movatecx(); // store
            } else if (*pc==RECORD) {
                // last value in stack is our BracketContext address
                pop(REGECX);  // address
                push(REGECX);// save
                // ecx= address of brcxt+4
                // eax == new value
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                push(REGECX);//  this
                RecordContext *brcxt1;
                char *brptr; // dummy ptr
                brcxt1=reinterpret_cast<RecordContext*>(brptr); // any valid address
                push(REGEAX);// our new byte value
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Update),REGEAX);
                calleax(); 
                pop(REGECX);
                jmov((int)&brcxt1->context,REGEAX);// mov eax,brcxt1->context
                *je++ =0x01,*je++ = 0xC1; //add    ecx,eax
                getatreg(REGECX);
                pop(REGECX);// our dest variable
                movatecx(); // store
            } else { 
                kprintf("TYPA bad! %d\n",*pc); exit(2); 
            }
            *pc++ = ((int)je << 8) | i;// skip type selector
        } else if (i==TYPR) {
            int operation=*pc&255;
            int context=*pc>>8;
            if (operation==DISTANCE) {
                push(REGEAX);// //save
                pop(REGECX);
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                char *brptr; // dummy ptr
                BracketContext<U8> *brcxt1;
                brcxt1=reinterpret_cast<BracketContext<U8>*>(brptr);
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Reset),REGEAX);
                calleax();
            } else if (operation==COLUMN) {
                /*    push(REGEAX);// //save
            pop(REGECX);
            *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
            char *colptr; // dummy ptr, not used
            ColummnContext *colcxt1;
            brcxt1=reinterpret_cast<BracketContext<U8>*>(colptr);
            jmov((int)reinterpret_cast<int>((void *) &colcxt1->Reset),REGEAX);
            calleax();*/
            } else if (operation==TABLE) {
                push(REGEAX);  //save
                pop(REGECX);
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                char *brptr; // dummy ptr
                TableContext *brcxt1;
                brcxt1=reinterpret_cast<TableContext*>(brptr);
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Reset),REGEAX);
                calleax();
            } else if (operation==RECORD) {
                push(REGEAX);  //save
                pop(REGECX);
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                char *brptr; // dummy ptr
                RecordContext *brcxt1;
                brcxt1=reinterpret_cast<RecordContext*>(brptr);
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->Reset),REGEAX);
                calleax();
            } else if (operation==DICTIONARY) {
                jmov((int)*(pc+1),REGEAX);
                *pc++ = ((int)je << 8) | i;// skip type selector
                push(REGEAX);  //dest
                push(REGEAX);  //dest
                getatreg(REGEAX);
                pop(REGECX);
                *je++ = 0x83; *je++ = 0xC1; *je++ = 0x04; // add    ecx,0x4
                push(REGEAX);  //save
                char *brptr; // dummy ptr
                DictionaryContext *brcxt1;
                brcxt1=reinterpret_cast<DictionaryContext*>(brptr);
                jmov((int)reinterpret_cast<int>((void *) &brcxt1->get),REGEAX);
                calleax();
                pop(REGECX);// our dest variable
                movatecx(); // store
            } else {
                kprintf("TYPR bad! %d\n",*pc); exit(2); 
            }
            *pc++ = ((int)je << 8) | i; // skip type reset selector
        } else if (i == PRTF || i == EXIT) {
            if (i == PRTF) { tmp = (int)printf;  }
            else if (i == EXIT) { continue;  }
            else   { kprintf("Illegal call %d\n",i); exit(2);  }
            
            u=i;
            if (*pc++ == ADJ) { i = *pc++; } else { kprintf("no ADJ after native proc!\n"); exit(2); }
            *je++ = 0xb9; *(int*)je = i << 2; je += 4; dprintf("\tmov ecx, 0x%x\n", i << 2 );  // movl $(4 * n), %ecx;
            *(int*)je = 0xce29e689; je += 4; dprintf("\tmov esi,esp\n\tsub esi,ecx\n"); // mov %esp, %esi; sub %ecx, %esi;  -- %esi will adjust the stack
            *(int*)je = 0x8302e9c1; je += 4; dprintf("\tshr ecx,2\n");// shr $2, %ecx; and                -- alignment of %esp for OS X
            *(int*)je = 0x895af0e6; je += 4; // $0xfffffff0, %esi; pop %edx; mov..
            *(int*)je = 0xe2fc8e54; je += 4; // ..%edx, -4(%esi,%ecx,4); loop..  -- reversing args order

            *(int*)je = 0xe8f487f9; je += 4; // ..<'pop' offset>; xchg %esi, %esp; call    -- saving old stack in %esi
            dprintf("\tand esi,0xfffffff0\n"); 
            dprintf("\tpop edx\n\tmov DWORD PTR [esi+ecx*4-0x4],edx\n"); 
            dprintf("\tloop 0x00000006\n\tcall "); 
            if (u == PRTF) {  dprintf("printf"); }
            *(int*)je = tmp - (int)(je + 4); je = je + 4; // <*tmp offset>;
            dprintf(" //  %x\n", *(int*)(je-4) );  //print running memory address
            dprintf("// %x: ",pc-2); 
            dprintf("ADJ 0x%x\n",i<<2);
            *(int*)je = 0xf487; je += 2;         // xchg %esi, %esp  -- ADJ, back to old stack without arguments
            dprintf("\txchg esp,esi\n");
        }
        else { kprintf("code generation failed for %d!\n", i); return -1; }
    }
    //int codesize=je -jitmem;

    // second pass, relocation
    pc = text + 1;
    int prev=0;
    while (pc <= e) {
        prev=i;
        i = *pc & 0xff;
        je = (char*)(((unsigned)*pc++ >> 8) | ((unsigned)jitmem & 0xff000000)); // MSB is restored from jitmem
        if (i == JSR || i == JMP || i == BZ || i == BNZ|| i == LBZ || i == LBNZ) {
            tmp = (*(unsigned*)(*pc++) >> 8) | ((unsigned)jitmem & 0xff000000); // extract address
            if      (i == JSR || i == JMP) { je += 1; *(int*)je = tmp - (int)(je + 4); }
            //else if (i == BZ && (prev==LT || prev==GT) ) { je += 2; *(int*)je = tmp - (int)(je + 4); }
            else if (i == BZ  || i == BNZ  ||i == LBZ  || i == LBNZ) { je += 4; *(int*)je = tmp - (int)(je + 4); }
        }
        else if (i < LEV) { ++pc; }
        else if (i>VMXD/*TYPR*/) {kprintf("code generation failed. relocation error\n"); return -1; }
    }
    /*printf(" Code size %d \n",codesize);
    FILE *jitfile=fopen("jit1.bin","wb+");
    fwrite (jitmem, 1, codesize,jitfile);
    fclose(jitfile);*/
    jitReadonly(jitmem,poolsz); //make jit read/execute only
    return 0;
}

int VM::decode(int info, int len) {
    if (doJIT){
        int (*jitmain)(int, int); // pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)(int, int) >(*(unsigned*)(iddecode->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        return jitmain(len, info);
    } else {
        // setup stack
       data =data0;
       bp=sp = (int *)((int)sp0 + poolsz);
       *--sp = EXIT; // call exit if main returns
       *--sp = PSH; t = sp;
       *--sp = info;
       *--sp = len; 
       *--sp = (int)t;
       return dovm((int *)iddecode->Val);
    }
}

int VM::encode(int info, int len) {
    if (doJIT){
        int (*jitmain)(int, int); // pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)(int, int) >(*(unsigned*)(idencode->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        return jitmain(len,info);
   } else {
        // setup stack
        data =data0;
        bp=sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; t = sp;
        *--sp = info;
        *--sp = len; 
        *--sp = (int)t;
        return dovm((int *)idencode->Val);
    }
}

int VM::detect(int c4, int pos) {
    if (doJIT){
        int (*jitmain)(int, int); // pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)(int, int) >(*(unsigned*)(iddetect->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        return jitmain(pos,c4);
   } else {
        // setup stack
        data =data0;
        bp=sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; t = sp;
        *--sp = c4;
        *--sp = pos; 
        *--sp = (int)t;
        return dovm((int *)iddetect->Val);
    }
}

int VM::block(int info1, int info2) {
    if (doJIT){
        int (*jitmain)(int, int); // pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)(int, int) >(*(unsigned*)(idp->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        return jitmain(info2,info1);
    }else{
        // setup stack
        data =data0;
        bp=sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; t = sp;
        *--sp = info1;
        *--sp = info2; 
        *--sp = (int)t;
        return dovm((int *)idp->Val);
    }
}

#define jitadr jep
// do byteupdate
void VM::getBytePredictionJIT() {
    jepbyte = (char*)mmap(0, poolsz, PROT_EXEC | PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (!jepbyte) { printf("could not mmap(%d) jepbyte executable memory\n", poolsz); quit(""); }
    char *jep;
    jep=jepbyte;
    //printf("Init JIT byte update\n");
    enters();
    if(x.cInputs>=0){
        // Do prediction into mixers
        for (int j=0; j<x.cInputs+1; j++) {  
            for (int i=0; i<totalc; i++){
                int inputIndex=mcomp[i].mixerIndex;    // index  
                int component=mcomp[i].component;
                if (j==inputIndex && (mcomp[i].prIndex)==0 && (component==vmCM)) {
                    int componentIndex=mcomp[i].compIndex;    // component index
                    //  printcomponent(component);printf("(%d).mix\n",componentIndex);
                    switch (component) { // select component and mix
                    /*case vmRCM:{
                        //rcmA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(rcmA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &rcmA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }*/
                    case  vmCM: {
                        //cmC[componentIndex].mix(inputIndex,x);
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov(cmC[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &cmC[componentIndex].mix),REGEAX);
                        calleax();
                        //printf("%d ",componentIndex);
                        break;
                    }
                    default:
                        quit("Bad component in VM:getBytePredictionJIT mix.\n");
                        break;
                    }
                }
            }
        }
    }
    leaves();
    //printf("\n");
    // JIT done
   // printf("Get Byte prediction JIT done. Code size %d bytes.\n",int(jep-jepbyte));
}

#define jitadr jep
// generate update
void VM::getPredictionJIT() {
    jitmempr = (char*)mmap(0, poolsz, PROT_EXEC | PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (!jitmempr) { printf("could not mmap(%d) jitmempr executable memory\n", poolsz); quit(""); }
    char *jep;
    jep=jitmempr;
    enters();
    if(x.cInputs>=0){
        // Do prediction into mixers
        for (int j=0; j<x.cInputs+1; j++) {  
            for (int i=0; i<totalc; i++){
                int inputIndex=mcomp[i].mixerIndex;    // index  
                int component=mcomp[i].component;
                if (j==inputIndex && (mcomp[i].prIndex)==0 && (component==vmST|| component==vmSMC || 
                                                               component==vmDHS || component==vmRCM ||
                                                               component==vmSCM || component==vmCM && isDualMode==false ||
                                                               component==vmSM ||component==vmSK )) {
                    int componentIndex=mcomp[i].compIndex;    // component index
                    //  printcomponent(component);printf("(%d).mix\n",componentIndex);
                    switch (component) { // select component and mix
                    case vmSMC: {
                        //x.mxInputs[inputIndex].add(smA[componentIndex].p()));
                        jmov(smA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &smA[componentIndex].p),REGEAX);
                        calleax();
                        
                        pusheax();
                        jmov(x.mxInputs[inputIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &x.mxInputs[inputIndex].add),REGEAX);
                        calleax();
                        break;
                    }
                    case vmRCM:{
                        //rcmA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(rcmA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &rcmA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case vmSCM: {
                        //scmA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(scmA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &scmA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case  vmCM: {
                        //cmC[componentIndex].mix(inputIndex,x);
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov(cmC[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &cmC[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case  vmST: {
                        //x.mxInputs[inputIndex].add(stA[componentIndex].pr);
                        jmov((int)&stA[componentIndex].pr,REGEAX);
                        getatreg(REGEAX); // get value from stA[componentIndex].pr
                        pusheax();
                        jmov(x.mxInputs[inputIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &x.mxInputs[inputIndex].add),REGEAX);
                        calleax();
                        break;
                    }
                    case vmSM: {
                        //smcA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov(smcA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &smcA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case vmSK:{
                        //skA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(skA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &skA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case vmDHS:{
                        //dhsA[componentIndex].mix(x,inputIndex);
                        jmov(inputIndex,REGEAX);
                        pusheax();
                        jmov((int)&x,REGEAX);
                        pusheax();
                        jmov(dhsA[componentIndex].This(),REGECX) ;
                        jmov((int)reinterpret_cast<int>((void *) &dhsA[componentIndex].mix),REGEAX);
                        calleax();
                        break;
                    }
                    case vmERR: 
                    case vmBYT: 
                        break;
                    default:
                        quit("Bad component in VM:getPredictionJIT mix.\n");
                        break;
                    }
                }
            }
        }
    }
    // Do global predictions
    for (int i=0; i<totalc; i++){
        int prindex=mcomp[i].prIndex;
        int index=mcomp[i].compIndex;
        // individual components
        if ((prindex>0 || mcomp[i].component==vmMM )){
          //if (compnr!=vmMM )printcomponent(compnr),printf("%d ",i);
         // printcomponent(compnr);printf("(%d)\n",index);
            prindex--;
            switch (mcomp[i].component) {
            case vmSMC: {
                //prSize[prindex]=smA[index].pr;
                jmov((int)&smA[index].pr,REGEAX);
                getatreg(REGEAX); // get value from smA[index].pr
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }  
            case vmAPM1: {
                //prSize[prindex]=apm1A[index].p(prSize[apm1A[index].p1],x.y);
                jmov((int)&x.y,REGEAX);
                getatreg(REGEAX); // get value from x.y
                pusheax();
                jmov((int)&prSize[apm1A[index].p1],REGEAX);
                getatreg(REGEAX); // get value from prSize[apm1A[index].p1
                pusheax();
                jmov(apm1A[index].This(),REGECX) ;
                jmov((int)reinterpret_cast<int>((void *) &apm1A[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmAPM2: {
                //prSize[prindex]=apm2A[index].p(prSize[apm2A[index].p1],x.y);
                jmov((int)&x.y,REGEAX);
                getatreg(REGEAX); // get value from x.y
                pusheax();
                jmov((int)&prSize[apm2A[index].p1],REGEAX);
                getatreg(REGEAX); // get value from prSize[apm2A[index].p1
                pusheax();
                jmov(apm2A[index].This(),REGECX) ;
                jmov((int)reinterpret_cast<int>((void *) &apm2A[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmUAS: {
                //prSize[prindex]=uasA[index].p(x.y);
                jmov((int)&x.y,REGEAX);
                getatreg(REGEAX); // get value from x.y
                pusheax();
                jmov(uasA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &uasA[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmDS: {
                // prSize[prindex]=dsA[index].p();
                jmov(dsA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &dsA[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
             case vmDHS: {
                //prSize[prindex]=dhsA[index].p();
                jmov(dhsA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &dhsA[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmAVG: {
                //prSize[prindex]=avA[index].average(&prSize[0]);
                jmov((int)&prSize[0],REGEAX);
                pusheax();
                jmov(avA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &avA[index].average),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmLMX: {
                //prSize[prindex]=lmxA[index].average(&prSize[0]);
                jmov((int)&prSize[0],REGEAX);
                pusheax();
                jmov(lmxA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &lmxA[index].average),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmMM: {
                int const mixnr=mcomp[i].mixerIndex;
                // mmA[index].p=prSize[ mmA[index].i1];
                // x.mxInputs[mixnr].add(mmA[index].pr());
                if (mmA[index].c1==0 ||  mmA[index].c1==1) { // fixme
                    // m.add(stretch(p)>>1)
                    jmov((int)&prSize[mmA[index].i1],REGEAX);
                    jmov((int)&strt[0],REGECX); // ecx=str address
                    getatreg(REGEAX); // eax=value of pr
                    *jep++ =0x0f,*jep++ = 0xbf,*jep++ = 0x04,*jep++ = 0x41; //movsx  eax,WORD PTR [ecx+eax*2]
                    if (mmA[index].c1==1) {
                        *jep++ = 0xD1, *jep++ =0xF8 ; //sar    eax,1 // m.add(stretch(p)>>1)
                    }
                } else {
                    jmov((int)&prSize[0],REGEAX);
                    pusheax();
                    jmov(mmA[index].This(),REGECX);
                    jmov((int)reinterpret_cast<int>((void *) &mmA[index].pr),REGEAX);
                    calleax();
                }
                pusheax();
                jmov(x.mxInputs[mixnr].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &x.mxInputs[mixnr].add),REGEAX);
                calleax();
                break; 
            }
            case vmMX: {
                //prSize[prindex]=mxA[index].p();
                jmov(mxA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &mxA[index].p),REGEAX);
                calleax();
                jmov((int)&prSize[prindex],REGECX);
                movatecx();
                break;
            }
            case vmST: {
                 //prSize[prindex]=stA[index].pr1;
                 jmov((int)&stA[index].pr1,REGEAX);
                 getatreg(REGEAX); // get value from stA[index].pr1
                 jmov((int)&prSize[prindex],REGECX);
                 movatecx();
                 break;
            }
            case vmSK: {
                // prSize[prindex]=skA[index].p();
                jmov(skA[index].This(),REGECX);
                jmov((int)reinterpret_cast<int>((void *) &skA[index].p),REGEAX);
                calleax();
                jmov((int)&this->prSize[prindex],REGECX);
                movatecx();
                break; 
            }
            default:{
                quit("Bad component in VM:getPredictionJIT predict.\n");
                break;}
            }
        }
    }  
    leaves();
    //printf("Get prediction JIT done. Code size %d bytes.\n",int(jep-jitmempr));
    jitmemup = (char*)mmap(0, poolsz, PROT_EXEC | PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (!jitmemup) { printf("could not mmap(%d) jitmempr executable memory\n", poolsz); quit(""); }

    jep=jitmemup;
    enters();
    if(x.cInputs>=0){
        // Do mixer updates if any
        for (int i=0; i<totalc; i++){
            if (mcomp[i].prIndex>0 && mcomp[i].component==vmMX){
                //mxA[mcomp[i].compIndex].update(x.y);
                jmov((int)&x.y,REGEAX);
                getatreg(REGEAX); // get value from x.y
                pusheax();
                jmov(mxA[mcomp[i].compIndex].This(),REGECX) ;
                jmov((int)reinterpret_cast<int>((void *) &mxA[mcomp[i].compIndex].update),REGEAX);
                calleax();
            }
        }
        xoreax();
        for (int j=0; j<x.cInputs+1; j++) {
           // x.mxInputs[j].ncount=0; //reset
            jmov((int)&x.mxInputs[j].ncount,REGECX);
            movatecx();
        }
    }
    leaves();
    // JIT done
    //printf("Get update JIT done. Code size %d bytes.\n",int(jep-jitmemup));
}

int VM::getPrediction() {
    if (doJIT){
        void (*func_ptr)(void) = (void (*)(void))jitmempr;
        func_ptr();
        return prSize[totalPR]; //final prediction
    } else {
    // Do prediction into mixers
    if(x.cInputs>=0){     
        for (int j=0; j<x.cInputs+1; j++) {  
            for (int i=0; i<totalc; i++){
                int inputIndex=mcomp[i].mixerIndex;    // index  
                int component=mcomp[i].component;
                if (j==inputIndex && (mcomp[i].prIndex)==0 && (component==vmST|| component==vmSMC || component==vmDHS || component==vmRCM ||component==vmSCM || component==vmCM && isDualMode==false ||component==vmSM ||component==vmSK )) {
                    int componentIndex=mcomp[i].compIndex;    // component index
                    //  printcomponent(component);printf("(%d).mix\n",componentIndex);
                    switch (component) { // select component and mix
                    case vmSMC: {
                        //if (doPRT==true) printf("x.mxInputs[%d].add(stretch(smA[%d].pr));\n",inputIndex,componentIndex);
                        x.mxInputs[inputIndex].add(stretch(smA[componentIndex].pr));
                        break;
                    }
                    case vmRCM:{
                        //if (doPRT==true) printf("rcmA[%d].mix(x,%d);\n",componentIndex,inputIndex);
                        rcmA[componentIndex].mix(x,inputIndex);

                        break;
                    }
                    case vmSCM: {
                        //if (doPRT==true) printf("scmA[%d].mix(x,%d);\n",componentIndex,inputIndex);
                        scmA[componentIndex].mix(x,inputIndex);
                        break;
                    }
                    case  vmCM: {
                        //if (doPRT==true) printf("cmC[%d]->mix(%d);\n",componentIndex,inputIndex);
                        cmC[componentIndex].mix(inputIndex,x);
                        break;
                    }
                    case  vmST: {
                        //if (doPRT==true) printf("x.mxInputs[%d].add(stA[%d].pr);\n",inputIndex,componentIndex);
                        x.mxInputs[inputIndex].add(stA[componentIndex].pr);
                        break;
                    }
                    case vmSM: {
                        //i/f (doPRT==true) printf("smcA[%d].mix(x,%d);\n",componentIndex,inputIndex);
                        smcA[componentIndex].mix(inputIndex);
                        break;
                    }
                    case vmSK:{
                        //if (doPRT==true) printf("skA[%d].mix(x,%d);\n",componentIndex,inputIndex);
                        skA[componentIndex].mix(x,inputIndex);
                        break;
                    }
                    case vmDHS:{
                        //if (doPRT==true) printf("skDHS[%d].mix(x,%d);\n",componentIndex,inputIndex);
                        dhsA[componentIndex].mix(x,inputIndex);
                        break;
                    }
                    case vmERR: 
                        break;
                    case vmBYT: 
                        break;
                    default:
                        quit("VM mxp error\n");
                        break;
                    }
                }
            }
        }
    }
    // Do global predictions
    for (int i=0; i<totalc; i++){
        int prindex=mcomp[i].prIndex;
        int index=mcomp[i].compIndex;
        // individual components
        if ((prindex>0 || mcomp[i].component==vmMM )){
          //if (compnr!=vmMM )printcomponent(compnr),printf("%d ",i);
         // printcomponent(compnr);printf("(%d)\n",index);
            prindex--;
            switch (mcomp[i].component) {
            case vmSMC: {//if (doPRT==true) printf("prSize[%d]=smA[%d].pr;\n",prindex,index);
                prSize[prindex]=smA[index].pr;
                break;
            }  
            case vmAPM1: {//if (doPRT==true) printf("prSize[%d]=apm1A[%d].p(prSize[%d],x.y);\n",prindex,index,apm1A[index].p1);
                prSize[prindex]=apm1A[index].p(prSize[apm1A[index].p1],x.y);
                break;
            }
            case vmAPM2: {//if (doPRT==true) printf("prSize[%d]=apm2A[%d].p(prSize[%d],x.y);\n",prindex,index,apm2A[index].p1);
                prSize[prindex]=apm2A[index].p(prSize[apm2A[index].p1],x.y);
                break;
            }
            case vmUAS: {//if (doPRT==true) printf("prSize[%d]=uasA[%d].p(x.y);\n",prindex,index);
                prSize[prindex]=uasA[index].p(x.y);
                break;
            }
            case vmDS: {//if (doPRT==true) printf("prSize[%d]=dsA[%d].p();\n",prindex,index);
                prSize[prindex]=dsA[index].p();
                break;
            }
             case vmDHS: {//if (doPRT==true) printf("prSize[%d]=dhsA[%d].p();\n",prindex,index);
                prSize[prindex]=dhsA[index].p();
                break;
            }
            case vmAVG: {//if (doPRT==true) printf("prSize[%d]=avA[%d].average(&prSize[0]);\n",prindex,index);
                prSize[prindex]=avA[index].average(&prSize[0]);
                break;
            }
            case vmLMX: {//if (doPRT==true) printf("prSize[%d]=lmxA[%d].average(&prSize[0]);\n",prindex,index);
                prSize[prindex]=lmxA[index].average(&prSize[0]);
                break;
            }
            case vmMM: {
                int const mixnr=mcomp[i].mixerIndex;
                /*if (doPRT==true) {
                    if (mmA[index].c1==1) printf("x.mxInputs[%d].add(stretch(prSize[%d])>>1  );\n",mixnr,mmA[index].i1);
                    else if (mmA[index].c1>1) printf("x.mxInputs[%d].add( ((prSize[%d]-2048)>>%d));\n",mixnr,mmA[index].i1,mmA[index].c1);
                    else printf("x.mxInputs[%d].add(stretch(prSize[%d])  );\n",mixnr,mmA[index].i1); 
                }*/
                //mmA[index].p=prSize[ mmA[index].i1];
                //x.mxInputs[mixnr].add(mmA[index].pr());
                x.mxInputs[mixnr].add(mmA[index].pr(&prSize[0]));
                
                break; 
            }
            case vmMX: {//if (doPRT==true) printf("prSize[%d]=mxA[%d].p();\n",prindex,index);
                prSize[prindex]=mxA[index].p();
                break;
            }
            case vmST: {//if (doPRT==true) printf("prSize[%d]=stA[%d].pr1;\n",prindex,index);
                 prSize[prindex]=stA[index].pr1;
                 break;
            }
            case vmSK: {//if (doPRT==true) printf("prSize[%d]=skA[%d].p();\n",prindex,index);
                 prSize[prindex]=skA[index].p();
                 break; 
            }
            default: {
                quit("VM vmi error\n");
                break;}
            }
            //if (compnr!=vmMM )printf("%d ",prSize[prindex-1]);
        }
    }
   // if (compnr!=vmMM )
   //printf("  \n" );
   /* if (doPRT==true) {
        printf("p=prSize[%d];\n",totalPR ); 
        printf("// end of getPrediction\n"); 
        doPRT=false; // disable until next stream
    }*/
    return prSize[totalPR]; //final prediction
    }
}

void VM::updateComponents1() {
        void (*func_ptr)(void) = (void (*)(void))jitmemup;
        func_ptr();
}
void VM::updateComponents() {

    if(x.cInputs>=0){
        for (int i=0; i<totalc; i++){
            // individual components
            if (mcomp[i].prIndex>0 && mcomp[i].component==vmMX){
                mxA[mcomp[i].compIndex].update(x.y);
                //if (doPRT==true) printf(" mxA[%d].update(x.y); \n",mcomp[i].compIndex ); 
            }
        }
        for (int j=0; j<x.cInputs+1; j++) {
            x.mxInputs[j].ncount=0; //reset
           // if (doPRT==true) printf("x.mxInputs[%d].ncount=0;  \n",j ); 
        }
       // if (doPRT==true)printf("// end of updateComponents\n"); 
    }

}

int VM::doupdate1(int y, int c0, int bpos, U32 c4) {
    if (isDualMode==true){
        if (bpos==0){
            data =data0;
            bp=sp = (int *)((int)sp0 + poolsz);
            *--sp = EXIT; // call exit if main returns
            *--sp = PSH; t = sp;
            *--sp = c4;
            *--sp = (int)t;
            dovm((int *)idbyteupdate->Val);
        }
        // Do prediction into mixers
        if(x.cInputs>=0){     
            for (int j=0; j<x.cInputs+1; j++) {  
                for (int i=0; i<totalc; i++){
                    int inputIndex=mcomp[i].mixerIndex;    // index  
                    int component=mcomp[i].component;
                    if (j==inputIndex && (mcomp[i].prIndex)==0 && (component==vmCM)) {
                        int componentIndex=mcomp[i].compIndex;    // component index
                        //  printcomponent(component);printf("(%d).mix\n",componentIndex);
                        switch (component) { // select component and mix
                        case  vmCM: {
                            //if (doPRT==true) printf("cmC[%d]->mix(%d);\n",componentIndex,inputIndex);
                            cmC[componentIndex].mix(inputIndex,x);
                            break;
                        }
                        default:
                            quit("VM mxp error\n");
                            break;
                        }
                    }
                }
            }
        }
        isBitUpdate=true;
        data =data0;
        bp=sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; t = sp;
        *--sp = y;
        *--sp = c0;
        *--sp = bpos;
        *--sp = c4;
        *--sp = (int)t;
        dovm((int *)idbitupdate->Val);
        isBitUpdate=false;
    } else {
        data =data0;
        bp=sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; t = sp;
        *--sp = y;
        *--sp = c0;
        *--sp = bpos;
        *--sp = c4;
        *--sp = (int)t;
        dovm((int *)idupdate->Val);
    }
    return 0;
}

int VM::doupdate2(int y, int c0, int bpos, U32 c4) {
    if (isDualMode==true){
        // First: update byte components and collect info
        // mix bytecomponents
        // then update bitcomponents
        if (bpos==0){
            int (*jitmain)(U32);
            jitmain = reinterpret_cast< int(*)(U32) >(*(unsigned*)( idbyteupdate->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
            jitmain(c4);            
        }
        // mix bytecomponents
        void (*func_ptr)(void) = (void (*)(void))jepbyte;
        func_ptr();
        isBitUpdate=true;
        // update bit
        int (*jitmain)(U32, int, int, int);
        jitmain = reinterpret_cast< int(*)(U32, int, int, int) >(*(unsigned*)( idbitupdate->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        jitmain(c4,bpos,c0,y);    
        isBitUpdate=false;
    } else {
        int (*jitmain)(U32, int, int, int); // pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)(U32, int, int, int) >(*(unsigned*)( idupdate->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
        jitmain(c4,bpos,c0,y);
    }
    return 0;
}

static char cVar[] = "char else enum if int short return for sizeof distance column table record dictionary while printf vmi vmx read write exit void block update main detect decode encode byteupdate bitupdate";
static char cEnums[] = "enum {SMC=1,APM1,DS,AVG,SCM,RCM,CM,MX,ST,MM,DHS,SM,SK,APM2,ERR,UAS,LMX,STA,BYT};enum {NONE=0,START,INFO,END,RESET=0xfffffffe,REQUEST=0xffffffff};enum {false=0,true=1};";
  
int VM::initvm() {
    int i;
    idCount=0;
  poolsz = 2024*1024; // arbitrary size
  alloc(sym,1024);
  if (!(text = le = e = (int *)malloc(poolsz))) { kprintf("could not malloc(%d) text area\n", poolsz); return -1; }
  if (!(sp =sp0= (int *)malloc(poolsz))) { kprintf("could not malloc(%d) stack area\n", poolsz); return -1; }
  data0 =data;
  memset(e,    0, poolsz);
  memset(sp, 0, poolsz);
  // parse main types and commands
  p = cVar;
  i = Char; while (i <= While) { next(); id->Tk = i++; } // add keywords to symbol table
  i = PRTF; while (i <= EXIT) { next(); id->Class = Sys; id->Type = iINT; id->Val = i++; } // add library to symbol table
  next(); id->Tk = Char; // handle void type
  next(); idp = id;      // keep track of block
  next(); idupdate = id; // keep track of updater
  next(); idmain = id;   // keep track of main
  next(); iddetect = id;
  next(); iddecode = id;
  next(); idencode = id;
  next(); idbyteupdate = id;
  next(); idbitupdate = id;

  // parse component enum, detect enum and bool enum
  // user can redefine them
  p = cEnums;
  next();
  while (tk == Enum) {
  
  if (tk == Enum) {
      next();
      if (tk != '{') next();
      if (tk == '{') {
        next();
        i = 0;
        while (tk != '}') {
          if (tk != Id) { kprintf("%d: bad enum identifier %d\n", line, tk); return -1; }
          next();
          if (tk == Assign) {
            next();
            if (tk != Num) { kprintf("%d: bad enum initializer\n", line); return -1; }
            i = ival;
            next();
          }
          id->Class = Num; id->Type = iINT; id->Val = i++;
          if (tk == Comma) next();
        }
        next();
      }
    }
    next();
    }
  // parse our user model
  idCount=0;
  p=mod.data; 
  line = 1;
  next();
  while (tk) {
    bt = iINT; // basetype
    if (tk == Int) next();
    else if (tk == Char) { next(); bt = rCHAR; }
    else if (tk == Short) { next(); bt = sSHORT; }
    else if (tk == Distance) { 
        // implement type:
        //     distance true/false array(char/short) brcxt;
        //
        // brcxt=brcxt+c;  // add new
        // brcxt=0; //reset
        // a=brcxt; // get current bracket and distance
        bool brB=false;
        int bType=0; // bracket type 0=char, 1=short
        int bSize=0; // 2,4,6,...
        U8 *bAddr8=0; // address
        //U16 *bAddr16=0; // address
        //printf("Distance ");
        next();
        if (id->Tk!=Id || id->Class!=Num) { kprintf("Error: bad distance boolean type at line %d\n", line); return -1; }
        i=id->Val; // boolean value
        if (i!=0  && i!=1) { kprintf("Error: bad distance boolean value at line %d\n", line); return -1; }
        if (i) brB=true;
        //if (i) printf("true ");  else printf("false ");
        next();
        if (id->Class==Glo && id->Type<(iINT+PTR)) {
            if (id->UBound&1) { kprintf("Error: bad distance element size at line %d\n", line); return -1; }
            bSize=id->UBound;
            
            int a=id->Type-PTR;
            //if (a==sSHORT) bType=1,bAddr16=(U16 *)id->Val+4; //data after pointerprintf("short ");
            //else 
            if (a==rCHAR) bAddr8=(U8 *)id->Val+4; else { kprintf("Error: bad distance char type %d\n", line); return -1; }
        }else { kprintf("Error: bad distance element type at line %d\n", line); return -1; }
        next();
        if (tk!=Id) { kprintf("Error: distance id  not declared or has bad type %d\n", line); return -1; }
        /*int j=0;
        char a;
        while (a=*(id->Name+j), a >='a' && a <='z' || a >='A' && a <='Z'|| a =='_'){
           printf("%c",*(id->Name+j) ) ;
           j++;
        }*/
        id->Class = Glo; // add here our new type
        id->HClass =Distance;
        id->Type = dDST;
        id->Val =(int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
        if (bType==0){
            BracketContext<U8> brcxt;
            brcxt.Init(bAddr8,bSize,brB);
            memcpy(data,&brcxt,sizeof(brcxt));
            data = data + sizeof(brcxt);
            
        } else {
             kprintf("Error: bad distance char type %d\n", line); return -1;
        }
        next();
        continue;
    }
    else if (tk == Dictionary) { 
        // implement type:
        //     dictionary "filename" dict;
        //
        // a=dict; // count of words
        //printf("Dictionary ");
        char *bAddr8=0; // address
        next();
        if (tk!='"') { kprintf("Error: bad Dictionary string %d\n", line); return -1; }
        expr(Assign);
        if (strlen((char*)*e)==0) { kprintf("Error: bad Dictionary string %d\n", line); return -1; }
       // printf("%s\n",*e);
        char *file=(char*)*e;
        e=e-2; // remove IMM val
        //
        
        if (id->Class==Glo && id->Type<(iINT+PTR)) {
            if (id->UBound!=64) { kprintf("Error: bad Dictionary element size at line %d\n", line); return -1; }
           // bSize=id->UBound;
            int a=id->Type-PTR;
            if (a==rCHAR) bAddr8=(char *)id->Val+4; else { kprintf("Error: bad Dictionary char type %d\n", line); return -1; }
            //printf("Dictionary input/output string size %d\n",id->UBound);
             //printf("%s\n",bAddr8);
        }else { kprintf("Error: bad Dictionary element type at line %d\n", line); return -1; }
        next();
        if (tk!=Id) { kprintf("Error: distance id  not declared or has bad type %d\n", line); return -1; }
        /*int j=0;
        char a;
        while (a=*(id->Name+j), a >='a' && a <='z' || a >='A' && a <='Z'|| a =='_'){
           printf("%c",*(id->Name+j) ) ;
           j++;
        }*/
        id->Class = Glo; // add here our new type
        id->HClass =Dictionary;
        id->Type = dDST;
        id->Val =(int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
        DictionaryContext brcxt;
        brcxt.Init(this,file,bAddr8);
        memcpy(data,&brcxt,sizeof(brcxt));
        *((int *)data-1)=brcxt.maxlines;
        data = data + sizeof(brcxt);
            
        next();
        continue;
    } else if (tk == Column) { 
        // implement type:
        //     column len macfc nlchar colcxt;
        //
        // colcxt=colcxt+c;  // add new
        // colcxt=0; //reset
        // a=colcxt; // get current fc isNL col
        bool brB=false;
        int clen=0;
        int maxfc=0;
        int nlchar=0;
        // clen value
        next();
        if (tk!=Num) { kprintf("Error: bad Column value len at line %d\n", line); return -1; }
        clen=ival;
        // maxfc value
        next();
        if (tk!=Num) { kprintf("Error: bad Column value fc at line %d\n", line); return -1; }
        maxfc=ival; // clen value
        // maxfc value
        next();
        if (tk!=Num) { kprintf("Error: bad Column value nl at line %d\n", line); return -1; }
        nlchar=ival; // clen value
        next();
        if (tk!=Id) { kprintf("Error: Column id  not declared or has bad type %d\n", line); return -1; }
        /*int j=0;
        char a;
        while (a=*(id->Name+j), a >='a' && a <='z' || a >='A' && a <='Z'|| a =='_'){
           printf("%c",*(id->Name+j) ) ;
           j++;
        }*/
        id->Class = Glo; // add here our new type
        id->HClass =Column;
        id->Type = dDST;
        id->Val =(int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
        ColumnContext colcxt;
        colcxt.Init(this,clen,maxfc,nlchar);
        memcpy(data,&colcxt,sizeof(colcxt));
        data = data + sizeof(colcxt);
        //printf("Column cxt size: %d\n",sizeof(colcxt));
        next();
        continue;
    }
    else if (tk == Table) { 
        // implement type:
        //     table startendchars nlchar rows columns tblcxt;
        //
        // tblcxt=tblcxt+c;  // add new
        // tblcxt=0; //reset
        // a=tblcxt; // get current bracket and distance

        int bType=0; // bracket type 0=char, 1=short
        int bSize=0; // 2,4,6,...
        U32 *bAddr1=0; // address1
        U32 *bAddr2=0; // address2
        next();
        if (id->Class==Glo && id->Type==(iINT+PTR)) {
            if (id->UBound&1) { kprintf("Error: bad table element size at line %d\n", line); return -1; }
            bSize=id->UBound;
            //printf("%d ",bSize);
            int a=id->Type-PTR;
            //if (a==sSHORT) bType=1,bAddr16=(U16 *)id->Val+4; //data after pointerprintf("short ");
            //else 
            if (a==iINT) bAddr1=(U32 *)id->Val+1; else { kprintf("Error: bad table int type %d\n", line); return -1; }
        }else { kprintf("Error: bad table element type at line %d\n", line); return -1; }
        
        next();
        if (id->Class==Glo && id->Type==(iINT+PTR)) {
            if (id->UBound&1) { kprintf("Error: bad table element size at line %d\n", line); return -1; }
            bSize=id->UBound;
            //printf("%d ",bSize);
            int a=id->Type-PTR;
            //if (a==sSHORT) bType=1,bAddr16=(U16 *)id->Val+4; //data after pointerprintf("short ");
            //else 
            if (a==iINT) bAddr2=(U32 *)id->Val+1; else { kprintf("Error: bad table int type %d\n", line); return -1; }
        }else { kprintf("Error: bad table element type at line %d\n", line); return -1; }
        //next();
        //if (tk!=Id) { kprintf("Error: distance id  not declared or has bad type %d\n", line); return -1; }
        /*int j=0;
        char a;
        while (a=*(id->Name+j), a >='a' && a <='z' || a >='A' && a <='Z'|| a =='_'){
           printf("%c",*(id->Name+j) ) ;
           j++;
        }*/
         next();
        if (tk!=Id) { kprintf("Error: table id  not declared or has bad type %d\n", line); return -1; }
        id->Class = Glo; // add here our new type
        id->HClass =Table;
        id->Type = dDST;
        id->Val =(int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
        if (bType==0){
            TableContext brcxt;
            brcxt.Init(&x,bAddr1,bAddr2,bSize);
            memcpy(data,&brcxt,sizeof(brcxt));
            data = data + sizeof(brcxt);
        } else {
            kprintf("Error: bad distance char type %d\n", line); return -1; 
        }
        next();
        continue;
    } else if (tk == Record) { 
        // implement type:
        //     table startendchars nlchar rows columns tblcxt;
        //
        // tblcxt=tblcxt+c;  // add new
        // tblcxt=0; //reset
        // a=tblcxt; // get current bracket and distance

        next();
        if (tk!=Id) { kprintf("Error: record id  not declared or has bad type %d\n", line); return -1; }
        id->Class = Glo; // add here our new type
        id->HClass =Record;
        id->Type = dDST;
        id->Val =(int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
        RecordContext brcxt;
        brcxt.Init(this);
        memcpy(data,&brcxt,sizeof(brcxt));
        data = data + sizeof(brcxt);
        next();
        continue;
    } else if (tk == Enum) {
      next();
      if (tk != '{') next();
      if (tk == '{') {
        next();
        i = 0;
        while (tk != '}') {
          if (tk != Id) { kprintf("%d: bad enum identifier %d\n", line, tk); return -1; }
          next();
          if (tk == Assign) {
            next();
            if (tk != Num) { kprintf("%d: bad enum initializer\n", line); return -1; }
            i = ival;
            next();
          }
          id->Class = Num; id->Type = iINT; id->Val = i++;
          if (tk == Comma) next();
        }
        next();
      }
    }
    while (tk != ';' && tk != '}') {
      ty = bt;
      //while (tk == Mul) { next(); ty = ty + PTR; } // global pointer char *name;
      if (tk != Id) { 
        if (tk== Assign) { // look for global variable assignment
            int negative=0;
            next();
            if (tk == Sub) negative=1,next(); // if negative number
            if (tk != Num) { kprintf("%d: bad glabal variable value. %d\n", line,tk); return -1; }
            *((int *)data)=negative?(int)-ival:(int)ival;
            id->Class = Glo;
            id->Val =(int)data;
            data = data + 4;
            next();
            if (tk == Comma){next();}
            continue;
        }
        kprintf("%d: bad global declaration\n", line); return -1;
      }
      if (id->Class) { kprintf("%d: duplicate global definition\n", line); return -1; }
      next();
      id->Type = ty;
      if (tk == '(') { // function
        /*nam=(char*)id[Name];
        printf("\n// FUNCTION: ");
        for (int y=0;y<id[IDLen];y++)printf("%c",nam[y]);
        printf("\n");*/
        id->Class = Fun;
        id->Val = (int)(e + 1);
        if (idbitupdate->Val!=0 && idbyteupdate->Val==0) { kprintf("%d: byteupdate must be first\n", line); return -1; }
        next(); i = 0;
        while (tk != ')') {
          ty = iINT;
          if (tk == Int) next();
          else if (tk == Char) { next(); ty = rCHAR; }
          else if (tk == Short) { next(); ty = sSHORT; }
          while (tk == Mul) { next(); ty = ty + PTR; }
          if (tk != Id) { kprintf("%d: bad parameter declaration\n", line); return -1; }
          if (id->Class == Loc) { kprintf("%d: duplicate parameter definition\n", line); return -1; }
          id->HClass = id->Class; id->Class = Loc;
          id->HType  = id->Type;  id->Type = ty;
          id->HVal   = id->Val;   id->Val = i++;
          next();
          if (tk == Comma) next();
        }
        next();
        if (tk != '{') { kprintf("%d: bad function definition\n", line); return -1; }
        loc = ++i;
        next();
        while (tk == Int || tk == Char || tk == Short) {
          bt = (tk == Int) ? iINT : (tk == Short) ? sSHORT : rCHAR;;
          next();
          while (tk != ';') {
            ty = bt;
            //while (tk == Mul) { next(); ty = ty + PTR; } // local pointer ex. char *name;
            if (tk != Id) { kprintf("%d: bad local declaration\n", line); return -1; }
            if (id->Class == Loc) { kprintf("%d: duplicate local definition\n", line); return -1; }
            id->HClass = id->Class; id->Class = Loc;
            id->HType  = id->Type;  id->Type = ty;
            id->HVal   = id->Val;   id->Val = ++i;
            next();
            if (tk == Comma) next();
          }
          next();
        }
         *++e=(  ENT);*++e=( i - loc);
        while (tk != '}') stmt();
        if (*e!=LEV)
        *++e=( LEV);
        id = sym; // unwind symbol table locals
        while (id->Tk) {
          if (id->Class == Loc) {
            id->Class = id->HClass;
            id->Type = id->HType;
            id->Val = id->HVal;
          }
          id = id + 1;
        }
      }
      else if (tk == Brak) { // global array
        id->Class = Glo;
        id->Val =(int)data;
        *((int *)data)=(int)(data+4); //store pointer, data starts after pointer
        data = data + 4;
        int idz=id->Type;
        id->Type=id->Type+ PTR;
        if (idz>iINT) { kprintf("%d: bad global declaration only char, short, int \n", line); return -1; }
        next(); i = ival; //array size
        if (tk==Id && id->Class == Num) { // array size is enum value
           i=id->Val,tk=Num;              // set to Number toke with enum value
        }
        if (tk != Num )          { kprintf("%d: bad global declaration: size=%d\n", line,i); return -1; }
        if (i==0 )               { kprintf("%d: array to small \n", line); return -1; }
        id->UBound=i;
        next();if (tk != ']')    { kprintf("%d: missing closing braket\n", line); return -1; }
        next();if (tk != Assign) { kprintf("%d: missing array assingn \n", line); return -1; }
        next();if (tk != '{')    { kprintf("%d: missing array { \n", line); return -1; }
        next();
        int count=0;
        if (tk == '}') { 
           // type name[size]={}
           // dynamic array, i=size idz=sizeof(type)
           data = data - 4; // move back to our data pointer
           *((int *)data)=(int)vmmalloc(this,i,(idz==0)?1:(idz==1)?2:4);
           //printf("Dynamic Alloc: %d(%d)\n",(idz==0)?1:(idz==1)?2:4,i);
           data = data + 4; // for next data pointer
           next();          // next token 
        }else{
            // static content type name[size]={0,...}
            // store data after pointer
            // allow add/or of array elements
            //int lastIval=0;
            bool doOr=false;
            bool doAdd=false;
            while (tk != '}') {
              if (tk != Num) { 
                 if (tk==Id && id->Class == Num) { //  is value a enum
                    ival=id->Val;//,tk=Num;              // get enum value
                    //kprintf("enumb value  %d, line %d \n",ival, line);
                 }else {
                    kprintf("%d: bad glabal array value\n", line); return -1;} 
                 }
              if (doOr==true){
                doOr=false;
                //do 
                //kprintf("enumb or value  %d, line %d \n",ival, line);
                if (idz==0) data--,*((char *)data) |=(char)ival,data = data + 1;
                if (idz==1) data = data - 2,*((short *)data) |=(short)ival,data = data + 2;
                if (idz==2) data = data - 4,*((int *)data)|=(int)ival,data = data + 4;
              }else if (doAdd==true){
                doAdd=false;
                //do unsigned add in array elements
                //kprintf("enumb or value  %d, line %d \n",ival, line);
                if (idz==0) data--,*((unsigned char *)data) +=(unsigned char)ival,data = data + 1;
                if (idz==1) data = data - 2,*((unsigned short *)data) +=(unsigned short)ival,data = data + 2;
                if (idz==2) data = data - 4,*((unsigned int *)data)+=(unsigned int)ival,data = data + 4;
              }else{
                if (idz==0) *((char *)data) =(char)ival,data = data + 1;
                if (idz==1) *((short *)data) =(short)ival,data = data + 2;
                if (idz==2) *((int *)data)=(int)ival,data = data + 4;
                count++;
              }
              if (count > i) { kprintf("%d: array out of bounds \n", line); return -1; }
              next();
              if (tk == Comma) next(); //else { kprintf("%d: comma expected \n", line); return -1; }
              if (tk == Or) doOr=true,next(); //else { kprintf("%d: comma expected \n", line); return -1; }
              if (tk == Add) doAdd=true,next(); //else { kprintf("%d: comma expected \n", line); return -1; }
           }
            if (count != i) { 
              for (int m=count;m<i;m++){
                if (idz==0) *((char *)data) =0,data = data + 1;
                if (idz==1) *((short *)data) =0,data = data + 2;
                if (idz==2) *((int *)data)=0,data = data + 4;
              }
            //kprintf("%d: too few array elements\n", line); return -1; 
            }
       }
      }
      else {
        id->Class = Glo;
        id->Val = (int)data;
        *((int *)data)=0;
        data = data + sizeof(int);
      }
      if (tk == Comma) next();
    }
    next();
  }
    if (vmMode==VMCOMPRESS && (idbyteupdate->Val && idbitupdate->Val) ) {
       isDualMode=true;
       //printf("Dual mode\n");
    }
    else {
    
       if (vmMode==VMDETECT && (iddetect->Val==0 ||  idmain->Val==0) ) quit("detect or main not defined");
       if (vmMode==VMENCODE && (idencode->Val==0 ||  idmain->Val==0) ) quit("encode or main not defined");
       if (vmMode==VMDECODE && (iddecode->Val==0 ||  idmain->Val==0) ) quit("decode or main not defined");
       if (vmMode==VMCOMPRESS && (idp->Val==0 || idupdate->Val==0 || idmain->Val==0) ) quit("block update or main not defined");
       
    }
       if (vmMode==VMCOMPRESS) doGetComponents=true; // get component info on first run
/*
    id=sym;
    int s=0;
    for (int i=0;i<idCount;i++)  {
        if (id[i].Tk==0) break;
        printf("%3d %3d %d 0x%08x 0x%08x ",id[i].Tk,id[i].Class ,id[i].Type,id[i].Val,id[i].ePTR);
        if (id[i].Class==Glo) {
            int a=0;
            if (id[i].Type>iINT) a=id[i].Type-PTR,printf(" *");
            if (a==iINT) printf("int ");
            else if (a==sSHORT) printf("short ");
            else if (a==rCHAR) printf("char "); 
        }
        int j=0;
        char a;
        ident_s *t=&id[i];
        while (a=*(t->Name+j), a >='a' && a <='z' || a >='A' && a <='Z'|| a =='_'){
           printf("%c",*(t->Name+j) ) ;
           j++;
        }
        printf("\n" ) ;
        //id = id + 1;
        s=i;
      }
      printf("Symbol count %d\n",s);*/

    if (doJIT){
        // setup stack and call stage 1
        bp = sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; 
        t = sp;
        *--sp = (int)t;
        int r = dovm((int *)idmain->Val);
        initComp();
        if (dojit()!=0) return -1;
        // run jit code and call stage 1
        int (*jitmain)(); // c4 vm pushes first argument first, unlike cdecl
        jitmain = reinterpret_cast< int(*)() >(*(unsigned*)(idmain->Val) >> 8 | ((unsigned)jitmem & 0xff000000));
         r = jitmain();
        // stage 2
        //r = jitmain();
        return r;
    } else {
        // setup stack and call stage 1
        bp = sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; 
        t = sp;
        *--sp = (int)t;
        int r = dovm((int *)idmain->Val);
        initComp();
        // stage 2
        bp = sp = (int *)((int)sp0 + poolsz);
        *--sp = EXIT; // call exit if main returns
        *--sp = PSH; 
        t = sp;
        *--sp = (int)t;
        r = dovm((int *)idmain->Val);
        return r;
    }
    //no code after this, JIT kills this
}
 



