#define RECBUFFER 0x10000*7
struct RecordContext {
    VM *x;
    U32 context; 
    int rlen[3];
    int rcount[2];
    int *cpos1, *cpos2, *cpos3, *cpos4;
    U8 *buffer;
    int pos;
    void Init(VM *bd) {
        x=bd;
        cpos1=(int*)vmmalloc(bd, 256, sizeof(int)) ;
        cpos2=(int*)vmmalloc(bd, 256, sizeof(int)) ;
        cpos3=(int*)vmmalloc(bd, 256, sizeof(int)) ;
        cpos4=(int*)vmmalloc(bd, 256, sizeof(int)) ;
        buffer=(U8*)vmmalloc(bd, RECBUFFER, sizeof(U8)) ;
        Reset();
    }
    void Reset(){
        context=pos=0;
        // run length and 2 candidates
        rlen[0] = 2; 
        rlen[1] = 0; 
        rlen[2] = 0; 
        // candidate counts
        rcount[0] = 0;
        rcount[1] = 0;
        memset(buffer,   0, RECBUFFER);
        memset(cpos1,    0, sizeof(cpos1));
        memset(cpos2,    0, sizeof(cpos2));
        memset(cpos3,    0, sizeof(cpos3));
        memset(cpos4,    0, sizeof(cpos4));
    }
    int buf(int i){
        return buffer[(pos-i)&(RECBUFFER-1)];
    }
    void __attribute__ ((noinline)) Update(U32 byte) {
        buffer[pos&(RECBUFFER-1)]=U8(byte);
        pos++;
        
        int r=0;
        r=min(0x10000-1,pos-cpos1[byte]);
      if ( r>1) {
        if ( r==cpos1[byte]-cpos2[byte]  && r==cpos2[byte]-cpos3[byte] && (r>32 || r==cpos3[byte]-cpos4[byte])
          && (r>10 || ((byte==buf(r*5+1)) && byte==buf(r*6+1)))) {      
          if (r==rlen[1]) ++rcount[0];
          else if (r==rlen[2]) ++rcount[1];
          else if (rcount[0]>rcount[1]) rlen[2]=r, rcount[1]=1;
          else rlen[1]=r, rcount[0]=1;
        }
      }    

    // check candidate lengths
    for (int i=0; i < 2; i++) {
      if (rcount[i] > max(0,12-(int)ilog2(rlen[i+1]))){
        if (rlen[0] != rlen[i+1]){
            if ( (rlen[i+1] > rlen[0]) && (rlen[i+1] % rlen[0] == 0) ){
            // maybe we found a multiple of the real record size..?
            // in that case, it is probably an immediate multiple (2x).
            // that is probably more likely the bigger the length, so
            // check for small lengths too
            if ((rlen[0] > 32) && (rlen[i+1] == rlen[0]*2)){
              rcount[0]>>=1;
              rcount[1]>>=1;
              continue;
            }
          }
          rlen[0] = rlen[i+1];
          rcount[i] = 0;
        }
        else
          // we found the same length again, that's positive reinforcement that
          // this really is the correct record size, so give it a little boost
          rcount[i]>>=2;

        // if the other candidate record length is orders of
        // magnitude larger, it will probably never have enough time
        // to increase its counter before it's reset again. and if
        // this length is not a multiple of the other, than it might
        // really be worthwhile to investigate it, so we won't set its
        // counter to 0
        if (rlen[i+1]<<4 > rlen[1+(i^1)])
          rcount[i^1] = 0;
      }
    }
    assert(rlen[0]>0);
    // context
    context=min(0x10000-1,rlen[0])+0x10000*min(255,  pos-cpos1[byte]);
    // update last context positions
    cpos4[byte]=cpos3[byte];
    cpos3[byte]=cpos2[byte];
    cpos2[byte]=cpos1[byte];
    cpos1[byte]=pos;
    }
};
