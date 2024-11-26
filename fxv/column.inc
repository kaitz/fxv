#define LF 10
// Table/row & column context
struct ColumnS {
    U32 linepos;
    U8 fc;
    vec<U8,1024*2> bytes; // max lenght 1024*2 chars
};
#define WIKIHEADER GREATERTHAN
#define WIKITABLE  '-'
struct ColumnContext {
    VM *x;
    ColumnS col[4];           // Content of last 3 + current row
    int rows;
    bool NL,isTemp;
    int limit;  // column lenght limit
    U8 nlChar;
    U8 maxFC;
    U32 context;
    void Init(VM *bd, int l,int mfc,int nlc) {
        x=bd;
        maxFC=mfc;
        rows=0;
        nlChar=nlc;
        limit=l;
        for (int i=0;i<4;i++) {
            vec_new(&col[i].bytes);
            col[i].linepos=0;
            col[i].fc=0;
        }
        NL=isTemp=false;
        context=0;
    }
   
    U8 lastfc(int i=0){
        return col[(rows-i)&3].fc;
    }
    bool isNewLine(){
        return NL;
    }
    int  __attribute__ ((noinline)) collen(int i=0,int l=0){
        return min((l?l:limit), vec_size(&col[(rows-i)&3].bytes)+1);
    }
    int nlpos(int i=0){
        return col[(rows-i)&3].linepos;
    }
    U8  __attribute__ ((noinline)) colb(int i=1,int j=0,int l=0){
        return  vec_at(&col[(rows-i)&3].bytes,collen()-(1+j));
    }
    void Update(int byte) {
        // Column
        NL=false;
        if (U8(byte)==nlChar){
            vec_push( &col[rows].bytes,U8(byte));
            rows++;
            rows=rows&3;
            vec_reset(&col[rows].bytes); // reset new line.
            col[rows].fc=0;
            col[rows].linepos=x->x.blpos-1;
        }  else{
            vec_push( &col[rows].bytes,U8(byte)); // set new byte to line
            if (collen()==2) {
                col[rows].fc=min(U8(byte),maxFC);
                NL=true;
            }
        }
        context=lastfc()+collen()*256+(colb()*256+NL*0x10000)*256;
    }
};

