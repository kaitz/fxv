

// TYPA commands
enum { DISTANCE=1, COLUMN=2, TABLE=3, RECORD=4, DICTIONARY=5
};
// tokens and classes (operators last and in precedence order)
enum { Num=128, Fun     , Sys   , Glo  , Loc   , Id   , Load , Enter , Char, Else, Enum, If , Int, Short, Return, For,
       Sizeof , Distance, Column, Table, Record, Dictionary, While, Comma, Assign, Cond, Lor , Lan , Or , Xor, And  , Eq    , Ne,
       Lt     , Gt      , Le    , Ge   , Shl   , Shr  , Add  , Sub   , Mul , Div , Mod , Inc, Dec, Brak
};
// opcodes
enum { LEA ,IMM ,JMP  ,JSR  ,BZ  ,BNZ ,LBZ ,LBNZ,ENT ,ADJ ,LEV ,LI  ,LS  ,LC  ,SI  ,SS  ,
       SC  ,PSH ,OR   ,XOR  ,AND ,EQ  ,NE  ,LT  ,GT  ,LE  ,GE  ,SHL ,SHR ,ADD ,SUB ,MUL ,
       DIV ,MOD ,BOUND,VTHIS,PRTF,VMI ,VMX ,READ,WRTE,EXIT,TYPA,TYPR,NOP ,VMXD
};

// types
enum { rCHAR, sSHORT, iINT, PTR, dDST };

enum { VMCOMPRESS=0,VMDETECT,VMENCODE,VMDECODE,VMTUNE};

#ifdef VMMSG
#define kprintf printf // prints error messages
//#define dprintf printf // prints x86 asm to console
#define dprintf(...)
#else
#define kprintf(...)    
#define dprintf(...)    
#endif

struct Component{
    // 0xppiiccmm pr index,component index, component, mixer index
    U8 prIndex;
    U8 compIndex;
    U8 component;
    U8 mixerIndex;
};


class VM {
private:
    
    Array<char> data1; // data
    char *p, *lp, // current position in source code
     *data,*data0,   // data/bss pointer,
     *jitmem, // executable memory for JIT-compiled native code
     *je;     // current position in emitted native code
    int *e, *le, *text,*codestart,  // current position in emitted code

    //*n,       // current node in abstract syntax tree
    tk,       // current token
    ival,     // current token value
    ty,       // current expression type
    loc,      // local variable offset
    line;     // current line number
    int  bt,   poolsz;
    
    int *pc, *sp,*sp0, *bp, cycle; // vm registers
    int *t,*pc0,tmp; // temps
    int a;
    int initvm(); 
    char *jitmempr; // get prediction JIT code memory
    char *jitmemup; // get update JIT code memory
    char *jepbyte; // byteupdate
    // identifier
    struct ident_s {
      int Tk;
      int Hash;
      char *Name;
      int Class;
      int Type;
      int Val;
      int Stype;
      int HClass;
      int HType;
      int HVal;
      int UBound;
      int ePTR;
    } *id,  // currently parsed identifier
      *sym; // symbol table (simple list of identifiers)
    int idCount;
    struct ident_s *idmain,  *idp, *idupdate, *iddetect, *iddecode, *idencode, *idbitupdate,* idbyteupdate;
    void getPredictionJIT();
    void initComp();
public:
    object &mod; // Our current model in plain text
    int debug;    // print executed instructions
    BlockData& x;
    int vmMode;
    Array<char*> mem; //array of allocated memory
    Array<int> memSize; // size in bytes of memory in mem[mindex]    
    Array<int*> membound; //array of allocated memory
    Array<int> prSize; // size in bytes of memory in mem[mindex]   
    int totalPR;
    Array<Component> mcomp;  //component list set in vmi
    int smc, apm1, rcm, scm, cm, mx,st,av,ds,mm,dhs,sm,sk,apm2,em,uas,lmx,statec,byt;
    ContextMap *cmC;
    MixMap1 *mmA;
    Mixer1 *mxA;
    StaticMap *stA;
    DynamicSMap *dsA;
    DynamicHSMap *dhsA;
    APM1 *apm1A;
    StateMapContext *smA;
    AvgMap  *avA;   
    SmallStationaryContextMap *scmA;    
    StationaryMap *smcA;    
    RunContextMap *rcmA;  
    SkMap *skA;  
    APM2 *apm2A;
    ErrMap *emA;
    UAS *uasA;
    LmxMap  *lmxA;   
    ByteMap *bmA;
    COMPONENT *comp;
    int compIdx;
    int totalc;  //total number of components
    int currentc; //current component, used in vmi
    int mindex;        // count fo memory allocations
    FILE *inFile, *outFile; // files for decoding and encoding
    int inpos;
    VMParam *parm;
    int posParamCount;
    StateTable  *vmstate;

    bool isDualMode;
    bool isBitUpdate;
    bool doGetComponents;
    VM(object &mod,BlockData& bd,int mode, VMParam *p=0);
    ~VM() ;
    void next();
    void expr(int lev);
    void stmt();
    int dovm(int *ttt);
    void gen(int *n);
    int dojit();
    int detect(int c4,int pos);
    int decode(int info,int len);
    int encode(int info,int len);
    int block(int info1,int info2);
    int doupdate1(int y, int c0, int bpos,U32 c4);
    int doupdate2(int y, int c0, int bpos,U32 c4);
    int getPrediction( );
    void updateComponents();
    inline void updateComponents1();
    void getBytePredictionJIT();

};

// alloc function in interpreted code
// keep track of pointers and sizes in bytes
// no bounds test
// Automatically free memory when VM is deleted 
char* vmmalloc(VM* v, size_t i, int w) {
  char*ptr= (char*)calloc(i*w,1);
  if (ptr==0) perror("mem error "),printf("%d ",i),quit("VM mem alloc fail");
  //printf("Alloc %d %d %x\n",i,w,ptr);
  v->mem.resize(v->mem.size()+1);
  v->mindex =v->mem.size();
  v->mem[v->mindex-1]=ptr;
  v->memSize.resize(v->mindex);
  v->memSize[v->mindex-1]=i*w;
  return ptr;
}
