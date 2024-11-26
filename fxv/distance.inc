// This part is based on cmix BracketContext
template <typename T  = U8>
struct BracketContext {
    U32 context;           // bracket byte and distance
    vec<int,512> active;   // vector for brackets, max 512 elements
    vec<int,512> distance; // vector for distance, max 512 elements
    const T *element;
    int elementCount;
    bool doPop;            // set true for quotes
    int limit;
    T cxt,dst;

    void Init(T *d,const int e,int pop=false,int l=(1<< (sizeof(T)*8))) {
        elementCount=e;
        element=d;
        context=cxt=dst=0;
        doPop=pop;
        limit=l;
        //printf("Bracket limit: %d Size:%d\n",l,1<< (sizeof(T)*8));
        vec_new(&active);
        vec_new(&distance);
    }
    void __attribute__ ((noinline)) Reset(){
        vec_reset(&active);
        vec_reset(&distance);
        context=cxt=dst=0;
    }
    bool Find(int b){
        bool found=false;
        for (int i=0;i<elementCount;i=i+2) if (element[i]==b) {
            found=true;
            break;
        }
        return found;
    }
    bool FindEnd(int b,int c){
        bool found=false;
        for (int i=0;i<elementCount;i=i+2) if (element[i]==b&&element[i+1]==c) found=true;
        return found;
    }
    int last(){
        return vec_prev(&active);
    }
    int getcxt(){
        return context;
    }
    void __attribute__ ((noinline)) Update(int byte) {
        bool pop=false;
        if (!vec_empty(&active)) {
            if (FindEnd(vec_at(&active,vec_size(&active)-1) , byte) || vec_at(&distance,vec_size(&distance)-1) >= limit) {
                vec_pop(&active);
                vec_pop(&distance);
                pop=doPop;
            } else {
                vec_i(&distance,vec_size(&distance)-1);
            }
        }
        if (pop==false && Find(byte)) {
            vec_push( &active,byte);
            vec_push( &distance,0);
        }
        if (!vec_empty(&active)) {
            cxt=vec_at(&active,vec_size(&active)-1);
            dst=min(vec_at(&distance,vec_size(&distance)-1),(1<< (sizeof(T)*8))-1);
            context = (1<< ((sizeof(T)*8))) * cxt+dst;
        } else {
            context=cxt=dst=0;
        }
    }
};
