// Vector for different contexts

template <typename T = int,const int S=256 >
struct vec {
    T cxt[S];
    static constexpr int capacity=S;
    int size=0;
};

template <typename T = int,const int S>
void vec_new(vec<T,S>* o){
    o->size=0;
}
template <typename T = int,const int S>
int vec_size(vec<T,S> *o){
    return o->size;
}
template <typename T  = int,const int S>
void vec_push( vec<T,S> *o, const T element){
    assert(o->size<S);
    o->cxt[o->size++]=element;
    o->size=o->size&(o->capacity-1); // roll over
}
template <typename T  = int,const int S>
int vec_at(vec<T,S> *o, const int index){
    assert(index<S);
    return o->cxt[index];
}
template <typename T  = int,const int S>
T &vec_ref(vec<T,S> *o, const int index){
    assert(index<S);
    return o->cxt[index];
}
template <typename T  = int,const int S>
void vec_i(vec<T,S> *o, const int index){
    o->cxt[index]++;
}
template <typename T  = int,const int S>
void vec_pop(vec<T,S> *o){
    assert(o->size<S);
    if (o->size>0) o->cxt[o->size]=0, o->size--; // no rollback
}
template <typename T  = int,const int S>
void vec_reset(vec<T,S> *o){
    o->cxt[0]=0;
    o->size=0;
}
template <typename T  = int,const int S>
bool vec_empty(vec<T,S> *o){
    return (o->size==0)?true:false;
}
template <typename T  = int,const int S>
int vec_prev(vec<T,S> *o){
    assert(o->size>=0);
    return (o->size>1)?(o->cxt[o->size-2]):0; // no rollback
}
