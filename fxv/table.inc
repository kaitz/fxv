struct TableContext {
    BlockData *x;
    vec<U32,16*2> row[4];   // Content of table row cell positions, max 4 rows, max 16*2 positions
    int cellCount,rows,abovecellpos,abovecellpos1;
    U8 nlChar;
    U32 *startend;
    U32 startendsize; // size of 6, startid,mask,endid,mask,nlchar,mask
    U32 *cellID;
    U32 cellIDcount;
    bool isActive;
    U32 context; 
    void Init(BlockData *bd, U32 *d, U32 *e, int idco) {
        x=bd;
        startend=d;
        startendsize=6;
        cellID=e;
        cellIDcount=idco;
        rows=abovecellpos=cellCount=abovecellpos1=context=0;
        isActive=false;
        for (int i=0;i<4;i++) vec_new(&row[i]);
       // printf("%x %x\n",startend[0],startend[1]);
       // printf("%x %x\n",startend[2],startend[3]);
        //printf("%x %x\n",startend[4],startend[5]);
    }
    void Reset(){
        resetCells();
        rows=abovecellpos=cellCount=abovecellpos1=context=0;
        isActive=false;
    }
    bool Find(U32 b){
        bool found=false;
        for (int i=0;i<cellIDcount;i=i+2) if (cellID[i]==(b&cellID[i+1])) {
            found=true;
            break;
        }
        return found;
    }
    void __attribute__ ((noinline)) Update(U32 byte) {
        // Start and end of table
        if ( (byte&startend[1])==startend[0] && isActive==false) isActive=true;
        else if ( (byte&startend[3])==startend[2] ) isActive=false,context=0,resetCells();
        if (isActive==true){
            if ((byte&startend[5])==startend[4]){
                rows++;
                rows=rows&3;
                vec_reset(&row[rows]); // reset new row.
                vec_push( &row[rows],U32(x->blpos));
                cellCount=abovecellpos=abovecellpos1=0;
            }
            bool newcell=false;
            // Cells
            if (Find(byte)) vec_push( &row[rows],U32(x->blpos)),cellCount++,newcell=true;
            // Advence above cell pos
            if (abovecellpos ) {
                abovecellpos++;
                // When above cell is shorter reset
                if (abovecellpos>abovecellpos1) abovecellpos=abovecellpos1=0;
            }
            // If more then one cell get above cell based on current row cell
            if(newcell==true && cellsCount() >0){
                // Get current above cell pos
                context=abovecellpos=cellPos(cellCount-1);
                abovecellpos1=cellPos(cellCount);
            }
        }
    }
    int cellsCount(int rowsel=1){
        return vec_size(&row[(rows-rowsel)&3]);
    }
    int cellPos(int cellID,int rowsel=1){
        int total=cellsCount(rowsel)-1;
        total=min(total,cellID);
        return vec_at(&row[(rows-rowsel)&3],total);
    }
    void resetCells(){
        for (int i=0;i<4;i++) vec_reset(&row[i]);
    }
};

