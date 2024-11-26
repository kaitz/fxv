
struct DictionaryContext{
    VM *x;
    char **dictW;           // String pointer array
    U8 *dictWLen;           // Dictionary word lengh table 
    U32 *dictHash;          // Dictionary word hash table 
    U32 *dictEncode;        // String to encoded index array
    int *codeword2sym;  // For decoding encoded dictionry index
    int dict1size;
    int dict2size;
    int dict12size;
    int sizeDict;           // Current size of words
    U32 lastCW=0;
    int maxlines;           // Maximum words defined in input file
    char *text;             // Input/output string in VM
    U32 context;            // Contex used by VM
    void Init(VM *bd, char *file, char *str){
        x=bd;
        text=str;
        // for decode
        dict1size=80;
        dict2size=32;
        dict12size=dict1size*dict2size;
        context=0;
        codeword2sym=(int*)vmmalloc(bd, 256, sizeof(int)) ;
        FILE *f=fopen(file,"rb");
        if (f==0) printf("Can not open dictionary file %s\n",file),exit(1);
        // printf("%s\n",file);
        loaddict(f);
        fclose(f);
        dosym();
    }
    int wfgets(char *str, int count, FILE  *fp) {
        int c, i = 0;
        while (i<count-1 && ((c=getc(fp))!=EOF)) {
            str[i++]=c; if (c=='\n')str[i-1]=0;
            if (c=='\n')
            break;
        }
        str[i]=0;
        return i;
    }
    bool compare(int j, int l){
        int i=0;
        for (i=0; i<l; i++) {
            char a=text[i];
            if (a!=dictW[j][i]) break;
        }
        return i==l;
    }
    int getIndex() {
        //printf("GetIndex 0\n");
        U32 hash=0;
        int len=0;
        for (int i=0; i<64; i++) {
            char a=text[i];
            if (a==0) break;
            hash=hash*773+(hash>>16)*11+(U8)a;
            len++;
        }
        //printf("Len %d %s\n",len,text);
        if (len) {
            for (int i=0; i<sizeDict; i++) {
                if (dictHash[i]==hash && dictWLen[i]==len) {
                    if (compare(i,len)) {
                        context=i;
                        // printf("Found string %s, len %d, %d\n",dictW[i],len,context);
                        return context;
                    }
                }
            }
        }
        context= -1;
        return context;
    }
    int  get(int a) {
        if (a==0) return getIndex();
        else if (a==1) return getEncodeIndex();
        else if (a==2) return getEncodeIndexSub();
        else if (a>127) return getWordFromIndex(a);
        else printf("Bad range in dectionary operation selection.\n"), exit(1);
        return -1;
    }
    int getEncodeIndex() {
        // printf("GetIndex 1\n");
        U32 hash=0;
        int len=0;
        for (int i=0; i<64; i++) {
            char a=text[i];
            if (a==0) break;
            hash=hash*773+(hash>>16)*11+(U8)a;
            len++;
        }
        //printf("Len %d\n",len);
        if (len) {
            for (int i=0; i<sizeDict; i++) {
                if (dictHash[i]==hash && dictWLen[i]==len) {
                    if (compare(i,len)) {
                        context=dictEncode[i];
                        //printf("Found string %s, len %d, %d\n",dictW[i],len,context);
                        return context;
                    }
                }
            }
        }
        context= -1;
        return context;
    }
    int getEncodeIndexSub() {
        // printf("GetIndex 2\n");
        U32 hash=0;
        int len=0;
        for (int i=0; i<64; i++) {
            char a=text[i];
            if (a==0) break;
            len++;
        }
        if (len<5){
            context= -1;
            return context;
        }
        int maxlen=len-1;
        int lastvalid=-1;
        int validID=0;
        //printf("Find: %s\n",text);
        for (int i=0; i<len; i++) {
            context= -1;
            hash=0;
            for (int j=0; j<maxlen; j++) {
                char a=text[j];
                if (a==0) break;
                hash=hash*773+(hash>>16)*11+(U8)a;
            }
            for (int k=0; k<sizeDict; k++) {
                if (dictHash[k]==hash && dictWLen[k]==maxlen) {
                    if (compare(k,maxlen)) {
                        lastvalid=context=dictEncode[k];
                        validID=k;
                        lastvalid=lastvalid+maxlen*0x1000000;
                    }
                }
            }
            if (context!=-1) maxlen++;
            else break;
        }
        if ((maxlen-1)<=(len-3)) lastvalid=-1;
        //if (lastvalid>0)printf("Found sub-string %s of %s, len %d of %d, %d\n",dictW[validID],text,maxlen,len,context);
        context= lastvalid;
        return context;
    }
    int getWordFromIndex(int a) {
        //printf("Decode to string\n");
        decodeWord(a);
        if (lastCW>=44880) return -1;
        //printf("Decoded into word:%d %s %d\n",a,dictW[lastCW],lastCW);
        memcpy(text, dictW[lastCW],dictWLen[lastCW]+1);
        return dictWLen[lastCW];
    }
    void buildHashTable() {
        for (int i=0; i<sizeDict; i++) {
            dictHash[i]=0;
            for (int j=0; j<dictWLen[i]; j++) {
                dictHash[i]=dictHash[i]*773+(dictHash[i]>>16)*11+dictW[i][j];
            }
        }
    }
    // load whole file and parse in memory
    void loaddict(FILE  *file){
        int len=0;
        int line_count=0;
        char *s=(char *)calloc(64,1);
        len=wfgets(s, 64, file);
        maxlines=atoi(s);
        // get data lenght
        int pos=ftell(file);
        fseek (file, 0 ,SEEK_END);
        int isize=ftell(file)-pos; 
        fseek (file, pos ,SEEK_SET);
        if (isize==0) printf("No words found\n"),exit(1);
        if (maxlines<=0) printf("Defined word count in dictionary is bad %d\n",maxlines),exit(1);
        //printf("Words in dictionary %d\n",maxlines);
        dictW=(char **)vmmalloc(x, maxlines, sizeof(char *));
        dictWLen=(U8 *)vmmalloc(x, maxlines, sizeof(char));
        dictHash=(U32 *)vmmalloc(x, maxlines, sizeof(U32));
        dictEncode=(U32 *)vmmalloc(x, maxlines, sizeof(U32));
        char *words=(char *)vmmalloc(x, isize+1, sizeof(char)); // Alloc whole file and load into memory
        fread(words, 1,isize,file); 
        int curpos=0;
        len=0;
        const int dictBound1 = 80, dictBound2 = dictBound1 + 3840, dictBound3 = dictBound2 + 40960;
        // Parse line by line
        while (curpos <isize && *words!=0)  {
            if (line_count==maxlines &&  line_count<44880 ) {
                break;
            }
            dictW[line_count]=words;
            while (*words!='\n' && *words!=0) {
                words++;
                curpos++;
                len++;
            }
            *words=0; // set string end
            words++; curpos++;
            dictWLen[line_count]=U8(len);
            if (len==0) printf("Warning: dictionary word lenght zero. Line %d\n",line_count); // Maybe fail here?
            if (len>64) printf("Warning: dictionary word lenght to large (max 64). Line %d, len %d\n",line_count,len);
            //printf("%d,%d,%s\n",len,line_count,dictW[line_count]);
            U32 codeword=0;
            if (line_count < dictBound1) {
                codeword = 0x80 + line_count;
            } else if (line_count < dictBound2) {
                codeword = 0xD0 + ((line_count-dictBound1) / 80);
                codeword += (0x80 + ((line_count-dictBound1) % 80)) << 8;
            } else if (line_count < dictBound3) {
                codeword = 0xF0 + (((line_count-dictBound2) / 80) / 32);
                codeword += (0xD0 + (((line_count-dictBound2) / 80) % 32)) << 8;
                codeword += (0x80 + ((line_count-dictBound2) % 80)) << 16;
            } else printf("Bad codeword.\n");
            dictEncode[line_count]=codeword;
            line_count++;
            len=0;
        }
        if (curpos<isize) {
            printf("To many words in dictionary\n");
        }
        if (line_count<maxlines) {
            printf("To few words in dictionary.\n");
        }
        free(s);
        // printf("Loaded %d words\n",line_count);
        context=sizeDict=line_count;
        buildHashTable();
    }
    inline int decodeCodeWord(int cw) {
        int i=0;
        int c=cw&255;
        if (codeword2sym[c]<dict1size) {
            i=codeword2sym[c];
            return i;
        }

        i=dict1size*(codeword2sym[c]-dict1size);
        c=(cw>>8)&255;

        if (codeword2sym[c]<dict1size) {
            i+=codeword2sym[c];
            return i+dict1size;
        }

        i=(i-dict12size)*dict2size;
        i+=dict1size*(codeword2sym[c]-dict1size);

        c=(cw>>16)&255;
        i+=codeword2sym[c];
        return i+80*49;
    }

    void dosym(){
        for (int c=0; c<256; c++){
            codeword2sym[c]=0;
        }
        int charsUsed=0;
        for (int c=128; c<256; c++) {
            codeword2sym[c]=charsUsed;
            charsUsed++;
        }
    }

    void decodeWord(int c){
        lastCW=decodeCodeWord(c);
        assert(lastCW>=0 && lastCW<44880);
    }
};
