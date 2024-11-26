# Contents
[Config files](#config-files)

[Config files structure](#config-files-structure)

[Example detection config](#example-detection-config)

[Example compression model](#compression-model)



## Config files
Configuration file ends with .pxv file extension.
* .det files are used for detecting some type of data
* .enc files are used to transform a detected file before compressing
* .dec files are used to transform a detected file after decompressing (reverse transform)
* .cfg files are models used for compression and decompression

By default .cfg and .dec files are compressed and added to the archive if they are used. Models can be excluded from the archive.

### Config files structure

Smallest config file we can have has one stream where all data is going.
Model we use is myalgo.cfg. 

Our data type is 0 and there are no detect, decode or encode config files. Data type 0 and its corresponding stream must be always present in the config file.
It is called the default data type.

```
// Streams
stream 0
model myalgo.cfg
// Default type
type 0
detect -1
decode -1
encode -1
compress 0
```

fxcm.pxv structure is below. Whe have 2 streams. One for unknown data and one for detected text data.
Compressor attempts to detect text using text.det file. If it is found then it uses text.enc (when compressing) for transform and text.dec for reverse transform (when decompressing). 
Final data is compressed by stream 1 model fxcm.cfg. cfg files can be in subdirectories.

```
// Streams
stream 0
model fxcm\null.cfg
stream 1
model fxcm\fxcm.cfg
// Default type
type 0
detect -1
decode -1
encode -1
compress 0
// Type 1 for text
type 1
detect fxcm\text.det
decode fxcm\text.dec
encode fxcm\text.enc
compress 1
```


### Example detection config
```c
// My custom X type detection
int buf0,buf1,mystart;
int type,state,jstart,jend;
enum {DEFAULT=1, YOURTYPE}; //internal enum
// function will report its state 
// or if i=-1 then state results otherwise i is pos
// c4 is last 4 bytes
void reset() {
    state=NONE, type=DEFAULT, jstart=jend=buf0=buf1=mystart=0;
}
int detect(int c4, int i) {
    // If detect state parameters requested
    if (i==REQUEST) {
        if (state==NONE)  return 0xffffffff;  // No state
        if (state==START) return jstart;      // Report data start
        if (state==END)   return jend;        // Report data end
        if (state==INFO)  return 0xffffffff;  // Report info if any
    }
    if (i==RESET) {
        reset();
        return 0xffffffff;
    }
    buf1=(buf1<<8)|(buf0>>24);
    buf0=c4;
    // Detect header - is first four bytes 0xFFFFFFFF
    if (buf1==0xFFFFFFFF && mystart==0){
        mystart=i;
    }
    // Found possible start, report
    if (type==DEFAULT && mystart) {
        type=YOURTYPE;
        state=START; 
        jstart=mystart-4;
        return state;
    }
    // Found end, report if our type
    if (i-mystart>0x100) {
        if (type==YOURTYPE){
            state=END;
            type=DEFAULT;
            jend=i;
            return state;
      }
      state=NONE;
      type=DEFAULT;
      mystart=0;
    }
    return NONE;
}

int main() {
    reset();
}
```

### Compression model

Below is a model used by archiver itself to compress .det, .enc, .dec and .cfg files.

```c
// Update is called by VM for every input bit
// y    - last bit
// c0   - last 0-7 bits of the partial byte with a leading 1 bit
// c4   - last 4 whole bytes, packed.
// bpos - bit pos in c0 (0 to 7)
// pr   - last prediction
int t[5]={};

int update(int y, int c0, int bpos, int c4) {
    int i;
    if (bpos==0) {
        for (i=4; i>0; --i) t[i]=h2(t[i-1], c4&0xff);
    }
    for (i=1; i<5; ++i) 
    vmx(DS, 0, c0 | (t[i]<<8));
    vmx(APM1, 0, c0);
    return 0;
}

// Called at the start of every new data type
// a - info
// b - reserved (not used and set to 0)
void block(int a, int b) {
}

// Called once at the start of compression
int main() { 
    vmi(DS,   0, 18, 1023, 4, 0);        // pr[0]..pr[3]
    vmi(AVG,  0,   0,   0, 1, 0);        // pr[4]=avg(pr[0],pr[1])
    vmi(AVG,  1,   0,   2, 3, 0);        // pr[5]=avg(pr[2],pr[3])
    vmi(AVG,  2,   0,   4, 5, 0);        // pr[6]=avg(pr[4],pr[5])
    vmi(APM1, 0, 256,   7, 6, 0);        // pr[7]=apm(pr[6]) rate 7 
                                         // pr[7] is final prediction
}

```
