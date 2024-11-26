# Components
Table of available components.
|Name| Short name| ID|Prediction|Mixer input|
| --- | --- | --- | --- | --- |  
|State Context Map| [SMC](#smc)|1|yes|yes|
|Adaptive Probability Map 1| [APM1](#apm1)|2|yes|no|
|Dynamic State Map| [DS](#ds) |3|yes|no|
|Average Map| [AVG](#avg) |4|yes|no|
|Small Stationary Context Map| [SCM](#scm)|5|no|yes|
|Run Context Map |[RCM](#rcm)|6|no|yes|
|Context Map| [CM](#cm)|7|no|yes|
|Mixer |[MX](#mx) |8|yes|no|
|Static Map |[ST](#st)|9|yes|yes|
|Mixer Map| [MM](#mm)|10|no|yes|
|Dynamic Hash State Map| [DHS](#dhs) |11|yes|yes|
|Stationary Map| [SM](#sm) |12|no|yes|
|Direct Map| [SK](#sk) |13|no|yes|
|Adaptive Probability Map 2| [APM2](#apm2)|14|yes|no|
|Error Map| [ERR](#err)|15|no|no|
|Unaligned Sparse Map| [UAS](#uas)|17|yes|no|
|Linear Mixer| [LMX](#lmx)|18|yes|no|
|State Table| [STA](#sta)|19|no|no|
|Byte Map| [BYT](#byt)|20|no|no|

Adaptive Probability Mapping (APM) is the same as Secondary Symbol Estimation (SSE).

__Yes__ in 'Prediction' means that component outputs prediction into internal array pr[]. Last value in the array will be the final prediction.

This array size depends on how many non-mixer prediction components there are.

__Yes__ in 'Mixer input' means that component outputs prediction into internal array inputs[], which is used by MX components as inputs.

If both are yes then only one output can be selected in the vmi function.

## vmi - initialize component
Initialize a given component used in prediction or as an helper.

This function is usable only in function main().

## vmx - set component context
In the update function we set or get data to components with function vmx.

vmx(Component,Index,Data);

set component context. Only usable in function update(...)
```c
// Set APM1(index) context = c0
// first parameter is component ID
// second parameter is component index
// third parameter is context
vmx(APM1, index, c0);
```
Some components output data. In this example BYT component outputs data that was given at initialization in function main().
This is useful in the development and tuning phase as BYT component value changes on tuning mode.
```c
// Get BYT(index) data
// first parameter is component ID
// second parameter is component index
// third parameter is not used
a = vmx(BYT, index, 0);
```

# Individual components

### SMC
```c
// Create SMC component (0) - State Context Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is memory size=x^2
// forth parameter is limit=1...1023, if 0 then defaults to 1023
// fifth parameter is output (pr index=-1,MX input=0...MX)
//                    output=-1 use pr[index] where index=0...lastComponent
//                    output>=0 select MX component as output
// sixth parameter is nil
// vmi(SMC,index,size,limit,output);
```
Prediction to mixer:
```c
// in update
vmx(SMC,0,val1);      //  set component SMC(0) context to val1
vmx(SMC,1,val2);      //  set component SMC(1) context to val2
// in main
vmi(SMC,0,0x100,1023,0,0);  //  mixer[0].add(smc(0).predict())
vmi(SMC,1,0x100,1023,0,0);  //  mixer[0].add(smc(1).predict())
```
Direct prediction:
```c
// in update
vmx(SMC,0,val1);      //  set component SMC(0) context to val1
vmx(SMC,1,val2);      //  set component SMC(1) context to val2
// in main
vmi(SMC,0,0x100,1023,-1,0);  //  pr[0]=smc(0).predict()
vmi(SMC,1,0x100,1023,-1,0);  //  pr[1]=smc(1).predict()
```

### APM1
```c
// Create APM1 component (0) - Adaptive Probability Map 1
// 
// first parameter is component ID
// second parameter is component
// third parameter is size=x^2
// forth parameter is rate
// fifth parameter is predictionIndex
// sixth parameter is nil
//
// vmi(APM1,index,size,rate,predictionIndex,0);

// in update
vmx(SMC,0,val1);          //  set component SMC(0) context to val1
vmx(APM1,0,val2);      //  set component APM1(0) context to val2
// in main
vmi(SMC,0,0x100,1023,-1,0);  //  pr[0]=smc(0).predict()
vmi(APM1,0,0x1000,7,0,0);   //  pr[1]=apm(pr[0])
```

### DS
```c
// Create DS component (0) - Dynamic State Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is number of memory bits x in lower 16 bits, and statetable index y in upper 16 bits
// forth parameter is limit for state map, default 1023
// fifth parameter is number of contexts N
// sixth parameter is nil
//
vmi(DS,0,x+y*0x10000,1023,N,0);
```

### AVG
```c
// Create AVG component (0) - Average Map
// Calculate average and output prediction - pr[3]=(pr[1]+pr[2]+1)>>1
// first parameter is component ID
// second parameter is component index
// third parameter is average parameters x and y where pr[next]=(pr[z]*x+pr[w]*y+1)>>(x+y)
// forth not used, set 0
// fifth parameter is index into pr[] array as z and w
// sixth parameter is nil
//
vmi(AVG,0,x+y*256,0,z+w*256,0);
```

### SCM
```c
// Create SCM component (0) - Small Stationary Context Map
// 
// first parameter is component ID
// second parameter is component
// third parameter is input size in bits
// forth parameter is nil
// fifth parameter is mixer index
// sixth parameter is nil
//
vmi(SCM,0,8,0,0,0);  // input is 8 bits, use mixer 0
```

### RCM
```c
// Create RCM component (0) - Run Context Map
// 
// first parameter is component ID
// second parameter is component
// third parameter is memory*4096, must be power of two
// forth parameter is unused
// fifth parameter is predictionIndex
vmi(RCM,0,memory,0,0,0);
```

### CM
```c
// Create CM component (0) - Context Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is:
//                   memory in lower 24 bits (must be power of two)
//                   and statetable index si in upper 8 bits (si=0 is default STA)
// forth parameter is:
//                   number of contexts x (max 64) 
//                   run mul y (default 4)
//                   mixer prediction mul z (default 32) and w (default 12)
// fifth parameter is: mixer index mi
// sixth parameter is:
//                    v (default 8) and u (default 32). 
//                    r is bit 28, set to 1 to disable random state updates
//
// Parameters y, z, w, v, u are tunable.

// update
cmx(CM,0,cxt);
//main
//
vmi(CM, 0, memory*4096+(si<<24), x+y*0x100+z*0x10000+w*0x1000000, mi, v*0x100+u*0x10000+|(r<<28));
```

CM components can be used in two modes. Single mode and Dual mode. 

### MX
```c
// Create MX component (0) - Mixer
// 
// first parameter is component ID
// second parameter is component index
// third parameter is shift (default 64), error (default 0), mul (default 28), all tunable
// forth parameter is context size
// fifth parameter is mixer index
// sixth parameter is nil
//
// -Internal overview-
// Update:
//  err=((y<<12)-pr)*mul/4;
//    if (err>=-error && err<=error) err=0;
//    train(..., err);
// Predict:
//  dot_product(...)*shift>>11;

vmi(MX,0,shift+256*error+0x1000000*mul,1,0,0,0);
```

### ST
```c
// Create ST component (0) - Static Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is value m where pr=((m-128)*16) if fift parameter is >=0, or pr=(m*16) if -1
// forth parameter is nil
// fifth parameter is output (pr index=-1, MX input=0...MX)
//                    output=-1 use pr[index] where index=0...lastComponent
//                    output>=0 select MX component as output
// sixth parameter is nil
//
//main
//
vmi(ST,0,144,0,0,0); // set probability that next bit is 1 to 56%
```

### MM
```c
// Create MM component (0) - Mixer Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is option:
//         0 adds stretch(pr) to mixer
//         1 adds stretch(pr) >> 1  to mixer
//         3...x  (pr-2048 >> 3...x) to mixer
// forth parameter is pr index
// fifth parameter is mixer index
// sixth parameter is nil
//
```
```c
// mixer[2].add(stretch(pr[1]))
vmi(MM,0,0,1,2,0);
```

### DHS
```c
// Create DHS component (0) - Dynamic Hash State Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is input bits in lower 16 bits, and state table index in upper 16 bits
// forth parameter is memory bits, memory usage is ((1<<bits)*(1<<memory))
// fifth parameter is m - mixer index or -1 if no mixer is used
// sixth parameter is n - number of contexts
//
//main
// creates DHS with 10 context using 256MB of memory where state count per context is 16 (1<<4)
vmi(DHS,0,4,24,m,n);

// update
// set DHS contexts at the start of state update
for ( i=0; i<10; i++) {  vmx(DHS,0,cxt[i]);}

// update DHS contexts states where j is in range 0...16
vmx(DHS,0,j);
```

### SM
```c
// Create SM component (0) - Stationary Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is memory_bits
// forth parameter is input_bits (low 8 bits), memory usage is N=((1<<memory_bits)*((1<<input_bits)-1)). pr mul value
// fifth parameter is n - number of contexts
// sixth parameter is nil
//

// update
// set SM contexts at the start of state update
vmx(SM, 0,val);

//main
// creates SM
vmi(SM,0,16,3,n,0);
```
### SK
```c
// Create SK component (0) - Direct Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is nil
// forth parameter is nil
// fifth parameter is mixer index
// sixth parameter is nil
//
//main
// creates SK
vmi(SK,0,0,0,0,0);

// update
// set SK value to be added to mixer, range -2047..2047
vmx(SK, 0,val);
```

### APM2
```c
// Create APM2 component (0) - Adaptive Probability Map 2
// 
// first parameter is component ID
// second parameter is component index
// third parameter is size
// forth parameter is x
// fifth parameter is predictionIndex
// sixth parameter is nil
//
//
// vmi(APM1,index,size,rate,predictionIndex);

// in update
vmx(SMC,0,val1);          //  set component SMC(0) context to val1
vmx(APM2,0,val2);      //  set component APM1(0) context to val2
// in main
vmi(SMC,0,0x10,1023,-1,0);  //  pr[0]=smc(0).predict()
vmi(APM2,0,0x1000,7,0,0);   //  pr[1]=apm(pr[0])
```

### ERR
```c
// Create ERR component (0) - Error Map
//     map last final prediction to 0,1,3 with user provided ranges x and y
// first parameter is component ID
// second parameter is component index
// third parameter is x and y, where x is first input (range 0-2047) and y is second input (range 0-4094).
//                    x and y parameters are tunable.
// forth parameter is nil
// fifth parameter is nil
// sixth parameter is nil
//
// output: based on input thresholds output values are 1 or 3. 
//         val=y?finalprediction^4095:finalprediction;
//         0 
//         1 - val>x (for low)
//         3 - val>y (for high)

// update
// get mapped ERR value at bpos
a=vmx(ERR, bpos, 0);

//main
// creates ERR map for every bit in byte
for (i=0; i<8; i++)  
   vmi(ERR,i,x+(y<<16),0,0,0);

```

### UAS
```c
// Create UAS component (0) - Unaligned Sparse Map
// 
// first parameter is component ID
// second parameter is component index
// third parameter is x, input size in bits. x parameter is tunable.
// forth parameter is y, input mask - ignored
// fifth parameter is z, update rate, default 5
// sixth parameter is nil
//
// output: prediction when last 8 bits of val are set, otherwise prediction is 2048

// update
// set UAS value, where val is shifted output of component ERR outputs
a=vmx(ERR, bpos, 0); // get ERR component value for every bit pos
erra=(erra<<1)|(a&1);
if (bpos==0){ 
    val=erra;
}
vmx(UAS, 0,val);

//main
// create UAS and ERR
vmi(UAS,0,x,y,z,0);
for (i=0; i<8; i++) vmi(ERR,i,e_low[i]+(e_high[i]<<16),0,0);

```

### LMX
```c
// Create LMX component (0) - Linear Mixer
// 
// first parameter is component ID
// second parameter is component index
// third parameter is x and y, where x is first input and y is second input 
// fourth parameter is weight z, if z==0 default value is 2048. z parameter is tunable.
// fifth parameter is nil
// sixth parameter is nil
//
// output: x+(((y-x)*z)>>12);
// if z > 2048 then we favor input y

// update
// :none

//main
// creates LMX
vmi(LMX,0,x+y*256,z,0,0);

```

### STA
```c
// Create STA component (0) - State Table
// 
// first parameter is component ID
// second parameter is component index
// third parameter is u and v (in range 1-63)
// forth parameter is w and x (in range 1-63)
// fifth parameter is y, z and u (y range 1-63) (z range 1-32) (u range 2-32)
// sixth parameter is nil
//
// output: generates statetable;
//         if parameters==0 use default values as 42,41,13,6,5,16,14
// In every component that use STA as input needs to add index+1 if user defined STA is used

// update
// :none

//main
// creates STA
vmi(STA,0,u+v*0x10000,w+x*0x10000,y+z*0x10000+u*0x1000000,0);
// use STA(0) (1<<24) in CM(0), with 2*4096 memory, one context, output to MX(0)
vmi(CM,0,2*4096+(1<<24),1,0,0);
// create MX(0) with default parameters, context size is 1
vmi(MX,0,0,1,0,0);
```

### BYT
```c
// Create BYT component (0) - Byte Map
// Map byte value to new value in tune mode, otherwise output is the same as input. 
// Useful in model creation phase. Grouping chars, automatic "best" context selection for CM.
//
// first parameter is component ID
// second parameter is component index
// third parameter is value x in range 0-255, tunable
// forth parameter is value y of output range, must be between value x
// fifth parameter is nil
// sixth parameter is nil
//

//main
// creates BYT with value x
vmi(BYT,0,x,y,0,0);

// update
// get BYT value based on input val x in range 0-y
a=vmx(BYT, 0, 0);
```
