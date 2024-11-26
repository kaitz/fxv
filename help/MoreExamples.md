# More Examples and Tuning
[Smallest compression configuration and model](#smallest-compression-configuration-and-model)

[fpaq0](#fpaq0)

[Using DS and STA component 1](#using-ds-and-sta-component-1)

[Using DS and STA component 2](#using-ds-and-sta-component-2)


## Smallest compression configuration and model
### Model

```c
int update(int y, int c0, int bpos, int c4) { 
}

void block(int a, int b) { 
}

int main() { 
    vmi(ST,0,128,0,-1,0);
}
```
### Config

```
stream 0
model mymodel.cfg
type 0
detect -1
decode -1
encode -1
compress 0
```
When compressing enwik8 we have a result: Data 100000000->100000004

So no compression. Parameter 128 means that our probability is 50% for the next bit being 1.

How to improve this result? We can tune our parameters.

Run our program with parameters:
```
fxv.exe -2 -j -o20 -r140 -f -k -cmymodel.pxv
```
Output from compressor is:
```
 SA InitState: 20000194 bytes cooling rate=0.90000, seed 1729969538, time: 52
 Max runs: 140
 Full tune
[8] [rate: 25.0%], time: 52, temperature: 2431410.246011907, proposal: 19925677
 best: 19925677 (radius: 0.05000)
[62] [rate: 14.7%], time: 52, temperature: 16087.567494501, proposal: 1992467088
 best: 19924670 (radius: 0.01000)
[140] [rate: 19.6%], time: 52, temperature: 8148.615455965, proposal: 3530820977

  best value found: 19924670 bytes


vmST(0): 118
Compressing mymodel.cfg   stream(0).  Total 100000000
S[0] compressed 100000249->99647509 bytes. Data 100000000->99647319, model 249->190.
```

Our new parameter is **118** and data was compressed from 100000000 to 99647319.

[Top](#more-examples-and-tuning)

## fpaq0 
Source: https://mattmahoney.net/dc/fpaq0.cpp

We can use **ST** component to predict the next bit at every update. For that lets use fpaq0 as an example.
```c
int cxt;  // Context: last 0-8 bits with a leading 1
int ct0[256]={};
int ct1[256]={};

int update(int y, int c0, int bpos, int c4) { 
    if (y) ct1[cxt]=ct1[cxt]+1;
    else ct0[cxt]=ct0[cxt]+1;
    
    if (ct0[cxt] > 65534 || ct1[cxt] > 65534) {
      ct0[cxt]=ct0[cxt]>>1;
      ct1[cxt]=ct1[cxt]>>1;
    }
    cxt=cxt+cxt+y;
    if (cxt>= 256)    cxt=1;
    vmx(ST,0, (4096*(ct1[cxt]+1))/(ct0[cxt]+ct1[cxt]+2));
}

void block(int a, int b) { 
}

int main() { 
    cxt=1;
    vmi(ST,0,128,0,-1,0);
}

```
Data 100000000->63375532

fpaq0.c is 63391013

Can we somehow tune something? Yes.
Let's add **BYT** component and set its current value and max count to **65534**.
 
```c
int cxt;  // Context: last 0-8 bits with a leading 1
int ct0[256]={};
int ct1[256]={};

int update(int y, int c0, int bpos, int c4) {
    int a;
    if (y) ct1[cxt]=ct1[cxt]+1;
    else ct0[cxt]=ct0[cxt]+1;
    a=vmx(BYT,0,0);
    if (ct0[cxt] > a || ct1[cxt] > a) {
      ct0[cxt]=ct0[cxt]>>1;
      ct1[cxt]=ct1[cxt]>>1;
    }
    cxt=cxt+cxt+y;
    if (cxt>= 256)    cxt=1;
    vmx(ST,0, (4096*(ct1[cxt]+1))/(ct0[cxt]+ct1[cxt]+2));
}

void block(int a, int b) { 
}

int main() { 
    cxt=1;
    vmi(ST,0,128,0,-1,0);
    vmi(BYT,0,65534,65534,0,0);
}
```

Run our program with parameters:
```
fxv.exe -2 -pr -o10 -r70 -j -k -cmymodel.pxv
```
Output from compressor is:
```
 SA InitState: 6367079 bytes cooling rate=0.90000, seed 1729981019, time: 45
 Max runs: 70
[1] [rate: 100.0%], time: 45, temperature: 1492.887281542, proposal: 6360204
 best: 6360204 (radius: 0.05000)
[2] [rate: 100.0%], time: 45, temperature: 1129.382800189, proposal: 6356677
 best: 6356677 (radius: 0.05000)
[3] [rate: 100.0%], time: 45, temperature: 1016.444520170, proposal: 6274692
 best: 6274692 (radius: 0.05000)
[11] [rate: 36.4%], time: 45, temperature: 437.546036718, proposal: 6226236
 best: 6226236 (radius: 0.05000)
[15] [rate: 100.0%], time: 45, temperature: 8379.712028323, proposal: 6187646
 best: 6187646 (radius: 0.01000)
[43] [rate: 100.0%], time: 45, temperature: 16.068895830, proposal: 618757240
 best: 6187572 (radius: 0.00500)
[70] [rate:  3.6%], time: 45, temperature: 1156.981970124, proposal: 63180594

  best value found: 6187572 bytes

BYT(0): 111
Compressing tests\test4.cfg   stream(0).  Total 100000000
S[0] compressed 100000666->61612525 bytes. Data 100000000->61612145, model 666->380.
```
So improvement was from  63375532 to 61612145. Our new parameter value is 111 as max count.
To replicate above use seed command line option -s1729981019

[Top](#more-examples-and-tuning)

## Using DS and STA component 1
At first we used only the DS component. It uses the default STA component.
```c
int update(int y, int c0, int bpos, int c4) {
    vmx(DS,0,c0);
    return 0;
}

void block(int a, int b) { 
}

int main() { 
    vmi(DS,0,8,1023,1,0);
}
```
Result: Data 100000000->63344331

Add STA component
```c
int update(int y, int c0, int bpos, int c4) {
    vmx(DS,0,c0);
    return 0;
}

void block(int a, int b) { 
}

int main() { 
    vmi(STA,0,0,0,0,0);
    vmi(DS,0,8+(1<<16),1023,1,0);
}
```
Tune:
```
fxv.exe     -2   -pq  -o10 -r70   -j      -k          -ctests.pxv
```
```
 SA InitState: 6353186 bytes cooling rate=0.90000, seed 1730587785, time: 40
 Max runs: 70
[3] [rate: 66.7%], time: 40, temperature: 9289.591539997, proposal: 62060756
 best: 6206075 (radius: 0.05000)
[5] [rate: 60.0%], time: 40, temperature: 7524.569147397, proposal: 6205450
 best: 6205450 (radius: 0.05000)
[12] [rate: 50.0%], time: 40, temperature: 3598.978097036, proposal: 6195365
 best: 6195365 (radius: 0.05000)
[13] [rate: 53.8%], time: 40, temperature: 3239.080287332, proposal: 6182213
 best: 6182213 (radius: 0.05000)
[20] [rate: 16.7%], time: 40, temperature: 3465.091517910, proposal: 61637497
 best: 6163749 (radius: 0.01000)
[29] [rate: 20.0%], time: 40, temperature: 1342.447450298, proposal: 6161757
 best: 6161757 (radius: 0.01000)
[30] [rate: 25.0%], time: 40, temperature: 1208.202705269, proposal: 6161457
 best: 6161457 (radius: 0.01000)
[32] [rate: 27.8%], time: 40, temperature: 978.644191268, proposal: 61605642
 best: 6160564 (radius: 0.01000)
[40] [rate: 26.9%], time: 40, temperature: 421.274234598, proposal: 6154750
 best: 6154750 (radius: 0.01000)
[43] [rate: 100.0%], time: 40, temperature: 211.718559928, proposal: 6153775
 best: 6153775 (radius: 0.00500)
[57] [rate: 46.7%], time: 40, temperature: 4058.804022530, proposal: 61529423
 best: 6152942 (radius: 0.00500)
[58] [rate: 50.0%], time: 40, temperature: 3652.923620277, proposal: 6151667
 best: 6151667 (radius: 0.00500)
[70] [rate: 32.1%], time: 40, temperature: 1031.693524875, proposal: 6320573

  best value found: 6151667 bytes


STA(0): 2+(48<<16),13|(6<<16),45+(10<<16)+(24<<24)  2 48 13 6 45 10 24
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000220->61295927 bytes. Data 100000000->61295771, model 220->156.
```
Lets change STA component parameters and tune again.
```c
vmi(STA,0,2+(48<<16),13|(6<<16),45+(10<<16)+(24<<24),0);
```
This time replace -k with -f option. So that all parameters are changed at every tune run.
```
fxv.exe     -2   -pq  -o10 -r70   -j      -f          -ctests.pxv
```
```
 SA InitState: 6151692 bytes cooling rate=0.90000, seed 1730588954, time: 40
 Max runs: 70
 Full tune
[4] [rate: 25.0%], time: 40, temperature: 1246.966945429, proposal: 6148191
 best: 6148191 (radius: 0.05000)
[16] [rate: 50.0%], time: 40, temperature: 7897.102285308, proposal: 61481249
 best: 6148124 (radius: 0.01000)
[22] [rate: 62.5%], time: 40, temperature: 4196.843935607, proposal: 6147244
 best: 6147244 (radius: 0.01000)
[26] [rate: 50.0%], time: 40, temperature: 2753.549306151, proposal: 6145087
 best: 6145087 (radius: 0.01000)
[43] [rate: 100.0%], time: 40, temperature: 0.000000000, proposal: 614508774

  best value found: 6145087 bytes

STA(0): 3+(51<<16),13|(5<<16),47+(8<<16)+(23<<24)  3 51 13 5 47 8 23
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000257->61240831 bytes. Data 100000000->61240650, model 257->181.
```
For the final STA tuning round lets remove -f option.
```
fxv.exe     -2   -pq  -o10 -r70   -j   -ctests.pxv
```
```
 SA InitState: 6145087 bytes cooling rate=0.90000, seed 1730589418, time: 41
 Max runs: 70
[4] [rate: 25.0%], time: 40, temperature: 62.264799870, proposal: 6145056
 best: 6145056 (radius: 0.05000)
[5] [rate: 40.0%], time: 40, temperature: 56.038319883, proposal: 6145019
 best: 6145019 (radius: 0.05000)
[6] [rate: 50.0%], time: 40, temperature: 50.434487895, proposal: 6143359
 best: 6143359 (radius: 0.05000)
[10] [rate: 50.0%], time: 40, temperature: 33.090067508, proposal: 6143015
 best: 6143015 (radius: 0.05000)
[13] [rate: 53.8%], time: 40, temperature: 24.122659213, proposal: 6142037
 best: 6142037 (radius: 0.05000)
[16] [rate: 100.0%], time: 40, temperature: 0.000000000, proposal: 6142037

  best value found: 6142037 bytes

STA(0): 7+(55<<16),13|(5<<16),45+(9<<16)+(23<<24)  7 55 13 5 45 9
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000256->61211869 bytes. Data 100000000->61211688, model 256->181.
```

Now let's tune DS parameters.
```
fxv.exe     -2   -pc  -o10 -r70   -j   -ctests.pxv
```
```
 SA InitState: 6142037 bytes cooling rate=0.90000, seed 1730589708, time: 40
 Max runs: 70
[1] [rate: 100.0%], time: 40, temperature: 142.448590064, proposal: 6141381
 best: 6141381 (radius: 0.05000)
[43] [rate: 100.0%], time: 40, temperature: 0.000000000, proposal: 61413817

  best value found: 6141381 bytes

vmDS(0): 1020
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000256->61205182 bytes. Data 100000000->61205001, model 256->181.
```

Let's do our final "final" STA tuning round.
```
fxv.exe     -2   -pq  -o10 -r70   -j   -ctests.pxv
```
```
 SA InitState: 6141381 bytes cooling rate=0.90000, seed 1730590283, time: 40
 Max runs: 70
[1] [rate: 100.0%], time: 40, temperature: 924.395804731, proposal: 6137124
 best: 6137124 (radius: 0.05000)
[3] [rate: 100.0%], time: 40, temperature: 415.978112129, proposal: 6134080
 best: 6134080 (radius: 0.05000)
[6] [rate: 83.3%], time: 40, temperature: 303.248043742, proposal: 61338759
 best: 6133875 (radius: 0.05000)
[11] [rate: 72.7%], time: 40, temperature: 179.064937349, proposal: 6131992
 best: 6131992 (radius: 0.05000)
[14] [rate: 71.4%], time: 40, temperature: 130.538339328, proposal: 6131706
 best: 6131706 (radius: 0.05000)
[19] [rate: 40.0%], time: 40, temperature: 44.324094823, proposal: 61315966
 best: 6131596 (radius: 0.01000)
[24] [rate: 60.0%], time: 40, temperature: 26.172934752, proposal: 6128987
 best: 6128987 (radius: 0.01000)
[43] [rate: 100.0%], time: 40, temperature: 0.000000000, proposal: 6128987

  best value found: 6128987 bytes

STA(0): 9+(58<<16),13|(5<<16),42+(13<<16)+(20<<24)  9 58 13 5 42 13 20
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000256->61096399 bytes. Data 100000000->61096218, model 256->181.
```
Let's do our super final "final" STA tuning round.
```
fxv.exe     -2   -pq  -o10 -r70   -j   -ctests.pxv
```
```
 SA InitState: 6128988 bytes cooling rate=0.90000, seed 1730590833, time: 41
 Max runs: 70
[2] [rate: 50.0%], time: 40, temperature: 191.741013760, proposal: 6128455
 best: 6128455 (radius: 0.05000)
[8] [rate: 37.5%], time: 40, temperature: 101.899036094, proposal: 6127871
 best: 6127871 (radius: 0.05000)
[10] [rate: 50.0%], time: 40, temperature: 82.538219236, proposal: 6127719
 best: 6127719 (radius: 0.05000)
[43] [rate: 100.0%], time: 40, temperature: 0.000000000, proposal: 61277199

  best value found: 6127719 bytes

STA(0): 11+(58<<16),13|(5<<16),39+(13<<16)+(20<<24)  11 58 13 5 39 13 20
Compressing tests\test8.cfg   stream(0).  Total 100000000
S[0] compressed 100000257->61079953 bytes. Data 100000000->61079771, model 257->182.
```

Initial result: Data 100000000->63344331

Final result:   Data 100000000->61079771


Our final model with updated parameters is:
```c
int update(int y, int c0, int bpos, int c4) {
    vmx(DS,0,c0);
    return 0;
}

void block(int a, int b) { 
}

int main() { 
    vmi(STA,0,11+(58<<16),13|(5<<16),39+(13<<16)+(20<<24),0);
    vmi(DS,0,8+(1<<16),1020,1,0);
}
```

This model is similar to [fpaqX](http://www.mattmahoney.net/dc/text.html#5586)

Our final State table (parameters 11 58 13 5 39 13 20):
```
{  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7,  9, 2, 0}, // 0-3
{  8, 11, 1, 1},{  8, 11, 1, 1},{ 10, 12, 0, 2},{ 13, 15, 3, 0}, // 4-7
{ 14, 17, 2, 1},{ 14, 17, 2, 1},{ 16, 19, 1, 2},{ 16, 19, 1, 2}, // 8-11
{ 18, 20, 0, 3},{ 21, 23, 4, 0},{ 22, 25, 3, 1},{ 22, 25, 3, 1}, // 12-15
{ 24, 27, 2, 2},{ 24, 27, 2, 2},{ 26, 29, 1, 3},{ 26, 29, 1, 3}, // 16-19
{ 28, 30, 0, 4},{ 31, 33, 5, 0},{ 32, 35, 4, 1},{ 32, 35, 4, 1}, // 20-23
{ 34, 37, 3, 2},{ 34, 37, 3, 2},{ 36, 39, 2, 3},{ 36, 39, 2, 3}, // 24-27
{ 38, 41, 1, 4},{ 38, 41, 1, 4},{ 40, 42, 0, 5},{ 43, 45, 6, 0}, // 28-31
{ 44, 47, 5, 1},{ 44, 47, 5, 1},{ 46, 49, 4, 2},{ 46, 49, 4, 2}, // 32-35
{ 48, 51, 3, 3},{ 48, 51, 3, 3},{ 50, 53, 2, 4},{ 50, 53, 2, 4}, // 36-39
{ 52, 55, 1, 5},{ 52, 55, 1, 5},{ 54, 56, 0, 6},{ 57, 59, 7, 0}, // 40-43
{ 58, 61, 6, 1},{ 58, 61, 6, 1},{ 60, 25, 5, 2},{ 60, 25, 5, 2}, // 44-47
{ 24, 63, 4, 3},{ 24, 63, 4, 3},{ 62, 27, 3, 4},{ 62, 27, 3, 4}, // 48-51
{ 26, 65, 2, 5},{ 26, 65, 2, 5},{ 64, 67, 1, 6},{ 64, 67, 1, 6}, // 52-55
{ 66, 68, 0, 7},{ 69, 71, 8, 0},{ 70, 73, 7, 1},{ 70, 73, 7, 1}, // 56-59
{ 72, 35, 6, 2},{ 72, 35, 6, 2},{ 74, 77, 4, 4},{ 74, 77, 4, 4}, // 60-63
{ 38, 79, 2, 6},{ 38, 79, 2, 6},{ 78, 81, 1, 7},{ 78, 81, 1, 7}, // 64-67
{ 80, 82, 0, 8},{ 83, 85, 9, 0},{ 84, 87, 8, 1},{ 84, 87, 8, 1}, // 68-71
{ 86, 47, 7, 2},{ 86, 47, 7, 2},{ 88, 63, 5, 4},{ 88, 63, 5, 4}, // 72-75
{ 62, 91, 4, 5},{ 62, 91, 4, 5},{ 52, 93, 2, 7},{ 52, 93, 2, 7}, // 76-79
{ 92, 95, 1, 8},{ 92, 95, 1, 8},{ 94, 96, 0, 9},{ 83, 98,10, 0}, // 80-83
{ 97,100, 9, 1},{ 97,100, 9, 1},{ 99, 47, 8, 2},{ 99, 47, 8, 2}, // 84-87
{101, 75, 6, 4},{101, 75, 6, 4},{ 76,104, 4, 6},{ 76,104, 4, 6}, // 88-91
{ 52,106, 2, 8},{ 52,106, 2, 8},{105,108, 1, 9},{105,108, 1, 9}, // 92-95
{107, 96, 0,10},{109,112,10, 1},{109,112,10, 1},{111, 61, 9, 2}, // 96-99
{111, 61, 9, 2},{113, 89, 7, 4},{113, 89, 7, 4},{ 90,116, 4, 7}, // 100-103
{ 90,116, 4, 7},{ 64,118, 2, 9},{ 64,118, 2, 9},{117,120, 1,10}, // 104-107
{117,120, 1,10},{121,122,11, 1},{121,122,11, 1},{122, 73,10, 2}, // 108-111
{122, 73,10, 2},{123, 89, 8, 4},{123, 89, 8, 4},{ 90,124, 4, 8}, // 112-115
{ 90,124, 4, 8},{ 78,125, 2,10},{ 78,125, 2,10},{125,126, 1,11}, // 116-119
{125,126, 1,11},{127,128,12, 1},{128, 73,11, 2},{129,102, 9, 4}, // 120-123
{103,130, 4, 9},{ 78,131, 2,11},{131,132, 1,12},{133, 59,13, 1}, // 124-127
{ 58, 87,12, 2},{134,114,10, 4},{115,135, 4,10},{ 92, 67, 2,12}, // 128-131
{ 66,136, 1,13},{137, 59,14, 1},{138,123,11, 4},{124,139, 4,11}, // 132-135
{ 66,140, 1,14},{141, 71,15, 1},{142,129,12, 4},{130,143, 4,12}, // 136-139
{ 80,144, 1,15},{145, 71,16, 1},{146,129,13, 4},{130,147, 4,13}, // 140-143
{ 80,148, 1,16},{149, 85,17, 1},{150,134,14, 4},{135,151, 4,14}, // 144-147
{ 94,152, 1,17},{153, 85,18, 1},{154,138,15, 4},{139,155, 4,15}, // 148-151
{ 94,156, 1,18},{157, 98,19, 1},{158,142,16, 4},{143,159, 4,16}, // 152-155
{107,160, 1,19},{161, 98,20, 1},{162,146,17, 4},{147,163, 4,17}, // 156-159
{107,164, 1,20},{165, 98,21, 1},{166,146,18, 4},{147,167, 4,18}, // 160-163
{107,168, 1,21},{169, 98,22, 1},{170,150,19, 4},{151,171, 4,19}, // 164-167
{107,172, 1,22},{173, 98,23, 1},{174,150,20, 4},{151,175, 4,20}, // 168-171
{107,176, 1,23},{177, 98,24, 1},{178,150,21, 4},{151,179, 4,21}, // 172-175
{107,180, 1,24},{181, 98,25, 1},{182,150,22, 4},{151,183, 4,22}, // 176-179
{107,184, 1,25},{185, 98,26, 1},{186,150,23, 4},{151,187, 4,23}, // 180-183
{107,188, 1,26},{189, 98,27, 1},{190,150,24, 4},{151,191, 4,24}, // 184-187
{107,192, 1,27},{193, 98,28, 1},{194,150,25, 4},{151,195, 4,25}, // 188-191
{107,196, 1,28},{197, 98,29, 1},{198,150,26, 4},{151,199, 4,26}, // 192-195
{107,200, 1,29},{201, 98,30, 1},{202,150,27, 4},{151,203, 4,27}, // 196-199
{107,204, 1,30},{205, 98,31, 1},{206,150,28, 4},{151,207, 4,28}, // 200-203
{107,208, 1,31},{209, 98,32, 1},{210,150,29, 4},{151,211, 4,29}, // 204-207
{107,212, 1,32},{213, 98,33, 1},{214,150,30, 4},{151,215, 4,30}, // 208-211
{107,216, 1,33},{217, 98,34, 1},{218,150,31, 4},{151,219, 4,31}, // 212-215
{107,220, 1,34},{221, 98,35, 1},{222,150,32, 4},{151,223, 4,32}, // 216-219
{107,224, 1,35},{225, 98,36, 1},{226,150,33, 4},{151,227, 4,33}, // 220-223
{107,228, 1,36},{229, 98,37, 1},{230,150,34, 4},{151,231, 4,34}, // 224-227
{107,232, 1,37},{233, 98,38, 1},{234,150,35, 4},{151,235, 4,35}, // 228-231
{107,236, 1,38},{237, 98,39, 1},{238,150,36, 4},{151,239, 4,36}, // 232-235
{107,240, 1,39},{241, 98,40, 1},{242,150,37, 4},{151,243, 4,37}, // 236-239
{107,244, 1,40},{245, 98,41, 1},{ 97,150,38, 4},{151,108, 4,38}, // 240-243
{107,246, 1,41},{247, 98,42, 1},{107,248, 1,42},{249, 98,43, 1}, // 244-247
{107,250, 1,43},{251, 98,44, 1},{107,252, 1,44},{253, 98,45, 1}, // 248-251
{107,254, 1,45},{255, 98,46, 1},{107,  0, 1,46},{  1, 98,47, 1}, // 252-255
```
[Top](#more-examples-and-tuning)

## Using DS and STA component 2
```c
int update(int y, int c0, int bpos, int c4) {
    vmx(DS,0,(c4&0xff)*256+c0);
    return 0;
}

void block(int a,int b) { 
}

int main() { 
    vmi(DS,0,16,1023,1,0);
}
```
Data 100000000->47904276
```c
int update(int y, int c0, int bpos ,int c4) {
    vmx(DS,0,(c4&0xff)*256+c0);
    return 0;
}

void block(int a, int b) { 
}

int main() { 
    vmi(STA,0,0,0,0,0);
    vmi(DS,0,16+(1<<16),1023,1,0);
}
```

```
 SA InitState: 9611994 bytes cooling rate=0.90000, seed 1730628911, time: 90
 Max runs: 70
[5] [rate: 40.0%], time: 85, temperature: 31806.337243663, proposal: 9577081
 best: 9577081 (radius: 0.05000)
[7] [rate: 57.1%], time: 85, temperature: 25763.133167367, proposal: 9395709
 best: 9395709 (radius: 0.05000)
[8] [rate: 62.5%], time: 85, temperature: 23186.819850630, proposal: 9375965
 best: 9375965 (radius: 0.05000)
[9] [rate: 66.7%], time: 85, temperature: 20868.137865567, proposal: 9302909
 best: 9302909 (radius: 0.05000)
[10] [rate: 70.0%], time: 85, temperature: 18781.324079011, proposal: 9302851
 best: 9302851 (radius: 0.05000)
[11] [rate: 72.7%], time: 85, temperature: 16903.191671110, proposal: 9248615
 best: 9248615 (radius: 0.05000)
[13] [rate: 69.2%], time: 85, temperature: 13691.585253599, proposal: 9248594
 best: 9248594 (radius: 0.05000)
[14] [rate: 71.4%], time: 85, temperature: 12322.426728239, proposal: 9239650
 best: 9239650 (radius: 0.05000)
[18] [rate: 50.0%], time: 83, temperature: 42317.022418182, proposal: 9224540
 best: 9224540 (radius: 0.01000)
[20] [rate: 50.0%], time: 83, temperature: 34276.788158727, proposal: 9223597
 best: 9223597 (radius: 0.01000)
[21] [rate: 57.1%], time: 83, temperature: 30849.109342854, proposal: 9217497
 best: 9217497 (radius: 0.01000)
[24] [rate: 70.0%], time: 83, temperature: 22489.000710941, proposal: 9210083
 best: 9210083 (radius: 0.01000)
[49] [rate: 57.1%], time: 83, temperature: 24641.564186982, proposal: 9203591
 best: 9203591 (radius: 0.00500)
[50] [rate: 62.5%], time: 83, temperature: 22177.407768284, proposal: 9202276
 best: 9202276 (radius: 0.00500)
[54] [rate: 58.3%], time: 82, temperature: 14550.597236771, proposal: 91898910
 best: 9189891 (radius: 0.00500)
[59] [rate: 64.7%], time: 82, temperature: 8591.982162341, proposal: 91856397
 best: 9185639 (radius: 0.00500)
[70] [rate: 42.9%], time: 82, temperature: 2696.255043959, proposal: 9393318

  best value found: 9185639 bytes

STA(0): 2+(36<<16),62|(6<<16),32+(4<<16)+(29<<24)  2 36 62 6 32 4 29
Compressing tests\test9.cfg   stream(0).  Total 100000000
S[0] compressed 100000235->45740697 bytes. Data 100000000->45740528, model 235->169.
```

```
 SA InitState: 9185664 bytes cooling rate=0.90000, seed 1730629656, time: 84
 Max runs: 70
 Full tune
[5] [rate: 20.0%], time: 84, temperature: 1630.730938641, proposal: 9178456
 best: 9178456 (radius: 0.05000)
[12] [rate: 16.7%], time: 84, temperature: 779.973552686, proposal: 9176687
 best: 9176687 (radius: 0.05000)
[19] [rate: 20.0%], time: 84, temperature: 681.957858920, proposal: 91764247
 best: 9176424 (radius: 0.01000)
[26] [rate: 16.7%], time: 84, temperature: 326.178329852, proposal: 9174840
 best: 9174840 (radius: 0.01000)
[32] [rate: 16.7%], time: 84, temperature: 173.344537795, proposal: 9174549
 best: 9174549 (radius: 0.01000)
[41] [rate: 14.8%], time: 83, temperature: 67.157225598, proposal: 91713523
 best: 9171352 (radius: 0.01000)
[43] [rate: 100.0%], time: 83, temperature: 0.000000000, proposal: 9171352

  best value found: 9171352 bytes

STA(0): 4+(35<<16),61|(5<<16),34+(4<<16)+(23<<24)  4 35 61 5 34 4 23
Compressing tests\test9.cfg   stream(0).  Total 100000000
S[0] compressed 100000271->45682906 bytes. Data 100000000->45682712, model 271->194.
```

```
 SA InitState: 9171351 bytes cooling rate=0.90000, seed 1730630201, time: 83
 Max runs: 70
[11] [rate: 72.7%], time: 83, temperature: 1578.900983871, proposal: 9170739
 best: 9170739 (radius: 0.05000)
[12] [rate: 75.0%], time: 83, temperature: 1421.010885484, proposal: 9168804
 best: 9168804 (radius: 0.05000)
[15] [rate: 100.0%], time: 83, temperature: 170.677731388, proposal: 9168793
 best: 9168793 (radius: 0.01000)
[20] [rate: 83.3%], time: 83, temperature: 55.990829782, proposal: 91687853
 best: 9168785 (radius: 0.01000)
[21] [rate: 85.7%], time: 83, temperature: 50.391746804, proposal: 9168411
 best: 9168411 (radius: 0.01000)
[33] [rate: 68.4%], time: 83, temperature: 14.232117692, proposal: 9167441
 best: 9167441 (radius: 0.01000)
[43] [rate: 100.0%], time: 83, temperature: 0.000000000, proposal: 9167441

  best value found: 9167441 bytes

STA(0): 4+(34<<16),58|(5<<16),39+(2<<16)+(24<<24)  4 34 58 5 39 2 24
Compressing tests\test9.cfg   stream(0).  Total 100000000
S[0] compressed 100000271->45663820 bytes. Data 100000000->45663627, model 271->193.
```
DS
```
 SA InitState: 9167442 bytes cooling rate=0.90000, seed 1730630759, time: 85
 Max runs: 70
[70] [rate: 10.7%], time: 85, temperature: 16.793983087, proposal: 91686399

  best value found: 9167442 bytes

Tune failed. No improvement. Exit.
Compressing tests\test9.cfg   stream(0).  Total 100000000
S[0] compressed 100000271->45663821 bytes. Data 100000000->45663627, model 271->194.
```

```
 SA InitState: 9167442 bytes cooling rate=0.90000, seed 1730631521, time: 84
 Max runs: 70
[4] [rate: 25.0%], time: 82, temperature: 346.325963121, proposal: 9167432
 best: 9167432 (radius: 0.05000)
[16] [rate: 100.0%], time: 82, temperature: 0.000000000, proposal: 91674322

  best value found: 9167432 bytes

STA(0): 4+(34<<16),58|(5<<16),39+(3<<16)+(24<<24)  4 34 58 5 39 3 24
Compressing tests\test9.cfg   stream(0).  Total 100000000
S[0] compressed 100000271->45663811 bytes. Data 100000000->45663617, model 271->194.
```
Final model
```c
int update(int y, int c0, int bpos, int c4) {
    vmx(DS,0,(c4&0xff)*256+c0);
    return 0;
}

void block(int a, int b) { 
}

int main() { 
    vmi(STA,0,4+(34<<16),58|(5<<16),39+(3<<16)+(24<<24),0);
    vmi(DS,0,16+(1<<16),1023,1,0);
}
```
47904276->45663617
```
{  1,  2, 0, 0},{  3,  5, 1, 0},{  4,  6, 0, 1},{  7,  8, 2, 0}, // 0-3
{  8,  9, 1, 1},{  8,  9, 1, 1},{  9, 10, 0, 2},{  7, 11, 3, 0}, // 4-7
{ 11, 12, 2, 1},{ 12, 13, 1, 2},{ 13, 10, 0, 3},{ 14, 15, 3, 1}, // 8-11
{ 15, 16, 2, 2},{ 16, 17, 1, 3},{ 18, 19, 4, 1},{ 19, 20, 3, 2}, // 12-15
{ 20, 21, 2, 3},{ 21, 22, 1, 4},{ 23, 24, 5, 1},{ 24, 25, 4, 2}, // 16-19
{ 25, 26, 3, 3},{ 26, 27, 2, 4},{ 27, 28, 1, 5},{ 29, 30, 6, 1}, // 20-23
{ 30, 15, 5, 2},{ 15, 31, 4, 3},{ 31, 16, 3, 4},{ 16, 32, 2, 5}, // 24-27
{ 32, 33, 1, 6},{ 34, 35, 7, 1},{ 35, 19, 6, 2},{ 36, 37, 4, 4}, // 28-31
{ 21, 38, 2, 6},{ 38, 39, 1, 7},{ 40, 41, 8, 1},{ 41, 24, 7, 2}, // 32-35
{ 42, 31, 5, 4},{ 31, 43, 4, 5},{ 27, 44, 2, 7},{ 44, 45, 1, 8}, // 36-39
{ 46, 47, 9, 1},{ 47, 24, 8, 2},{ 48, 36, 6, 4},{ 37, 49, 4, 6}, // 40-43
{ 27, 50, 2, 8},{ 50, 51, 1, 9},{ 52, 53,10, 1},{ 53, 30, 9, 2}, // 44-47
{ 54, 42, 7, 4},{ 43, 55, 4, 7},{ 32, 56, 2, 9},{ 56, 57, 1,10}, // 48-51
{ 58, 59,11, 1},{ 59, 35,10, 2},{ 60, 42, 8, 4},{ 43, 61, 4, 8}, // 52-55
{ 38, 62, 2,10},{ 62, 63, 1,11},{ 64, 65,12, 1},{ 65, 35,11, 2}, // 56-59
{ 66, 48, 9, 4},{ 49, 67, 4, 9},{ 38, 68, 2,11},{ 68, 69, 1,12}, // 60-63
{ 70, 71,13, 1},{ 71, 41,12, 2},{ 72, 54,10, 4},{ 55, 73, 4,10}, // 64-67
{ 44, 74, 2,12},{ 74, 75, 1,13},{ 76, 77,14, 1},{ 77, 47,13, 2}, // 68-71
{ 78, 60,11, 4},{ 61, 79, 4,11},{ 50, 80, 2,13},{ 80, 81, 1,14}, // 72-75
{ 82, 83,15, 1},{ 83, 47,14, 2},{ 84, 66,12, 4},{ 67, 85, 4,12}, // 76-79
{ 50, 86, 2,14},{ 86, 87, 1,15},{ 88, 89,16, 1},{ 89, 53,15, 2}, // 80-83
{ 90, 66,13, 4},{ 67, 91, 4,13},{ 56, 92, 2,15},{ 92, 93, 1,16}, // 84-87
{ 94, 95,17, 1},{ 95, 59,16, 2},{ 96, 72,14, 4},{ 73, 97, 4,14}, // 88-91
{ 62, 98, 2,16},{ 98, 99, 1,17},{100,101,18, 1},{101, 59,17, 2}, // 92-95
{102, 78,15, 4},{ 79,103, 4,15},{ 62,104, 2,17},{104,105, 1,18}, // 96-99
{106,107,19, 1},{107, 65,18, 2},{108, 84,16, 4},{ 85,109, 4,16}, // 100-103
{ 68,110, 2,18},{110,111, 1,19},{112,113,20, 1},{113, 71,19, 2}, // 104-107
{114, 90,17, 4},{ 91,115, 4,17},{ 74,116, 2,19},{116,117, 1,20}, // 108-111
{118,119,21, 1},{119, 71,20, 2},{120, 90,18, 4},{ 91,121, 4,18}, // 112-115
{ 74,122, 2,20},{122,123, 1,21},{124,125,22, 1},{125, 77,21, 2}, // 116-119
{126, 96,19, 4},{ 97,127, 4,19},{ 80,128, 2,21},{128,129, 1,22}, // 120-123
{130,131,23, 1},{131, 83,22, 2},{132,102,20, 4},{103,133, 4,20}, // 124-127
{ 86,134, 2,22},{134,135, 1,23},{136,131,24, 1},{137, 83,23, 2}, // 128-131
{138,108,21, 4},{109,139, 4,21},{ 86,140, 2,23},{134,141, 1,24}, // 132-135
{142,131,25, 1},{143, 83,24, 2},{144,114,22, 4},{115,145, 4,22}, // 136-139
{ 86,146, 2,24},{134,147, 1,25},{148,131,26, 1},{149, 83,25, 2}, // 140-143
{150,114,23, 4},{115,151, 4,23},{ 86,152, 2,25},{134,153, 1,26}, // 144-147
{154,131,27, 1},{155, 83,26, 2},{156,114,24, 4},{115,157, 4,24}, // 148-151
{ 86,158, 2,26},{134,159, 1,27},{160,131,28, 1},{161, 83,27, 2}, // 152-155
{162,114,25, 4},{115,163, 4,25},{ 86,164, 2,27},{134,165, 1,28}, // 156-159
{166,131,29, 1},{167, 83,28, 2},{168,114,26, 4},{115,169, 4,26}, // 160-163
{ 86,170, 2,28},{134,171, 1,29},{172,131,30, 1},{173, 83,29, 2}, // 164-167
{174,114,27, 4},{115,175, 4,27},{ 86,176, 2,29},{134,177, 1,30}, // 168-171
{178,131,31, 1},{179, 83,30, 2},{180,114,28, 4},{115,181, 4,28}, // 172-175
{ 86,182, 2,30},{134,183, 1,31},{184,131,32, 1},{185, 83,31, 2}, // 176-179
{186,114,29, 4},{115,187, 4,29},{ 86,188, 2,31},{134,189, 1,32}, // 180-183
{184,131,33, 1},{190, 83,32, 2},{191,114,30, 4},{115,192, 4,30}, // 184-187
{ 86,193, 2,32},{134,189, 1,33},{194, 83,33, 2},{195,114,31, 4}, // 188-191
{115,196, 4,31},{ 86,197, 2,33},{198, 83,34, 2},{199,114,32, 4}, // 192-195
{115,200, 4,32},{ 86,201, 2,34},{202, 83,35, 2},{203,114,33, 4}, // 196-199
{115,204, 4,33},{ 86,205, 2,35},{206, 83,36, 2},{207,114,34, 4}, // 200-203
{115,208, 4,34},{ 86,209, 2,36},{210, 83,37, 2},{211,114,35, 4}, // 204-207
{115,212, 4,35},{ 86,213, 2,37},{214, 83,38, 2},{215,114,36, 4}, // 208-211
{115,216, 4,36},{ 86,217, 2,38},{218, 83,39, 2},{219,114,37, 4}, // 212-215
{115,220, 4,37},{ 86,221, 2,39},{222, 83,40, 2},{107,114,38, 4}, // 216-219
{115,110, 4,38},{ 86,223, 2,40},{224, 83,41, 2},{ 86,225, 2,41}, // 220-223
{226, 83,42, 2},{ 86,227, 2,42},{228, 83,43, 2},{ 86,229, 2,43}, // 224-227
{230, 83,44, 2},{ 86,231, 2,44},{232, 83,45, 2},{ 86,233, 2,45}, // 228-231
{234, 83,46, 2},{ 86,235, 2,46},{236, 83,47, 2},{ 86,237, 2,47}, // 232-235
{238, 83,48, 2},{ 86,239, 2,48},{240, 83,49, 2},{ 86,241, 2,49}, // 236-239
{242, 83,50, 2},{ 86,243, 2,50},{244, 83,51, 2},{ 86,245, 2,51}, // 240-243
{246, 83,52, 2},{ 86,247, 2,52},{248, 83,53, 2},{ 86,249, 2,53}, // 244-247
{250, 83,54, 2},{ 86,251, 2,54},{252, 83,55, 2},{ 86,253, 2,55}, // 248-251
{254, 83,56, 2},{ 86,255, 2,56},{160, 83,57, 2},{ 86,165, 2,57}, // 252-255
```

[Top](#more-examples-and-tuning)
