# FX language

Language used by the compressor is a subset of C language.

## Basic data types
* [char](#char-short-int)
* [short](#char-short-int)
* [int](#char-short-int)
* [enum](#enum)
## Extended data types
* [distance](#distance)
* [column](#column)
* [table](#table)
* [dictionary](#dictionary)
* [record](#record)

Types with * are not yet implemented.

## Flow control
* if, '?'
* else
* return
* for
* while
* void, int

## Operators
* sizeof
* '//' - comment

### Math
* '+', '-', '/', '*', '%','~'
* '<<' - unsigned
* '>>' - unsigned
* '|' - unsigned
* '&' - unsigned
* '^' - unsigned

### Compare, assign
* '>', '<', '<=', '>=', '!=', '=='
* '=','!'

## Functions
* printf
* read
* write
* exit
* vmi
* vmx

## Required function prototypes depending on what type of config file is used:
* block
* update
* bitupdate
* byteupdate
* main
* detect
* decode
* encode


## Examples
```c
int a;               // Global variable a
int b[5]={};         // Global array where all values are zero, local arrays are not allowed
int c[5]={1,2};      // Global array where first two values are set
enum {VAL, VAL2=3};  // Enum
// Add two values and return value
int add(int g, int f) {
    return g+f;
}
void change() {
    int b;          // Local variable, must be defined at the beginning. Uninitialized.
    b=0;
    b=5+VAL2;       // Add to local variable sum of 5 and VAL2
    a=a+b;          // Add to global variable
    b[3]=add(a,b);  // Assign to global array
}
```

[Top](#fx-language)

# char, short, int
[Top](#fx-language)
# enum
Default enum values set by the compiler itself.
```c
// Used in models. (.cfg)
enum {SMC=1,APM1,DS,AVG,SCM,RCM,CM,MX,ST,MM,DHS,SM,SK,APM2,ERR,UAS,LMX,STA,BYT};
enum {false=0,true=1};
// Used in detection. (.det)
enum {NONE=0,START,INFO,END,RESET=0xfffffffe,REQUEST=0xffffffff};

```
[Top](#fx-language)
# distance
```
char charArray[2]={'(',')'};
distance true/false charArray brcxt;   // Create.
                                       //   true - when updating, remove last inserted byte found in charArray.
                                       //   false - allow nested bytes found in charArray
                                       //   charArray contains a pair of bytes. Array must be multiple of 2
                                       //     First byte is used to start a new context.
                                       //     Second byte removes last context if it matches to our current context
brcxt=brcxt+c;                         // Update
brcxt=0;                               // Reset
a=brcxt;                               // Get context. Return value is in format 0x0000XXYY
                                       //   YY is last byte detected in charArray
                                       //   XX is distance in bytes to the last byte YY, max lenght 255 bytes.
```
[Top](#fx-language)
# column
```
column maxColLen maxFirsChar nlChar colcxt;    // Create
colcxt=0;                                      // Reset
colcxt=colcxt+val;                             // Update, val is last byte
a=colcxt;                                      // Get context. Return value is in format 0xQQWWXXYY
                                               //   YY is first char in column (byte)
                                               //   XX is column length (max colLen) 
                                               //   WW is byte in above column
                                               //   QQ is 1 for new line, 0 for none
```
[Top](#fx-language)
# table
```
table startendchars celld tblcxt;             // Create
tblcxt=0;                                     // Reset
tblcxt=tblcxt+val;                            // Update, val is last 4 bytes (c4)
a=colcxt;                                     // Get context. Column position in above row.
```
[Top](#fx-language)

# dictionary
```
char string[64]={' '};
dictionary "filename" string dict;     // Create a dictionary from file "filename" and use string as input/output buffer.
                                       // String must be 64 bytes.
a=dict;                                // Return number of words loaded

dict=val;                              // val=0   Get dictionary index of string in char array
                                       // val=1   Get encoded dictionary index of string in char array
                                       // val=2   Get encoded dictionary index of substring (prefix) in char array 
                                       // val>127 Get string from encoded dictionary id into char array
a=dict;                                // Get value
                                       // val=0   Index
                                       // val=1   Encoded index
                                       // val=2   Encoded index in 0x00XXXXXX and length of in 0xYY000000
                                       // val>127 Length of decoded string, string contains decoded word
```
First line in the dictionary must contain a number of words. Each line contains one word. \n terminated.

Example:
```
4
dog
cat
line
box
```
Maximum allowed word count is 44880 words.

[Top](#fx-language)

# record
```
// Record context
record reccxt;        // Create 
reccxt=0;             // Reset
reccxt=reccxt+val;    // Update, val must be byte
rlen=reccxt&0xffff;   // Record length max 65535 bytes.
i=reccxt>>16;         // Distance to the c1 in bytes. Max length is 255 bytes.
```
[Top](#fx-language)
