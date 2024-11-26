
#ifdef WINDOWS
#define PROT_NONE       0
#define PROT_READ       1
#define PROT_WRITE      2
#define PROT_EXEC       4

#define MAP_FILE        0
#define MAP_SHARED      1
#define MAP_PRIVATE     2
#define MAP_TYPE        0xf
#define MAP_FIXED       0x10
#define MAP_ANONYMOUS   0x20
#define MAP_ANON        MAP_ANONYMOUS
#define MAP_FAILED      ((void *)-1)

void*   mmap(void *addr, size_t len, int prot, int flags, int fildes, off_t off);
#else
#include <sys/mman.h>
#endif
void jitReadonly(void *p, size_t size){
#ifdef WINDOWS
  DWORD oldProtect;
    VirtualProtect(p, size, PAGE_EXECUTE_READ, &oldProtect);
#else
    mprotect(p, size, PROT_READ | PROT_EXEC);
#endif
}

#ifdef WINDOWS
void* mmap(void *addr, size_t len, int prot, int flags, int fildes, off_t off){
    HANDLE fm, h;
    void * map = MAP_FAILED;
    const off_t maxSize = off + (off_t)len;
    h = (HANDLE)_get_osfhandle(fildes);
    fm = CreateFileMapping(h, NULL, PAGE_EXECUTE_READWRITE, 0, maxSize, NULL);
    map = MapViewOfFile(fm, FILE_MAP_READ | FILE_MAP_WRITE | FILE_MAP_EXECUTE, 0, off, len);
    CloseHandle(fm);
    return map;
}
#endif
