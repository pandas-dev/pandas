/* The include and declarations in this file exists because when cython translates to c it cannot translate an IF UNAME_SYSNAME=="Windows" to
 * a c #ifdef. That means the resulting c file will be different if produced on windows or on unix. We want our c files to be portable and 
 * therefore do the right includes for windows here, and just declare dummy symbols for unix. 
 */
 
#ifdef _WIN32

    #include <fcntl.h>
    #include <share.h>
    #include <io.h>
    #include <sys/stat.h>
    
    // Stuff for handling paths with international characters on windows
    void assign_fd(void *io_ctx, int fd) { ((unistd_io_ctx_t*)io_ctx)->fd = fd; }
    //ssize_t write(int fd, const void *buf, size_t nbyte){return 0;};
    //int close(int fd);
        
#else

    #include<sys/stat.h>
    #include<unistd.h>
    #include<fcntl.h>

    //int open(const char *path, int oflag, int mode);
    //int close(int fd);
    int _wsopen(const wchar_t *filename, int oflag, int shflag, int pmode){ return 0; };
    int _O_RDONLY;
    int _O_BINARY;
    int _O_WRONLY;
    int _O_TRUNC;
    int _O_CREAT;
    int _SH_DENYRW;  
    int _SH_DENYWR; 
    int _SH_DENYRD;
    int _SH_DENYNO;
    int _S_IWRITE;
    int _S_IREAD;
    void assign_fd(void *io_ctx, int fd){};
    int _close(int fd){ return 0; };
    ssize_t _write(int fd, const void *buf, size_t nbyte){return 0;};
    
#endif
