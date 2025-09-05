/*
 * Code for simulating pthreads API on Windows.  This is Git-specific,
 * but it is enough for Numexpr needs too.
 *
 * Copyright (C) 2009 Andrzej K. Haczewski <ahaczewski@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * DISCLAIMER: The implementation is Git-specific, it is subset of original
 * Pthreads API, without lots of other features that Git doesn't use.
 * Git also makes sure that the passed arguments are valid, so there's
 * no need for double-checking.
 */

#ifndef PTHREAD_H
#define PTHREAD_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Defines that adapt Windows API threads to pthreads API
 */
#define pthread_mutex_t CRITICAL_SECTION

#define pthread_mutex_init(a,b) InitializeCriticalSection((a))
#define pthread_mutex_destroy(a) DeleteCriticalSection((a))
#define pthread_mutex_lock EnterCriticalSection
#define pthread_mutex_unlock LeaveCriticalSection

/*
 * Implement simple condition variable for Windows threads, based on ACE
 * implementation.
 *
 * See original implementation: http://bit.ly/1vkDjo
 * ACE homepage: http://www.cse.wustl.edu/~schmidt/ACE.html
 * See also: http://www.cse.wustl.edu/~schmidt/win32-cv-1.html
 */
typedef struct {
	LONG waiters;
	int was_broadcast;
	CRITICAL_SECTION waiters_lock;
	HANDLE sema;
	HANDLE continue_broadcast;
} pthread_cond_t;

extern int pthread_cond_init(pthread_cond_t *cond, const void *unused);
extern int pthread_cond_destroy(pthread_cond_t *cond);
extern int pthread_cond_wait(pthread_cond_t *cond, CRITICAL_SECTION *mutex);
extern int pthread_cond_signal(pthread_cond_t *cond);
extern int pthread_cond_broadcast(pthread_cond_t *cond);

/*
 * Simple thread creation implementation using pthread API
 */
typedef struct {
	HANDLE handle;
	void *(*start_routine)(void*);
	void *arg;
} pthread_t;

extern int pthread_create(pthread_t *thread, const void *unused,
			  void *(*start_routine)(void*), void *arg);

/*
 * To avoid the need of copying a struct, we use small macro wrapper to pass
 * pointer to win32_pthread_join instead.
 */
#define pthread_join(a, b) win32_pthread_join(&(a), (b))

extern int win32_pthread_join(pthread_t *thread, void **value_ptr);

/*
 * The POSIX signal system has a more developed interface than what's in
 * Windows. We create a no-op shim layer to proivde enough of the API to
 * pretend to support what's used when creating threads on POSIX systems.
 */
typedef int sigset_t;
enum sigop {
    SIG_BLOCK,
    SIG_UNBLOCK,
    SIG_SETMASK
};

static inline int sigemptyset(sigset_t *sigs) { return 0; }
static inline int sigfillset(sigset_t *sigs) { return 0; }
static inline int sigaddset(sigset_t *sigs, int sig) { return 0; }
static inline int sigdelset(sigset_t *sigs, int sig) { return 0; }
static inline int pthread_sigmask(int how, sigset_t *newmask,
                                  sigset_t *oldmask) { return 0; }

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* PTHREAD_H */
