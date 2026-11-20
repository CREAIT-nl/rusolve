// nodeprof sampler: LD_PRELOAD CPU-time profiler for hosts where perf_event_open is
// not permitted (kernel.perf_event_paranoid=4). ITIMER_PROF fires in process CPU time,
// so the sample count is insensitive to machine load. The handler records the
// interrupted ip and walks the frame-pointer chain (build the target with frame
// pointers: -fno-omit-frame-pointer / -C force-frame-pointers=yes). Stack words are
// read with process_vm_readv, which fails with EFAULT instead of crashing the handler.
// At exit it writes <NODEPROF_OUT>: header, raw samples, then /proc/self/maps.
//   NODEPROF_OUT=prof.bin NODEPROF_US=2000 LD_PRELOAD=./sampler.so cmd ...
#define _GNU_SOURCE
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <ucontext.h>
#include <unistd.h>

#define DEPTH 48
#define MAXS 400000
static uint64_t (*buf)[DEPTH];
static volatile uint32_t nsamp;
static uint32_t nlost;

static int rd(uint64_t addr, uint64_t *out) {
  struct iovec l = {out, 8}, r = {(void *)addr, 8};
  return process_vm_readv(getpid(), &l, 1, &r, 1, 0) == 8;
}

static void handler(int sig, siginfo_t *si, void *uc_) {
  (void)sig; (void)si;
  ucontext_t *uc = uc_;
  uint32_t i = __sync_fetch_and_add(&nsamp, 1);
  if (i >= MAXS) { nlost++; return; }
  uint64_t *s = buf[i];
  uint64_t ip = uc->uc_mcontext.gregs[REG_RIP];
  uint64_t bp = uc->uc_mcontext.gregs[REG_RBP];
  uint64_t sp = uc->uc_mcontext.gregs[REG_RSP];
  int d = 0;
  s[d++] = ip;
  // leaf code without a frame (libc memcpy/memset): [rsp] is its return address, which the
  // frame-pointer chain skips; report.py keeps it only when ip is outside the main exe.
  { uint64_t w; s[d++] = rd(sp, &w) ? w : 0; }
  while (d < DEPTH) {
    uint64_t next, ret;
    if (bp < sp || (bp & 7) || bp - sp > (1ull << 30)) break;
    if (!rd(bp, &next) || !rd(bp + 8, &ret)) break;
    if (ret == 0) break;
    s[d++] = ret;
    if (next <= bp) break;
    bp = next;
  }
  while (d < DEPTH) s[d++] = 0;
}

__attribute__((destructor)) static void dump(void) {
  struct itimerval z = {{0, 0}, {0, 0}};
  setitimer(ITIMER_PROF, &z, NULL);
  const char *out = getenv("NODEPROF_OUT");
  if (!out) return;
  FILE *f = fopen(out, "wb");
  if (!f) return;
  uint32_t n = nsamp > MAXS ? MAXS : nsamp;
  uint32_t hdr[4] = {0x4e50524f, DEPTH, n, nlost};
  fwrite(hdr, 4, 4, f);
  fwrite(buf, sizeof(uint64_t) * DEPTH, n, f);
  FILE *m = fopen("/proc/self/maps", "r");
  char line[512];
  while (m && fgets(line, sizeof line, m)) fputs(line, f);
  if (m) fclose(m);
  fclose(f);
}

__attribute__((constructor)) static void init(void) {
  buf = malloc(sizeof(uint64_t) * DEPTH * MAXS);
  if (!buf) return;
  struct sigaction sa;
  memset(&sa, 0, sizeof sa);
  sa.sa_sigaction = handler;
  sa.sa_flags = SA_SIGINFO | SA_RESTART;
  sigaction(SIGPROF, &sa, NULL);
  const char *us = getenv("NODEPROF_US");
  long u = us ? atol(us) : 2000;
  struct itimerval t = {{0, u}, {0, u}};
  setitimer(ITIMER_PROF, &t, NULL);
}
