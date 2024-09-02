#line 0 "reactor-cpp.c"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
#line 1 "reactor-cpp.c"
#if 700 < 700
  #undef 700
  #define 700 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/Min/basilisk/src/common.h"
#line 1 "/home/Min/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/Min/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif _MPI

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 49 "/home/Min/basilisk/src/common.h"
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 81

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 105 "/home/Min/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/Min/basilisk/src/common.h", 180, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/Min/basilisk/src/common.h", 452, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define tracing(func, file, line) trace_push (&trace_func, func)
# define end_tracing(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/Min/basilisk/src/common.h", 556, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define tracing(...)
# define end_tracing(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/Min/basilisk/src/common.h", 666, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 668

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/Min/basilisk/src/common.h", 670, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 673


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/Min/basilisk/src/common.h",680);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/Min/basilisk/src/common.h",683);return _ret;}
end_tracing("mpi_all_reduce0","/home/Min/basilisk/src/common.h",684);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 692

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type * global = malloc ((elem)*sizeof(type)), * tmp = malloc ((elem)*sizeof(type));\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else if (!strcmp(#type, "unsigned char")) datatype = MPI_UNSIGNED_CHAR;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  free (global), free (tmp);\
  prof_stop();\
}\

#line 715


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 826 "/home/Min/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/Min/basilisk/src/common.h", 836, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 920 "/home/Min/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 924
norm += sq(n->y);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 927
n->y /= norm;
}

void origin (double x, double y, double z) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define _dirichlet(expr, ...) (2.*(expr) - val(_s,0,0,0))
#define _dirichlet_homogeneous(...) (- val(_s,0,0,0))
#define _dirichlet_face(expr,...) (expr)
#define _dirichlet_face_homogeneous(...) (0.)
#define _neumann(expr,...) (Delta*(expr) + val(_s,0,0,0))
#define _neumann_homogeneous(...) (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/home/Min/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/Min/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 965 "/home/Min/basilisk/src/common.h"



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/Min/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/Min/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 178 "/home/Min/basilisk/src/embed.h"
bool third;
  
#line 28 "/home/Min/basilisk/src/vof.h"
scalar * tracers, c;
  bool inverse;
  
#line 456 "/home/Min/basilisk/src/heights.h"
vector height;
  
#line 21 "/home/Min/basilisk/src/iforce.h"
scalar phi;
  
#line 22 "/home/Min/basilisk/src/tension.h"
double sigma;

#line 988 "/home/Min/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 1089
if (w.y.i != v.y.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/Min/basilisk/src/common.h", 1112, "s->i >= 0");
      v.x = *s++;
    } 
#line 1111
{
      if (!(s->i >= 0)) qassert ("/home/Min/basilisk/src/common.h", 1112, "s->i >= 0");
      v.y = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/Min/basilisk/src/common.h", 1143, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1142
{
      if (!(v->y.i >= 0)) qassert ("/home/Min/basilisk/src/common.h", 1143, "v->x.i >= 0");
      t.y = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 1154
if (_attribute[s.i].d.y != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
vector (* init_face_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
void (* scalar_clone) (scalar, scalar);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
#if _CADNA
  cadna_init (-1);
#endif
#if _MPI
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1}};
const vector unityf = {{_NVARMAX+2},{_NVARMAX+3}};
const scalar unity = {_NVARMAX+4};
const scalar zeroc = {_NVARMAX+5};



const vector unityf0 = {{_NVARMAX+6},{_NVARMAX+7}};
const scalar unity0 = {_NVARMAX+8};
        vector fm = {{_NVARMAX+6},{_NVARMAX+7}};
        scalar cm = {_NVARMAX+8};
#line 1281 "/home/Min/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = 1e30;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (const char * commands, bool overwrite)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, commands);
  }
}



typedef struct {
  double x;

  double y;




} _coord;





#line 1 "/home/Min/basilisk/src/grid/stencils.h"
#line 17 "/home/Min/basilisk/src/grid/stencils.h"










typedef struct {
  char * name;
  char * type;
  void * pointer;
  int * dimensions;
  int nd;
  char reduct;
  scalar data;
} NonLocal;

typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
  int parallel;
  scalar * listc;
  vectorl listf;
  scalar * dirty;
  void * data;
} ForeachData;


#define foreach_stencil(...) {\
  static ForeachData _loop = {\
    .fname = __FILE__, .line = __LINE__, .first = 1\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 62


#define end_foreach_stencil()\
  check_stencil (&_loop);\
  boundary_stencil (&_loop);\
  _loop.first = 0;\
}\

#line 69


#define foreach_vertex_stencil(...) foreach_stencil(__VA_ARGS__) _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_visible_stencil() end_foreach_stencil()

#define foreach_point_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_point_stencil() end_foreach_stencil()

#define foreach_region_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_region_stencil() end_foreach_stencil()

#define _stencil_is_face_x() { _loop.face |= (1 << 0);
#define end__stencil_is_face_x() }
#define _stencil_is_face_y() { _loop.face |= (1 << 1);
#define end__stencil_is_face_y() }
#define _stencil_is_face_z() { _loop.face |= (1 << 2);
#define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 100

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 103

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 106

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 109


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void check_stencil (ForeachData * loop)
{
  loop->listf = (vectorl){NULL};




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     loop->listc = list_append (loop->listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.x = list_append (loop->listf.x, s);
       }
       
#line 213
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.y = list_append (loop->listf.y, s);
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   loop->listc = list_append (loop->listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 242
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 260
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 290
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 296
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 308
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 315
_attribute[s.i].v.y.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 loop->dirty = list_append (loop->dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     loop->dirty = list_append (loop->dirty, d);}}
      }
    }
  }}}
}




void boundary_stencil (ForeachData * loop)
{
  bool flux = false;
  
    if (loop->listf.x)
      flux = true;
    
#line 344
if (loop->listf.y)
      flux = true;
  if (flux) {
#line 359 "/home/Min/basilisk/src/grid/stencils.h"
    boundary_face (loop->listf);
    
      pfree (loop->listf.x,__func__,__FILE__,__LINE__), loop->listf.x = NULL;
      
#line 361
pfree (loop->listf.y,__func__,__FILE__,__LINE__), loop->listf.y = NULL;
  }




  if (loop->listc) {






    boundary_internal (loop->listc, loop->fname, loop->line);
    pfree (loop->listc,__func__,__FILE__,__LINE__), loop->listc = NULL;
  }





  if (loop->dirty) {






    {scalar*_i=(scalar*)( loop->dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (loop->dirty,__func__,__FILE__,__LINE__), loop->dirty = NULL;
  }
}
#line 1462 "/home/Min/basilisk/src/common.h"
#line 14 "reactor-cpp.c"
#line 1 "reactor.c"
#line 1 "grid/multigrid.h"
#line 1 "/home/Min/basilisk/src/grid/multigrid.h"
#line 16 "/home/Min/basilisk/src/grid/multigrid.h"
typedef struct {
  Grid g;
  char ** d;
} Multigrid;

struct _Point {
  int i;

  int j;




  int level, n;
#ifdef foreach_block
  int l;
  #define _BLOCK_INDEX , point.l
#else
  #define _BLOCK_INDEX
#endif
};
static Point last_point;
#line 49 "/home/Min/basilisk/src/grid/multigrid.h"
static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*2;
  return sq(n);
}
#line 62 "/home/Min/basilisk/src/grid/multigrid.h"
#define data(k,l,m)\
  ((double *)&((Multigrid *)grid)->d[point.level][((point.i + k)*((1 << point.level) +\
       2*2) +\
      (point.j + l))*datasize]) 
#line 64

#line 89 "/home/Min/basilisk/src/grid/multigrid.h"
#define allocated(k,l,m) (point.i+k >= 0 && point.i+k < (1 << point.level) + 2*2 &&\
         point.j+l >= 0 && point.j+l < (1 << point.level) + 2*2)\

#line 91


#define allocated_child(k,l,m) (level < depth() &&\
         point.i > 0 && point.i <= (1 << point.level) + 2 &&\
         point.j > 0 && point.j <= (1 << point.level) + 2)\

#line 96

#line 117 "/home/Min/basilisk/src/grid/multigrid.h"
#define depth() (grid->depth)
#line 136 "/home/Min/basilisk/src/grid/multigrid.h"
#define fine(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level+1][((2*point.i-2 +k)*2*((1 << point.level) +\
        2) +\
     (2*point.j-2 +l))*datasize])[_index(a,m)]\

#line 141

#define coarse(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level-1][(((point.i+2)/2+k)*((1 << point.level)/2 +\
        2*2) +\
     (point.j+2)/2+l)*datasize])[_index(a,m)]\

#line 147

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  }; NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2; parent.j = (point.j + 2)/2;\

#line 157

#line 191 "/home/Min/basilisk/src/grid/multigrid.h"
#define foreach_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 208

#define end_foreach_level()\
\
 }\
\
  }\
}\

#line 215


#define foreach()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 234

#define end_foreach()\
\
 }\
\
  }\
}\

#line 241


#define is_active(cell) (true)
#define is_leaf(cell) (level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (stderr, "grid depths do not match. Aborting.\n");\
  if (!(0)) qassert ("/home/Min/basilisk/src/grid/multigrid.h", 249, "0");\
} while (0)\

#line 251

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/Min/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/Min/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 254 "/home/Min/basilisk/src/grid/multigrid.h"

#define foreach_face_generic()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k <= point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j <= point.n + 2; point.j++)\
\
\
\
        {\
\
   POINT_VARIABLES\

#line 272

#define end_foreach_face_generic()\
\
 }\
\
  }\
}\

#line 279


#define foreach_vertex()\
foreach_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 290

#define end_foreach_vertex() } end_foreach_face_generic()

#define is_coarse() (point.level < depth())
#line 321 "/home/Min/basilisk/src/grid/multigrid.h"
#define is_face_x() { int ig = -1; VARIABLES; if (point.j < point.n + 2) {
#define end_is_face_x() }}
#define is_face_y() { int jg = -1; VARIABLES; if (point.i < point.n + 2) {
#define end_is_face_y() }}

#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  point.n *= 2;\
  for (int _k = 0; _k < 2; _k++)\
    for (int _l = 0; _l < 2; _l++) {\
      point.i = _i + _k; point.j = _j + _l;\
      POINT_VARIABLES;\

#line 334

#define end_foreach_child()\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
  point.n /= 2;\
}\

#line 341

#define foreach_child_break() _k = _l = 2
#line 387 "/home/Min/basilisk/src/grid/multigrid.h"
#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/Min/basilisk/src/grid/neighbors.h"
#line 17 "/home/Min/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 393 "/home/Min/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Point p;
  p.level = depth(); p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < sq(p.n + 2*2); i++)
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     ((double *)(&((Multigrid *)grid)->d[p.level][i*datasize]))[s.i + b] = val;
      }}}
}
#line 433 "/home/Min/basilisk/src/grid/multigrid.h"
#define foreach_boundary_dir(l,d)\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l < 0 ? depth() : l;\
  point.n = 1 << point.level;\
  int * _i = &point.j;\
  if (d == left) {\
    point.i = 2;\
    ig = -1;\
  }\
  else if (d == right) {\
    point.i = point.n + 2 - 1;\
    ig = 1;\
  }\
  else if (d == bottom) {\
    point.j = 2;\
    _i = &point.i;\
    jg = -1;\
  }\
  else if (d == top) {\
    point.j = point.n + 2 - 1;\
    _i = &point.i;\
    jg = 1;\
  }\
  int _l;\
  OMP(omp for schedule(static))\
  for (_l = 0; _l < point.n + 2*2; _l++) {\
    *_i = _l;\
    {\
      POINT_VARIABLES\

#line 464

#define end_foreach_boundary_dir()\
    }\
  }\
}\

#line 469


#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.level, point.n _BLOCK_INDEX})\

#line 473

#define is_boundary(point) (point.i < 2 || point.i >= point.n + 2 ||\
    point.j < 2 || point.j >= point.n + 2)\

#line 476

#line 538 "/home/Min/basilisk/src/grid/multigrid.h"
#define foreach_boundary(b)\
  if (default_scalar_bc[b] != periodic_bc)\
    foreach_boundary_dir (depth(), b)\
      if (!is_boundary(point)) {\

#line 542

#define end_foreach_boundary() } end_foreach_boundary_dir()

#define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, void * data);

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  extern double (* default_scalar_bc[]) (Point, Point, scalar, void *);
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int d = 0; d < 2*2; d++)
    if (default_scalar_bc[d] == periodic_bc)
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   if (is_vertex_scalar (s))
     _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
   else if (_attribute[s.i].face) {
     vector v = _attribute[s.i].v;
     _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
   }
 }}}
  for (int bghost = 1; bghost <= 2; bghost++)
    for (int d = 0; d < 2*2; d++) {

      scalar * list = NULL, * listb = NULL;
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 2) % 2];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }}}

      if (list) {
 {foreach_boundary_dir (l, d) {
   scalar s, sb;
   {scalar*_i0=listb;scalar*_i1= list;if(_i0)for(sb=*_i0,s=*_i1;_i0->i>= 0;sb=*++_i0,s=*++_i1){ {
     if ((_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) || is_vertex_scalar (s)) {

       if (bghost == 1)
 
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
     }
     else

      
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL);
   }}}
 }end_foreach_boundary_dir();}
 pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 654 "/home/Min/basilisk/src/grid/multigrid.h"
#define VT _attribute[s.i].v.y


static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face) {
 scalar vt = VT;
 if (_attribute[vt.i].boundary[right] == periodic_bc)
   list1 = list_add (list1, s);
      }
      else if (_attribute[s.i].boundary[right] == periodic_bc)
 list1 = list_add (list1, s);
    }}}
  if (!list1)
    return;

  if (l == 0) {
    {foreach_level(0)
      {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double * v = &val(s,0,0,0);
 {foreach_neighbor()
   memcpy (&val(s,0,0,0), v, _attribute[s.i].block*sizeof(double));end_foreach_neighbor()}
      }}}end_foreach_level();}
    pfree (list1,__func__,__FILE__,__LINE__);
    return;
  }

  OMP_PARALLEL() {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;

    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*2; j++) {
 for (int i = 0; i < 2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,i,j,0), &val(s,i + point.n,j,0), _attribute[s.i].block*sizeof(double));}}
 for (int i = point.n + 2; i < point.n + 2*2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,i,j,0), &val(s,i - point.n,j,0), _attribute[s.i].block*sizeof(double));}}
      }
#line 711 "/home/Min/basilisk/src/grid/multigrid.h"
  }
  pfree (list1,__func__,__FILE__,__LINE__);
}

#line 657
static void periodic_boundary_level_y (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face) {
 scalar vt = VT;
 if (_attribute[vt.i].boundary[top] == periodic_bc)
   list1 = list_add (list1, s);
      }
      else if (_attribute[s.i].boundary[top] == periodic_bc)
 list1 = list_add (list1, s);
    }}}
  if (!list1)
    return;

  if (l == 0) {
    {foreach_level(0)
      {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double * v = &val(s,0,0,0);
 {foreach_neighbor()
   memcpy (&val(s,0,0,0), v, _attribute[s.i].block*sizeof(double));end_foreach_neighbor()}
      }}}end_foreach_level();}
    pfree (list1,__func__,__FILE__,__LINE__);
    return;
  }

  OMP_PARALLEL() {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;

    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*2; j++) {
 for (int i = 0; i < 2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,j,i,0), &val(s,j,i + point.n,0), _attribute[s.i].block*sizeof(double));}}
 for (int i = point.n + 2; i < point.n + 2*2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,j,i,0), &val(s,j,i - point.n,0), _attribute[s.i].block*sizeof(double));}}
      }
#line 711 "/home/Min/basilisk/src/grid/multigrid.h"
  }
  pfree (list1,__func__,__FILE__,__LINE__);
}

#undef VT





void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++)
    pfree (m->d[l],__func__,__FILE__,__LINE__);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  N = 1 << depth();

  grid->n = grid->tn = 1 << 2*depth();

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);





   {
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_x;
    add_boundary (b);
  } 
#line 759
{
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_y;
    add_boundary (b);
  }


  m->d = (char **) pmalloc(sizeof(Point *)*(depth() + 1),__func__,__FILE__,__LINE__);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l)*datasize;
    m->d[l] = (char *) pmalloc (len,__func__,__FILE__,__LINE__);


    double * v = (double *) m->d[l];
    for (int i = 0; i < len/sizeof(double); i++)
      v[i] = undefined;
  }
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l);
    p->d[l] = (char *) prealloc (p->d[l], (len*(datasize + size))*sizeof(char),__func__,__FILE__,__LINE__);
    char * data = p->d[l] + (len - 1)*datasize;
    for (int i = len - 1; i > 0; i--, data -= datasize)
      memmove (data + i*size, data, datasize);
  }
  datasize += size;
}
#line 806 "/home/Min/basilisk/src/grid/multigrid.h"
Point locate (double xp, double yp, double zp)
{
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  point.level = -1, point.n = 1 << depth();
#line 825 "/home/Min/basilisk/src/grid/multigrid.h"
  point.i = (xp - X0)/L0*point.n + 2;
  if (point.i < 2 || point.i >= point.n + 2)
    return point;

  point.j = (yp - Y0)/L0*point.n + 2;
  if (point.j < 2 || point.j >= point.n + 2)
    return point;







  point.level = depth();
  return point;
}

#line 1 "grid/multigrid-common.h"
#line 1 "/home/Min/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/Min/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/Min/basilisk/src/grid/events.h"




static int END_EVENT = 1234567890;
static double TEND_EVENT = 1234567890;
static double TEPS = 1e-9;

static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = - TEND_EVENT;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == - TEND_EVENT) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == TEND_EVENT) {
 ev->i = END_EVENT; ev->t = - TEND_EVENT;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != - TEND_EVENT)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->i = -1; ev->t = - TEND_EVENT;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/Min/basilisk/src/grid/events.h", 88, "Events");
  if (!(!event.last)) qassert ("/home/Min/basilisk/src/grid/events.h", 89, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/Min/basilisk/src/grid/events.h", 93, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 136 "/home/Min/basilisk/src/grid/events.h"
static bool overload_event() { return true; }

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (!overload_event() || iter == ev->i || fabs (t - ev->t) <= TEPS*t) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = END_EVENT; tnext = 1e30;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if (overload_event() && (!cond || cond1) && (tnext != 1e30 || inext != END_EVENT)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != 1e30 && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/Min/basilisk/src/grid/events.h", 252, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt*(1. + TEPS))
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/Min/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  double Delta_x = Delta;\
\
  double Delta_y = Delta;\
\
\
\
\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/home/Min/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 47 "/home/Min/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

#define foreach_point(...)\
{\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  coord _p = { __VA_ARGS__ };\
  Point point = locate (_p.x, _p.y, _p.z);\
  if (point.level >= 0) {\
    POINT_VARIABLES\

#line 57

#define end_foreach_point() }}

#define foreach_region(p, box, n)\
  OMP_PARALLEL() { NOT_UNUSED (p);\
    coord p = {0, 0, box[0].z};\
    OMP(omp for schedule(static))\
      for (int _i = 0; _i < (int) n.x; _i++) {\
 p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);\
 for (int _j = 0; _j < (int) n.y; _j++) {\
   p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);\
   Point point = locate (p.x, p.y, p.z);\
   if (point.level >= 0) {\
     int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
     POINT_VARIABLES\

#line 72

#define end_foreach_region() }}}}





static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
  }
  _attribute[sb.i].name = pstrdup (bname,__func__,__FILE__,__LINE__);
  all = list_append (all, sb);
}

#define interpreter_set_int(...)
#define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
 init_block_scalar (sb, name, ext, n, block);
 interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/Min/basilisk/src/grid/cartesian-common.h", 124, "nvar + block <= _NVARMAX");

  if (_attribute == NULL)
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  else
    _attribute = (_Attributes *)
      prealloc (_attribute, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}

scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = alloc_block_scalar (name, ext.x, block);
    
#line 166
v.y = alloc_block_scalar (name, ext.y, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 190
vb.y.i = v.y.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 193
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 196
_attribute[v.y.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 206
vb.y.i = v.y.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 209
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 212
_attribute[v.y.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = alloc_block_vector (cname, 1);
  } 
#line 221
{
    sprintf (cname, ext.y, name);
    t.y = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  
    t.x.x = alloc_block_scalar (name, ext.x, 1);
    
#line 234
t.y.y = alloc_block_scalar (name, ext.y, 1);

    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;
#line 248 "/home/Min/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 273
init_const_scalar (v.y, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 280
v.y.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = _attribute[clone.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[clone.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[clone.i].boundary_homogeneous;
  if (!(_attribute[src.i].block > 0 && _attribute[clone.i].block == _attribute[src.i].block)) qassert ("/home/Min/basilisk/src/grid/cartesian-common.h", 291, "src.block > 0 && clone.block == src.block");
  pfree (_attribute[clone.i].depends,__func__,__FILE__,__LINE__);
  _attribute[clone.i] = _attribute[src.i];
  _attribute[clone.i].name = cname;
  _attribute[clone.i].boundary = boundary;
  _attribute[clone.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[clone.i].boundary[i] = _attribute[src.i].boundary[i];
    _attribute[clone.i].boundary_homogeneous[i] = _attribute[src.i].boundary_homogeneous[i];
  }
  _attribute[clone.i].depends = list_copy (_attribute[src.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 318
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/Min/basilisk/src/grid/cartesian-common.h", 369, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 421
list1.y = list_append (list1.y, v.y);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);
    
#line 424
pfree (list1.y,__func__,__FILE__,__LINE__);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/Min/basilisk/src/grid/cartesian-common.h",440);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/Min/basilisk/src/grid/cartesian-common.h",443);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 452
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 467
pfree (listf.y,__func__,__FILE__,__LINE__);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/Min/basilisk/src/grid/cartesian-common.h",475);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 485
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;
  _attribute[s.i].name = pname;
  if (block < 0)
    _attribute[s.i].block = block;
  else
    _attribute[s.i].block = block > 0 ? block : 1;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 533
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 545
_attribute[s.i].d.y = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 560
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_scalar (v.y, cname);
    }
    else
      cartesian_init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 580
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  } 
#line 592
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_vector (t.y, cname);
    }
    else
      cartesian_init_vector (t.y, NULL);
  }






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

void output_cells (FILE * fp, coord c, double size)
{
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && size > 0. &&
   (o.x > c.x + size || o.x < c.x - size))
 inside = false;
      
#line 627
if (inside && size > 0. &&
   (o.y > c.y + size || o.y < c.y - size))
 inside = false;
    if (inside) {
      Delta /= 2.;



      fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 655 "/home/Min/basilisk/src/grid/cartesian-common.h"
    }
  }end_foreach();}
  fflush (fp);
}
#line 667 "/home/Min/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){



    fprintf (fp, "x y %s ", _attribute[v.i].name);}}



  fputc ('\n', fp);
#line 734 "/home/Min/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }}}
 fputc ('\n', fp);
      }
#line 764 "/home/Min/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_face_vector = cartesian_init_face_vector;
  init_tensor = cartesian_init_tensor;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  scalar_clone = cartesian_scalar_clone;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
      double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;







  x = (xp - x)/Delta - _attribute[v.i].d.x/2.;
  y = (yp - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 834 "/home/Min/basilisk/src/grid/cartesian-common.h"
}


#line 805
static void _stencil_interpolate_linear (Point point, scalar v,
_stencil_undefined * xp,_stencil_undefined * yp,_stencil_undefined * zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;      







        
        
  
       

_stencil_val(v,0,0,0);_stencil_val(v, o_stencil,0,0);
_stencil_val(v,0,o_stencil,0); _stencil_val(v,o_stencil,o_stencil,0);

  
#line 820
return         
    ;
#line 834 "/home/Min/basilisk/src/grid/cartesian-common.h"
}

     
double interpolate (scalar v, double xp, double yp, double zp,
      bool linear)
{tracing("interpolate","/home/Min/basilisk/src/grid/cartesian-common.h",837);
  double val = 1e30;
  foreach_point_stencil (1,{(NonLocal[]){{"v","scalar",(void *)&v,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"val","double",(void *)&val,NULL,0,'m'},{"zp","double",(void *)&zp,NULL,0},{"yp","double",(void *)&yp,NULL,0},{"xp","double",(void *)&xp,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\n}","\n    \n// #line 842 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\nval = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);"})
    { _stencil_interpolate_linear (point, v, NULL, NULL, NULL); _stencil_val(v,0,0,0);    }end_foreach_point_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (min:val)){
#line 841
foreach_point (xp, yp, zp)
    val = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);end_foreach_point();mpi_all_reduce_array(&val,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 843
{end_tracing("interpolate","/home/Min/basilisk/src/grid/cartesian-common.h",843);return val;}
end_tracing("interpolate","/home/Min/basilisk/src/grid/cartesian-common.h",844);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v,
   bool linear)
{tracing("interpolate_array","/home/Min/basilisk/src/grid/cartesian-common.h",847);
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len++;}}
  for (int i = 0; i < n; i++) {
    double * w = v;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      *(w++) = 1e30;}}
    foreach_point_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"v","double",(void *)v,NULL,1},{"list","scalar",(void *)list,NULL,1},{"len","int",(void *)&len,NULL,0},{"i","int",(void *)&i,NULL,0},{"a","coord",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\n}"," \n// #line 857 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\n{\n      int j = 0;\n      {forin (scalar, s , list)\n v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z); endforin()}\n    }"}) {   
      
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 { _stencil_val(s,0,0,0); _stencil_interpolate_linear (point, s, NULL, NULL, NULL);    }}}
    }end_foreach_point_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:v[:len])){
#line 857
foreach_point (a[i].x, a[i].y, a[i].z) {
      int j = 0;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);}}
    }end_foreach_point();mpi_all_reduce_array(v,double,MPI_MIN,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 862
v = w;
  }
end_tracing("interpolate_array","/home/Min/basilisk/src/grid/cartesian-common.h",864);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 885
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return 1e30;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/Min/basilisk/src/grid/cartesian-common.h", 924, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 2; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 2; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 986
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 2; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}




#define dimensional(...)

#define show_dimension_internal(...)
#define display_value(...)
#define interpreter_verbosity(...)
#line 4 "/home/Min/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  } 
#line 44
{




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = 1e30;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    {foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  {foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  {foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    {foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/home/Min/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/Min/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/Min/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static void multigrid_setup_vector (vector v)
{
   {
    _attribute[v.x.i].prolongation = refine_bilinear;
    _attribute[v.x.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.y.i].prolongation = refine_bilinear;
    _attribute[v.y.i].restriction = restriction_average;
  }
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 260
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  
    multigrid_setup_vector (t.x);
    
#line 269
multigrid_setup_vector (t.y);
  return t;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 294 "/home/Min/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 325 "/home/Min/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 347 "/home/Min/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 385 "/home/Min/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector = multigrid_init_vector;
  init_face_vector = multigrid_init_face_vector;
  init_tensor = multigrid_init_tensor;
  restriction = multigrid_restriction;
  debug = multigrid_debug;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil(1,{(NonLocal[]){{"size","scalar",(void *)&size,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 454 \"/home/Min/basilisk/src/grid/multigrid-common.h\"\nval_out_(size,0,0,0) = 1;"})
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 453
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), l); };
  }
}
#line 844 "/home/Min/basilisk/src/grid/multigrid.h"

void dimensions (int nx, int ny, int nz)
{





}
#line 2 "reactor.c"
#line 1 "embed.h"
#line 1 "/home/Min/basilisk/src/embed.h"
#line 12 "/home/Min/basilisk/src/embed.h"
#line 1 "fractions.h"
#line 1 "/home/Min/basilisk/src/fractions.h"
#line 12 "/home/Min/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/Min/basilisk/src/geometry.h"
#line 35 "/home/Min/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double __tmp = n1; n1 = n2; n2 = __tmp; } while(0);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 163 "/home/Min/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 267 "/home/Min/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 270
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
  return line_area(n1.x, n1.y, alpha);
}
#line 292 "/home/Min/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
      
#line 297
if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 382 "/home/Min/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 388
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
    
#line 399
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 424
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}
#line 512 "/home/Min/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 518
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
    
#line 535
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  } 
#line 543
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }

   {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 551
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
#line 13 "/home/Min/basilisk/src/fractions.h"





#line 1 "myc2d.h"
#line 1 "/home/Min/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l - c_r);
  my0 = 0.5*(c_b - c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;
}
#line 13 "/home/Min/basilisk/src/fractions.h"





#line 1 "myc2d.h"
#line 1 "/home/Min/basilisk/src/myc2d.h"





static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
  
   


_stencil_val(c,-1,1,0); _stencil_val(c,0,1,0); _stencil_val(c,1,1,0); 


     
#line 14
_stencil_val(c,-1,-1,0); _stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0); 
     _stencil_val(c,1,-1,0); _stencil_val(c,1,0,0); _stencil_val(c,1,1,0); 
     _stencil_val(c,-1,-1,0); _stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0);
            
      
   
            
     
       
    



     
      





_stencil_val(c,-1,-1,0);_stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0); 


     
  


      
#line 35
_stencil_val(c,1,-1,0);_stencil_val(c,1,0,0); _stencil_val(c,1,1,0);  
     
        _stencil_val(c,-1,-1,0);_stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0);
       _stencil_val(c,-1,1,0);_stencil_val(c,0,1,0); _stencil_val(c,1,1,0);    
        
        
    
      
       
       
   
    
      
       


  
        
        
    
      
       
       
   



      
  

  
#line 64
return ;
}
#line 19 "/home/Min/basilisk/src/fractions.h"
#line 120 "/home/Min/basilisk/src/fractions.h"
     
void fractions (scalar Phi, scalar c,
  vector s, double val)
{tracing("fractions","/home/Min/basilisk/src/fractions.h",121);

  vector   as=(s).x.i>0?(s):new_face_vector("as");
#line 136 "/home/Min/basilisk/src/fractions.h"
  vector p;
  p.x = as.y; p.y = as.x;
#line 146 "/home/Min/basilisk/src/fractions.h"
  foreach_face_stencil(1,{(NonLocal[]){{"p","vector",(void *)&p,NULL,0},{"Phi","scalar",(void *)&Phi,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 146 \"/home/Min/basilisk/src/fractions.h\"\n{is_face_y(){ {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.x,0,0,0) = 1. - val(p.x,0,0,0);\n    }\n// # 171 \"/home/Min/basilisk/src/fractions.h\"\n    else\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);\n  }}end_is_face_y()\n// #line 146\nis_face_x(){ {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.y,0,0,0) = 1. - val(p.y,0,0,0);\n    }\n// # 171 \"/home/Min/basilisk/src/fractions.h\"\n    else\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);\n  }}end_is_face_x()}"}){_stencil_is_face_y(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);






      
#line 159
_stencil_val_a(p.x,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.x,0,0,0);_stencil_val_a(p.x,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);_stencil_val_a(p.x,0,0,0);       }}





           
#line 171 "/home/Min/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_y()
#line 146
_stencil_is_face_x(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);






      
#line 159
_stencil_val_a(p.y,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.y,0,0,0);_stencil_val_a(p.y,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);_stencil_val_a(p.y,0,0,0);       }}





           
#line 171 "/home/Min/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_x()}end_foreach_face_stencil();
#line 146 "/home/Min/basilisk/src/fractions.h"
  {foreach_face_generic(){is_face_y(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 171 "/home/Min/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }}end_is_face_y()
#line 146
is_face_x(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 171 "/home/Min/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }}end_is_face_x()}end_foreach_face_generic();}
#line 196 "/home/Min/basilisk/src/fractions.h"
  scalar s_z = c;
  foreach_stencil(1,{(NonLocal[]){{"Phi","scalar",(void *)&Phi,NULL,0},{"s_z","scalar",(void *)&s_z,NULL,0},{"p","vector",(void *)&p,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}","\n\n  \n// #line 199 \"/home/Min/basilisk/src/fractions.h\"\n{\n// # 231 \"/home/Min/basilisk/src/fractions.h\"\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 233\n{\n      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);\n      nn += fabs(n.y);\n    }\n\n\n\n\n\n    if (nn == 0.)\n      val_out_(s_z,0,0,0) = val(p.x,0,0,0);\n    else {\n\n\n\n\n\n      \n n.x /= nn;\n \n// #line 251\nn.y /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {\n     real a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);\n     alpha += n.x*a + n.y*(i - 0.5);\n     ni++;\n   }\n   \n// #line 261\nif (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {\n     real a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);\n     alpha += n.y*a + n.x*(i - 0.5);\n     ni++;\n   }}\n// # 274 \"/home/Min/basilisk/src/fractions.h\"\n      if (ni == 0)\n val_out_(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));\n      else if (ni != 4)\n val_out_(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);\n      else {\n\n\n\n val_out_(s_z,0,0,0) = 0.;\n\n      }\n    }\n  }"})

  {    
#line 231 "/home/Min/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.x,0,0,0);_stencil_val_a(s_z,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);_stencil_val_a(s_z,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{



 _stencil_val_a(s_z,0,0,0);  

      }}}
#line 274 "/home/Min/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}end_foreach_stencil();
  {
#line 197
foreach()

  {
#line 231 "/home/Min/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 233
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 251
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 274 "/home/Min/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  }end_foreach();}if(!(s).x.i)delete((scalar*)((vector[]){as,{{-1},{-1}}}));
#line 351 "/home/Min/basilisk/src/fractions.h"
end_tracing("fractions","/home/Min/basilisk/src/fractions.h",351);}
#line 395 "/home/Min/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/home/Min/basilisk/src/fractions.h", 399, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 400
{
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 408
n.y /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 423
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 429
n.y /= nn;}
    else
      {
 n.x = 1./2;
 
#line 432
n.y = 1./2;}
    return n;
  }
  return mycs (point, c);
}






#line 418
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 426
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 435
return;
}
#line 445 "/home/Min/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/Min/basilisk/src/fractions.h",446);
  foreach_stencil(1,{(NonLocal[]){{"n","vector",(void *)&n,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}"," \n// #line 448 \"/home/Min/basilisk/src/fractions.h\"\n{\n\n\n\n\n\n    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {\n      val_out_(alpha,0,0,0) = 0.;\n      \n val_out_(n.x,0,0,0) = 0.;\n \n// #line 457\nval_out_(n.y,0,0,0) = 0.;\n    }\n    else {\n\n\n\n\n\n\n      coord m = mycs (point, c);\n      \n val_out_(n.x,0,0,0) = m.x;\n \n// #line 468\nval_out_(n.y,0,0,0) = m.y;\n      val_out_(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);\n    }\n  }"}) {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 457
{_stencil_val_a(n.y,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 468
{_stencil_val_a(n.y,0,0,0);  }
_stencil_val(c,0,0,0);
      
#line 469
_stencil_val_a(alpha,0,0,0);    
    }}





          
    
  
#line 471
}end_foreach_stencil();
  {
#line 448
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 457
val(n.y,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 468
val(n.y,0,0,0) = m.y;
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}
#line 489 "/home/Min/basilisk/src/fractions.h"
end_tracing("reconstruction","/home/Min/basilisk/src/fractions.h",489);}
#line 509 "/home/Min/basilisk/src/fractions.h"
     
void output_facets (scalar c, FILE * fp, vector s)
{tracing("output_facets","/home/Min/basilisk/src/fractions.h",510);
  foreach_stencil(1,{(NonLocal[]){{"fp","not implemented yet",(void *)fp,NULL,1},{"s","vector",(void *)&s,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 292 \"/home/Min/basilisk/src/geometry.h\"\nint facets (coord n, real alpha, coord p[2])\n{\n  int i = 0;\n  for (real s = -0.5; s <= 0.5; s += 1.)\n    {\n      if (fabs (n.y) > 1e-4 && i < 2) {\n real a = (alpha - s*n.x)/n.y;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].x = s;\n   p[i++].y = a;\n }\n      }\n      \n// #line 297\nif (fabs (n.x) > 1e-4 && i < 2) {\n real a = (alpha - s*n.y)/n.x;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].y = s;\n   p[i++].x = a;\n }\n      }}\n  return i;\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/Min/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}","\n    \n// #line 513 \"/home/Min/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = facet_normal (point, c, s);\n      real alpha = line_alpha (val(c,0,0,0), n);\n\n\n\n      coord segment[2];\n      if (facets (n, alpha, segment) == 2)\n fprintf (fp, \"%g %g\\n%g %g\\n\\n\",\n   x + segment[0].x*Delta, y + segment[0].y*Delta,\n   x + segment[1].x*Delta, y + segment[1].y*Delta);\n// # 533 \"/home/Min/basilisk/src/fractions.h\"\n    }"})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0); 



            
           
        
     
 
#line 533 "/home/Min/basilisk/src/fractions.h"
    }        }end_foreach_stencil();
  {
#line 512
foreach()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);



      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 533 "/home/Min/basilisk/src/fractions.h"
    }end_foreach();}

  fflush (fp);
end_tracing("output_facets","/home/Min/basilisk/src/fractions.h",536);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/Min/basilisk/src/fractions.h",545);
  double area = 0.;
  foreach_stencil (1,{(NonLocal[]){{"area","double",(void *)&area,NULL,0,'+'},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}","\n    \n// #line 549 \"/home/Min/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = mycs (point, c), p;\n      real alpha = line_alpha (val(c,0,0,0), n);\n      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);\n    }"})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
          
    }        }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)){
#line 548
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 554
{end_tracing("interface_area","/home/Min/basilisk/src/fractions.h",554);return area;}
end_tracing("interface_area","/home/Min/basilisk/src/fractions.h",555);}
#line 13 "/home/Min/basilisk/src/embed.h"






scalar  cs={0};
vector  fs={{1},{2}};

double (* metric_embed_factor) (Point, coord) = NULL;
#line 70 "/home/Min/basilisk/src/embed.h"

static inline double embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));
  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert ("/home/Min/basilisk/src/embed.h", 74, "cs[i] && cs[i-1]");
  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))
    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +
     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);
  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;
}

#line 71
static inline double embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));
  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert ("/home/Min/basilisk/src/embed.h", 74, "cs[i] && cs[i-1]");
  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))
    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +
     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);
  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;
}
#line 70 "/home/Min/basilisk/src/embed.h"

static void _stencil_embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.x,i,-1,0);_stencil_val(fs.x,i,1,0);
_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);
_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.y,i,o_stencil,0    ); _stencil_val(fs.y,i-1,o_stencil,0    ); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);
    { _stencil_val(fs.x,i,0,0);_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0); 
_stencil_val(fs.x,i,0,0);_stencil_val(a,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0);    
       
#line 77
}
_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0);  
      
         
  
#line 78
return  ;
}

#line 71
static void _stencil_embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.y,-1,i,0);_stencil_val(fs.y,1,i,0);
_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);
_stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0);
    { _stencil_val(fs.y,0,i,0);_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0); 
_stencil_val(fs.y,0,i,0);_stencil_val(a,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0);    
       
#line 77
}
_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0);  
      
         
  
#line 78
return  ;
}


static inline double embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));
  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?
    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :
    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));
}

#line 82
static inline double embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));
  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?
    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :
    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));
}



#line 82
static void _stencil_embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.x,i,-1,0);_stencil_val(fs.x,i,1,0);
_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.y,i,o_stencil,0    ); _stencil_val(fs.y,i-1,o_stencil,0    ); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0); 
_stencil_val(fs.x,i,0,0);_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0); _stencil_val(fs.x,i,0,0);_stencil_val(a,i,o_stencil,0); _stencil_val(cs,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);_stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);
_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);
  
#line 85
return                  
                    
    ;
}

#line 82
static void _stencil_embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.y,-1,i,0);_stencil_val(fs.y,1,i,0);
_stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0); 
_stencil_val(fs.y,0,i,0);_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0); _stencil_val(fs.y,0,i,0);_stencil_val(a,o_stencil,i,0); _stencil_val(cs,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0); _stencil_val(cs,o_stencil,i-1,0);_stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0);
_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);
  
#line 85
return                  
                    
    ;
}
#line 177 "/home/Min/basilisk/src/embed.h"

#line 222 "/home/Min/basilisk/src/embed.h"
static inline
double embed_geometry (Point point, coord * p, coord * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  *n = facet_normal (point, cs, fs);
  double alpha = line_alpha (val(cs,0,0,0), *n);
  double area = line_length_center(*n,alpha,p);
  normalize (n);
  return area;
}
#line 222 "/home/Min/basilisk/src/embed.h"
static void 
_stencil_embed_geometry (Point point,_stencil_undefined * p,_stencil_undefined * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
_stencil_facet_normal (point, cs, fs);     
   
  
#line 226
_stencil_val(cs,0,0,0);   
   
  
  return ;
}





static inline
double embed_area_center (Point point, double * x1, double * y1, double * z1)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double area = 0.;
  if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {
    coord n, p;
    area = embed_geometry (point, &p, &n);
    *x1 += p.x*Delta, *y1 += p.y*Delta, *z1 += p.z*Delta;
  }
  return area;
}
#line 255 "/home/Min/basilisk/src/embed.h"
double embed_interpolate (Point point, scalar s, coord p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!(2 == 2)) qassert ("/home/Min/basilisk/src/embed.h", 257, "dimension == 2");
  int i = sign(p.x), j = sign(p.y);
  if (val(cs,i,0,0) && val(cs,0,j,0) && val(cs,i,j,0))

    return ((val(s,0,0,0)*(1. - fabs(p.x)) + val(s,i,0,0)*fabs(p.x))*(1. - fabs(p.y)) +
     (val(s,0,j,0)*(1. - fabs(p.x)) + val(s,i,j,0)*fabs(p.x))*fabs(p.y));
  else {


    double val = val(s,0,0,0);
     {
      int i = sign(p.x);
      if (val(cs,i,0,0))
 val += fabs(p.x)*(val(s,i,0,0) - val(s,0,0,0));
      else if (val(cs,-i,0,0))
 val += fabs(p.x)*(val(s,0,0,0) - val(s,-i,0,0));
    } 
#line 267
{
      int i = sign(p.y);
      if (val(cs,0,i,0))
 val += fabs(p.y)*(val(s,0,i,0) - val(s,0,0,0));
      else if (val(cs,0,-i,0))
 val += fabs(p.y)*(val(s,0,0,0) - val(s,0,-i,0));
    }
    return val;
  }
}
#line 255 "/home/Min/basilisk/src/embed.h"
static void _stencil_embed_interpolate (Point point, scalar s,_stencil_undefined * p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;         
      
  
_stencil_val(cs,o_stencil,0,0); _stencil_val(cs,0,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);{

    {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);
_stencil_val(s,0,o_stencil,0); _stencil_val(s,o_stencil,o_stencil,0);         
      
#line 262
} 
{  


     _stencil_val(s,0,0,0);
     {   
      
_stencil_val(cs,o_stencil,0,0);{
 {_stencil_val(s,o_stencil,0,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,o_stencil,0,0);
 {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);   } }}
       
      
    
#line 273
} 
#line 267
{   
      
_stencil_val(cs,0,o_stencil,0);{
 {_stencil_val(s,0,o_stencil,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,0,o_stencil,0);
 {_stencil_val(s,0,0,0);_stencil_val(s,0, o_stencil,0);   } }}
       
      
    
#line 273
} 
    
  }}
     
  

#line 276
}
#line 285 "/home/Min/basilisk/src/embed.h"
struct Cleanup {
  scalar c;
  vector s;
  double smin;
  bool opposite;
};

     
int fractions_cleanup (scalar c, vector s,
         double smin, bool opposite)
{tracing("fractions_cleanup","/home/Min/basilisk/src/embed.h",293);







  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {




    foreach_face_stencil(1,{(NonLocal[]){{"smin","double",(void *)&smin,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"s","vector",(void *)&s,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 309 \"/home/Min/basilisk/src/embed.h\"\n{is_face_x(){\n      if (val(s.x,0,0,0) && ((!val(c,0,0,0) || !val(c,-1,0,0)) || val(s.x,0,0,0) < smin))\n val_out_(s.x,0,0,0) = 0.;}end_is_face_x()\n// #line 309\nis_face_y(){\n      if (val(s.y,0,0,0) && ((!val(c,0,0,0) || !val(c,0,-1,0)) || val(s.y,0,0,0) < smin))\n val_out_(s.y,0,0,0) = 0.;}end_is_face_y()}"}){_stencil_is_face_x(){
      {_stencil_val(s.x,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,-1,0,0); _stencil_val(s.x,0,0,0);
 {_stencil_val_a(s.x,0,0,0);  }        }}end__stencil_is_face_x()
#line 309
_stencil_is_face_y(){
      {_stencil_val(s.y,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,0,-1,0); _stencil_val(s.y,0,0,0);
 {_stencil_val_a(s.y,0,0,0);  }        }}end__stencil_is_face_y()}end_foreach_face_stencil();




    {
#line 309
foreach_face_generic(){is_face_x(){
      if (val(s.x,0,0,0) && ((!val(c,0,0,0) || !val(c,-1,0,0)) || val(s.x,0,0,0) < smin))
 val(s.x,0,0,0) = 0.;}end_is_face_x()
#line 309
is_face_y(){
      if (val(s.y,0,0,0) && ((!val(c,0,0,0) || !val(c,0,-1,0)) || val(s.y,0,0,0) < smin))
 val(s.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}

    changed = 0;
    foreach_stencil(1,{(NonLocal[]){{"changed","int",(void *)&changed,NULL,0,'+'},{"opposite","bool",(void *)&opposite,NULL,0},{"s","vector",(void *)&s,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 315 \"/home/Min/basilisk/src/embed.h\"\nif (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {\n int n = 0;\n  {\n   for (int i = 0; i <= 1; i++)\n     if (val(s.x,i,0,0) > 0.)\n       n++;\n// # 330 \"/home/Min/basilisk/src/embed.h\"\n   if (opposite && val(s.x,0,0,0) == 0. && val(s.x,1,0,0) == 0.)\n     val_out_(c,0,0,0) = 0., changed++;\n } \n// #line 317\n{\n   for (int i = 0; i <= 1; i++)\n     if (val(s.y,0,i,0) > 0.)\n       n++;\n// # 330 \"/home/Min/basilisk/src/embed.h\"\n   if (opposite && val(s.y,0,0,0) == 0. && val(s.y,0,1,0) == 0.)\n     val_out_(c,0,0,0) = 0., changed++;\n }\n\n\n\n\n\n\n\n if (n < 2)\n   val_out_(c,0,0,0) = 0., changed++;\n      }"})
      {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
 
  {
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.x,i,0,0);
          } 









_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);
     {_stencil_val_a(c,0,0,0);   }
#line 330 "/home/Min/basilisk/src/embed.h"
          
 
} 
#line 317
{
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.y,0,i,0);
          } 









_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);
     {_stencil_val_a(c,0,0,0);   }
#line 330 "/home/Min/basilisk/src/embed.h"
          
 
}
   







{_stencil_val_a(c,0,0,0);   }







    
      
#line 342
}      }end_foreach_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:changed)){
#line 314
foreach()
      if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 int n = 0;
  {
   for (int i = 0; i <= 1; i++)
     if (val(s.x,i,0,0) > 0.)
       n++;
#line 330 "/home/Min/basilisk/src/embed.h"
   if (opposite && val(s.x,0,0,0) == 0. && val(s.x,1,0,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 } 
#line 317
{
   for (int i = 0; i <= 1; i++)
     if (val(s.y,0,i,0) > 0.)
       n++;
#line 330 "/home/Min/basilisk/src/embed.h"
   if (opposite && val(s.y,0,0,0) == 0. && val(s.y,0,1,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 }







 if (n < 2)
   val(c,0,0,0) = 0., changed++;
      }end_foreach();mpi_all_reduce_array(&changed,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

    
#line 344
schanged += changed;
  }
  if (changed)
    fprintf (ferr, "WARNING: fractions_cleanup() did not converge after "
      "%d iterations\n", i);
  {end_tracing("fractions_cleanup","/home/Min/basilisk/src/embed.h",349);return schanged;}
end_tracing("fractions_cleanup","/home/Min/basilisk/src/embed.h",350);}
#line 374 "/home/Min/basilisk/src/embed.h"

static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  
    n.x = - n.x;
    
#line 380
n.y = - n.y;
  double d[2], v[2] = {1e30,1e30};
  bool defined = true;
  
    if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;
    
#line 384
if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;

      if (val(fs.x,i + (i < 0),j,0) && val(fs.y,i,j,0) && val(fs.y,i,j+1,0) &&
   val(cs,i,j-1,0) && val(cs,i,j,0) && val(cs,i,j+1,0))
 v[l] = ((((val(s,i,j-1,0)))*((y1) - 1.) + ((val(s,i,j+1,0)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,0)))*((y1) - 1.)*((y1) + 1.));
#line 418 "/home/Min/basilisk/src/embed.h"
      else
 break;
    }
  if (v[0] == 1e30) {





    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != 1e30)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}

#line 375
static inline double dirichlet_gradient_y (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  
    n.y = - n.y;
    
#line 380
n.x = - n.x;
  double d[2], v[2] = {1e30,1e30};
  bool defined = true;
  
    if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;
    
#line 384
if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.y);
      d[l] = (i - p.y)/n.y;
      double y1 = p.x + d[l]*n.x;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;

      if (val(fs.y,j,i + (i < 0),0) && val(fs.x,j,i,0) && val(fs.x,j+1,i,0) &&
   val(cs,j-1,i,0) && val(cs,j,i,0) && val(cs,j+1,i,0))
 v[l] = ((((val(s,j-1,i,0)))*((y1) - 1.) + ((val(s,j+1,i,0)))*((y1) + 1.))*(y1)/2. - ((val(s,j,i,0)))*((y1) - 1.)*((y1) + 1.));
#line 418 "/home/Min/basilisk/src/embed.h"
      else
 break;
    }
  if (v[0] == 1e30) {





    d[0] = max(1e-3, fabs(p.y/n.y));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != 1e30)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}
#line 374 "/home/Min/basilisk/src/embed.h"

static void _stencil_dirichlet_gradient_x (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;       
  
       
  
  
  
    {_stencil_val(fs.x,o_stencil,0,0  ); 
          }
    
#line 384
{_stencil_val(fs.y,0,o_stencil,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                         
       
         
      
      
        

_stencil_val(fs.x,    o_stencil,o_stencil,0); _stencil_val(fs.y,o_stencil,o_stencil,0); _stencil_val(fs.y,o_stencil,o_stencil,0);
   _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);
#line 394
{
 
{_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);              }
 
#line 419
}

            
#line 418 "/home/Min/basilisk/src/embed.h"
      
    
}         
     
   
     





   
     
  







return   ;
}

#line 375
static void _stencil_dirichlet_gradient_y (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;       
  
       
  
  
  
    {_stencil_val(fs.y,0,o_stencil,0  ); 
          }
    
#line 384
{_stencil_val(fs.x,o_stencil,0,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                         
       
         
      
      
        

_stencil_val(fs.y,o_stencil,    o_stencil,0); _stencil_val(fs.x,o_stencil,o_stencil,0); _stencil_val(fs.x,o_stencil,o_stencil,0);
   _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);
#line 394
{
 
{_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);              }
 
#line 419
}

            
#line 418 "/home/Min/basilisk/src/embed.h"
      
    
}         
     
   
     





   
     
  







return   ;
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
      coord n, coord p, double bc, double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
    
#line 447
if (fabs(n.y) >= fabs(n.x))
      return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
#line 458 "/home/Min/basilisk/src/embed.h"
  return 1e30;
}


#line 442
static void _stencil_dirichlet_gradient (Point point, scalar s, scalar cs,
_stencil_undefined * n,_stencil_undefined * p,_stencil_undefined * bc,_stencil_undefined * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  
      
{ _stencil_dirichlet_gradient_x (point, s, cs,NULL ,NULL ,NULL ,NULL );}
      
#line 448
{ _stencil_dirichlet_gradient_y (point, s, cs,NULL ,NULL ,NULL ,NULL );}
       
#line 458 "/home/Min/basilisk/src/embed.h"
  return ;
}

bid embed;
#line 470 "/home/Min/basilisk/src/embed.h"
static inline
coord embed_gradient (Point point, vector u, coord p, coord n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord dudn;
   {
    bool dirichlet = false;
    double vb = _attribute[u.x.i].boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
    }
    else
      dudn.x = vb;
    if (dudn.x == 1e30)
      dudn.x = 0.;
  } 
#line 474
{
    bool dirichlet = false;
    double vb = _attribute[u.y.i].boundary[embed] (point, point, u.y, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.y = dirichlet_gradient (point, u.y, cs, n, p, vb, &val);
    }
    else
      dudn.y = vb;
    if (dudn.y == 1e30)
      dudn.y = 0.;
  }
  return dudn;
}
#line 470 "/home/Min/basilisk/src/embed.h"
static void 
_stencil_embed_gradient (Point point, vector u,_stencil_undefined * p,_stencil_undefined * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
  
   {     
    
    default_stencil ( point,((scalar[]){ u.x,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.x, cs,NULL ,NULL ,NULL ,NULL ); 
      
    
#line 480
}   
     
    
        
     
       
  
#line 485
} 
#line 474
{     
    
    default_stencil ( point,((scalar[]){ u.y,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.y, cs,NULL ,NULL ,NULL ,NULL ); 
      
    
#line 480
}   
     
    
        
     
       
  
#line 485
}
  return ;
}
#line 508 "/home/Min/basilisk/src/embed.h"
     
void embed_force (scalar p, vector u, vector mu, coord * Fp, coord * Fmu)
{tracing("embed_force","/home/Min/basilisk/src/embed.h",509);
  coord Fps = {0}, Fmus = {0};
  if(!is_constant(fm.x)){
  
#line 512
foreach_stencil (1,{(NonLocal[]){{"Fmus","coord",(void *)&Fmus,NULL,0,'+'},{"u","vector",(void *)&u,NULL,0},{"val","double",(void *)&val,NULL,0},{"embed","bid",(void *)&embed,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"Fps","coord",(void *)&Fps,NULL,0,'+'},{"p","scalar",(void *)&p,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"fs","vector",(void *)&fs,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 375 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real dirichlet_gradient_y (Point point, scalar s, scalar cs,\n        coord n, coord p, real bc,\n        real * coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  \n    n.y = - n.y;\n    \n// #line 380\nn.x = - n.x;\n  real d[2], v[2] = {1e30,1e30};\n  bool defined = true;\n  \n    if (defined && !val(fs.y,0,(n.y > 0.),0))\n      defined = false;\n    \n// #line 384\nif (defined && !val(fs.x,(n.x > 0.),0,0))\n      defined = false;\n  if (defined)\n    for (int l = 0; l <= 1; l++) {\n      int i = (l + 1)*sign(n.y);\n      d[l] = (i - p.y)/n.y;\n      real y1 = p.x + d[l]*n.x;\n      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;\n      y1 -= j;\n\n      if (val(fs.y,j,i + (i < 0),0) && val(fs.x,j,i,0) && val(fs.x,j+1,i,0) &&\n   val(cs,j-1,i,0) && val(cs,j,i,0) && val(cs,j+1,i,0))\n v[l] = ((((val(s,j-1,i,0)))*((y1) - 1.) + ((val(s,j+1,i,0)))*((y1) + 1.))*(y1)/2. - ((val(s,j,i,0)))*((y1) - 1.)*((y1) + 1.));\n// # 418 \"/home/Min/basilisk/src/embed.h\"\n      else\n break;\n    }\n  if (v[0] == 1e30) {\n\n\n\n\n\n    d[0] = max(1e-3, fabs(p.y/n.y));\n    *coef = - 1./(d[0]*Delta);\n    return bc/(d[0]*Delta);\n  }\n\n\n\n\n\n  *coef = 0.;\n  if (v[1] != 1e30)\n    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);\n  return (bc - v[0])/(d[0]*Delta);\n}\n// # 374 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real dirichlet_gradient_x (Point point, scalar s, scalar cs,\n        coord n, coord p, real bc,\n        real * coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  \n    n.x = - n.x;\n    \n// #line 380\nn.y = - n.y;\n  real d[2], v[2] = {1e30,1e30};\n  bool defined = true;\n  \n    if (defined && !val(fs.x,(n.x > 0.),0,0))\n      defined = false;\n    \n// #line 384\nif (defined && !val(fs.y,0,(n.y > 0.),0))\n      defined = false;\n  if (defined)\n    for (int l = 0; l <= 1; l++) {\n      int i = (l + 1)*sign(n.x);\n      d[l] = (i - p.x)/n.x;\n      real y1 = p.y + d[l]*n.y;\n      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;\n      y1 -= j;\n\n      if (val(fs.x,i + (i < 0),j,0) && val(fs.y,i,j,0) && val(fs.y,i,j+1,0) &&\n   val(cs,i,j-1,0) && val(cs,i,j,0) && val(cs,i,j+1,0))\n v[l] = ((((val(s,i,j-1,0)))*((y1) - 1.) + ((val(s,i,j+1,0)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,0)))*((y1) - 1.)*((y1) + 1.));\n// # 418 \"/home/Min/basilisk/src/embed.h\"\n      else\n break;\n    }\n  if (v[0] == 1e30) {\n\n\n\n\n\n    d[0] = max(1e-3, fabs(p.x/n.x));\n    *coef = - 1./(d[0]*Delta);\n    return bc/(d[0]*Delta);\n  }\n\n\n\n\n\n  *coef = 0.;\n  if (v[1] != 1e30)\n    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);\n  return (bc - v[0])/(d[0]*Delta);\n}\n\n\n// #line 442 \"/home/Min/basilisk/src/embed.h\"\nreal dirichlet_gradient (Point point, scalar s, scalar cs,\n      coord n, coord p, real bc, inout real  coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n  \n    if (fabs(n.x) >= fabs(n.y))\n      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);\n    \n// #line 447\nif (fabs(n.y) >= fabs(n.x))\n      return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);\n// # 458 \"/home/Min/basilisk/src/embed.h\"\n  return 1e30;\n}\n// # 470 \"/home/Min/basilisk/src/embed.h\"\nstatic inline\ncoord embed_gradient (Point point, vector u, coord p, coord n)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  coord dudn;\n   {\n    bool dirichlet = false;\n    real vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);\n    if (dirichlet) {\n      real val;\n      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, val);\n    }\n    else\n      dudn.x = vb;\n    if (dudn.x == 1e30)\n      dudn.x = 0.;\n  } \n// #line 474\n{\n    bool dirichlet = false;\n    real vb = u.y.boundary[embed] (point, point, u.y, &dirichlet);\n    if (dirichlet) {\n      real val;\n      dudn.y = dirichlet_gradient (point, u.y, cs, n, p, vb, val);\n    }\n    else\n      dudn.y = vb;\n    if (dudn.y == 1e30)\n      dudn.y = 0.;\n  }\n  return dudn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 255 \"/home/Min/basilisk/src/embed.h\"\nreal embed_interpolate (Point point, scalar s, coord p)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (!(2 == 2)) qassert (\"/home/Min/basilisk/src/embed.h\", 257, \"dimension == 2\");\n  int i = sign(p.x), j = sign(p.y);\n  if (val(cs,i,0,0) && val(cs,0,j,0) && val(cs,i,j,0))\n\n    return ((val(s,0,0,0)*(1. - fabs(p.x)) + val(s,i,0,0)*fabs(p.x))*(1. - fabs(p.y)) +\n     (val(s,0,j,0)*(1. - fabs(p.x)) + val(s,i,j,0)*fabs(p.x))*fabs(p.y));\n  else {\n\n\n    real val = val(s,0,0,0);\n     {\n      int i = sign(p.x);\n      if (val(cs,i,0,0))\n val += fabs(p.x)*(val(s,i,0,0) - val(s,0,0,0));\n      else if (val(cs,-i,0,0))\n val += fabs(p.x)*(val(s,0,0,0) - val(s,-i,0,0));\n    } \n// #line 267\n{\n      int i = sign(p.y);\n      if (val(cs,0,i,0))\n val += fabs(p.y)*(val(s,0,i,0) - val(s,0,0,0));\n      else if (val(cs,0,-i,0))\n val += fabs(p.y)*(val(s,0,0,0) - val(s,0,-i,0));\n    }\n    return val;\n  }\n}\n\nOMP(omp declare reduction (+ : coord :\n      omp_out.x += omp_in.x,\n      omp_out.y += omp_in.y,\n      omp_out.z += omp_in.z))\n// # 920 \"/home/Min/basilisk/src/common.h\"\nvoid normalize (coord * n)\n{\n  real norm = 0.;\n  \n    norm += sq(n->x);\n    \n// #line 924\nnorm += sq(n->y);\n  norm = sqrt(norm);\n  \n    n->x /= norm;\n    \n// #line 927\nn->y /= norm;\n}\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/Min/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}\n// # 177 \"/home/Min/basilisk/src/embed.h\"\n\n// # 222 \"/home/Min/basilisk/src/embed.h\"\nstatic inline\nreal embed_geometry (Point point, coord * p, coord * n)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  *n = facet_normal (point, cs, fs);\n  real alpha = line_alpha (val(cs,0,0,0), *n);\n  real area = line_length_center(*n,alpha,p);\n  normalize (n);\n  return area;\n}","\n    \n// #line 513 \"/home/Min/basilisk/src/embed.h\"\nif (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {\n\n\n\n\n\n\n      coord n, b;\n      real area = embed_geometry (point, &b, &n);\n      area *= pow (Delta, 2 - 1);\n      real Fn = area*embed_interpolate (point, p, b);\n      \n Fps.x += Fn*n.x;\n \n// #line 525\nFps.y += Fn*n.y;\n// # 534 \"/home/Min/basilisk/src/embed.h\"\n      if (constant(mu.x) != 0.) {\n real mua = 0., fa = 0.;\n  {\n   mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);\n   fa += val(fm.x,0,0,0) + val(fm.x,1,0,0);\n } \n// #line 536\n{\n   mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);\n   fa += val(fm.y,0,0,0) + val(fm.y,0,1,0);\n }\n mua /= fa;\n// # 599 \"/home/Min/basilisk/src/embed.h\"\n if (!(2 == 2)) qassert (\"/home/Min/basilisk/src/embed.h\", 599, \"dimension == 2\");\n coord dudn = embed_gradient (point, u, b, n);\n \n   Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);\n   \n// #line 602\nFmus.y -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.y*n.x);\n      }\n    }"})
    {_stencil_val(cs,0,0,0); _stencil_val(cs,0,0,0); {    






      
       _stencil_embed_geometry (point,NULL ,NULL );   
            
      _stencil_embed_interpolate (point, p,NULL );
       
  
#line 534 "/home/Min/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {      
 
  { 
_stencil_val(mu.x,0,0,0); _stencil_val(mu.x,1,0,0); 
     _stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0); 
    
 
#line 539
} 
#line 536
{ 
_stencil_val(mu.y,0,0,0); _stencil_val(mu.y,0,1,0); 
     _stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0); 
    
 
#line 539
}  
      
#line 599 "/home/Min/basilisk/src/embed.h"
     
  _stencil_embed_gradient (point, u,NULL ,NULL );
  
        
      }
    }      }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:Fmus)reduction(+:Fps)){
#line 512
foreach ()
    if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {






      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, 2 - 1);
      double Fn = area*embed_interpolate (point, p, b);
      
 Fps.x += Fn*n.x;
 
#line 525
Fps.y += Fn*n.y;
#line 534 "/home/Min/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {
 double mua = 0., fa = 0.;
  {
   mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
   fa += val(fm.x,0,0,0) + val(fm.x,1,0,0);
 } 
#line 536
{
   mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
   fa += val(fm.y,0,0,0) + val(fm.y,0,1,0);
 }
 mua /= fa;
#line 599 "/home/Min/basilisk/src/embed.h"
 if (!(2 == 2)) qassert ("/home/Min/basilisk/src/embed.h", 599, "dimension == 2");
 coord dudn = embed_gradient (point, u, b, n);
 
   Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
   
#line 602
Fmus.y -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.y*n.x);
      }
    }end_foreach();mpi_all_reduce_array(&Fmus.x,double,MPI_SUM,2);mpi_all_reduce_array(&Fps.x,double,MPI_SUM,2);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 604
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 512
foreach_stencil (1,{(NonLocal[]){{"Fmus","coord",(void *)&Fmus,NULL,0,'+'},{"u","vector",(void *)&u,NULL,0},{"val","double",(void *)&val,NULL,0},{"embed","bid",(void *)&embed,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"Fps","coord",(void *)&Fps,NULL,0,'+'},{"p","scalar",(void *)&p,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"fs","vector",(void *)&fs,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 375 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real dirichlet_gradient_y (Point point, scalar s, scalar cs,\n        coord n, coord p, real bc,\n        real * coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  \n    n.y = - n.y;\n    \n// #line 380\nn.x = - n.x;\n  real d[2], v[2] = {1e30,1e30};\n  bool defined = true;\n  \n    if (defined && !val(fs.y,0,(n.y > 0.),0))\n      defined = false;\n    \n// #line 384\nif (defined && !val(fs.x,(n.x > 0.),0,0))\n      defined = false;\n  if (defined)\n    for (int l = 0; l <= 1; l++) {\n      int i = (l + 1)*sign(n.y);\n      d[l] = (i - p.y)/n.y;\n      real y1 = p.x + d[l]*n.x;\n      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;\n      y1 -= j;\n\n      if (val(fs.y,j,i + (i < 0),0) && val(fs.x,j,i,0) && val(fs.x,j+1,i,0) &&\n   val(cs,j-1,i,0) && val(cs,j,i,0) && val(cs,j+1,i,0))\n v[l] = ((((val(s,j-1,i,0)))*((y1) - 1.) + ((val(s,j+1,i,0)))*((y1) + 1.))*(y1)/2. - ((val(s,j,i,0)))*((y1) - 1.)*((y1) + 1.));\n// # 418 \"/home/Min/basilisk/src/embed.h\"\n      else\n break;\n    }\n  if (v[0] == 1e30) {\n\n\n\n\n\n    d[0] = max(1e-3, fabs(p.y/n.y));\n    *coef = - 1./(d[0]*Delta);\n    return bc/(d[0]*Delta);\n  }\n\n\n\n\n\n  *coef = 0.;\n  if (v[1] != 1e30)\n    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);\n  return (bc - v[0])/(d[0]*Delta);\n}\n// # 374 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real dirichlet_gradient_x (Point point, scalar s, scalar cs,\n        coord n, coord p, real bc,\n        real * coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  \n    n.x = - n.x;\n    \n// #line 380\nn.y = - n.y;\n  real d[2], v[2] = {1e30,1e30};\n  bool defined = true;\n  \n    if (defined && !val(fs.x,(n.x > 0.),0,0))\n      defined = false;\n    \n// #line 384\nif (defined && !val(fs.y,0,(n.y > 0.),0))\n      defined = false;\n  if (defined)\n    for (int l = 0; l <= 1; l++) {\n      int i = (l + 1)*sign(n.x);\n      d[l] = (i - p.x)/n.x;\n      real y1 = p.y + d[l]*n.y;\n      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;\n      y1 -= j;\n\n      if (val(fs.x,i + (i < 0),j,0) && val(fs.y,i,j,0) && val(fs.y,i,j+1,0) &&\n   val(cs,i,j-1,0) && val(cs,i,j,0) && val(cs,i,j+1,0))\n v[l] = ((((val(s,i,j-1,0)))*((y1) - 1.) + ((val(s,i,j+1,0)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,0)))*((y1) - 1.)*((y1) + 1.));\n// # 418 \"/home/Min/basilisk/src/embed.h\"\n      else\n break;\n    }\n  if (v[0] == 1e30) {\n\n\n\n\n\n    d[0] = max(1e-3, fabs(p.x/n.x));\n    *coef = - 1./(d[0]*Delta);\n    return bc/(d[0]*Delta);\n  }\n\n\n\n\n\n  *coef = 0.;\n  if (v[1] != 1e30)\n    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);\n  return (bc - v[0])/(d[0]*Delta);\n}\n\n\n// #line 442 \"/home/Min/basilisk/src/embed.h\"\nreal dirichlet_gradient (Point point, scalar s, scalar cs,\n      coord n, coord p, real bc, inout real  coef)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n  \n    if (fabs(n.x) >= fabs(n.y))\n      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);\n    \n// #line 447\nif (fabs(n.y) >= fabs(n.x))\n      return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);\n// # 458 \"/home/Min/basilisk/src/embed.h\"\n  return 1e30;\n}\n// # 470 \"/home/Min/basilisk/src/embed.h\"\nstatic inline\ncoord embed_gradient (Point point, vector u, coord p, coord n)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  coord dudn;\n   {\n    bool dirichlet = false;\n    real vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);\n    if (dirichlet) {\n      real val;\n      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, val);\n    }\n    else\n      dudn.x = vb;\n    if (dudn.x == 1e30)\n      dudn.x = 0.;\n  } \n// #line 474\n{\n    bool dirichlet = false;\n    real vb = u.y.boundary[embed] (point, point, u.y, &dirichlet);\n    if (dirichlet) {\n      real val;\n      dudn.y = dirichlet_gradient (point, u.y, cs, n, p, vb, val);\n    }\n    else\n      dudn.y = vb;\n    if (dudn.y == 1e30)\n      dudn.y = 0.;\n  }\n  return dudn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 255 \"/home/Min/basilisk/src/embed.h\"\nreal embed_interpolate (Point point, scalar s, coord p)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (!(2 == 2)) qassert (\"/home/Min/basilisk/src/embed.h\", 257, \"dimension == 2\");\n  int i = sign(p.x), j = sign(p.y);\n  if (val(cs,i,0,0) && val(cs,0,j,0) && val(cs,i,j,0))\n\n    return ((val(s,0,0,0)*(1. - fabs(p.x)) + val(s,i,0,0)*fabs(p.x))*(1. - fabs(p.y)) +\n     (val(s,0,j,0)*(1. - fabs(p.x)) + val(s,i,j,0)*fabs(p.x))*fabs(p.y));\n  else {\n\n\n    real val = val(s,0,0,0);\n     {\n      int i = sign(p.x);\n      if (val(cs,i,0,0))\n val += fabs(p.x)*(val(s,i,0,0) - val(s,0,0,0));\n      else if (val(cs,-i,0,0))\n val += fabs(p.x)*(val(s,0,0,0) - val(s,-i,0,0));\n    } \n// #line 267\n{\n      int i = sign(p.y);\n      if (val(cs,0,i,0))\n val += fabs(p.y)*(val(s,0,i,0) - val(s,0,0,0));\n      else if (val(cs,0,-i,0))\n val += fabs(p.y)*(val(s,0,0,0) - val(s,0,-i,0));\n    }\n    return val;\n  }\n}\n\nOMP(omp declare reduction (+ : coord :\n      omp_out.x += omp_in.x,\n      omp_out.y += omp_in.y,\n      omp_out.z += omp_in.z))\n// # 920 \"/home/Min/basilisk/src/common.h\"\nvoid normalize (coord * n)\n{\n  real norm = 0.;\n  \n    norm += sq(n->x);\n    \n// #line 924\nnorm += sq(n->y);\n  norm = sqrt(norm);\n  \n    n->x /= norm;\n    \n// #line 927\nn->y /= norm;\n}\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/Min/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}\n// # 177 \"/home/Min/basilisk/src/embed.h\"\n\n// # 222 \"/home/Min/basilisk/src/embed.h\"\nstatic inline\nreal embed_geometry (Point point, coord * p, coord * n)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  *n = facet_normal (point, cs, fs);\n  real alpha = line_alpha (val(cs,0,0,0), *n);\n  real area = line_length_center(*n,alpha,p);\n  normalize (n);\n  return area;\n}","\n    \n// #line 513 \"/home/Min/basilisk/src/embed.h\"\nif (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {\n\n\n\n\n\n\n      coord n, b;\n      real area = embed_geometry (point, &b, &n);\n      area *= pow (Delta, 2 - 1);\n      real Fn = area*embed_interpolate (point, p, b);\n      \n Fps.x += Fn*n.x;\n \n// #line 525\nFps.y += Fn*n.y;\n// # 534 \"/home/Min/basilisk/src/embed.h\"\n      if (constant(mu.x) != 0.) {\n real mua = 0., fa = 0.;\n  {\n   mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);\n   fa += _const_fm.x + _const_fm.x;\n } \n// #line 536\n{\n   mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);\n   fa += _const_fm.y + _const_fm.y;\n }\n mua /= fa;\n// # 599 \"/home/Min/basilisk/src/embed.h\"\n if (!(2 == 2)) qassert (\"/home/Min/basilisk/src/embed.h\", 599, \"dimension == 2\");\n coord dudn = embed_gradient (point, u, b, n);\n \n   Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);\n   \n// #line 602\nFmus.y -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.y*n.x);\n      }\n    }"})
    {_stencil_val(cs,0,0,0); _stencil_val(cs,0,0,0); {    






      
       _stencil_embed_geometry (point,NULL ,NULL );   
            
      _stencil_embed_interpolate (point, p,NULL );
       
  
#line 534 "/home/Min/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {      
 
  { 
_stencil_val(mu.x,0,0,0); _stencil_val(mu.x,1,0,0);
;; 
    
 
#line 539
} 
#line 536
{ 
_stencil_val(mu.y,0,0,0); _stencil_val(mu.y,0,1,0);
;; 
    
 
#line 539
}  
      
#line 599 "/home/Min/basilisk/src/embed.h"
     
  _stencil_embed_gradient (point, u,NULL ,NULL );
  
        
      }
    }      }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:Fmus)reduction(+:Fps)){
#line 512
foreach ()
    if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {






      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, 2 - 1);
      double Fn = area*embed_interpolate (point, p, b);
      
 Fps.x += Fn*n.x;
 
#line 525
Fps.y += Fn*n.y;
#line 534 "/home/Min/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {
 double mua = 0., fa = 0.;
  {
   mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
   fa += _const_fm.x + _const_fm.x;
 } 
#line 536
{
   mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
   fa += _const_fm.y + _const_fm.y;
 }
 mua /= fa;
#line 599 "/home/Min/basilisk/src/embed.h"
 if (!(2 == 2)) qassert ("/home/Min/basilisk/src/embed.h", 599, "dimension == 2");
 coord dudn = embed_gradient (point, u, b, n);
 
   Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
   
#line 602
Fmus.y -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.y*n.x);
      }
    }end_foreach();mpi_all_reduce_array(&Fmus.x,double,MPI_SUM,2);mpi_all_reduce_array(&Fps.x,double,MPI_SUM,2);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 604
}

  *Fp = Fps; *Fmu = Fmus;
end_tracing("embed_force","/home/Min/basilisk/src/embed.h",607);}
#line 616 "/home/Min/basilisk/src/embed.h"
double embed_vorticity (Point point, vector u, coord p, coord n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;




  coord dudn = embed_gradient (point, u, p, n);
#line 633 "/home/Min/basilisk/src/embed.h"
  return dudn.y*n.x - dudn.x*n.y;
}
#line 655 "/home/Min/basilisk/src/embed.h"
double embed_flux (Point point, scalar s, vector mu, double * val)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 656
if(!is_constant(fm.x)){{





  *val = 0.;
  if (val(cs,0,0,0) >= 1. || val(cs,0,0,0) <= 0.)
    return 0.;





  bool dirichlet = false;
  double grad = _attribute[s.i].boundary[embed] (point, point, s, &dirichlet);
  if (!grad && !dirichlet)
    return 0.;





  coord n = facet_normal (point, cs, fs), p;
  double alpha = line_alpha (val(cs,0,0,0), n);
  double area = line_length_center(n,alpha,&p);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, p);





  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    grad = dirichlet_gradient (point, s, cs, n, p, grad, &coef);
  }




  double mua = 0., fa = 0.;
   {
    mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
    fa += val(fm.x,0,0,0) + val(fm.x,1,0,0);
  } 
#line 699
{
    mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    fa += val(fm.y,0,0,0) + val(fm.y,0,1,0);
  }
  *val = - mua/(fa + 1e-30)*grad*area/Delta;
  return - mua/(fa + 1e-30)*coef*area/Delta;
}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);

#line 656
{





  *val = 0.;
  if (val(cs,0,0,0) >= 1. || val(cs,0,0,0) <= 0.)
    return 0.;





  bool dirichlet = false;
  double grad = _attribute[s.i].boundary[embed] (point, point, s, &dirichlet);
  if (!grad && !dirichlet)
    return 0.;





  coord n = facet_normal (point, cs, fs), p;
  double alpha = line_alpha (val(cs,0,0,0), n);
  double area = line_length_center(n,alpha,&p);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, p);





  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    grad = dirichlet_gradient (point, s, cs, n, p, grad, &coef);
  }




  double mua = 0., fa = 0.;
   {
    mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
    fa += _const_fm.x + _const_fm.x;
  } 
#line 699
{
    mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    fa += _const_fm.y + _const_fm.y;
  }
  *val = - mua/(fa + 1e-30)*grad*area/Delta;
  return - mua/(fa + 1e-30)*coef*area/Delta;
}}

#line 705
}
#line 714 "/home/Min/basilisk/src/embed.h"
#undef _neumann
#define _neumann(expr, ...) (data ? embed_area_center (point, &x, &y, &z),\
       *((bool *)data) = false, (expr) :\
       Delta*(expr) + val(_s,0,0,0))\

#line 718

#undef _neumann_homogeneous
#define _neumann_homogeneous(...) (data ? *((bool *)data) = false, (0) :\
    val(_s,0,0,0))\

#line 722

#undef _dirichlet
#define _dirichlet(expr, ...) (data ? embed_area_center (point, &x, &y, &z),\
       *((bool *)data) = true, (expr) :\
       2.*(expr) - val(_s,0,0,0))\

#line 727

#undef _dirichlet_homogeneous
#define _dirichlet_homogeneous(...) (data ? *((bool *)data) = true, (0) :\
      - val(_s,0,0,0))\

#line 731

#line 742 "/home/Min/basilisk/src/embed.h"
static inline double bilinear_embed (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!coarse(cs,0,0,0) || !coarse(cs,child.x,0,0))
    return coarse(s,0,0,0);

  if (!coarse(cs,0,child.y,0) || !coarse(cs,child.x,child.y,0))
    return coarse(s,0,0,0);







  return bilinear (point, s);
}
#line 789 "/home/Min/basilisk/src/embed.h"
     
void update_tracer (scalar f, vector uf, vector flux, double dt)
{tracing("update_tracer","/home/Min/basilisk/src/embed.h",790);
#line 801 "/home/Min/basilisk/src/embed.h"
  scalar  e=new_scalar("e");
  if(!is_constant(cm)){
  
#line 802
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"e","scalar",(void *)&e,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 802 \"/home/Min/basilisk/src/embed.h\"\n{\n\n\n\n\n    if (val(cs,0,0,0) <= 0.)\n      val_out_(e,0,0,0) = 0.;\n\n\n\n\n\n\n    else if (val(cs,0,0,0) >= 1.) {\n      \n val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/Delta;\n \n// #line 817\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/Delta;\n      val_out_(e,0,0,0) = 0.;\n    }\n// # 830 \"/home/Min/basilisk/src/embed.h\"\n    else {\n      real umax = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (fabs(val(uf.x,i,0,0)) > umax)\n     umax = fabs(val(uf.x,i,0,0));\n   \n// #line 834\nif (fabs(val(uf.y,0,i,0)) > umax)\n     umax = fabs(val(uf.y,0,i,0));}\n      real dtmax = Delta*val(cm,0,0,0)/(umax + 1e-30);\n\n\n\n\n      real F = 0.;\n      \n F += val(flux.x,0,0,0) - val(flux.x,1,0,0);\n \n// #line 843\nF += val(flux.y,0,0,0) - val(flux.y,0,1,0);\n      F /= Delta*val(cm,0,0,0);\n\n\n\n\n\n\n      if (dt <= dtmax) {\n val_out_(f,0,0,0) += dt*F;\n val_out_(e,0,0,0) = 0.;\n      }\n\n\n\n\n\n\n\n      else {\n val_out_(f,0,0,0) += dtmax*F;\n real scs = 0.;\n {foreach_neighbor(1)\n   scs += sq(val(cm,0,0,0));end_foreach_neighbor()}\n val_out_(e,0,0,0) = (dt - dtmax)*F*val(cm,0,0,0)/scs;\n      }\n    }\n  }"}) {




_stencil_val(cs,0,0,0);{
      {_stencil_val_a(e,0,0,0);  } 






{_stencil_val(cs,0,0,0);{ {
      
 {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val_r(f,0,0,0);   }
 
#line 817
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val_r(f,0,0,0);   }
      _stencil_val_a(e,0,0,0);  
    } 
#line 830
{   
      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(uf.x,i,0,0);
     {_stencil_val(uf.x,i,0,0);  }   }
   
#line 834
{_stencil_val(uf.y,0,i,0);
     {_stencil_val(uf.y,0,i,0);  }   }}     
      _stencil_val(cm,0,0,0);   




      
      
 { _stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);  }
 
#line 843
{ _stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);  }
_stencil_val(cm,0,0,0);






{ {
 _stencil_val_r(f,0,0,0);  
 _stencil_val_a(e,0,0,0);  
      } 







{
 _stencil_val_r(f,0,0,0);     
 
 {foreach_neighbor(1)
   {_stencil_val(cm,0,0,0);  }end_foreach_neighbor()}
_stencil_val(cm,0,0,0);
 
#line 867
_stencil_val_a(e,0,0,0);    
      }}
        






         







      
    
#line 869
}}   
#line 830 "/home/Min/basilisk/src/embed.h"
    
#line 869
}}




       






    
  
#line 870
}end_foreach_stencil();{
#line 802
foreach() {




    if (val(cs,0,0,0) <= 0.)
      val(e,0,0,0) = 0.;






    else if (val(cs,0,0,0) >= 1.) {
      
 val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/Delta;
 
#line 817
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/Delta;
      val(e,0,0,0) = 0.;
    }
#line 830 "/home/Min/basilisk/src/embed.h"
    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (fabs(val(uf.x,i,0,0)) > umax)
     umax = fabs(val(uf.x,i,0,0));
   
#line 834
if (fabs(val(uf.y,0,i,0)) > umax)
     umax = fabs(val(uf.y,0,i,0));}
      double dtmax = Delta*val(cm,0,0,0)/(umax + 1e-30);




      double F = 0.;
      
 F += val(flux.x,0,0,0) - val(flux.x,1,0,0);
 
#line 843
F += val(flux.y,0,0,0) - val(flux.y,0,1,0);
      F /= Delta*val(cm,0,0,0);






      if (dt <= dtmax) {
 val(f,0,0,0) += dt*F;
 val(e,0,0,0) = 0.;
      }







      else {
 val(f,0,0,0) += dtmax*F;
 double scs = 0.;
 {foreach_neighbor(1)
   scs += sq(val(cm,0,0,0));end_foreach_neighbor()}
 val(e,0,0,0) = (dt - dtmax)*F*val(cm,0,0,0)/scs;
      }
    }
  }end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 802
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"e","scalar",(void *)&e,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 802 \"/home/Min/basilisk/src/embed.h\"\n{\n\n\n\n\n    if (val(cs,0,0,0) <= 0.)\n      val_out_(e,0,0,0) = 0.;\n\n\n\n\n\n\n    else if (val(cs,0,0,0) >= 1.) {\n      \n val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/Delta;\n \n// #line 817\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/Delta;\n      val_out_(e,0,0,0) = 0.;\n    }\n// # 830 \"/home/Min/basilisk/src/embed.h\"\n    else {\n      real umax = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (fabs(val(uf.x,i,0,0)) > umax)\n     umax = fabs(val(uf.x,i,0,0));\n   \n// #line 834\nif (fabs(val(uf.y,0,i,0)) > umax)\n     umax = fabs(val(uf.y,0,i,0));}\n      real dtmax = Delta*_const_cm/(umax + 1e-30);\n\n\n\n\n      real F = 0.;\n      \n F += val(flux.x,0,0,0) - val(flux.x,1,0,0);\n \n// #line 843\nF += val(flux.y,0,0,0) - val(flux.y,0,1,0);\n      F /= Delta*_const_cm;\n\n\n\n\n\n\n      if (dt <= dtmax) {\n val_out_(f,0,0,0) += dt*F;\n val_out_(e,0,0,0) = 0.;\n      }\n\n\n\n\n\n\n\n      else {\n val_out_(f,0,0,0) += dtmax*F;\n real scs = 0.;\n {foreach_neighbor(1)\n   scs += sq(_const_cm);end_foreach_neighbor()}\n val_out_(e,0,0,0) = (dt - dtmax)*F*_const_cm/scs;\n      }\n    }\n  }"}) {




_stencil_val(cs,0,0,0);{
      {_stencil_val_a(e,0,0,0);  } 






{_stencil_val(cs,0,0,0);{ {
      
 {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val_r(f,0,0,0);   }
 
#line 817
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val_r(f,0,0,0);   }
      _stencil_val_a(e,0,0,0);  
    } 
#line 830
{   
      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(uf.x,i,0,0);
     {_stencil_val(uf.x,i,0,0);  }   }
   
#line 834
{_stencil_val(uf.y,0,i,0);
     {_stencil_val(uf.y,0,i,0);  }   }}
;   




      
      
 { _stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);  }
 
#line 843
{ _stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);  }
;






{ {
 _stencil_val_r(f,0,0,0);  
 _stencil_val_a(e,0,0,0);  
      } 







{
 _stencil_val_r(f,0,0,0);     
 
 {foreach_neighbor(1)
   {;  }end_foreach_neighbor()}
;
 
#line 867
_stencil_val_a(e,0,0,0);    
      }}
        






         







      
    
#line 869
}}   
#line 830 "/home/Min/basilisk/src/embed.h"
    
#line 869
}}




       






    
  
#line 870
}end_foreach_stencil();
  {
#line 802
foreach() {




    if (val(cs,0,0,0) <= 0.)
      val(e,0,0,0) = 0.;






    else if (val(cs,0,0,0) >= 1.) {
      
 val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/Delta;
 
#line 817
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/Delta;
      val(e,0,0,0) = 0.;
    }
#line 830 "/home/Min/basilisk/src/embed.h"
    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (fabs(val(uf.x,i,0,0)) > umax)
     umax = fabs(val(uf.x,i,0,0));
   
#line 834
if (fabs(val(uf.y,0,i,0)) > umax)
     umax = fabs(val(uf.y,0,i,0));}
      double dtmax = Delta*_const_cm/(umax + 1e-30);




      double F = 0.;
      
 F += val(flux.x,0,0,0) - val(flux.x,1,0,0);
 
#line 843
F += val(flux.y,0,0,0) - val(flux.y,0,1,0);
      F /= Delta*_const_cm;






      if (dt <= dtmax) {
 val(f,0,0,0) += dt*F;
 val(e,0,0,0) = 0.;
      }







      else {
 val(f,0,0,0) += dtmax*F;
 double scs = 0.;
 {foreach_neighbor(1)
   scs += sq(_const_cm);end_foreach_neighbor()}
 val(e,0,0,0) = (dt - dtmax)*F*_const_cm/scs;
      }
    }
  }end_foreach();}}





  foreach_stencil(1,{(NonLocal[]){{"cs","scalar",(void *)&cs,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"e","scalar",(void *)&e,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 876 \"/home/Min/basilisk/src/embed.h\"\n{\n    real se = 0.;\n    {foreach_neighbor(1)\n      se += val(e,0,0,0);end_foreach_neighbor()}\n    val_out_(f,0,0,0) += val(cs,0,0,0)*se;\n  }"}) {   
    
    {foreach_neighbor(1)
      { _stencil_val(e,0,0,0); }end_foreach_neighbor()} 
_stencil_val(cs,0,0,0);
    
#line 880
_stencil_val_r(f,0,0,0); 
  }end_foreach_stencil();





  {
#line 876
foreach() {
    double se = 0.;
    {foreach_neighbor(1)
      se += val(e,0,0,0);end_foreach_neighbor()}
    val(f,0,0,0) += val(cs,0,0,0)*se;
  }end_foreach();}delete((scalar*)((scalar[]){e,{-1}}));
end_tracing("update_tracer","/home/Min/basilisk/src/embed.h",882);}
#line 895
static int metric_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 895 "/home/Min/basilisk/src/embed.h"
      static int metric(const int i,const double t,Event *_ev){tracing("metric","/home/Min/basilisk/src/embed.h",895);
{
  if (is_constant (fm.x)) {
    
      if (!(constant (fm.x) == 1.)) qassert ("/home/Min/basilisk/src/embed.h", 899, "constant (fm.x) == 1.");
      
#line 899
if (!(constant (fm.y) == 1.)) qassert ("/home/Min/basilisk/src/embed.h", 899, "constant (fm.x) == 1.");
    fm = fs;
  }
  foreach_face_stencil(1,{(NonLocal[]){{"fs","vector",(void *)&fs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 902 \"/home/Min/basilisk/src/embed.h\"\n{is_face_x(){\n    val_out_(fs.x,0,0,0) = 1.;}end_is_face_x()\n// #line 902\nis_face_y(){\n    val_out_(fs.y,0,0,0) = 1.;}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val_a(fs.x,0,0,0);  }}end__stencil_is_face_x()
#line 902
_stencil_is_face_y(){
    {_stencil_val_a(fs.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 902
foreach_face_generic(){is_face_x(){
    val(fs.x,0,0,0) = 1.;}end_is_face_x()
#line 902
is_face_y(){
    val(fs.y,0,0,0) = 1.;}end_is_face_y()}end_foreach_face_generic();}
  if (is_constant (cm)) {
    if (!(constant (cm) == 1.)) qassert ("/home/Min/basilisk/src/embed.h", 905, "constant (cm) == 1.");
    cm = cs;
  }
  foreach_stencil(1,{(NonLocal[]){{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 909 \"/home/Min/basilisk/src/embed.h\"\nval_out_(cs,0,0,0) = 1.;"})
    {_stencil_val_a(cs,0,0,0);  }end_foreach_stencil();
  {
#line 908
foreach()
    val(cs,0,0,0) = 1.;end_foreach();}
#line 932 "/home/Min/basilisk/src/embed.h"
  restriction (((scalar[]){cs, fs.x, fs.y,{-1}}));
}{end_tracing("metric","/home/Min/basilisk/src/embed.h",933);return 0;}end_tracing("metric","/home/Min/basilisk/src/embed.h",933);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}





#line 938
      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/Min/basilisk/src/embed.h",938); {
  display ("draw_vof (c = 'cs', s = 'fs', filled = -1, "
    "fc = {0.5,0.5,0.5}, order = 2);"
#line 1425 "/home/Min/basilisk/src/common.h"
, false
#line 940 "/home/Min/basilisk/src/embed.h"
);
}{end_tracing("defaults","/home/Min/basilisk/src/embed.h",941);return 0;}end_tracing("defaults","/home/Min/basilisk/src/embed.h",941);}
#line 3 "reactor.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/Min/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/Min/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/Min/basilisk/src/run.h"
#line 9 "/home/Min/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/Min/basilisk/src/utils.h"







double DT = 1e30, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf = {0};





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if _MPI
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Multigrid"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 137 \"/home/Min/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (sq(Delta)*val(cm,0,0,0)) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (sq(Delta)*val(cm,0,0,0));\n      avg += (sq(Delta)*val(cm,0,0,0))*v;\n      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);\n    }"}
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (sq(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val(cm,0,0,0));
      avg += (sq(Delta)*val(cm,0,0,0))*v;
      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 137 \"/home/Min/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (sq(Delta)*_const_cm) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (sq(Delta)*_const_cm);\n      avg += (sq(Delta)*_const_cm)*v;\n      rms += (sq(Delta)*_const_cm)*sq(v);\n    }"}
)
    {_stencil_val(f,0,0,0);; {   
      _stencil_val(f,0,0,0);

;
;
; 
       
    
#line 143
}       }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (sq(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*_const_cm);
      avg += (sq(Delta)*_const_cm)*v;
      rms += (sq(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 165 \"/home/Min/basilisk/src/utils.h\"\nif ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (sq(Delta)*val(cm,0,0,0));\n      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);\n      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }"}
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {
      volume += (sq(Delta)*val(cm,0,0,0));
      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 165 \"/home/Min/basilisk/src/utils.h\"\nif ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (sq(Delta)*_const_cm);\n      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);\n      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }"}
)
    {; _stencil_val(f,0,0,0); {
;
;_stencil_val(f,0,0,0);
;_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {
      volume += (sq(Delta)*_const_cm);
      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/Min/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/Min/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/Min/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/Min/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
  foreach_stencil(1,{(NonLocal[]){{"fs","vector",(void *)&fs,NULL,0},{"g","vector",(void *)g,NULL,1},{"f","scalar",(void *)f,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 240 \"/home/Min/basilisk/src/utils.h\"\n{\n    scalar s; vector v;\n    {forin2 (s,v , f,g) {\n      if (s.gradient)\n { {\n\n   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))\n     val_out_(v.x,0,0,0) = 0.;\n   else\n\n     val_out_(v.x,0,0,0) = s.gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;\n } \n// #line 244\n{\n\n   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))\n     val_out_(v.y,0,0,0) = 0.;\n   else\n\n     val_out_(v.y,0,0,0) = s.gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;\n }}\n      else\n { {\n\n   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))\n     val_out_(v.x,0,0,0) = 0.;\n   else\n\n     val_out_(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);\n } \n// #line 253\n{\n\n   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))\n     val_out_(v.y,0,0,0) = 0.;\n   else\n\n     val_out_(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);\n }}\n    } endforin2()}\n  }"}) {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
     {_stencil_val_a(v.x,0,0,0);  }

     
{_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);_stencil_val_a(v.x,0,0,0);   }}

      
   
 
#line 251
} 
#line 244
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
     {_stencil_val_a(v.y,0,0,0);  }

     
{_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);_stencil_val_a(v.y,0,0,0);   }}

      
   
 
#line 251
}}
      else
 { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
     {_stencil_val_a(v.x,0,0,0);  }

     
{_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);_stencil_val_a(v.x,0,0,0);   }}

      
   
 
#line 260
} 
#line 253
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
     {_stencil_val_a(v.y,0,0,0);  }

     
{_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);_stencil_val_a(v.y,0,0,0);   }}

      
   
 
#line 260
}}
    }}}
  }end_foreach_stencil();
  {
#line 240
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {

   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
     val(v.x,0,0,0) = 0.;
   else

     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{

   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
     val(v.y,0,0,0) = 0.;
   else

     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }}
      else
 { {

   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
     val(v.x,0,0,0) = 0.;
   else

     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{

   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
     val(v.y,0,0,0) = 0.;
   else

     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
#line 280 "/home/Min/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/Min/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 1e-30)*Delta);"})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();{
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 1e-30)*Delta);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/Min/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 1e-30)*Delta);"})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 1e-30)*Delta);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/Min/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 1e-30)*Delta);"})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 1e-30)*Delta);end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/Min/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 1e-30)*Delta);"})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 1e-30)*Delta);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  if(!is_constant(cm)){
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 296 \"/home/Min/basilisk/src/utils.h\"\n{\n    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }"}) {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 296 \"/home/Min/basilisk/src/utils.h\"\n{\n    if ((sq(Delta)*_const_cm) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }"}) {
; {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
#line 340 "/home/Min/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/Min/basilisk/src/utils.h", 343, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }\
\
      if (_n == 2) {\

#line 380

#line 410 "/home/Min/basilisk/src/utils.h"
#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

#line 1 "./output.h"
#line 1 "/home/Min/basilisk/src/output.h"
#line 37 "/home/Min/basilisk/src/output.h"
     
void output_field (scalar * list,
     FILE * fp,
     int n,
     bool linear,
     coord box[2])
{tracing("output_field","/home/Min/basilisk/src/output.h",38);
  n++;
  int len = list_len (list);
  double Delta = 0.999999*(box[1].x - box[0].x)/(n - 1);
  int ny = (box[1].y - box[0].y)/Delta + 1;
  double ** field = (double **) matrix_new (n, ny, len*sizeof(double)), * v = field[0];
  for (int i = 0; i < n*ny*len; i++, v++)
    *v = 1e30;
  coord box1[2] = {{box[0].x - Delta/2., box[0].y - Delta/2.},
     {box[0].x + (n - 0.5)*Delta, box[0].y + (ny - 0.5)*Delta}};
  coord cn = {n, ny}, p;




  foreach_region_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"len","int",(void *)&len,NULL,0},{"list","scalar",(void *)list,NULL,1},{"field","double",(void *)field,NULL,2},{"cn","coord",(void *)&cn,NULL,0},{"box1","coord",(void *)box1,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\n}","\n\n  \n// #line 60 \"/home/Min/basilisk/src/output.h\"\n{\n    real ** alias = field;\n    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;\n    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;\n    int k = 0;\n    {forin (scalar, s , list)\n      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0); endforin()}\n  }"})

  {                     
    
    
    
    
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_interpolate_linear (point, s, NULL, NULL, NULL); _stencil_val(s,0,0,0);      }}}
  }end_foreach_region_stencil();




  {
#line 58
foreach_region (p, box1, cn)

  {
    double ** alias = field;
    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;
    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;
    int k = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0);}}
  }end_foreach_region();}

  if (pid() == 0) {
    fprintf (fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + box[0].x;
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + box[0].y;

 fprintf (fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', fp);
      }
      fputc ('\n', fp);
    }
    fflush (fp);
  }

  matrix_free (field);
end_tracing("output_field","/home/Min/basilisk/src/output.h",92);}
#line 120 "/home/Min/basilisk/src/output.h"
     
void output_matrix (scalar f, FILE * fp, int n, bool linear)
{tracing("output_matrix","/home/Min/basilisk/src/output.h",121);
  float fn = n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      v = interpolate (f, xp, yp
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0.
#line 135 "/home/Min/basilisk/src/output.h"
, linear);
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
end_tracing("output_matrix","/home/Min/basilisk/src/output.h",140);}
#line 149 "/home/Min/basilisk/src/output.h"
typedef void (* Colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} Color;

Color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  Color c;
  if (val == 1e30) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i >= 0 && i < 127 - 1)) qassert ("/home/Min/basilisk/src/output.h", 281, "i >= 0 && i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 300 "/home/Min/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);



    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/Min/basilisk/src/output.h", 382, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/Min/basilisk/src/output.h", 456, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 534 "/home/Min/basilisk/src/output.h"
     
void output_ppm (scalar f,
   FILE * fp,
   int n,
   char * file,
   double min, double max, double spread,
   double z,
   bool linear,
   coord box[2],
   scalar mask,
   Colormap map,
   char * opt)
{tracing("output_ppm","/home/Min/basilisk/src/output.h",535);

  if (!min && !max) {
    stats s = statsf (f);
    if (spread < 0.)
      min = s.min, max = s.max;
    else {
      double avg = s.sum/s.volume;
      min = avg - spread*s.stddev; max = avg + spread*s.stddev;
    }
  }
  box[0].z = z, box[1].z = z;

  coord cn = {n}, p;
  double delta = (box[1].x - box[0].x)/n;
  cn.y = (int)((box[1].y - box[0].y)/delta);
  if (((int)cn.y) % 2) cn.y++;

  Color ** ppm = (Color **) matrix_new (cn.y, cn.x, sizeof(Color));
  unsigned char * ppm0 = &ppm[0][0].r;
  int len = 3*cn.x*cn.y;
  memset (ppm0, 0, len*sizeof (unsigned char));
  double cmap[127][3];
  (* map) (cmap);




  foreach_region_stencil (1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0},{"min","double",(void *)&min,NULL,0},{"cmap","double",(void *)cmap,(int[]){1273,0},0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"ppm","Color",(void *)ppm,NULL,2},{"f","scalar",(void *)&f,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"mask","scalar",(void *)&mask,NULL,0},{"cn","coord",(void *)&cn,NULL,0},{"box","coord",(void *)box,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 261 \"/home/Min/basilisk/src/output.h\"\nColor colormap_color (real cmap[127][3],\n        real val, real min, real max)\n{\n  Color c;\n  if (val == 1e30) {\n    c.r = c.g = c.b = 0;\n    return c;\n  }\n  int i;\n  real coef;\n  if (max != min)\n    val = (val - min)/(max - min);\n  else\n    val = 0.;\n  if (val <= 0.) i = 0, coef = 0.;\n  else if (val >= 1.) i = 127 - 2, coef = 1.;\n  else {\n    i = val*(127 - 1);\n    coef = val*(127 - 1) - i;\n  }\n  if (!(i >= 0 && i < 127 - 1)) qassert (\"/home/Min/basilisk/src/output.h\", 281, \"i >= 0 && i < NCMAP - 1\");\n  unsigned char * c1 = (unsigned char *) &c;\n  for (int j = 0; j < 3; j++)\n    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);\n  return c;\n}\n\n\n// #line 805 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/Min/basilisk/src/grid/cartesian-common.h\"\n}","\n\n  \n// #line 576 \"/home/Min/basilisk/src/output.h\"\n{\n    real v;\n    if (mask.i >= 0) {\n      if (linear) {\n real m = interpolate_linear (point, mask, p.x, p.y, p.z);\n if (m < 0.)\n   v = 1e30;\n else\n   v = interpolate_linear (point, f, p.x, p.y, p.z);\n      }\n      else {\n if (val(mask,0,0,0) < 0.)\n   v = 1e30;\n else\n   v = val(f,0,0,0);\n      }\n    }\n    else if (linear)\n      v = interpolate_linear (point, f, p.x, p.y, p.z);\n    else\n      v = val(f,0,0,0);\n    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;\n    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;\n    Color ** alias = ppm;\n    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);\n  }"})

  { 
    
    if (mask.i >= 0) {
      if (linear) {  
  _stencil_interpolate_linear (point, mask, NULL, NULL, NULL);
{ 
    
   
{ _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }}
    
 
      
#line 585
}
      else {
_stencil_val(mask,0,0,0);{
     
   
{ _stencil_val(f,0,0,0); }}
    
 
      
#line 591
}
    }
    else if (linear)
      { _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }
    else
      { _stencil_val(f,0,0,0); }                  
    
    
         
         
  }end_foreach_region_stencil();




  {
#line 574
foreach_region (p, box, cn)

  {
    double v;
    if (mask.i >= 0) {
      if (linear) {
 double m = interpolate_linear (point, mask, p.x, p.y, p.z);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate_linear (point, f, p.x, p.y, p.z);
      }
      else {
 if (val(mask,0,0,0) < 0.)
   v = 1e30;
 else
   v = val(f,0,0,0);
      }
    }
    else if (linear)
      v = interpolate_linear (point, f, p.x, p.y, p.z);
    else
      v = val(f,0,0,0);
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;
    Color ** alias = ppm;
    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);
  }end_foreach_region();}

  if (pid() == 0) {
    if (file)
      fp = open_image (file, opt);

    fprintf (fp, "P6\n%g %g 255\n", cn.x, cn.y);
    fwrite (ppm0, sizeof(unsigned char), 3*cn.x*cn.y, fp);

    if (file)
      close_image (file, fp);
    else
      fflush (fp);
  }

  matrix_free (ppm);
end_tracing("output_ppm","/home/Min/basilisk/src/output.h",617);}
#line 649 "/home/Min/basilisk/src/output.h"
     
void output_grd (scalar f,
   FILE * fp,
   double Delta,
   bool linear,
   double box[2][2],
   scalar mask)
{tracing("output_grd","/home/Min/basilisk/src/output.h",650);
  int nx = (box[1][0] - box[0][0])/Delta;
  int ny = (box[1][1] - box[0][1])/Delta;


  fprintf (fp, "ncols          %d\n", nx);
  fprintf (fp, "nrows          %d\n", ny);
  fprintf (fp, "xllcorner      %g\n", box[0][0]);
  fprintf (fp, "yllcorner      %g\n", box[0][1]);
  fprintf (fp, "cellsize       %g\n", Delta);
  fprintf (fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + box[0][0] + Delta/2., v;
      if (mask.i >= 0) {
 double m = interpolate (mask, xp, yp
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0.
#line 674 "/home/Min/basilisk/src/output.h"
, linear);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate (f, xp, yp
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0.
#line 678 "/home/Min/basilisk/src/output.h"
, linear);
      }
      else
 v = interpolate (f, xp, yp
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0.
#line 681 "/home/Min/basilisk/src/output.h"
, linear);
      if (v == 1e30)
 fprintf (fp, "-9999 ");
      else
 fprintf (fp, "%f ", v);
    }
    fprintf (fp, "\n");
  }

  fflush (fp);
end_tracing("output_grd","/home/Min/basilisk/src/output.h",691);}
#line 718 "/home/Min/basilisk/src/output.h"
static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (FILE * fp,
   scalar * list,
   char * file,
   bool translate)
{tracing("output_gfs","/home/Min/basilisk/src/output.h",739);
  char * fname = file;

#if _MPI



  FILE * sfp = fp;
  if (file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    fp = NULL;
  }
#endif

  bool opened = false;
  if (fp == NULL) {
    if (fname == NULL)
      fp = fout;
    else if (!(fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * slist = list ? list : list_copy (all);

  restriction (slist);
  fprintf (fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (slist != NULL && slist[0].i != -1) {
    scalar s = slist[0];
    char * name = replace (_attribute[s.i].name, '.', '_', translate);
    fprintf (fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(slist); i++) {
      scalar s = slist[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', translate);
 fprintf (fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (fp, " ");
  }
  fprintf (fp, "} {\n");
  fprintf (fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (fp, "  VariableTracerVOF f\n");
  fprintf (fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if _MPI
  long header;
  if ((header = ftell (fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



  {foreach_cell() {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 851 "/home/Min/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, fp);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

#if _MPI
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", fp);
  fflush (fp);

  if (!list)
    pfree (slist,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (fp);

#if _MPI
  if (file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (sfp == NULL)
 sfp = fout;
      fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, fp)) > 0)
 fwrite (buffer, 1, l, sfp);
      fflush (sfp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
end_tracing("output_gfs","/home/Min/basilisk/src/output.h",919);}
#line 943 "/home/Min/basilisk/src/output.h"
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !_MPI
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/Min/basilisk/src/output.h",989);
  char * name = NULL;
  if (!fp) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/Min/basilisk/src/output.h", 1005, "fp");

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, slist);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (slist,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/Min/basilisk/src/output.h",1038);}
#else
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/Min/basilisk/src/output.h",1041);
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, slist);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

  {foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (slist,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/Min/basilisk/src/output.h",1112);}
#endif

     
bool restore (const char * file,
       scalar * list,
       FILE * fp)
{tracing("restore","/home/Min/basilisk/src/output.h",1116);
  if (!fp && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/Min/basilisk/src/output.h",1121);return false;}
  if (!(fp)) qassert ("/home/Min/basilisk/src/output.h", 1122, "fp");

  struct DumpHeader header = {0};
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }
#line 1153 "/home/Min/basilisk/src/output.h"
  init_grid (1 << header.depth);



  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (slist));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (slist,__func__,__FILE__,__LINE__);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }
#line 1232 "/home/Min/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});



  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }}}
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}


  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (slist,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/Min/basilisk/src/output.h",1282);return true;}
end_tracing("restore","/home/Min/basilisk/src/output.h",1283);}
#line 431 "/home/Min/basilisk/src/utils.h"
#line 12 "/home/Min/basilisk/src/run.h"

     
void run (void)
{tracing("run","/home/Min/basilisk/src/run.h",14);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/home/Min/basilisk/src/run.h",37);}




static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}





#line 42
      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/home/Min/basilisk/src/run.h",42); {
  display ("box();"
#line 1425 "/home/Min/basilisk/src/common.h"
, false
#line 43 "/home/Min/basilisk/src/run.h"
);
}{end_tracing("defaults_0","/home/Min/basilisk/src/run.h",44);return 0;}end_tracing("defaults_0","/home/Min/basilisk/src/run.h",44);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = TEND_EVENT)!=0;*ip=i;*tp=t;return ret;}






#line 50
      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/home/Min/basilisk/src/run.h",50); {
  display ("", true);
}{end_tracing("cleanup","/home/Min/basilisk/src/run.h",52);return 0;}end_tracing("cleanup","/home/Min/basilisk/src/run.h",52);}
#line 28 "/home/Min/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/Min/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  if(!is_constant(fm.x)){
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"fm","vector",(void *)&fm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}","\n// #line 6 \"/home/Min/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n      if (!(val(fm.x,0,0,0))) qassert (\"/home/Min/basilisk/src/timestep.h\", 10, \"fm.x[]\");\n      dt *= val(fm.x,0,0,0);\n\n\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n      if (!(val(fm.y,0,0,0))) qassert (\"/home/Min/basilisk/src/timestep.h\", 10, \"fm.x[]\");\n      dt *= val(fm.y,0,0,0);\n\n\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);

_stencil_val(fm.x,0,0,0); 
_stencil_val(fm.x,0,0,0);    

          
       



         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);

_stencil_val(fm.y,0,0,0); 
_stencil_val(fm.y,0,0,0);    

          
       



         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));

      if (!(val(fm.x,0,0,0))) qassert ("/home/Min/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= val(fm.x,0,0,0);



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));

      if (!(val(fm.y,0,0,0))) qassert ("/home/Min/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= val(fm.y,0,0,0);



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}","\n// #line 6 \"/home/Min/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n      if (!(_const_fm.x)) qassert (\"/home/Min/basilisk/src/timestep.h\", 10, \"fm.x[]\");\n      dt *= _const_fm.x;\n\n\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n      if (!(_const_fm.y)) qassert (\"/home/Min/basilisk/src/timestep.h\", 10, \"fm.x[]\");\n      dt *= _const_fm.y;\n\n\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);

;
;    

          
       



         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);

;
;    

          
       



         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));

      if (!(_const_fm.x)) qassert ("/home/Min/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= _const_fm.x;



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));

      if (!(_const_fm.y)) qassert ("/home/Min/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= _const_fm.y;



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/home/Min/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/Min/basilisk/src/bcg.h"
#line 11 "/home/Min/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}),((vector[]) {g,{{-1},{-1}}}));




  if(!is_constant(fm.x) && !is_constant(src)){




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/Min/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(src)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/Min/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/Min/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/Min/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/Min/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
#line 58 "/home/Min/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/Min/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}






void advection (scalar * tracers, vector u, double dt,
  scalar * src)
{




  scalar * psrc = src;
  if (!src)
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      const scalar zero = new_const_scalar("zero",9, 0.);
      src = list_append (src, zero);
    }}}
  if (!(list_len (tracers) == list_len (src))) qassert ("/home/Min/basilisk/src/bcg.h", 80, "list_len (tracers) == list_len (src)");

  scalar f, source;
  {scalar*_i0=src;scalar*_i1= tracers;if(_i0)for(source=*_i0,f=*_i1;_i0->i>= 0;source=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, u, flux, dt, source);





    update_tracer (f, u, flux, dt);delete((scalar*)((vector[]){flux,{{-1},{-1}}}));

  }}}

  if (!psrc)
    pfree (src,__func__,__FILE__,__LINE__);
}
#line 30 "/home/Min/basilisk/src/navier-stokes/centered.h"

#line 1 "./viscosity-embed.h"
#line 1 "/home/Min/basilisk/src/viscosity-embed.h"
#line 31 "/home/Min/basilisk/src/viscosity-embed.h"
#line 1 "./poisson.h"
#line 1 "/home/Min/basilisk/src/poisson.h"
#line 32 "/home/Min/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}





    else
      {foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear_embed(point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil(1,{(NonLocal[]){{"da","scalar",(void *)da,NULL,1},{"a","scalar",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 84 \"/home/Min/basilisk/src/poisson.h\"\n{\n    scalar s, ds;\n    {forin2 (s, ds , a, da)\n     \n val_out_(s,0,0,0) += val(ds,0,0,0); endforin2()}\n  }"}) {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 { _stencil_val(ds,0,0,0);_stencil_val_r(s,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
#line 102 "/home/Min/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/Min/basilisk/src/poisson.h"
mgstats mg_solve (scalar * a, scalar * b,
    double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data),
    void (* relax) (scalar * da, scalar * res, int depth,
      void * data),
    void * data,
    int nrelax,
    scalar * res,
    int minlevel,
    double tolerance)
{





  scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  scalar rhs = b[0];
  foreach_stencil (1,{(NonLocal[]){{"rhs","scalar",(void *)&rhs,NULL,0},{"sum","double",(void *)&sum,NULL,0,'+'},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 161 \"/home/Min/basilisk/src/poisson.h\"\nsum += val(rhs,0,0,0);"})
    { _stencil_val(rhs,0,0,0); }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 160
foreach ()
    sum += val(rhs,0,0,0);end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 162
s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;




  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, relax, data,
       s.nrelax,
       minlevel,
       grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);
#line 192 "/home/Min/basilisk/src/poisson.h"
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = minlevel;




  if (s.resa > tolerance) {
    scalar v = a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d tolerance: %g\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }




  if (!pres)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 251 "/home/Min/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;

  double (* embed_flux) (Point, scalar, vector, double *);

};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 289 "/home/Min/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = 0., val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = 0., val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = 0., val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = 0., val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
#line 331 "/home/Min/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;
#line 365 "/home/Min/basilisk/src/poisson.h"
  if(!is_constant(lambda) && !is_constant(alpha.x)){
#line 365 "/home/Min/basilisk/src/poisson.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"p","not implemented yet",(void *)p,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"lambda","scalar",(void *)&lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}"," \n// #line 365 \"/home/Min/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) += (val(alpha.x,0,0,0)*(a.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -\n  val(alpha.x,1,0,0)*(a.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;\n      \n// #line 368\nval_out_(res,0,0,0) += (val(alpha.y,0,0,0)*(a.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -\n  val(alpha.y,0,1,0)*(a.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;\n\n    if (p->embed_flux) {\n      real c, e = p->embed_flux (point, a, alpha, &c);\n      val_out_(res,0,0,0) += c - e*val(a,0,0,0);\n    }\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 366
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(alpha.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
  _stencil_val(alpha.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, a, 1);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}
      
#line 368
{_stencil_val(alpha.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
  _stencil_val(alpha.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, a, 1);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}

    if (p->embed_flux) {   
      default_stencil ( point,((scalar[]){alpha.x , alpha.y, a,{-1}}) );
_stencil_val(a,0,0,0);
      
#line 373
_stencil_val_r(res,0,0,0);    
    }

_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }

        
  
#line 378
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 365
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*(_attribute[a.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 368
val(res,0,0,0) += (val(alpha.y,0,0,0)*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*(_attribute[a.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 378
}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
#line 365 "/home/Min/basilisk/src/poisson.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"p","not implemented yet",(void *)p,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"_const_lambda","double",(void *)&_const_lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}"," \n// #line 365 \"/home/Min/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) += (val(alpha.x,0,0,0)*(a.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -\n  val(alpha.x,1,0,0)*(a.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;\n      \n// #line 368\nval_out_(res,0,0,0) += (val(alpha.y,0,0,0)*(a.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -\n  val(alpha.y,0,1,0)*(a.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;\n\n    if (p->embed_flux) {\n      real c, e = p->embed_flux (point, a, alpha, &c);\n      val_out_(res,0,0,0) += c - e*val(a,0,0,0);\n    }\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 366
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(alpha.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
  _stencil_val(alpha.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, a, 1);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}
      
#line 368
{_stencil_val(alpha.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
  _stencil_val(alpha.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, a, 1);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}

    if (p->embed_flux) {   
      default_stencil ( point,((scalar[]){alpha.x , alpha.y, a,{-1}}) );
_stencil_val(a,0,0,0);
      
#line 373
_stencil_val_r(res,0,0,0);    
    }

_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }

        
  
#line 378
}end_foreach_stencil();
#line 365 "/home/Min/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 365
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*(_attribute[a.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 368
val(res,0,0,0) += (val(alpha.y,0,0,0)*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*(_attribute[a.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 378
}else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 365 "/home/Min/basilisk/src/poisson.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"alpha","vector",(void *)&alpha,NULL,0},{"p","not implemented yet",(void *)p,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"lambda","scalar",(void *)&lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}"," \n// #line 365 \"/home/Min/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) += (_const_alpha.x*(a.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -\n  _const_alpha.x*(a.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;\n      \n// #line 368\nval_out_(res,0,0,0) += (_const_alpha.y*(a.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -\n  _const_alpha.y*(a.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;\n\n    if (p->embed_flux) {\n      real c, e = p->embed_flux (point, a, alpha, &c);\n      val_out_(res,0,0,0) += c - e*val(a,0,0,0);\n    }\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 366
_stencil_val_a(res,0,0,0);  
    
      {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, a, 1);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}
      
#line 368
{; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, a, 1);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}

    if (p->embed_flux) {   
      default_stencil ( point,((scalar[]){alpha.x , alpha.y, a,{-1}}) );
_stencil_val(a,0,0,0);
      
#line 373
_stencil_val_r(res,0,0,0);    
    }

_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }

        
  
#line 378
}end_foreach_stencil();
#line 365 "/home/Min/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 365
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*(_attribute[a.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 368
val(res,0,0,0) += (_const_alpha.y*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*(_attribute[a.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 378
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 365 "/home/Min/basilisk/src/poisson.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"alpha","vector",(void *)&alpha,NULL,0},{"p","not implemented yet",(void *)p,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"_const_lambda","double",(void *)&_const_lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}"," \n// #line 365 \"/home/Min/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) += (_const_alpha.x*(a.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -\n  _const_alpha.x*(a.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;\n      \n// #line 368\nval_out_(res,0,0,0) += (_const_alpha.y*(a.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -\n  _const_alpha.y*(a.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;\n\n    if (p->embed_flux) {\n      real c, e = p->embed_flux (point, a, alpha, &c);\n      val_out_(res,0,0,0) += c - e*val(a,0,0,0);\n    }\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 366
_stencil_val_a(res,0,0,0);  
    
      {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, a, 1);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}
      
#line 368
{; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, a, 1);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 368
_stencil_val_r(res,0,0,0);                       
}

    if (p->embed_flux) {   
      default_stencil ( point,((scalar[]){alpha.x , alpha.y, a,{-1}}) );
_stencil_val(a,0,0,0);
      
#line 373
_stencil_val_r(res,0,0,0);    
    }

_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }

        
  
#line 378
}end_foreach_stencil();
#line 365 "/home/Min/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 365
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*(_attribute[a.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, a, 1) : (val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 368
val(res,0,0,0) += (_const_alpha.y*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*(_attribute[a.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, a, 1) : (val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 378
}

  return maxres;
}
#line 392 "/home/Min/basilisk/src/poisson.h"
mgstats poisson (scalar a, scalar b,
           vector alpha,
           scalar lambda,
   double tolerance,
   int nrelax,
   int minlevel,
   scalar * res,
   double (* flux) (Point, scalar, vector, double *))
{






  if (alpha.x.i < 0)
    alpha = unityf;
  if (lambda.i < 0) {
    const scalar zeroc = new_const_scalar("zeroc",10, 0.);
    lambda = zeroc;
  }




  restriction (((scalar[]){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {a, b, alpha, lambda, tolerance, nrelax, minlevel, res };

  if (!flux && _attribute[a.i].boundary[embed] != symmetry)
    p.embed_flux = embed_flux;
  else
    p.embed_flux = flux;

  mgstats s = mg_solve ((
#line 125
scalar *
#line 434
)((scalar[]){a,{-1}}),( 
#line 125
scalar *
#line 434
)((scalar[]) {b,{-1}}), residual, relax, &p
,
   
#line 435
nrelax, res, max(1, minlevel)
#line 133
, 
TOLERANCE
#line 435
);




  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 463 "/home/Min/basilisk/src/poisson.h"
     
mgstats project (vector uf, scalar p,
           vector alpha,
   double dt,
   int nrelax)
{tracing("project","/home/Min/basilisk/src/poisson.h",464);






  scalar  div=new_scalar("div");
  foreach_stencil(1,{(NonLocal[]){{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"div","scalar",(void *)&div,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 476 \"/home/Min/basilisk/src/poisson.h\"\n{\n    val_out_(div,0,0,0) = 0.;\n    \n      val_out_(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);\n      \n// #line 479\nval_out_(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);\n    val_out_(div,0,0,0) /= dt*Delta;\n  }"}) {
    _stencil_val_a(div,0,0,0);  
    
      { _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 479
{ _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val_r(div,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 476
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 479
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
#line 492 "/home/Min/basilisk/src/poisson.h"
  mgstats mgp = poisson (p, div, alpha
#line 393
,
( scalar) {-1}
#line 493
, TOLERANCE/sq(dt), nrelax
#line 396
, 
0, 
NULL, 
NULL
#line 493
);




  if(!is_constant(alpha.x)){




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n// #line 498 \"/home/Min/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*(p.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*(p.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val(alpha.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);            }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {_stencil_val(alpha.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);            }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*(_attribute[p.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*(_attribute[p.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n// #line 498 \"/home/Min/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*_const_alpha.x*(p.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*_const_alpha.y*(p.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}"}){_stencil_is_face_x(){
    {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);            }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);            }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*(_attribute[p.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*(_attribute[p.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/Min/basilisk/src/poisson.h",501);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/Min/basilisk/src/poisson.h",502);}
#line 32 "/home/Min/basilisk/src/viscosity-embed.h"

struct Viscosity {
  vector mu;
  scalar rho;
  double dt;
  double (* embed_flux) (Point, scalar, vector, double *);
};
#line 58 "/home/Min/basilisk/src/viscosity-embed.h"
static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));

  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  if(!is_constant(mu.x) && !is_constant(rho)){{foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += val(mu.x,0,0,0) + val(mu.x,1,0,0);
      
#line 70
avgmu += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
 
#line 78
a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0) + ((coord){0.,0.}).x + dt*d) + avgmu);
    } 
#line 72
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
 
#line 78
a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0) + ((coord){0.,0.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else if(is_constant(mu.x) && !is_constant(rho)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
  {
#line 67
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += _const_mu.x + _const_mu.x;
      
#line 70
avgmu += _const_mu.y + _const_mu.y;
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
 
#line 78
a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0) + ((coord){0.,0.}).x + dt*d) + avgmu);
    } 
#line 72
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
 
#line 78
a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0) + ((coord){0.,0.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else if(!is_constant(mu.x) && is_constant(rho)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
  {
#line 67
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += val(mu.x,0,0,0) + val(mu.x,1,0,0);
      
#line 70
avgmu += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
 
#line 78
a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho + ((coord){0.,0.}).x + dt*d) + avgmu);
    } 
#line 72
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
 
#line 78
a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho + ((coord){0.,0.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
  {
#line 67
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += _const_mu.x + _const_mu.x;
      
#line 70
avgmu += _const_mu.y + _const_mu.y;
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
 
#line 78
a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho + ((coord){0.,0.}).x + dt*d) + avgmu);
    } 
#line 72
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
 
#line 78
a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho + ((coord){0.,0.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}
#line 94 "/home/Min/basilisk/src/viscosity-embed.h"
}

static double residual_diffusion (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  if(!is_constant(mu.x) && !is_constant(rho)){
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"dt","double",(void *)&dt,NULL,0},{"rho","scalar",(void *)&rho,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n    \n// #line 129 \"/home/Min/basilisk/src/viscosity-embed.h\"\n{ {\n      scalar s = u.x;\n      real a = 0.;\n      \n a += val(mu.x,0,0,0)*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n \n// #line 133\na += val(mu.y,0,0,0)*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.x, mu, &c);\n val_out_(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));\n      }\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    } \n// #line 129\n{\n      scalar s = u.y;\n      real a = 0.;\n      \n a += val(mu.y,0,0,0)*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n \n// #line 133\na += val(mu.x,0,0,0)*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.y, mu, &c);\n val_out_(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));\n      }\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }}"})
    { {
      scalar s = u.x;   
      
      
 { _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0); _stencil_val(mu.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      }
 
#line 133
{ _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0); _stencil_val(mu.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      } 
_stencil_val(r.x,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(u.x,0,0,0);
      
#line 134
_stencil_val_a(res.x,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.x,{-1}}) );
_stencil_val(u.x,0,0,0);
 
#line 137
_stencil_val_r(res.x,0,0,0);    
      }
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }
          
    
#line 141
} 
#line 129
{
      scalar s = u.y;   
      
      
 { _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0); _stencil_val(mu.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      }
 
#line 133
{ _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0); _stencil_val(mu.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      } 
_stencil_val(r.y,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(u.y,0,0,0);
      
#line 134
_stencil_val_a(res.y,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.y,{-1}}) );
_stencil_val(u.y,0,0,0);
 
#line 137
_stencil_val_r(res.y,0,0,0);    
      }
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }
          
    
#line 141
}}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 128
foreach ()
    { {
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
 
#line 133
a += val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
      val(res.x,0,0,0) = val(r.x,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 129
{
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
 
#line 133
a += val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
      val(res.y,0,0,0) = val(r.y,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 141
}else if(is_constant(mu.x) && !is_constant(rho)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"mu","vector",(void *)&mu,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"rho","scalar",(void *)&rho,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n    \n// #line 129 \"/home/Min/basilisk/src/viscosity-embed.h\"\n{ {\n      scalar s = u.x;\n      real a = 0.;\n      \n a += _const_mu.x*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n \n// #line 133\na += _const_mu.y*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.x, mu, &c);\n val_out_(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));\n      }\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    } \n// #line 129\n{\n      scalar s = u.y;\n      real a = 0.;\n      \n a += _const_mu.y*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n \n// #line 133\na += _const_mu.x*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.y, mu, &c);\n val_out_(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));\n      }\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }}"})
    { {
      scalar s = u.x;   
      
      
 {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      }
 
#line 133
{; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      } 
_stencil_val(r.x,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(u.x,0,0,0);
      
#line 134
_stencil_val_a(res.x,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.x,{-1}}) );
_stencil_val(u.x,0,0,0);
 
#line 137
_stencil_val_r(res.x,0,0,0);    
      }
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }
          
    
#line 141
} 
#line 129
{
      scalar s = u.y;   
      
      
 {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      }
 
#line 133
{; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      } 
_stencil_val(r.y,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(u.y,0,0,0);
      
#line 134
_stencil_val_a(res.y,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.y,{-1}}) );
_stencil_val(u.y,0,0,0);
 
#line 137
_stencil_val_r(res.y,0,0,0);    
      }
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }
          
    
#line 141
}}end_foreach_stencil();
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 128
foreach ()
    { {
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
 
#line 133
a += _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
      val(res.x,0,0,0) = val(r.x,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 129
{
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
 
#line 133
a += _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
      val(res.y,0,0,0) = val(r.y,0,0,0) - (val(rho,0,0,0) + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 141
}else if(!is_constant(mu.x) && is_constant(rho)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"dt","double",(void *)&dt,NULL,0},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n    \n// #line 129 \"/home/Min/basilisk/src/viscosity-embed.h\"\n{ {\n      scalar s = u.x;\n      real a = 0.;\n      \n a += val(mu.x,0,0,0)*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n \n// #line 133\na += val(mu.y,0,0,0)*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - (_const_rho + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.x, mu, &c);\n val_out_(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));\n      }\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    } \n// #line 129\n{\n      scalar s = u.y;\n      real a = 0.;\n      \n a += val(mu.y,0,0,0)*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n \n// #line 133\na += val(mu.x,0,0,0)*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - (_const_rho + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.y, mu, &c);\n val_out_(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));\n      }\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }}"})
    { {
      scalar s = u.x;   
      
      
 { _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0); _stencil_val(mu.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      }
 
#line 133
{ _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0); _stencil_val(mu.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      } 
_stencil_val(r.x,0,0,0);;_stencil_val(u.x,0,0,0);
      
#line 134
_stencil_val_a(res.x,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.x,{-1}}) );
_stencil_val(u.x,0,0,0);
 
#line 137
_stencil_val_r(res.x,0,0,0);    
      }
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }
          
    
#line 141
} 
#line 129
{
      scalar s = u.y;   
      
      
 { _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0); _stencil_val(mu.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      }
 
#line 133
{ _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0); _stencil_val(mu.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      } 
_stencil_val(r.y,0,0,0);;_stencil_val(u.y,0,0,0);
      
#line 134
_stencil_val_a(res.y,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.y,{-1}}) );
_stencil_val(u.y,0,0,0);
 
#line 137
_stencil_val_r(res.y,0,0,0);    
      }
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }
          
    
#line 141
}}end_foreach_stencil();
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 128
foreach ()
    { {
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
 
#line 133
a += val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
      val(res.x,0,0,0) = val(r.x,0,0,0) - (_const_rho + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 129
{
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - val(mu.y,0,1,0)*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
 
#line 133
a += val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - val(mu.x,1,0,0)*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
      val(res.y,0,0,0) = val(r.y,0,0,0) - (_const_rho + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 141
}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"mu","vector",(void *)&mu,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 71 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_gradient_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))\n    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +\n     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);\n  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 70 \"/home/Min/basilisk/src/embed.h\"\n\nstatic inline real embed_face_gradient_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert (\"/home/Min/basilisk/src/embed.h\", 74, \"cs[i] && cs[i-1]\");\n  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))\n    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +\n     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);\n  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;\n}","\n    \n// #line 129 \"/home/Min/basilisk/src/viscosity-embed.h\"\n{ {\n      scalar s = u.x;\n      real a = 0.;\n      \n a += _const_mu.x*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n \n// #line 133\na += _const_mu.y*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - (_const_rho + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.x, mu, &c);\n val_out_(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));\n      }\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    } \n// #line 129\n{\n      scalar s = u.y;\n      real a = 0.;\n      \n a += _const_mu.y*(s.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(s.third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);\n \n// #line 133\na += _const_mu.x*(s.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(s.third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - (_const_rho + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;\n      if (embed_flux) {\n real c, d = embed_flux (point, u.y, mu, &c);\n val_out_(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));\n      }\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }}"})
    { {
      scalar s = u.x;   
      
      
 {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      }
 
#line 133
{; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      } 
_stencil_val(r.x,0,0,0);;_stencil_val(u.x,0,0,0);
      
#line 134
_stencil_val_a(res.x,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.x,{-1}}) );
_stencil_val(u.x,0,0,0);
 
#line 137
_stencil_val_r(res.x,0,0,0);    
      }
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }
          
    
#line 141
} 
#line 129
{
      scalar s = u.y;   
      
      
 {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);; _stencil_val(fs.y,0,1,0); _stencil_val(fs.y,0,1,0); _stencil_embed_face_gradient_y (point, s, 1);_stencil_val(s,0,1,0); _stencil_val(s,0,1 -1,0);                      }
 
#line 133
{; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);; _stencil_val(fs.x,1,0,0); _stencil_val(fs.x,1,0,0); _stencil_embed_face_gradient_x (point, s, 1);_stencil_val(s,1,0,0); _stencil_val(s,1 -1,0,0);                      } 
_stencil_val(r.y,0,0,0);;_stencil_val(u.y,0,0,0);
      
#line 134
_stencil_val_a(res.y,0,0,0);       
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){mu.x , mu.y, u.y,{-1}}) );
_stencil_val(u.y,0,0,0);
 
#line 137
_stencil_val_r(res.y,0,0,0);    
      }
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }
          
    
#line 141
}}end_foreach_stencil();
#line 128 "/home/Min/basilisk/src/viscosity-embed.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 128
foreach ()
    { {
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
 
#line 133
a += _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
      val(res.x,0,0,0) = val(r.x,0,0,0) - (_const_rho + ((coord){0.,0.}).x)*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 129
{
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta) - _const_mu.y*(_attribute[s.i].third && val(fs.y,0,1,0) < 1. && val(fs.y,0,1,0) > 0. ? embed_face_gradient_y (point, s, 1) : (val(s,0,1,0) - val(s,0,1 -1,0))/Delta);
 
#line 133
a += _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta) - _const_mu.x*(_attribute[s.i].third && val(fs.x,1,0,0) < 1. && val(fs.x,1,0,0) > 0. ? embed_face_gradient_x (point, s, 1) : (val(s,1,0,0) - val(s,1 -1,0,0))/Delta);
      val(res.y,0,0,0) = val(r.y,0,0,0) - (_const_rho + ((coord){0.,0.}).y)*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 141
}

  return maxres;
}



double TOLERANCE_MU = 0.;

     
mgstats viscosity (vector u, vector mu, scalar rho, double dt,
     int nrelax, scalar * res)
{tracing("viscosity","/home/Min/basilisk/src/viscosity-embed.h",151);
  vector  r=new_vector("r");
  foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"rho","scalar",(void *)&rho,NULL,0},{"r","vector",(void *)&r,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 156 \"/home/Min/basilisk/src/viscosity-embed.h\"\n{\n      val_out_(r.x,0,0,0) = val(rho,0,0,0)*val(u.x,0,0,0);\n      \n// #line 157\nval_out_(r.y,0,0,0) = val(rho,0,0,0)*val(u.y,0,0,0);}"})
    {
      { _stencil_val(rho,0,0,0);_stencil_val(u.x,0,0,0);_stencil_val_a(r.x,0,0,0); }
      
#line 157
{ _stencil_val(rho,0,0,0);_stencil_val(u.y,0,0,0);_stencil_val_a(r.y,0,0,0); }}end_foreach_stencil();
  {
#line 155
foreach()
    {
      val(r.x,0,0,0) = val(rho,0,0,0)*val(u.x,0,0,0);
      
#line 157
val(r.y,0,0,0) = val(rho,0,0,0)*val(u.y,0,0,0);}end_foreach();}

  restriction (((scalar[]){mu.x,mu.y, rho,{-1}}));
  struct Viscosity p = { mu, rho, dt };
  p.embed_flux = _attribute[u.x.i].boundary[embed] != antisymmetry ? embed_flux : NULL;
  { mgstats _ret= mg_solve ((scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}})
,
     
#line 163
residual_diffusion, relax_diffusion, &p, nrelax, res
, 1

, TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity","/home/Min/basilisk/src/viscosity-embed.h",166);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity","/home/Min/basilisk/src/viscosity-embed.h",167);}
#line 32 "/home/Min/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/Min/basilisk/src/navier-stokes/centered.h"
scalar  p={3};
vector  u={{4},{5}},  g={{6},{7}};
scalar  pf={8};
vector  uf={{9},{10}};
#line 70 "/home/Min/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+0},{_NVARMAX+1}}, a = {{_NVARMAX+0},{_NVARMAX+1}}, alpha = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar rho = {_NVARMAX+4};
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}
static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}








static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}
static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)), point, neighbor, _s, data));}}}}}
#line 126
static int defaults_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 126 "/home/Min/basilisk/src/navier-stokes/centered.h"
      static int defaults_1(const int i,const double t,Event *_ev){tracing("defaults_1","/home/Min/basilisk/src/navier-stokes/centered.h",126);
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;





  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    if(!is_constant(fm.x)){
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 146 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}"}){_stencil_is_face_x(){
      { _stencil_val(fm.x,0,0,0);_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      { _stencil_val(fm.y,0,0,0);_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 146 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}"}){_stencil_is_face_x(){
      {;_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      {;_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}end_foreach_face_generic();}}
  }
#line 178 "/home/Min/basilisk/src/navier-stokes/centered.h"
  foreach_stencil(1,{(NonLocal[]){{"t","double",(void *)&t,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 179 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{\n      dimensional (val(u.x,0,0,0) == Delta/t);\n      \n// #line 180\ndimensional (val(u.y,0,0,0) == Delta/t);}"})
    {
      {_stencil_val(u.x,0,0,0);   }
      
#line 180
{_stencil_val(u.y,0,0,0);   }}end_foreach_stencil();
#line 178 "/home/Min/basilisk/src/navier-stokes/centered.h"
  {foreach()
    {
      dimensional (val(u.x,0,0,0) == Delta/t);
      
#line 180
dimensional (val(u.y,0,0,0) == Delta/t);}end_foreach();}
}{end_tracing("defaults_1","/home/Min/basilisk/src/navier-stokes/centered.h",181);return 0;}end_tracing("defaults_1","/home/Min/basilisk/src/navier-stokes/centered.h",181);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}






#line 187
      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/home/Min/basilisk/src/navier-stokes/centered.h",187);
  display ("squares (color = 'u.x', spread = -1);"
#line 1425 "/home/Min/basilisk/src/common.h"
, false
#line 188 "/home/Min/basilisk/src/navier-stokes/centered.h"
);{end_tracing("default_display","/home/Min/basilisk/src/navier-stokes/centered.h",188);return 0;}end_tracing("default_display","/home/Min/basilisk/src/navier-stokes/centered.h",188);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 196
      static int init(const int i,const double t,Event *_ev){tracing("init","/home/Min/basilisk/src/navier-stokes/centered.h",196);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 199 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*(u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*(u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);                 }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);                 }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 199 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*(u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*(u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}"}){_stencil_is_face_x(){
    {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);                 }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);                 }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}end_foreach_face_generic();}}




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/home/Min/basilisk/src/navier-stokes/centered.h",213);return 0;}end_tracing("init","/home/Min/basilisk/src/navier-stokes/centered.h",213);}








static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 222 "/home/Min/basilisk/src/navier-stokes/centered.h"
      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/home/Min/basilisk/src/navier-stokes/centered.h",222); dtmax = DT;{end_tracing("set_dtmax","/home/Min/basilisk/src/navier-stokes/centered.h",222);return 0;}end_tracing("set_dtmax","/home/Min/basilisk/src/navier-stokes/centered.h",222);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 224
      static int stability(const int i,const double t,Event *_ev){tracing("stability","/home/Min/basilisk/src/navier-stokes/centered.h",224); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/home/Min/basilisk/src/navier-stokes/centered.h",226);return 0;}end_tracing("stability","/home/Min/basilisk/src/navier-stokes/centered.h",226);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 234
static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 235
static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 236
static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 243
static int properties(const int i,const double t,Event *_ev){;return 0;}
#line 255 "/home/Min/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 258
{
    scalar s = new_scalar("s");
    du.y = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 264
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 265 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))\n   val_out_(du.x,0,0,0) = 0.;\n else\n\n   val_out_(du.x,0,0,0) = u.x.gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;\n      } \n// #line 265\n{\n\n        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))\n   val_out_(du.y,0,0,0) = 0.;\n else\n\n   val_out_(du.y,0,0,0) = u.y.gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;\n      }}"})
      { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
   {_stencil_val_a(du.x,0,0,0);  }

   
{_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);_stencil_val_a(du.x,0,0,0);   }}

           
 
      
#line 272
} 
#line 265
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
   {_stencil_val_a(du.y,0,0,0);  }

   
{_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);_stencil_val_a(du.y,0,0,0);   }}

           
 
      
#line 272
}}end_foreach_stencil();{
#line 264
foreach()
      { {

        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
   val(du.x,0,0,0) = 0.;
 else

   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 265
{

        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
   val(du.y,0,0,0) = 0.;
 else

   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }}end_foreach();}}
  else
    {
    
#line 274
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 275 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))\n   val_out_(du.x,0,0,0) = 0.;\n else\n\n   val_out_(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);\n    } \n// #line 275\n{\n\n        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))\n   val_out_(du.y,0,0,0) = 0.;\n else\n\n   val_out_(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);\n    }}"})
      { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
   {_stencil_val_a(du.x,0,0,0);  }

   
{_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);_stencil_val_a(du.x,0,0,0);   }}

           
 
    
#line 282
} 
#line 275
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
   {_stencil_val_a(du.y,0,0,0);  }

   
{_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);_stencil_val_a(du.y,0,0,0);   }}

           
 
    
#line 282
}}end_foreach_stencil();{
#line 274
foreach()
      { {

        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
   val(du.x,0,0,0) = 0.;
 else

   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 275
{

        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
   val(du.y,0,0,0) = 0.;
 else

   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }}end_foreach();}}

  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 285 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.x,0,0,0) *= val(fm.x,0,0,0);\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.y,0,0,0) *= val(fm.y,0,0,0);\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    } 







_stencil_val(fm.x,0,0,0);        

      







    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    } 







_stencil_val(fm.y,0,0,0);        

      







    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 285 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (_const_fm.y && _const_fm.y) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.x,0,0,0) *= _const_fm.x;\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (_const_fm.x && _const_fm.x) {\n      real fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.y,0,0,0) *= _const_fm.y;\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

;; {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    }







;        

      







    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

;; {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    }







;        

      







    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()}end_foreach_face_generic();}}

  delete ((scalar *)((vector[]){du,{{-1},{-1}}}));
}
#line 316
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 316 "/home/Min/basilisk/src/navier-stokes/centered.h"
      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/home/Min/basilisk/src/navier-stokes/centered.h",316);
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    advection ((scalar *)((vector[]){u,{{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1}}}));
  }
}{end_tracing("advection_term","/home/Min/basilisk/src/navier-stokes/centered.h",323);return 0;}end_tracing("advection_term","/home/Min/basilisk/src/navier-stokes/centered.h",323);}







static void correction (double dt)
{
  foreach_stencil(1,{(NonLocal[]){{"g","vector",(void *)&g,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 334 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(u.x,0,0,0) += dt*val(g.x,0,0,0);\n      \n// #line 335\nval_out_(u.y,0,0,0) += dt*val(g.y,0,0,0);}"})
    {
      {_stencil_val(g.x,0,0,0);_stencil_val_r(u.x,0,0,0);  }
      
#line 335
{_stencil_val(g.y,0,0,0);_stencil_val_r(u.y,0,0,0);  }}end_foreach_stencil();
  {
#line 333
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 335
val(u.y,0,0,0) += dt*val(g.y,0,0,0);}end_foreach();}
}








static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 345 "/home/Min/basilisk/src/navier-stokes/centered.h"
      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/home/Min/basilisk/src/navier-stokes/centered.h",345);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax
#line 152 "/home/Min/basilisk/src/viscosity-embed.h"
, NULL
#line 349 "/home/Min/basilisk/src/navier-stokes/centered.h"
);
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1}}}));
    foreach_face_stencil(1,{(NonLocal[]){{"af","vector",(void *)&af,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 359 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(af.x,0,0,0) = 0.;}end_is_face_x()\n// #line 359\nis_face_y(){\n      val_out_(af.y,0,0,0) = 0.;}end_is_face_y()}"}){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 359
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 359
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 359
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}
  }
}{end_tracing("viscous_term","/home/Min/basilisk/src/navier-stokes/centered.h",362);return 0;}end_tracing("viscous_term","/home/Min/basilisk/src/navier-stokes/centered.h",362);}
#line 381
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 381 "/home/Min/basilisk/src/navier-stokes/centered.h"
      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/home/Min/basilisk/src/navier-stokes/centered.h",381);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x) && !is_constant(a.x)){
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 384 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*((u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*((u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);                   }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);                   }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 384 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*((u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*((u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}"}){_stencil_is_face_x(){
    {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);                   }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);                   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 384 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*((u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*((u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);                   }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);                   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_y (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));\n  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?\n    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :\n    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));\n}\n\n\n\n// #line 82 \"/home/Min/basilisk/src/embed.h\"\nstatic inline real embed_face_value_x (Point point, scalar a, int i)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));\n  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?\n    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :\n    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));\n}","\n// #line 384 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*((u.x.third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*((u.y.third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}"}){_stencil_is_face_x(){
    {; _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);                   }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {; _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);                   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}
}{end_tracing("acceleration","/home/Min/basilisk/src/navier-stokes/centered.h",386);return 0;}end_tracing("acceleration","/home/Min/basilisk/src/navier-stokes/centered.h",386);}
#line 395 "/home/Min/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    {;; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}"}){_stencil_is_face_x(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}





  trash (((vector[]){g,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 412 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 1e-30);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 1e-30);}"})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);_stencil_val_a(g.y,0,0,0);      }}end_foreach_stencil();{
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 1e-30);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 1e-30);}end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 412 \"/home/Min/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 1e-30);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 1e-30);}"})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);;;_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);;;_stencil_val_a(g.y,0,0,0);      }}end_foreach_stencil();
  {
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 1e-30);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 1e-30);}end_foreach();}}delete((scalar*)((vector[]){gf,{{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 421
      static int projection(const int i,const double t,Event *_ev){tracing("projection","/home/Min/basilisk/src/navier-stokes/centered.h",421);
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/home/Min/basilisk/src/navier-stokes/centered.h",430);return 0;}end_tracing("projection","/home/Min/basilisk/src/navier-stokes/centered.h",430);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 436
static int end_timestep(const int i,const double t,Event *_ev){;return 0;}
#line 4 "reactor.c"

#line 1 "two-phase.h"
#line 1 "/home/Min/basilisk/src/two-phase.h"
#line 13 "/home/Min/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/Min/basilisk/src/vof.h"
#line 27 "/home/Min/basilisk/src/vof.h"

#line 44 "/home/Min/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;








static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
 else
   return (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
  }
  return 0.;
}

#line 55
static double vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
 else
   return (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
  }
  return 0.;
}









#line 55
static void _stencil_vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;           
  
   _stencil_val(c,1,0,0); _stencil_val(c,0,0,0); _stencil_val(c,-1,0,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,-1,0,0); _stencil_val(t,0,0,0); _stencil_val(t,1,0,0);  }
 else
   {_stencil_val(t,1,0,0); _stencil_val(t,-1,0,0);  }
      }
 
{_stencil_val(t,1,0,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,-1,0,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}

#line 55
static void _stencil_vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;           
  
   _stencil_val(c,0,1,0); _stencil_val(c,0,0,0); _stencil_val(c,0,-1,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,0,-1,0); _stencil_val(t,0,0,0); _stencil_val(t,0,1,0);  }
 else
   {_stencil_val(t,0,1,0); _stencil_val(t,0,-1,0);  }
      }
 
{_stencil_val(t,0,1,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,0,-1,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}
#line 127
static int defaults_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 127 "/home/Min/basilisk/src/vof.h"
      static int defaults_2(const int i,const double t,Event *_ev){tracing("defaults_2","/home/Min/basilisk/src/vof.h",127);
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
    scalar * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
      _attribute[t.i].depends = list_add (_attribute[t.i].depends, c);}}
  }}}
}{end_tracing("defaults_2","/home/Min/basilisk/src/vof.h",134);return 0;}end_tracing("defaults_2","/home/Min/basilisk/src/vof.h",134);}





static int stability_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 140
      static int stability_0(const int i,const double t,Event *_ev){tracing("stability_0","/home/Min/basilisk/src/vof.h",140); {
  if (CFL > 0.5)
    CFL = 0.5;
}{end_tracing("stability_0","/home/Min/basilisk/src/vof.h",143);return 0;}end_tracing("stability_0","/home/Min/basilisk/src/vof.h",143);}
#line 157 "/home/Min/basilisk/src/vof.h"

static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/Min/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil(1,{(NonLocal[]){{"c","scalar",(void *)&c,NULL,0},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n\n\n\n\n\n\n// #line 55 \"/home/Min/basilisk/src/vof.h\"\nstatic real vof_concentration_gradient_x (Point point, scalar c, scalar t)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  static const real cmin = 0.5;\n  real cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);\n  if (t.inverse)\n    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;\n  if (cc >= cmin && t.gradient != zero) {\n    if (cr >= cmin) {\n      if (cl >= cmin) {\n if (t.gradient)\n   return t.gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;\n else\n   return (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);\n      }\n      else\n return (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;\n    }\n    else if (cl >= cmin)\n      return (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;\n  }\n  return 0.;\n}"," \n// #line 182 \"/home/Min/basilisk/src/vof.h\"\n{\n      scalar t, gf;\n      {forin2 (t,gf , tracers,gfl)\n val_out_(gf,0,0,0) = vof_concentration_gradient_x (point, c, t); endforin2()}\n    }"}) {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 { _stencil_vof_concentration_gradient_x (point, c, t);_stencil_val_a(gf,0,0,0); }}}
    }end_foreach_stencil();




    {
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_x (point, c, t);}}
    }end_foreach();}
  }






  reconstruction (c, n, alpha);
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_x(){ {\n\n\n\n\n\n\n    real un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 1e-30) > cfl)\n      cfl = un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :\n      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.x,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,i,0,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_x()"})_stencil_is_face_x(){ {       






    _stencil_val(fm.x,0,0,0); _stencil_val(uf.x,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);;_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 1e-30) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_x(){ {\n\n\n\n\n\n\n    real un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*_const_fm.x*s/(val(cm,0,0,0) + 1e-30) > cfl)\n      cfl = un*_const_fm.x*s/(val(cm,0,0,0) + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :\n      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.x,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,i,0,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_x()"})_stencil_is_face_x(){ {






; _stencil_val(uf.x,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {;_stencil_val(cm,0,0,0);
      {;_stencil_val(cm,0,0,0);    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);;_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*_const_fm.x*s/(val(cm,0,0,0) + 1e-30) > cfl)
      cfl = un*_const_fm.x*s/(val(cm,0,0,0) + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_x(){ {\n\n\n\n\n\n\n    real un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*val(fm.x,0,0,0)*s/(_const_cm + 1e-30) > cfl)\n      cfl = un*val(fm.x,0,0,0)*s/(_const_cm + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :\n      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.x,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,i,0,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_x()"})_stencil_is_face_x(){ {       






    _stencil_val(fm.x,0,0,0); _stencil_val(uf.x,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {_stencil_val(fm.x,0,0,0);;
      {_stencil_val(fm.x,0,0,0);;    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);;_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*val(fm.x,0,0,0)*s/(_const_cm + 1e-30) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(_const_cm + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_x(){ {\n\n\n\n\n\n\n    real un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*_const_fm.x*s/(_const_cm + 1e-30) > cfl)\n      cfl = un*_const_fm.x*s/(_const_cm + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :\n      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.x,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,i,0,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_x()"})_stencil_is_face_x(){ {






; _stencil_val(uf.x,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {;;
      {;;    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);;_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*_const_fm.x*s/(_const_cm + 1e-30) > cfl)
      cfl = un*_const_fm.x*s/(_const_cm + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 303 "/home/Min/basilisk/src/vof.h"
  foreach_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"tcl","scalar",(void *)tcl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"uf","vector",(void *)&uf,NULL,0},{"cc","scalar",(void *)&cc,NULL,0},{"flux","scalar",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 304 \"/home/Min/basilisk/src/vof.h\"\nif (val(cs,0,0,0) > 0.) {\n      val_out_(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/Delta;\n\n\n\n\n      scalar t, tc, tflux;\n      {forin3 (t, tc, tflux , tracers, tcl, tfluxl)\n val_out_(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/Delta; endforin3()}\n\n    }"})
    {_stencil_val(cs,0,0,0); {
_stencil_val(flux,0,0,0); _stencil_val(flux,1,0,0); _stencil_val(cc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);
      
#line 305
_stencil_val_r(c,0,0,0);     




      scalar t, tc, tflux;
      {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
 {_stencil_val(tflux,0,0,0); _stencil_val(tflux,1,0,0); _stencil_val(tc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val_r(t,0,0,0);     }}}

    }   }end_foreach_stencil();
#line 303 "/home/Min/basilisk/src/vof.h"
  {foreach()
    if (val(cs,0,0,0) > 0.) {
      val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/Delta;




      scalar t, tc, tflux;
      {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
 val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/Delta;}}

    }end_foreach();}


  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,{-1}}));
}

#line 158
static void sweep_y (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/Min/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil(1,{(NonLocal[]){{"c","scalar",(void *)&c,NULL,0},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n// #line 55 \"/home/Min/basilisk/src/vof.h\"\nstatic real vof_concentration_gradient_y (Point point, scalar c, scalar t)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  static const real cmin = 0.5;\n  real cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);\n  if (t.inverse)\n    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;\n  if (cc >= cmin && t.gradient != zero) {\n    if (cr >= cmin) {\n      if (cl >= cmin) {\n if (t.gradient)\n   return t.gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;\n else\n   return (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);\n      }\n      else\n return (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;\n    }\n    else if (cl >= cmin)\n      return (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;\n  }\n  return 0.;\n}"," \n// #line 182 \"/home/Min/basilisk/src/vof.h\"\n{\n      scalar t, gf;\n      {forin2 (t,gf , tracers,gfl)\n val_out_(gf,0,0,0) = vof_concentration_gradient_y (point, c, t); endforin2()}\n    }"}) {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 { _stencil_vof_concentration_gradient_y (point, c, t);_stencil_val_a(gf,0,0,0); }}}
    }end_foreach_stencil();




    {
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_y (point, c, t);}}
    }end_foreach();}
  }






  reconstruction (c, n, alpha);
  if(!is_constant(fm.y) && !is_constant(cm)){
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_y(){ {\n\n\n\n\n\n\n    real un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 1e-30) > cfl)\n      cfl = un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :\n      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.y,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,0,i,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_y()"})_stencil_is_face_y(){ {       






    _stencil_val(fm.y,0,0,0); _stencil_val(uf.y,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);;_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 1e-30) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(is_constant(fm.y) && !is_constant(cm)){_coord _const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_y(){ {\n\n\n\n\n\n\n    real un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*_const_fm.y*s/(val(cm,0,0,0) + 1e-30) > cfl)\n      cfl = un*_const_fm.y*s/(val(cm,0,0,0) + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :\n      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.y,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,0,i,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_y()"})_stencil_is_face_y(){ {






; _stencil_val(uf.y,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {;_stencil_val(cm,0,0,0);
      {;_stencil_val(cm,0,0,0);    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);;_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*_const_fm.y*s/(val(cm,0,0,0) + 1e-30) > cfl)
      cfl = un*_const_fm.y*s/(val(cm,0,0,0) + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(!is_constant(fm.y) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_y(){ {\n\n\n\n\n\n\n    real un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*val(fm.y,0,0,0)*s/(_const_cm + 1e-30) > cfl)\n      cfl = un*val(fm.y,0,0,0)*s/(_const_cm + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :\n      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.y,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,0,i,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_y()"})_stencil_is_face_y(){ {       






    _stencil_val(fm.y,0,0,0); _stencil_val(uf.y,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {_stencil_val(fm.y,0,0,0);;
      {_stencil_val(fm.y,0,0,0);;    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);;_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*val(fm.y,0,0,0)*s/(_const_cm + 1e-30) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(_const_cm + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else {_coord _const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"gfl","scalar",(void *)gfl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"flux","scalar",(void *)&flux,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"_val_higher_dimension","double",(void *)&_val_higher_dimension,NULL,0},{"n","vector",(void *)&n,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cfl","double",(void *)&cfl,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/Min/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}\n// # 267 \"/home/Min/basilisk/src/geometry.h\"\nreal rectangle_fraction (coord n, real alpha, coord a, coord b)\n{\n  coord n1;\n   {\n    alpha -= n.x*(b.x + a.x)/2.;\n    n1.x = n.x*(b.x - a.x);\n  } \n// #line 270\n{\n    alpha -= n.y*(b.y + a.y)/2.;\n    n1.y = n.y*(b.y - a.y);\n  }\n  return line_area(n1.x, n1.y, alpha);\n}","\n// #line 195 \"/home/Min/basilisk/src/vof.h\"\nis_face_y(){ {\n\n\n\n\n\n\n    real un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 1e-30), s = sign(un);\n    int i = -(s + 1.)/2.;\n\n\n\n\n\n    if (val(cs,0,0,0) >= 1.)\n\n    if (un*_const_fm.y*s/(_const_cm + 1e-30) > cfl)\n      cfl = un*_const_fm.y*s/(_const_cm + 1e-30);\n// # 225 \"/home/Min/basilisk/src/vof.h\"\n    real cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :\n      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),\n     (coord){-0.5, -0.5, -0.5},\n     (coord){s*un - 0.5, 0.5, 0.5});\n\n\n\n\n\n    val_out_(flux,0,0,0) = cf*val(uf.y,0,0,0);\n\n\n\n\n\n\n    scalar t, gf, tflux;\n    {forin3 (t,gf,tflux , tracers,gfl,tfluxl) {\n      real cf1 = cf, ci = val(c,0,i,0);\n      if (t.inverse)\n cf1 = 1. - cf1, ci = 1. - ci;\n      if (ci > 1e-10) {\n real ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;\n val_out_(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);\n      }\n      else\n val_out_(tflux,0,0,0) = 0.;\n    } endforin3()}\n  }}end_is_face_y()"})_stencil_is_face_y(){ {






; _stencil_val(uf.y,0,0,0);     
    





_stencil_val(cs,0,0,0);

    {;;
      {;;    }     }              
       
      
      





       
#line 225 "/home/Min/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);;_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;





    if (val(cs,0,0,0) >= 1.)

    if (un*_const_fm.y*s/(_const_cm + 1e-30) > cfl)
      cfl = un*_const_fm.y*s/(_const_cm + 1e-30);
#line 225 "/home/Min/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 303 "/home/Min/basilisk/src/vof.h"
  foreach_stencil(1,{(NonLocal[]){{"tfluxl","scalar",(void *)tfluxl,NULL,1},{"tcl","scalar",(void *)tcl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"uf","vector",(void *)&uf,NULL,0},{"cc","scalar",(void *)&cc,NULL,0},{"flux","scalar",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"cs","scalar",(void *)&cs,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 304 \"/home/Min/basilisk/src/vof.h\"\nif (val(cs,0,0,0) > 0.) {\n      val_out_(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/Delta;\n\n\n\n\n      scalar t, tc, tflux;\n      {forin3 (t, tc, tflux , tracers, tcl, tfluxl)\n val_out_(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/Delta; endforin3()}\n\n    }"})
    {_stencil_val(cs,0,0,0); {
_stencil_val(flux,0,0,0); _stencil_val(flux,0,1,0); _stencil_val(cc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);
      
#line 305
_stencil_val_r(c,0,0,0);     




      scalar t, tc, tflux;
      {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
 {_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,1,0); _stencil_val(tc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val_r(t,0,0,0);     }}}

    }   }end_foreach_stencil();
#line 303 "/home/Min/basilisk/src/vof.h"
  {foreach()
    if (val(cs,0,0,0) > 0.) {
      val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/Delta;




      scalar t, tc, tflux;
      {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
 val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/Delta;}}

    }end_foreach();}


  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,{-1}}));
}






void vof_advection (scalar * interfaces, int i)
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
#line 337 "/home/Min/basilisk/src/vof.h"
    scalar  cc=new_scalar("cc"), * tcl = NULL, * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {

      scalar tc = new_scalar("tc");
      tcl = list_append (tcl, tc);
#line 351 "/home/Min/basilisk/src/vof.h"
    }}}
    foreach_stencil(1,{(NonLocal[]){{"tcl","scalar",(void *)tcl,NULL,1},{"tracers","scalar",(void *)tracers,NULL,1},{"c","scalar",(void *)&c,NULL,0},{"cc","scalar",(void *)&cc,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 352 \"/home/Min/basilisk/src/vof.h\"\n{\n      val_out_(cc,0,0,0) = (val(c,0,0,0) > 0.5);\n\n      scalar t, tc;\n      {forin2 (t, tc , tracers, tcl) {\n if (t.inverse)\n   val_out_(tc,0,0,0) = val(c,0,0,0) < 0.5 ? val(t,0,0,0)/(1. - val(c,0,0,0)) : 0.;\n else\n   val_out_(tc,0,0,0) = val(c,0,0,0) > 0.5 ? val(t,0,0,0)/val(c,0,0,0) : 0.;\n      } endforin2()}\n\n    }"}) {
_stencil_val(c,0,0,0);
      
#line 353
_stencil_val_a(cc,0,0,0);    

      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   { _stencil_val(c,0,0,0); _stencil_val(t,0,0,0); _stencil_val(c,0,0,0);_stencil_val_a(tc,0,0,0);       }
 else
   { _stencil_val(c,0,0,0); _stencil_val(t,0,0,0);_stencil_val(c,0,0,0);_stencil_val_a(tc,0,0,0);      }
      }}}

    }end_foreach_stencil();
    {
#line 352
foreach() {
      val(cc,0,0,0) = (val(c,0,0,0) > 0.5);

      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   val(tc,0,0,0) = val(c,0,0,0) < 0.5 ? val(t,0,0,0)/(1. - val(c,0,0,0)) : 0.;
 else
   val(tc,0,0,0) = val(c,0,0,0) > 0.5 ? val(t,0,0,0)/val(c,0,0,0) : 0.;
      }}}

    }end_foreach();}






    void (* sweep[2]) (scalar, scalar, scalar *);
    int d = 0;
    
      sweep[d++] = sweep_x;
      
#line 373
sweep[d++] = sweep_y;
    for (d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc, tcl);
    delete (tcl), pfree (tcl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){cc,{-1}}));
  }}}
}

static int vof_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 380
      static int vof_0(const int i,const double t,Event *_ev){tracing("vof_0","/home/Min/basilisk/src/vof.h",380);
  vof_advection (interfaces, i);{end_tracing("vof_0","/home/Min/basilisk/src/vof.h",381);return 0;}end_tracing("vof_0","/home/Min/basilisk/src/vof.h",381);}
#line 14 "/home/Min/basilisk/src/two-phase.h"

scalar  f={11}, * interfaces =((scalar[]) {{11},{-1}});

#line 1 "two-phase-generic.h"
#line 1 "/home/Min/basilisk/src/two-phase-generic.h"
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector  alphav={{12},{13}};
scalar  rhov={14};

static int defaults_3_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 10
      static int defaults_3(const int i,const double t,Event *_ev){tracing("defaults_3","/home/Min/basilisk/src/two-phase-generic.h",10);
{
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");




  display ("draw_vof (c = 'f');"
#line 1425 "/home/Min/basilisk/src/common.h"
, false
#line 25 "/home/Min/basilisk/src/two-phase-generic.h"
);
}{end_tracing("defaults_3","/home/Min/basilisk/src/two-phase-generic.h",26);return 0;}end_tracing("defaults_3","/home/Min/basilisk/src/two-phase-generic.h",26);}
#line 50
static int tracer_advection_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 50 "/home/Min/basilisk/src/two-phase-generic.h"
static int tracer_advection_0(const int i,const double t,Event *_ev){
{
#line 79 "/home/Min/basilisk/src/two-phase-generic.h"
}return 0;}



static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}




#line 83
      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","/home/Min/basilisk/src/two-phase-generic.h",83);
{
  if(!is_constant(fm.x)){
  
#line 85
foreach_face_stencil(1,{(NonLocal[]){{"mu","vector",(void *)&mu,NULL,0},{"mu2","double",(void *)&mu2,NULL,0},{"mu1","double",(void *)&mu1,NULL,0},{"rho2","double",(void *)&rho2,NULL,0},{"rho1","double",(void *)&rho1,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 85 \"/home/Min/basilisk/src/two-phase-generic.h\"\n{is_face_x(){ {\n    real ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;\n    val_out_(alphav.x,0,0,0) = val(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);\n    if (mu1 || mu2) {\n      vector muv = mu;\n      val_out_(muv.x,0,0,0) = val(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);\n    }\n  }}end_is_face_x()\n// #line 85\nis_face_y(){ {\n    real ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;\n    val_out_(alphav.y,0,0,0) = val(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);\n    if (mu1 || mu2) {\n      vector muv = mu;\n      val_out_(muv.y,0,0,0) = val(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);\n    }\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {    
     _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0); 
_stencil_val(fm.x,0,0,0);
    
#line 87
_stencil_val_a(alphav.x,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu; 
_stencil_val(fm.x,0,0,0);
      
#line 90
_stencil_val_a(muv.x,0,0,0);     
    }
  }}end__stencil_is_face_x()
#line 85
_stencil_is_face_y(){ {    
     _stencil_val(f,0,-1,0);_stencil_val(f,0,0,0); 
_stencil_val(fm.y,0,0,0);
    
#line 87
_stencil_val_a(alphav.y,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu; 
_stencil_val(fm.y,0,0,0);
      
#line 90
_stencil_val_a(muv.y,0,0,0);     
    }
  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 85
foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 85
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 85
foreach_face_stencil(1,{(NonLocal[]){{"mu","vector",(void *)&mu,NULL,0},{"mu2","double",(void *)&mu2,NULL,0},{"mu1","double",(void *)&mu1,NULL,0},{"rho2","double",(void *)&rho2,NULL,0},{"rho1","double",(void *)&rho1,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 85 \"/home/Min/basilisk/src/two-phase-generic.h\"\n{is_face_x(){ {\n    real ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;\n    val_out_(alphav.x,0,0,0) = _const_fm.x/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);\n    if (mu1 || mu2) {\n      vector muv = mu;\n      val_out_(muv.x,0,0,0) = _const_fm.x*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);\n    }\n  }}end_is_face_x()\n// #line 85\nis_face_y(){ {\n    real ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;\n    val_out_(alphav.y,0,0,0) = _const_fm.y/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);\n    if (mu1 || mu2) {\n      vector muv = mu;\n      val_out_(muv.y,0,0,0) = _const_fm.y*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);\n    }\n  }}end_is_face_y()}"}){_stencil_is_face_x(){ {    
     _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);
;
    
#line 87
_stencil_val_a(alphav.x,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
;
      
#line 90
_stencil_val_a(muv.x,0,0,0);     
    }
  }}end__stencil_is_face_x()
#line 85
_stencil_is_face_y(){ {    
     _stencil_val(f,0,-1,0);_stencil_val(f,0,0,0);
;
    
#line 87
_stencil_val_a(alphav.y,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
;
      
#line 90
_stencil_val_a(muv.y,0,0,0);     
    }
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 85
foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = _const_fm.x/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = _const_fm.x*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 85
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = _const_fm.y/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = _const_fm.y*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()}end_foreach_face_generic();}}

  if(!is_constant(cm)){

  
#line 94
foreach_stencil(1,{(NonLocal[]){{"rho2","double",(void *)&rho2,NULL,0},{"rho1","double",(void *)&rho1,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"rhov","scalar",(void *)&rhov,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 95 \"/home/Min/basilisk/src/two-phase-generic.h\"\nval_out_(rhov,0,0,0) = val(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);"})
    { _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);_stencil_val_a(rhov,0,0,0);     }end_foreach_stencil();{
#line 94
foreach()
    val(rhov,0,0,0) = val(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

  
#line 94
foreach_stencil(1,{(NonLocal[]){{"rho2","double",(void *)&rho2,NULL,0},{"rho1","double",(void *)&rho1,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"rhov","scalar",(void *)&rhov,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 95 \"/home/Min/basilisk/src/two-phase-generic.h\"\nval_out_(rhov,0,0,0) = _const_cm*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);"})
    {;_stencil_val(f,0,0,0);_stencil_val_a(rhov,0,0,0);     }end_foreach_stencil();

  {
#line 94
foreach()
    val(rhov,0,0,0) = _const_cm*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}}





}{end_tracing("properties_0","/home/Min/basilisk/src/two-phase-generic.h",101);return 0;}end_tracing("properties_0","/home/Min/basilisk/src/two-phase-generic.h",101);}
#line 18 "/home/Min/basilisk/src/two-phase.h"
#line 6 "reactor.c"

#line 1 "curvature.h"
#line 1 "/home/Min/basilisk/src/curvature.h"
#line 68 "/home/Min/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/Min/basilisk/src/heights.h"
#line 29 "/home/Min/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/Min/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  const int complete = -1;

   {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = 1e30;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/Min/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/Min/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/Min/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/Min/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = 1e30;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);
    }
  } 
#line 59
{







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = 1e30;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/Min/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/Min/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/Min/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/Min/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = 1e30;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);
    }
  }
}
#line 49 "/home/Min/basilisk/src/heights.h"
static void _stencil_half_column (Point point, scalar c, vector h, vector cs, int j)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    






  

   {       







     _stencil_val_o(c,0,0,0);            







    
    
    if (j == 1) {




_stencil_val_o(h.x,0,0,0);{ 
      




{     
 _stencil_val_o(h.x,0,0,0); 
_stencil_val_o(h.x,0,0,0);
    
     
      
#line 92
}}   




         




      





      
      
    
#line 100
}
#line 109 "/home/Min/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) { 
_stencil_val_o(c,i*j,0,0); _stencil_val_o(cs.x,(i - 2)*j,0,0); 
           
  
 
       
         
  
 
        







     
   
  
   
        
        
           
        




             
#line 138 "/home/Min/basilisk/src/heights.h"
              
              
#line 156 "/home/Min/basilisk/src/heights.h"
      
        
    }





    if (j == -1) {







_stencil_val_o(c,0,0,0); _stencil_val_o(c,0,0,0);{
 
{_stencil_val_a(h.x,0,0,0);  }
{
 {_stencil_val_a(h.x,0,0,0);  }





 
{_stencil_val_a(h.x,0,0,0);        }}    
      
#line 183
}







                
              
      
    
#line 184
}
    else {
#line 203
{
 {_stencil_val_a(h.x,0,0,0);  }
 
{_stencil_val_a(h.x,0,0,0);        }}
      
#line 195 "/home/Min/basilisk/src/heights.h"
                
   





         
      
    


}
  } 
#line 59
{       







     _stencil_val_o(c,0,0,0);            







    
    
    if (j == 1) {




_stencil_val_o(h.y,0,0,0);{ 
      




{     
 _stencil_val_o(h.y,0,0,0); 
_stencil_val_o(h.y,0,0,0);
    
     
      
#line 92
}}   




         




      





      
      
    
#line 100
}
#line 109 "/home/Min/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) { 
_stencil_val_o(c,0,i*j,0); _stencil_val_o(cs.y,0,(i - 2)*j,0); 
           
  
 
       
         
  
 
        







     
   
  
   
        
        
           
        




             
#line 138 "/home/Min/basilisk/src/heights.h"
              
              
#line 156 "/home/Min/basilisk/src/heights.h"
      
        
    }





    if (j == -1) {







_stencil_val_o(c,0,0,0); _stencil_val_o(c,0,0,0);{
 
{_stencil_val_a(h.y,0,0,0);  }
{
 {_stencil_val_a(h.y,0,0,0);  }





 
{_stencil_val_a(h.y,0,0,0);        }}    
      
#line 183
}







                
              
      
    
#line 184
}
    else {
#line 203
{
 {_stencil_val_a(h.y,0,0,0);  }
 
{_stencil_val_a(h.y,0,0,0);        }}
      
#line 195 "/home/Min/basilisk/src/heights.h"
                
   





         
      
    


}
  }
}
#line 222 "/home/Min/basilisk/src/heights.h"
static void column_propagation (vector h)
{
  foreach_stencil (0,{(NonLocal[]){{"h","vector",(void *)&h,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},NULL})
    for (int i = -2; i <= 2; i++)
      {
 {_stencil_val(h.x,i,0,0);
_stencil_val(h.x,i,0,0);_stencil_val(h.x,0,0,0);
   { _stencil_val(h.x,i,0,0);_stencil_val_a(h.x,0,0,0);   }      
       
#line 229
}
 
#line 227
{_stencil_val(h.y,0,i,0);
_stencil_val(h.y,0,i,0);_stencil_val(h.y,0,0,0);
   { _stencil_val(h.y,0,i,0);_stencil_val_a(h.y,0,0,0);   }      
       
#line 229
}}end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 224
foreach ()
    for (int i = -2; i <= 2; i++)
      {
 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
 
#line 227
if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;}end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


#line 230
}
#line 239 "/home/Min/basilisk/src/heights.h"
     
void heights (scalar c, vector h)
{tracing("heights","/home/Min/basilisk/src/heights.h",240);







  vector  s=new_vector("s");
  
    for (int i = 0; i < nboundary; i++)
      _attribute[s.x.i].boundary[i] = _attribute[c.i].boundary[i];
    
#line 251
for (int i = 0; i < nboundary; i++)
      _attribute[s.y.i].boundary[i] = _attribute[c.i].boundary[i];






  for (int j = -1; j <= 1; j += 2) {





    foreach_stencil(1,{(NonLocal[]){{"j","int",(void *)&j,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"s","vector",(void *)&s,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 266 \"/home/Min/basilisk/src/heights.h\"\n{\n        val_out_(s.x,0,0,0) = val(c,2*j,0,0);\n        \n// #line 267\nval_out_(s.y,0,0,0) = val(c,0,2*j,0);}"})
      {
        { _stencil_val(c,2*j,0,0);_stencil_val_a(s.x,0,0,0); }
        
#line 267
{ _stencil_val(c,0,2*j,0);_stencil_val_a(s.y,0,0,0); }}end_foreach_stencil();





    {
#line 265
foreach()
      {
        val(s.x,0,0,0) = val(c,2*j,0,0);
        
#line 267
val(s.y,0,0,0) = val(c,0,2*j,0);}end_foreach();}





    foreach_stencil (1,{(NonLocal[]){{"j","int",(void *)&j,NULL,0},{"s","vector",(void *)&s,NULL,0},{"h","vector",(void *)&h,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 18 \"/home/Min/basilisk/src/two-phase.h\" 2\n// # 6 \"reactor.c\" 2\n\n// # 1 \"curvature.h\" 1\n// # 1 \"/home/Min/basilisk/src/curvature.h\"\n// # 68 \"/home/Min/basilisk/src/curvature.h\"\n// # 1 \"heights.h\" 1\n// # 1 \"/home/Min/basilisk/src/heights.h\"\n// # 29 \"/home/Min/basilisk/src/heights.h\"\nstatic inline real height (real H) {\n  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;\n}\n// # 49 \"/home/Min/basilisk/src/heights.h\"\nstatic void half_column (Point point, scalar c, vector h, vector cs, int j)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n  const int complete = -1;\n\n   {\n\n\n\n\n\n\n\n    real S = val(c,0,0,0), H = S, ci, a;\n\n\n\n\n\n\n\n    typedef struct { int s; real h; } HState;\n    HState state = {0, 0};\n    if (j == 1) {\n\n\n\n\n      if (val(h.x,0,0,0) == 300.)\n state.s = complete, state.h = 1e30;\n\n\n\n\n      else {\n int s = (val(h.x,0,0,0) + 20./2.)/100.;\n state.h = val(h.x,0,0,0) - 100.*s;\n state.s = s - 1;\n      }\n\n\n\n\n\n      if (state.s != complete)\n S = state.s, H = state.h;\n    }\n// # 109 \"/home/Min/basilisk/src/heights.h\"\n    for (int i = 1; i <= 4; i++) {\n      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);\n      H += ci;\n\n\n\n\n      if (S > 0. && S < 1.) {\n S = ci;\n if (ci <= 0. || ci >= 1.) {\n\n\n\n\n\n\n\n   H -= i*ci;\n   break;\n }\n      }\n// # 138 \"/home/Min/basilisk/src/heights.h\"\n      else if (S >= 1. && ci <= 0.) {\n H = (H - 0.5)*j + (j == -1)*20.;\n S = complete;\n break;\n      }\n      else if (S <= 0. && ci >= 1.) {\n H = (i + 0.5 - H)*j + (j == 1)*20.;\n S = complete;\n break;\n      }\n// # 156 \"/home/Min/basilisk/src/heights.h\"\n      else if (S == ci && modf(H, &a))\n break;\n    }\n\n\n\n\n\n    if (j == -1) {\n\n\n\n\n\n\n\n      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||\n       (S > 0. && S < 1.)))\n val_out_(h.x,0,0,0) = 300.;\n      else if (S == complete)\n val_out_(h.x,0,0,0) = H;\n      else\n\n\n\n\n\n val_out_(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));\n    }\n    else {\n// # 195 \"/home/Min/basilisk/src/heights.h\"\n      if (state.s != complete ||\n   (S == complete && fabs(height(H)) < fabs(height(state.h))))\n state.s = S, state.h = H;\n\n\n\n\n\n      if (state.s != complete)\n val_out_(h.x,0,0,0) = 1e30;\n      else\n val_out_(h.x,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);\n    }\n  } \n// #line 59\n{\n\n\n\n\n\n\n\n    real S = val(c,0,0,0), H = S, ci, a;\n\n\n\n\n\n\n\n    typedef struct { int s; real h; } HState;\n    HState state = {0, 0};\n    if (j == 1) {\n\n\n\n\n      if (val(h.y,0,0,0) == 300.)\n state.s = complete, state.h = 1e30;\n\n\n\n\n      else {\n int s = (val(h.y,0,0,0) + 20./2.)/100.;\n state.h = val(h.y,0,0,0) - 100.*s;\n state.s = s - 1;\n      }\n\n\n\n\n\n      if (state.s != complete)\n S = state.s, H = state.h;\n    }\n// # 109 \"/home/Min/basilisk/src/heights.h\"\n    for (int i = 1; i <= 4; i++) {\n      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);\n      H += ci;\n\n\n\n\n      if (S > 0. && S < 1.) {\n S = ci;\n if (ci <= 0. || ci >= 1.) {\n\n\n\n\n\n\n\n   H -= i*ci;\n   break;\n }\n      }\n// # 138 \"/home/Min/basilisk/src/heights.h\"\n      else if (S >= 1. && ci <= 0.) {\n H = (H - 0.5)*j + (j == -1)*20.;\n S = complete;\n break;\n      }\n      else if (S <= 0. && ci >= 1.) {\n H = (i + 0.5 - H)*j + (j == 1)*20.;\n S = complete;\n break;\n      }\n// # 156 \"/home/Min/basilisk/src/heights.h\"\n      else if (S == ci && modf(H, &a))\n break;\n    }\n\n\n\n\n\n    if (j == -1) {\n\n\n\n\n\n\n\n      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||\n       (S > 0. && S < 1.)))\n val_out_(h.y,0,0,0) = 300.;\n      else if (S == complete)\n val_out_(h.y,0,0,0) = H;\n      else\n\n\n\n\n\n val_out_(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));\n    }\n    else {\n// # 195 \"/home/Min/basilisk/src/heights.h\"\n      if (state.s != complete ||\n   (S == complete && fabs(height(H)) < fabs(height(state.h))))\n state.s = S, state.h = H;\n\n\n\n\n\n      if (state.s != complete)\n val_out_(h.y,0,0,0) = 1e30;\n      else\n val_out_(h.y,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);\n    }\n  }\n}","\n      \n// #line 274 \"/home/Min/basilisk/src/heights.h\"\nhalf_column (point, c, h, s, j);"})
      _stencil_half_column (point, c, h, s, j);end_foreach_stencil();





    {
#line 273
foreach ()
      half_column (point, c, h, s, j);end_foreach();}
  }




  column_propagation (h);delete((scalar*)((vector[]){s,{{-1},{-1}}}));
end_tracing("heights","/home/Min/basilisk/src/heights.h",281);}
#line 455 "/home/Min/basilisk/src/heights.h"

#line 69 "/home/Min/basilisk/src/curvature.h"



static double kappa_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.y,i,0,0) == 1e30 || orientation(val(h.y,i,0,0)) != ori)
      return 1e30;
  double hx = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}

#line 72
static double kappa_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.x,0,i,0) == 1e30 || orientation(val(h.x,0,i,0)) != ori)
      return 1e30;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}


static coord normal_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n = {1e30, 1e30, 1e30};
  if (val(h.y,0,0,0) == 1e30)
    return n;
  int ori = orientation(val(h.y,0,0,0));
  if (val(h.y,-1,0,0) != 1e30 && orientation(val(h.y,-1,0,0)) == ori) {
    if (val(h.y,1,0,0) != 1e30 && orientation(val(h.y,1,0,0)) == ori)
      n.x = (val(h.y,-1,0,0) - val(h.y,1,0,0))/2.;
    else
      n.x = val(h.y,-1,0,0) - val(h.y,0,0,0);
  }
  else if (val(h.y,1,0,0) != 1e30 && orientation(val(h.y,1,0,0)) == ori)
    n.x = val(h.y,0,0,0) - val(h.y,1,0,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;
}

#line 84
static coord normal_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n = {1e30, 1e30, 1e30};
  if (val(h.x,0,0,0) == 1e30)
    return n;
  int ori = orientation(val(h.x,0,0,0));
  if (val(h.x,0,-1,0) != 1e30 && orientation(val(h.x,0,-1,0)) == ori) {
    if (val(h.x,0,1,0) != 1e30 && orientation(val(h.x,0,1,0)) == ori)
      n.y = (val(h.x,0,-1,0) - val(h.x,0,1,0))/2.;
    else
      n.y = val(h.x,0,-1,0) - val(h.x,0,0,0);
  }
  else if (val(h.x,0,1,0) != 1e30 && orientation(val(h.x,0,1,0)) == ori)
    n.y = val(h.x,0,0,0) - val(h.x,0,1,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.y));
  n.y /= nn;
  n.x = 1./nn;
  return n;
}
#line 181 "/home/Min/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
    
#line 195
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormKappa __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 213 "/home/Min/basilisk/src/curvature.h"
  double kappa = 1e30;
  
    if (kappa == 1e30) {
      kappa = n.x.kappa (point, h);
      if (kappa != 1e30) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
    
#line 215
if (kappa == 1e30) {
      kappa = n.y.kappa (point, h);
      if (kappa != 1e30) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }

  if (kappa != 1e30) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 249 "/home/Min/basilisk/src/curvature.h"
  }

  return kappa;
}
#line 181 "/home/Min/basilisk/src/curvature.h"
static void _stencil_height_curvature (Point point, scalar c,_stencil_undefined * h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;        
      
     
      






  
  
  
    { _stencil_val(c,1,0,0); _stencil_val(c,-1,0,0);     }
    
#line 195
{ _stencil_val(c,0,1,0); _stencil_val(c,0,-1,0);     }                      
  
      




     
#line 213 "/home/Min/basilisk/src/curvature.h"
     
     
           
      
          
  
 
      
       
     

  
        




       
#line 249 "/home/Min/basilisk/src/curvature.h"
   

  return ;
}






coord height_normal (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.normal = normal_x;
    
#line 273
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.normal = normal_y;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormNormal __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 290 "/home/Min/basilisk/src/curvature.h"
  coord normal = {1e30, 1e30, 1e30};
  
    if (normal.x == 1e30)
      normal = n.x.normal (point, h);
    
#line 292
if (normal.y == 1e30)
      normal = n.y.normal (point, h);

  return normal;
}
#line 332 "/home/Min/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/Min/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/Min/basilisk/src/parabola.h"






typedef struct {
  coord o;

  coord m;
  double ** M, rhs[3], a[3];
#line 21 "/home/Min/basilisk/src/parabola.h"
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  
    p->o.x = o.x;
    
#line 26
p->o.y = o.y;

  
    p->m.x = m.x;
    
#line 29
p->m.y = m.y;
  normalize (&p->m);
  int n = 3;
#line 65 "/home/Min/basilisk/src/parabola.h"
  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{

  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#line 111 "/home/Min/basilisk/src/parabola.h"
}

static double parabola_fit_solve (ParabolaFit * p)
{

  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  double pivmin = matrix_inverse (p->M, 3, 1e-10);
  if (pivmin) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else
    p->a[0] = p->a[1] = 0.;
#line 158 "/home/Min/basilisk/src/parabola.h"
  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;

  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#line 190 "/home/Min/basilisk/src/parabola.h"
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 333 "/home/Min/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      
 d2 += sq(p[i].x - p[j].x);
 
#line 349
d2 += sq(p[i].y - p[j].y);
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  coord ip[2 == 2 ? 6 : 27];
  int n = 0;




   {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != 1e30) {
 if (orientation(val(h.y,i,0,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != 1e30 && orientation(val(h.y,i,0,0)) == ori)
 ip[n].x = i, ip[n++].y = height(val(h.y,i,0,0));






  } 
#line 375
{





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != 1e30) {
 if (orientation(val(h.x,0,i,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != 1e30 && orientation(val(h.x,0,i,0)) == ori)
 ip[n].y = i, ip[n++].x = height(val(h.x,0,i,0));






  }





  if (independents (ip, n) < (2 == 2 ? 3 : 9))
    return 1e30;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  double area = line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, .1);
#line 440 "/home/Min/basilisk/src/curvature.h"
  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 362
static void _stencil_height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;          





  
  




   {      





    

    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.y,i,0,0); {
_stencil_val(h.y,i,0,0);  
   
      
#line 386
}   }     







    







    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.y,i,0,0);_stencil_val(h.y,i,0,0);
 {_stencil_val(h.y,i,0,0);     }       }






  } 
#line 375
{      





    

    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.x,0,i,0); {
_stencil_val(h.x,0,i,0);  
   
      
#line 386
}   }     







    







    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.x,0,i,0);_stencil_val(h.x,0,i,0);
 {_stencil_val(h.x,0,i,0);     }       }






  }    
    





             





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);          
  
                 
  

  
  
#line 440 "/home/Min/basilisk/src/curvature.h"
     
    
  
  



  return ;
}






static double centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
  {foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = line_alpha (val(c,0,0,0), m);
      double area = line_length_center(m,alpha,&fc);
      coord rn = {x,y,z};
      
 fc.x += (rn.x - r.x)/Delta;
 
#line 480
fc.y += (rn.y - r.y)/Delta;
      parabola_fit_add (&fit, fc, area);
    }end_foreach_neighbor()}
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 455
static void _stencil_centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);    
  
     
  





  
  {foreach_neighbor(1)
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0);      
      
      
       
       
      
    }      }end_foreach_neighbor()}       
  
  



  return ;
}
#line 504 "/home/Min/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) <= 0.)
   return true;
 
#line 509
if (val(c,0,i,0) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) >= 1.)
   return true;
 
#line 515
if (val(c,0,i,0) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}
#line 504 "/home/Min/basilisk/src/curvature.h"
static void _stencil_interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 509
{_stencil_val(c,0,i,0); 
      }}
  } 
{_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 515
{_stencil_val(c,0,i,0); 
      }}
  } 
    
}   
  
#line 519
}}
     
  
  
#line 520
return ;
}
#line 533 "/home/Min/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

     
cstats curvature (scalar c, scalar kappa,
    double sigma, bool add)
{tracing("curvature","/home/Min/basilisk/src/curvature.h",541);
  int sh = 0, f = 0, sa = 0, sc = 0;
#line 557 "/home/Min/basilisk/src/curvature.h"
  vector ch = _attribute[c.i].height,   h=(ch).x.i>0?(ch):new_vector("h");
  if (!ch.x.i)
    heights (c, h);





  scalar  k=new_scalar("k");
  scalar_clone (k, kappa);

  foreach_stencil(1,{(NonLocal[]){{"f","int",(void *)&f,NULL,0,'+'},{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"sh","int",(void *)&sh,NULL,0,'+'},{"h","vector",(void *)&h,NULL,0},{"k","scalar",(void *)&k,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 162 \"/home/Min/basilisk/src/parabola.h\"\nstatic real parabola_fit_curvature (ParabolaFit * p,\n          real kappamax, real * kmax)\n{\n  real kappa;\n\n  real dnm = 1. + sq(p->a[1]);\n  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);\n  if (kmax)\n    *kmax = fabs (kappa);\n// # 190 \"/home/Min/basilisk/src/parabola.h\"\n  if (fabs (kappa) > kappamax) {\n    if (kmax)\n      *kmax = kappamax;\n    return kappa > 0. ? kappamax : - kappamax;\n  }\n  return kappa;\n}\n\n\n// #line 200 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 264 \"/home/Min/basilisk/src/common.h\"\nstatic void pfree (void * ptr,\n     const char * func, const char * file, int line)\n{\n  sysfree (pmfunc_free (ptr, '-'));\n}\n\n\n// #line 1398 \"/home/Min/basilisk/src/common.h\"\nvoid matrix_free (void * m)\n{\n  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);\n  pfree (m,__func__,__FILE__,__LINE__);\n}\n\n\n// #line 1347 \"/home/Min/basilisk/src/common.h\"\nreal matrix_inverse (real ** m, int n, real pivmin)\n{\n  int indxc[n], indxr[n], ipiv[n];\n  int i, icol = 0, irow = 0, j, k, l, ll;\n  real big, dum, pivinv, minpiv = 1e30;\n\n  for (j = 0; j < n; j++)\n    ipiv[j] = -1;\n\n  for (i = 0; i < n; i++) {\n    big = 0.0;\n    for (j = 0; j < n; j++)\n      if (ipiv[j] != 0)\n for (k = 0; k < n; k++) {\n   if (ipiv[k] == -1) {\n     if (fabs (m[j][k]) >= big) {\n       big = fabs (m[j][k]);\n       irow = j;\n       icol = k;\n     }\n   }\n }\n    ipiv[icol]++;\n    if (irow != icol)\n      for (l = 0; l < n; l++)\n do { real __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);\n    indxr[i] = irow;\n    indxc[i] = icol;\n    if (fabs (m[icol][icol]) <= pivmin)\n      return 0.;\n    if (fabs (m[icol][icol]) < minpiv)\n      minpiv = fabs (m[icol][icol]);\n    pivinv = 1.0/m[icol][icol];\n    m[icol][icol] = 1.0;\n    for (l = 0; l < n; l++) m[icol][l] *= pivinv;\n    for (ll = 0; ll < n; ll++)\n      if (ll != icol) {\n dum = m[ll][icol];\n m[ll][icol] = 0.0;\n for (l = 0; l < n; l++)\n   m[ll][l] -= m[icol][l]*dum;\n      }\n  }\n  for (l = n - 1; l >= 0; l--) {\n    if (indxr[l] != indxc[l])\n      for (k = 0; k < n; k++)\n do { real __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);\n  }\n  return minpiv;\n}\n\n\n// #line 113 \"/home/Min/basilisk/src/parabola.h\"\nstatic real parabola_fit_solve (ParabolaFit * p)\n{\n\n  p->M[0][1] = p->M[1][0];\n  p->M[0][2] = p->M[2][0] = p->M[1][1];\n  p->M[1][2] = p->M[2][1];\n  real pivmin = matrix_inverse (p->M, 3, 1e-10);\n  if (pivmin) {\n    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];\n    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];\n  }\n  else\n    p->a[0] = p->a[1] = 0.;\n// # 158 \"/home/Min/basilisk/src/parabola.h\"\n  matrix_free (p->M);\n  return pivmin;\n}\n\n\n// #line 73 \"/home/Min/basilisk/src/parabola.h\"\nstatic void parabola_fit_add (ParabolaFit * p, coord m, real w)\n{\n\n  real x1 = m.x - p->o.x, y1 = m.y - p->o.y;\n  real x = p->m.y*x1 - p->m.x*y1;\n  real y = p->m.x*x1 + p->m.y*y1;\n  real x2 = w*x*x, x3 = x2*x, x4 = x3*x;\n  p->M[0][0] += x4;\n  p->M[1][0] += x3; p->M[1][1] += x2;\n  p->M[2][1] += w*x; p->M[2][2] += w;\n  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;\n// # 111 \"/home/Min/basilisk/src/parabola.h\"\n}\n\n\n// #line 154 \"/home/Min/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if _GNU_SOURCE\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 135 \"/home/Min/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/Min/basilisk/src/common.h\", 180, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 242 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmalloc (size_t size,\n         const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),\n         size, func, file, line, '+');\n}\n\n\n\n\n// #line 1338 \"/home/Min/basilisk/src/common.h\"\nvoid * matrix_new (int n, int p, size_t size)\n{\n  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));\n  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));\n  for (int i = 0; i < n; i++)\n    m[i] = a + i*p*size;\n  return m;\n}\n\nOMP(omp declare reduction (+ : coord :\n      omp_out.x += omp_in.x,\n      omp_out.y += omp_in.y,\n      omp_out.z += omp_in.z))\n// # 920 \"/home/Min/basilisk/src/common.h\"\nvoid normalize (coord * n)\n{\n  real norm = 0.;\n  \n    norm += sq(n->x);\n    \n// #line 924\nnorm += sq(n->y);\n  norm = sqrt(norm);\n  \n    n->x /= norm;\n    \n// #line 927\nn->y /= norm;\n}\n\n\n// #line 23 \"/home/Min/basilisk/src/parabola.h\"\nstatic void parabola_fit_init (ParabolaFit * p, coord o, coord m)\n{\n  \n    p->o.x = o.x;\n    \n// #line 26\np->o.y = o.y;\n\n  \n    p->m.x = m.x;\n    \n// #line 29\np->m.y = m.y;\n  normalize (&p->m);\n  int n = 3;\n// # 65 \"/home/Min/basilisk/src/parabola.h\"\n  p->M = (real **) matrix_new (n, n, sizeof(real));\n  for (int i = 0; i < n; i++) {\n    for (int j = 0; j < n; j++)\n      p->M[i][j] = 0.;\n    p->rhs[i] = 0.;\n  }\n}\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n// # 333 \"/home/Min/basilisk/src/curvature.h\" 2\n\n\n\n\n\n\nstatic int independents (coord * p, int n)\n{\n  if (n < 2)\n    return n;\n  int ni = 1;\n  for (int j = 1; j < n; j++) {\n    bool depends = false;\n    for (int i = 0; i < j && !depends; i++) {\n      real d2 = 0.;\n      \n d2 += sq(p[i].x - p[j].x);\n \n// #line 349\nd2 += sq(p[i].y - p[j].y);\n      depends = (d2 < sq(0.5));\n    }\n    ni += !depends;\n  }\n  return ni;\n}\n// # 18 \"/home/Min/basilisk/src/two-phase.h\" 2\n// # 6 \"reactor.c\" 2\n\n// # 1 \"curvature.h\" 1\n// # 1 \"/home/Min/basilisk/src/curvature.h\"\n// # 68 \"/home/Min/basilisk/src/curvature.h\"\n// # 1 \"heights.h\" 1\n// # 1 \"/home/Min/basilisk/src/heights.h\"\n// # 29 \"/home/Min/basilisk/src/heights.h\"\nstatic inline real height (real H) {\n  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;\n}\n\n\n// #line 33 \"/home/Min/basilisk/src/heights.h\"\nstatic inline int orientation (real H) {\n  return fabs(H) > 20./2.;\n}\n\n\n\n\n\n\n\n// #line 362 \"/home/Min/basilisk/src/curvature.h\"\nstatic real height_curvature_fit (Point point, scalar c, vector h)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n  coord ip[2 == 2 ? 6 : 27];\n  int n = 0;\n\n\n\n\n   {\n\n\n\n\n\n    int n1 = 0, n2 = 0;\n\n    for (int i = -1; i <= 1; i++)\n      if (val(h.y,i,0,0) != 1e30) {\n if (orientation(val(h.y,i,0,0))) n1++; else n2++;\n      }\n\n\n\n\n\n\n\n    int ori = (n1 > n2);\n\n\n\n\n\n\n\n    for (int i = -1; i <= 1; i++)\n      if (val(h.y,i,0,0) != 1e30 && orientation(val(h.y,i,0,0)) == ori)\n ip[n].x = i, ip[n++].y = height(val(h.y,i,0,0));\n\n\n\n\n\n\n  } \n// #line 375\n{\n\n\n\n\n\n    int n1 = 0, n2 = 0;\n\n    for (int i = -1; i <= 1; i++)\n      if (val(h.x,0,i,0) != 1e30) {\n if (orientation(val(h.x,0,i,0))) n1++; else n2++;\n      }\n\n\n\n\n\n\n\n    int ori = (n1 > n2);\n\n\n\n\n\n\n\n    for (int i = -1; i <= 1; i++)\n      if (val(h.x,0,i,0) != 1e30 && orientation(val(h.x,0,i,0)) == ori)\n ip[n].y = i, ip[n++].x = height(val(h.x,0,i,0));\n\n\n\n\n\n\n  }\n\n\n\n\n\n  if (independents (ip, n) < (2 == 2 ? 3 : 9))\n    return 1e30;\n\n\n\n\n\n  coord m = mycs (point, c), fc;\n  real alpha = line_alpha (val(c,0,0,0), m);\n  real area = line_length_center(m,alpha,&fc);\n  ParabolaFit fit;\n  parabola_fit_init (&fit, fc, m);\n\n  NOT_UNUSED(area);\n  parabola_fit_add (&fit, fc, .1);\n// # 440 \"/home/Min/basilisk/src/curvature.h\"\n  for (int i = 0; i < n; i++)\n    parabola_fit_add (&fit, ip[i], 1.);\n  parabola_fit_solve (&fit);\n  real kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;\n\n\n\n  return kappa;\n}\n// # 181 \"/home/Min/basilisk/src/curvature.h\"\nstatic real height_curvature (Point point, scalar c, vector h)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n  typedef struct {\n    real n;\n    real (* kappa) (Point, vector);\n  } NormKappa;\n  struct { NormKappa x, y, z; } n;\n  \n    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;\n    \n// #line 195\nn.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;\n  real (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);\n\n\n\n\n  if (fabs(n.x.n) < fabs(n.y.n))\n    do { NormKappa __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);\n// # 213 \"/home/Min/basilisk/src/curvature.h\"\n  real kappa = 1e30;\n  \n    if (kappa == 1e30) {\n      kappa = n.x.kappa (point, h);\n      if (kappa != 1e30) {\n kappaf = n.x.kappa;\n if (n.x.n < 0.)\n   kappa = - kappa;\n      }\n    }\n    \n// #line 215\nif (kappa == 1e30) {\n      kappa = n.y.kappa (point, h);\n      if (kappa != 1e30) {\n kappaf = n.y.kappa;\n if (n.y.n < 0.)\n   kappa = - kappa;\n      }\n    }\n\n  if (kappa != 1e30) {\n\n\n\n\n    if (fabs(kappa) > 1./Delta)\n      kappa = sign(kappa)/Delta;\n// # 249 \"/home/Min/basilisk/src/curvature.h\"\n  }\n\n  return kappa;\n}\n// # 504 \"/home/Min/basilisk/src/curvature.h\"\nstatic inline bool interfacial (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (val(c,0,0,0) >= 1.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) <= 0.)\n   return true;\n \n// #line 509\nif (val(c,0,i,0) <= 0.)\n   return true;}\n  }\n  else if (val(c,0,0,0) <= 0.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) >= 1.)\n   return true;\n \n// #line 515\nif (val(c,0,i,0) >= 1.)\n   return true;}\n  }\n  else\n    return true;\n  return false;\n}"," \n// #line 568 \"/home/Min/basilisk/src/curvature.h\"\n{\n\n\n\n\n    if (!interfacial (point, c))\n      val_out_(k,0,0,0) = 1e30;\n\n\n\n\n\n    else if ((val_out_(k,0,0,0) = height_curvature (point, c, h)) != 1e30)\n      sh++;\n    else if ((val_out_(k,0,0,0) = height_curvature_fit (point, c, h)) != 1e30)\n      f++;\n  }"}) {




_stencil_interfacial (point, c);{
      {_stencil_val_a(k,0,0,0);  } 





{ _stencil_height_curvature (point, c, NULL);_stencil_val_a(k,0,0,0);{
       
{ _stencil_height_curvature_fit (point, c, h);_stencil_val_a(k,0,0,0);
          }}    
    
#line 583
}}




     





    
  
#line 584
}end_foreach_stencil();

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:f)reduction(+:sh)){
#line 568
foreach() {




    if (!interfacial (point, c))
      val(k,0,0,0) = 1e30;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != 1e30)
      sh++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != 1e30)
      f++;
  }end_foreach();mpi_all_reduce_array(&f,int,MPI_SUM,1);mpi_all_reduce_array(&sh,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 586
foreach_stencil (1,{(NonLocal[]){{"sigma","double",(void *)&sigma,NULL,0},{"add","bool",(void *)&add,NULL,0},{"kappa","scalar",(void *)&kappa,NULL,0},{"sc","int",(void *)&sc,NULL,0,'+'},{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"sa","int",(void *)&sa,NULL,0,'+'},{"c","scalar",(void *)&c,NULL,0},{"k","scalar",(void *)&k,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 162 \"/home/Min/basilisk/src/parabola.h\"\nstatic real parabola_fit_curvature (ParabolaFit * p,\n          real kappamax, real * kmax)\n{\n  real kappa;\n\n  real dnm = 1. + sq(p->a[1]);\n  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);\n  if (kmax)\n    *kmax = fabs (kappa);\n// # 190 \"/home/Min/basilisk/src/parabola.h\"\n  if (fabs (kappa) > kappamax) {\n    if (kmax)\n      *kmax = kappamax;\n    return kappa > 0. ? kappamax : - kappamax;\n  }\n  return kappa;\n}\n\n\n// #line 200 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 264 \"/home/Min/basilisk/src/common.h\"\nstatic void pfree (void * ptr,\n     const char * func, const char * file, int line)\n{\n  sysfree (pmfunc_free (ptr, '-'));\n}\n\n\n// #line 1398 \"/home/Min/basilisk/src/common.h\"\nvoid matrix_free (void * m)\n{\n  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);\n  pfree (m,__func__,__FILE__,__LINE__);\n}\n\n\n// #line 1347 \"/home/Min/basilisk/src/common.h\"\nreal matrix_inverse (real ** m, int n, real pivmin)\n{\n  int indxc[n], indxr[n], ipiv[n];\n  int i, icol = 0, irow = 0, j, k, l, ll;\n  real big, dum, pivinv, minpiv = 1e30;\n\n  for (j = 0; j < n; j++)\n    ipiv[j] = -1;\n\n  for (i = 0; i < n; i++) {\n    big = 0.0;\n    for (j = 0; j < n; j++)\n      if (ipiv[j] != 0)\n for (k = 0; k < n; k++) {\n   if (ipiv[k] == -1) {\n     if (fabs (m[j][k]) >= big) {\n       big = fabs (m[j][k]);\n       irow = j;\n       icol = k;\n     }\n   }\n }\n    ipiv[icol]++;\n    if (irow != icol)\n      for (l = 0; l < n; l++)\n do { real __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);\n    indxr[i] = irow;\n    indxc[i] = icol;\n    if (fabs (m[icol][icol]) <= pivmin)\n      return 0.;\n    if (fabs (m[icol][icol]) < minpiv)\n      minpiv = fabs (m[icol][icol]);\n    pivinv = 1.0/m[icol][icol];\n    m[icol][icol] = 1.0;\n    for (l = 0; l < n; l++) m[icol][l] *= pivinv;\n    for (ll = 0; ll < n; ll++)\n      if (ll != icol) {\n dum = m[ll][icol];\n m[ll][icol] = 0.0;\n for (l = 0; l < n; l++)\n   m[ll][l] -= m[icol][l]*dum;\n      }\n  }\n  for (l = n - 1; l >= 0; l--) {\n    if (indxr[l] != indxc[l])\n      for (k = 0; k < n; k++)\n do { real __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);\n  }\n  return minpiv;\n}\n\n\n// #line 113 \"/home/Min/basilisk/src/parabola.h\"\nstatic real parabola_fit_solve (ParabolaFit * p)\n{\n\n  p->M[0][1] = p->M[1][0];\n  p->M[0][2] = p->M[2][0] = p->M[1][1];\n  p->M[1][2] = p->M[2][1];\n  real pivmin = matrix_inverse (p->M, 3, 1e-10);\n  if (pivmin) {\n    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];\n    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];\n  }\n  else\n    p->a[0] = p->a[1] = 0.;\n// # 158 \"/home/Min/basilisk/src/parabola.h\"\n  matrix_free (p->M);\n  return pivmin;\n}\n\n\n// #line 73 \"/home/Min/basilisk/src/parabola.h\"\nstatic void parabola_fit_add (ParabolaFit * p, coord m, real w)\n{\n\n  real x1 = m.x - p->o.x, y1 = m.y - p->o.y;\n  real x = p->m.y*x1 - p->m.x*y1;\n  real y = p->m.x*x1 + p->m.y*y1;\n  real x2 = w*x*x, x3 = x2*x, x4 = x3*x;\n  p->M[0][0] += x4;\n  p->M[1][0] += x3; p->M[1][1] += x2;\n  p->M[2][1] += w*x; p->M[2][2] += w;\n  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;\n// # 111 \"/home/Min/basilisk/src/parabola.h\"\n}\n\n\n// #line 154 \"/home/Min/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if _GNU_SOURCE\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 135 \"/home/Min/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/Min/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/Min/basilisk/src/common.h\", 180, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 242 \"/home/Min/basilisk/src/common.h\"\nstatic void * pmalloc (size_t size,\n         const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),\n         size, func, file, line, '+');\n}\n\n\n\n\n// #line 1338 \"/home/Min/basilisk/src/common.h\"\nvoid * matrix_new (int n, int p, size_t size)\n{\n  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));\n  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));\n  for (int i = 0; i < n; i++)\n    m[i] = a + i*p*size;\n  return m;\n}\n\nOMP(omp declare reduction (+ : coord :\n      omp_out.x += omp_in.x,\n      omp_out.y += omp_in.y,\n      omp_out.z += omp_in.z))\n// # 920 \"/home/Min/basilisk/src/common.h\"\nvoid normalize (coord * n)\n{\n  real norm = 0.;\n  \n    norm += sq(n->x);\n    \n// #line 924\nnorm += sq(n->y);\n  norm = sqrt(norm);\n  \n    n->x /= norm;\n    \n// #line 927\nn->y /= norm;\n}\n\n\n// #line 23 \"/home/Min/basilisk/src/parabola.h\"\nstatic void parabola_fit_init (ParabolaFit * p, coord o, coord m)\n{\n  \n    p->o.x = o.x;\n    \n// #line 26\np->o.y = o.y;\n\n  \n    p->m.x = m.x;\n    \n// #line 29\np->m.y = m.y;\n  normalize (&p->m);\n  int n = 3;\n// # 65 \"/home/Min/basilisk/src/parabola.h\"\n  p->M = (real **) matrix_new (n, n, sizeof(real));\n  for (int i = 0; i < n; i++) {\n    for (int j = 0; j < n; j++)\n      p->M[i][j] = 0.;\n    p->rhs[i] = 0.;\n  }\n}\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n\n// #line 455 \"/home/Min/basilisk/src/curvature.h\"\nstatic real centroids_curvature_fit (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n  coord m = mycs (point, c), fc;\n  real alpha = line_alpha (val(c,0,0,0), m);\n  line_length_center(m,alpha,&fc);\n  ParabolaFit fit;\n  parabola_fit_init (&fit, fc, m);\n\n\n\n\n\n  coord r = {x,y,z};\n  {foreach_neighbor(1)\n    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {\n      coord m = mycs (point, c), fc;\n      real alpha = line_alpha (val(c,0,0,0), m);\n      real area = line_length_center(m,alpha,&fc);\n      coord rn = {x,y,z};\n      \n fc.x += (rn.x - r.x)/Delta;\n \n// #line 480\nfc.y += (rn.y - r.y)/Delta;\n      parabola_fit_add (&fit, fc, area);\n    }end_foreach_neighbor()}\n  parabola_fit_solve (&fit);\n  real kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;\n\n\n\n  return kappa;\n}\n// # 504 \"/home/Min/basilisk/src/curvature.h\"\nstatic inline bool interfacial (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (val(c,0,0,0) >= 1.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) <= 0.)\n   return true;\n \n// #line 509\nif (val(c,0,i,0) <= 0.)\n   return true;}\n  }\n  else if (val(c,0,0,0) <= 0.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) >= 1.)\n   return true;\n \n// #line 515\nif (val(c,0,i,0) >= 1.)\n   return true;}\n  }\n  else\n    return true;\n  return false;\n}"," \n// #line 586 \"/home/Min/basilisk/src/curvature.h\"\n{\n\n\n\n\n\n    real kf;\n    if (val(k,0,0,0) < 1e30)\n      kf = val(k,0,0,0);\n    else if (interfacial (point, c)) {\n\n\n\n\n\n      real sk = 0., a = 0.;\n      {foreach_neighbor(1)\n if (val(k,0,0,0) < 1e30)\n   sk += val(k,0,0,0), a++;end_foreach_neighbor()}\n      if (a > 0.)\n kf = sk/a, sa++;\n      else\n\n\n\n\n kf = centroids_curvature_fit (point, c), sc++;\n    }\n    else\n      kf = 1e30;\n\n\n\n\n    if (kf == 1e30)\n      val_out_(kappa,0,0,0) = 1e30;\n    else if (add)\n      val_out_(kappa,0,0,0) += sigma*kf;\n    else\n      val_out_(kappa,0,0,0) = sigma*kf;\n  }"}) { 





    
_stencil_val(k,0,0,0);{
      { _stencil_val(k,0,0,0); } 
{_stencil_interfacial (point, c);{ {      





      
      {foreach_neighbor(1)
 {_stencil_val(k,0,0,0);
   { _stencil_val(k,0,0,0);  }   }end_foreach_neighbor()}




 


{ _stencil_centroids_curvature_fit (point, c);  }
         
    
      
    
#line 613
}
        
} 
    
#line 615
}}




{
      {_stencil_val_a(kappa,0,0,0);  } 
if (add)
      {_stencil_val_r(kappa,0,0,0);  }
    else
      {_stencil_val_a(kappa,0,0,0);  }}
       
    




       
    
  
#line 626
}end_foreach_stencil();

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:sc)reduction(+:sa)){
#line 586
foreach () {





    double kf;
    if (val(k,0,0,0) < 1e30)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
      {foreach_neighbor(1)
 if (val(k,0,0,0) < 1e30)
   sk += val(k,0,0,0), a++;end_foreach_neighbor()}
      if (a > 0.)
 kf = sk/a, sa++;
      else




 kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = 1e30;




    if (kf == 1e30)
      val(kappa,0,0,0) = 1e30;
    else if (add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  }end_foreach();mpi_all_reduce_array(&sc,int,MPI_SUM,1);mpi_all_reduce_array(&sa,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 643 "/home/Min/basilisk/src/curvature.h"
  { cstats _ret= (cstats){sh, f, sa, sc};delete((scalar*)((scalar[]){k,{-1}}));if(!(ch).x.i)delete((scalar*)((vector[]){h,{{-1},{-1}}}));{end_tracing("curvature","/home/Min/basilisk/src/curvature.h",643);return _ret;}}delete((scalar*)((scalar[]){k,{-1}}));
end_tracing("curvature","/home/Min/basilisk/src/curvature.h",644);}
#line 665 "/home/Min/basilisk/src/curvature.h"

static double pos_x (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (fabs(height(val(h.x,0,0,0))) > 1.)
    return 1e30;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.x - Z->x)*G->x;
    
#line 674
pos += (o.y - Z->y)*G->y;
  return pos;
}

#line 666
static double pos_y (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (fabs(height(val(h.y,0,0,0))) > 1.)
    return 1e30;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.y - Z->y)*G->y;
    
#line 674
pos += (o.x - Z->x)*G->x;
  return pos;
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  
    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
    
#line 699
n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormPos __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 716 "/home/Min/basilisk/src/curvature.h"
  double pos = 1e30;
  
    if (pos == 1e30)
      pos = n.x.pos (point, h, G, Z);
    
#line 718
if (pos == 1e30)
      pos = n.y.pos (point, h, G, Z);

  return pos;
}








#line 684
static void _stencil_height_position (Point point, scalar f,_stencil_undefined * h,
_stencil_undefined * G,_stencil_undefined * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;        
          
     
      






  
  
  
    { _stencil_val(f,1,0,0); _stencil_val(f,-1,0,0);     }
    
#line 699
{ _stencil_val(f,0,1,0); _stencil_val(f,0,-1,0);     }                
    




     
#line 716 "/home/Min/basilisk/src/curvature.h"
  
     
          
      

  return ;
}
#line 735 "/home/Min/basilisk/src/curvature.h"
void position (scalar f, scalar pos,
        coord G, coord Z, bool add)
{
#line 749 "/home/Min/basilisk/src/curvature.h"
  vector fh = _attribute[f.i].height,   h=(fh).x.i>0?(fh):new_vector("h");
  if (!fh.x.i)
    heights (f, h);
  foreach_stencil(1,{(NonLocal[]){{"pos","scalar",(void *)&pos,NULL,0},{"add","bool",(void *)&add,NULL,0},{"Z","coord",(void *)&Z,NULL,0},{"G","coord",(void *)&G,NULL,0},{"h","vector",(void *)&h,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 382 \"/home/Min/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 2 \"reactor.c\" 2\n// # 1 \"embed.h\" 1\n// # 1 \"/home/Min/basilisk/src/embed.h\"\n// # 12 \"/home/Min/basilisk/src/embed.h\"\n// # 1 \"fractions.h\" 1\n// # 1 \"/home/Min/basilisk/src/fractions.h\"\n// # 12 \"/home/Min/basilisk/src/fractions.h\"\n// # 1 \"geometry.h\" 1\n// # 1 \"/home/Min/basilisk/src/geometry.h\"\n// # 35 \"/home/Min/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/Min/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"myc2d.h\" 1\n// # 1 \"/home/Min/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n\n\n// #line 684 \"/home/Min/basilisk/src/curvature.h\"\nstatic real height_position (Point point, scalar f, vector h,\n          coord * G, coord * Z)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n  typedef struct {\n    real n;\n    real (* pos) (Point, vector, coord *, coord *);\n  } NormPos;\n  struct { NormPos x, y, z; } n;\n  \n    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;\n    \n// #line 699\nn.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;\n\n\n\n\n  if (fabs(n.x.n) < fabs(n.y.n))\n    do { NormPos __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);\n// # 716 \"/home/Min/basilisk/src/curvature.h\"\n  real pos = 1e30;\n  \n    if (pos == 1e30)\n      pos = n.x.pos (point, h, G, Z);\n    \n// #line 718\nif (pos == 1e30)\n      pos = n.y.pos (point, h, G, Z);\n\n  return pos;\n}\n// # 504 \"/home/Min/basilisk/src/curvature.h\"\nstatic inline bool interfacial (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (val(c,0,0,0) >= 1.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) <= 0.)\n   return true;\n \n// #line 509\nif (val(c,0,i,0) <= 0.)\n   return true;}\n  }\n  else if (val(c,0,0,0) <= 0.) {\n    for (int i = -1; i <= 1; i += 2)\n      {\n if (val(c,i,0,0) >= 1.)\n   return true;\n \n// #line 515\nif (val(c,0,i,0) >= 1.)\n   return true;}\n  }\n  else\n    return true;\n  return false;\n}"," \n// #line 752 \"/home/Min/basilisk/src/curvature.h\"\n{\n    if (interfacial (point, f)) {\n      real hp = height_position (point, f, h, &G, &Z);\n      if (hp == 1e30) {\n\n\n\n\n\n coord n = mycs (point, f), o = {x,y,z}, c;\n real alpha = line_alpha (val(f,0,0,0), n);\n line_length_center(n,alpha,&c);\n hp = 0.;\n \n   hp += (o.x + Delta*c.x - Z.x)*G.x;\n   \n// #line 766\nhp += (o.y + Delta*c.y - Z.y)*G.y;\n      }\n      if (add)\n val_out_(pos,0,0,0) += hp;\n      else\n val_out_(pos,0,0,0) = hp;\n    }\n    else\n      val_out_(pos,0,0,0) = 1e30;\n  }"}) {
_stencil_interfacial (point, f);{ {  
       _stencil_height_position (point, f, NULL,NULL ,NULL ); 
{      





  _stencil_mycs (point, f);     
 _stencil_val(f,0,0,0); 
 
  
 
         
      }
         
      
#line 768
if (add)
 {_stencil_val_r(pos,0,0,0);  }
      else
 {_stencil_val_a(pos,0,0,0);  }
    }
      
{_stencil_val_a(pos,0,0,0);  }}
     
    
  
#line 775
}end_foreach_stencil();
  {
#line 752
foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, &G, &Z);
      if (hp == 1e30) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = line_alpha (val(f,0,0,0), n);
 line_length_center(n,alpha,&c);
 hp = 0.;
 
   hp += (o.x + Delta*c.x - Z.x)*G.x;
   
#line 766
hp += (o.y + Delta*c.y - Z.y)*G.y;
      }
      if (add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = 1e30;
  }end_foreach();}if(!(fh).x.i)delete((scalar*)((vector[]){h,{{-1},{-1}}}));
#line 790 "/home/Min/basilisk/src/curvature.h"
}
#line 30 "/home/Min/basilisk/src/iforce.h"
static int defaults_4_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 8 "reactor.c"
#line 1 "tension.h"
#line 1 "/home/Min/basilisk/src/tension.h"
#line 15 "/home/Min/basilisk/src/tension.h"
#line 1 "iforce.h"
#line 1 "/home/Min/basilisk/src/iforce.h"
#line 20 "/home/Min/basilisk/src/iforce.h"










      static int defaults_4(const int i,const double t,Event *_ev){tracing("defaults_4","/home/Min/basilisk/src/iforce.h",30); {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
    if(!is_constant(a.x)){
    
#line 33
foreach_face_stencil(1,{(NonLocal[]){{"DT","double",(void *)&DT,NULL,0},{"a","vector",(void *)&a,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 33 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){ {\n      val_out_(a.x,0,0,0) = 0.;\n      dimensional (val(a.x,0,0,0) == Delta/sq(DT));\n    }}end_is_face_x()\n// #line 33\nis_face_y(){ {\n      val_out_(a.y,0,0,0) = 0.;\n      dimensional (val(a.y,0,0,0) == Delta/sq(DT));\n    }}end_is_face_y()}"}){_stencil_is_face_x(){ {
      _stencil_val_a(a.x,0,0,0);
_stencil_val(a.x,0,0,0);     
      
    
#line 36
}}end__stencil_is_face_x()
#line 33
_stencil_is_face_y(){ {
      _stencil_val_a(a.y,0,0,0);
_stencil_val(a.y,0,0,0);     
      
    
#line 36
}}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 33
foreach_face_generic(){is_face_x(){ {
      val(a.x,0,0,0) = 0.;
      dimensional (val(a.x,0,0,0) == Delta/sq(DT));
    }}end_is_face_x()
#line 33
is_face_y(){ {
      val(a.y,0,0,0) = 0.;
      dimensional (val(a.y,0,0,0) == Delta/sq(DT));
    }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
    
#line 33
foreach_face_stencil(1,{(NonLocal[]){{"DT","double",(void *)&DT,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 33 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){ {\n      _const_a.x = 0.;\n      dimensional (_const_a.x == Delta/sq(DT));\n    }}end_is_face_x()\n// #line 33\nis_face_y(){ {\n      _const_a.y = 0.;\n      dimensional (_const_a.y == Delta/sq(DT));\n    }}end_is_face_y()}"}){_stencil_is_face_x(){ {
      _stencil_val_a(a.x,0,0,0);
;     
      
    
#line 36
}}end__stencil_is_face_x()
#line 33
_stencil_is_face_y(){ {
      _stencil_val_a(a.y,0,0,0);
;     
      
    
#line 36
}}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 33
foreach_face_generic(){is_face_x(){ {
      _const_a.x = 0.;
      dimensional (_const_a.x == Delta/sq(DT));
    }}end_is_face_x()
#line 33
is_face_y(){ {
      _const_a.y = 0.;
      dimensional (_const_a.y == Delta/sq(DT));
    }}end_is_face_y()}end_foreach_face_generic();}}
  }
}{end_tracing("defaults_4","/home/Min/basilisk/src/iforce.h",38);return 0;}end_tracing("defaults_4","/home/Min/basilisk/src/iforce.h",38);}






static int acceleration_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 45
      static int acceleration_0(const int i,const double t,Event *_ev){tracing("acceleration_0","/home/Min/basilisk/src/iforce.h",45);
{





  scalar * list = NULL;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






      foreach_stencil(1,{(NonLocal[]){{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n \n// #line 63 \"/home/Min/basilisk/src/iforce.h\"\nval_out_(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.);"})
 {_stencil_val(f,0,0,0);_stencil_val_a(f,0,0,0);     }end_foreach_stencil();






      {
#line 62
foreach()
 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.);end_foreach();}
    }}}
#line 88 "/home/Min/basilisk/src/iforce.h"
  vector ia = a;
  if(!is_constant(fm.x) && !is_constant(alpha.x)){
  
#line 89
foreach_face_stencil(1,{(NonLocal[]){{"alpha","vector",(void *)&alpha,NULL,0},{"ia","vector",(void *)&ia,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"list","scalar",(void *)list,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 89 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :\n   0.;\n\n val_out_(ia.x,0,0,0) += val(alpha.x,0,0,0)/(val(fm.x,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;\n      } endforin()}}end_is_face_x()\n// #line 89\nis_face_y(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :\n   0.;\n\n val_out_(ia.y,0,0,0) += val(alpha.y,0,0,0)/(val(fm.y,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;\n      } endforin()}}end_is_face_y()}"}){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0); _stencil_val(fm.x,0,0,0); {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0); _stencil_val(fm.y,0,0,0); {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/(val(fm.x,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/(val(fm.y,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 89
foreach_face_stencil(1,{(NonLocal[]){{"alpha","vector",(void *)&alpha,NULL,0},{"ia","vector",(void *)&ia,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"list","scalar",(void *)list,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 89 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :\n   0.;\n\n val_out_(ia.x,0,0,0) += val(alpha.x,0,0,0)/(_const_fm.x + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;\n      } endforin()}}end_is_face_x()\n// #line 89\nis_face_y(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :\n   0.;\n\n val_out_(ia.y,0,0,0) += val(alpha.y,0,0,0)/(_const_fm.y + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;\n      } endforin()}}end_is_face_y()}"}){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);; {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.x,0,0,0);;_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);; {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.y,0,0,0);;_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/(_const_fm.x + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/(_const_fm.y + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 89
foreach_face_stencil(1,{(NonLocal[]){{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"ia","vector",(void *)&ia,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"list","scalar",(void *)list,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 89 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :\n   0.;\n\n val_out_(ia.x,0,0,0) += _const_alpha.x/(val(fm.x,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;\n      } endforin()}}end_is_face_x()\n// #line 89\nis_face_y(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :\n   0.;\n\n val_out_(ia.y,0,0,0) += _const_alpha.y/(val(fm.y,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;\n      } endforin()}}end_is_face_y()}"}){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0); _stencil_val(fm.x,0,0,0); {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0);





;_stencil_val(fm.x,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0); _stencil_val(fm.y,0,0,0); {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0);





;_stencil_val(fm.y,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/(val(fm.x,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/(val(fm.y,0,0,0) + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 89
foreach_face_stencil(1,{(NonLocal[]){{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"ia","vector",(void *)&ia,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"list","scalar",(void *)list,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 89 \"/home/Min/basilisk/src/iforce.h\"\n{is_face_x(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :\n   0.;\n\n val_out_(ia.x,0,0,0) += _const_alpha.x/(_const_fm.x + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;\n      } endforin()}}end_is_face_x()\n// #line 89\nis_face_y(){\n    {forin (scalar, f , list)\n      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {\n// # 101 \"/home/Min/basilisk/src/iforce.h\"\n scalar phi = f.phi;\n real phif =\n   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?\n   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :\n   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :\n   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :\n   0.;\n\n val_out_(ia.y,0,0,0) += _const_alpha.y/(_const_fm.y + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;\n      } endforin()}}end_is_face_y()}"}){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);; {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0);





;;_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);; {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0);





;;_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/(_const_fm.x + 1e-30)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 101 "/home/Min/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/(_const_fm.y + 1e-30)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}
#line 127 "/home/Min/basilisk/src/iforce.h"
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar[]){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }}}
  pfree (list,__func__,__FILE__,__LINE__);
}{end_tracing("acceleration_0","/home/Min/basilisk/src/iforce.h",133);return 0;}end_tracing("acceleration_0","/home/Min/basilisk/src/iforce.h",133);}
#line 36 "/home/Min/basilisk/src/tension.h"
static int stability_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 16 "/home/Min/basilisk/src/tension.h"






#line 36 "/home/Min/basilisk/src/tension.h"
      static int stability_1(const int i,const double t,Event *_ev){tracing("stability_1","/home/Min/basilisk/src/tension.h",36);
{





  double amin = 1e30, amax = -1e30, dmin = 1e30;
  if(!is_constant(fm.x) && !is_constant(alpha.x)){
  
#line 44
foreach_face_stencil (1,{(NonLocal[]){{"dmin","double",(void *)&dmin,NULL,0,'m'},{"amin","double",(void *)&amin,NULL,0,'m'},{"amax","double",(void *)&amax,NULL,0,'M'},{"alpha","vector",(void *)&alpha,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 44 \"/home/Min/basilisk/src/tension.h\"\n{is_face_x(){\n    if (val(fm.x,0,0,0) > 0.) {\n      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) > amax) amax = val(alpha.x,0,0,0)/val(fm.x,0,0,0);\n      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) < amin) amin = val(alpha.x,0,0,0)/val(fm.x,0,0,0);\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_x()\n// #line 44\nis_face_y(){\n    if (val(fm.y,0,0,0) > 0.) {\n      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) > amax) amax = val(alpha.y,0,0,0)/val(fm.y,0,0,0);\n      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) < amin) amin = val(alpha.y,0,0,0)/val(fm.y,0,0,0);\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val(fm.x,0,0,0); {
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {_stencil_val(fm.y,0,0,0); {
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) > amax) amax = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) < amin) amin = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) > amax) amax = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) < amin) amin = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else if(is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 44
foreach_face_stencil (1,{(NonLocal[]){{"dmin","double",(void *)&dmin,NULL,0,'m'},{"amin","double",(void *)&amin,NULL,0,'m'},{"amax","double",(void *)&amax,NULL,0,'M'},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 44 \"/home/Min/basilisk/src/tension.h\"\n{is_face_x(){\n    if (_const_fm.x > 0.) {\n      if (val(alpha.x,0,0,0)/_const_fm.x > amax) amax = val(alpha.x,0,0,0)/_const_fm.x;\n      if (val(alpha.x,0,0,0)/_const_fm.x < amin) amin = val(alpha.x,0,0,0)/_const_fm.x;\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_x()\n// #line 44\nis_face_y(){\n    if (_const_fm.y > 0.) {\n      if (val(alpha.y,0,0,0)/_const_fm.y > amax) amax = val(alpha.y,0,0,0)/_const_fm.y;\n      if (val(alpha.y,0,0,0)/_const_fm.y < amin) amin = val(alpha.y,0,0,0)/_const_fm.y;\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {; {
_stencil_val(alpha.x,0,0,0);; { _stencil_val(alpha.x,0,0,0);; }
_stencil_val(alpha.x,0,0,0);; { _stencil_val(alpha.x,0,0,0);; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {; {
_stencil_val(alpha.y,0,0,0);; { _stencil_val(alpha.y,0,0,0);; }
_stencil_val(alpha.y,0,0,0);; { _stencil_val(alpha.y,0,0,0);; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (val(alpha.x,0,0,0)/_const_fm.x > amax) amax = val(alpha.x,0,0,0)/_const_fm.x;
      if (val(alpha.x,0,0,0)/_const_fm.x < amin) amin = val(alpha.x,0,0,0)/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (val(alpha.y,0,0,0)/_const_fm.y > amax) amax = val(alpha.y,0,0,0)/_const_fm.y;
      if (val(alpha.y,0,0,0)/_const_fm.y < amin) amin = val(alpha.y,0,0,0)/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else if(!is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 44
foreach_face_stencil (1,{(NonLocal[]){{"dmin","double",(void *)&dmin,NULL,0,'m'},{"amin","double",(void *)&amin,NULL,0,'m'},{"amax","double",(void *)&amax,NULL,0,'M'},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 44 \"/home/Min/basilisk/src/tension.h\"\n{is_face_x(){\n    if (val(fm.x,0,0,0) > 0.) {\n      if (_const_alpha.x/val(fm.x,0,0,0) > amax) amax = _const_alpha.x/val(fm.x,0,0,0);\n      if (_const_alpha.x/val(fm.x,0,0,0) < amin) amin = _const_alpha.x/val(fm.x,0,0,0);\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_x()\n// #line 44\nis_face_y(){\n    if (val(fm.y,0,0,0) > 0.) {\n      if (_const_alpha.y/val(fm.y,0,0,0) > amax) amax = _const_alpha.y/val(fm.y,0,0,0);\n      if (_const_alpha.y/val(fm.y,0,0,0) < amin) amin = _const_alpha.y/val(fm.y,0,0,0);\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {_stencil_val(fm.x,0,0,0); {
;_stencil_val(fm.x,0,0,0); {;_stencil_val(fm.x,0,0,0); }
;_stencil_val(fm.x,0,0,0); {;_stencil_val(fm.x,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {_stencil_val(fm.y,0,0,0); {
;_stencil_val(fm.y,0,0,0); {;_stencil_val(fm.y,0,0,0); }
;_stencil_val(fm.y,0,0,0); {;_stencil_val(fm.y,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (_const_alpha.x/val(fm.x,0,0,0) > amax) amax = _const_alpha.x/val(fm.x,0,0,0);
      if (_const_alpha.x/val(fm.x,0,0,0) < amin) amin = _const_alpha.x/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (_const_alpha.y/val(fm.y,0,0,0) > amax) amax = _const_alpha.y/val(fm.y,0,0,0);
      if (_const_alpha.y/val(fm.y,0,0,0) < amin) amin = _const_alpha.y/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 44
foreach_face_stencil (1,{(NonLocal[]){{"dmin","double",(void *)&dmin,NULL,0,'m'},{"amin","double",(void *)&amin,NULL,0,'m'},{"amax","double",(void *)&amax,NULL,0,'M'},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 44 \"/home/Min/basilisk/src/tension.h\"\n{is_face_x(){\n    if (_const_fm.x > 0.) {\n      if (_const_alpha.x/_const_fm.x > amax) amax = _const_alpha.x/_const_fm.x;\n      if (_const_alpha.x/_const_fm.x < amin) amin = _const_alpha.x/_const_fm.x;\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_x()\n// #line 44\nis_face_y(){\n    if (_const_fm.y > 0.) {\n      if (_const_alpha.y/_const_fm.y > amax) amax = _const_alpha.y/_const_fm.y;\n      if (_const_alpha.y/_const_fm.y < amin) amin = _const_alpha.y/_const_fm.y;\n      if (Delta < dmin) dmin = Delta;\n    }}end_is_face_y()}"}){_stencil_is_face_x(){
    {; {
;; {;; }
;; {;; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {; {
;; {;; }
;; {;; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (_const_alpha.x/_const_fm.x > amax) amax = _const_alpha.x/_const_fm.x;
      if (_const_alpha.x/_const_fm.x < amin) amin = _const_alpha.x/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (_const_alpha.y/_const_fm.y > amax) amax = _const_alpha.y/_const_fm.y;
      if (_const_alpha.y/_const_fm.y < amin) amin = _const_alpha.y/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}
  double rhom = (1./amin + 1./amax)/2.;





  double sigma = 0.;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){
    sigma += _attribute[c.i].sigma;}}
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(3.14159265358979*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}{end_tracing("stability_1","/home/Min/basilisk/src/tension.h",64);return 0;}end_tracing("stability_1","/home/Min/basilisk/src/tension.h",64);}







static int acceleration_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 72
      static int acceleration_1(const int i,const double t,Event *_ev){tracing("acceleration_1","/home/Min/basilisk/src/tension.h",72);
{




  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature (f, phi, _attribute[f.i].sigma, true);
      else {
 phi = new_scalar("phi");
 curvature (f, phi, _attribute[f.i].sigma, false);
 _attribute[f.i].phi = phi;
      }
    }}}
}{end_tracing("acceleration_1","/home/Min/basilisk/src/tension.h",94);return 0;}end_tracing("acceleration_1","/home/Min/basilisk/src/tension.h",94);}
#line 9 "reactor.c"

#line 1 "reduced.h"
#line 1 "/home/Min/basilisk/src/reduced.h"
#line 19 "/home/Min/basilisk/src/reduced.h"
coord G = {0.,0.,0.}, Z = {0.,0.,0.};
#line 36
static int acceleration_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 36 "/home/Min/basilisk/src/reduced.h"
      static int acceleration_2(const int i,const double t,Event *_ev){tracing("acceleration_2","/home/Min/basilisk/src/reduced.h",36);
{
  scalar phi = _attribute[f.i].phi;
  coord G1;
  
    G1.x = (rho2 - rho1)*G.x;
    
#line 41
G1.y = (rho2 - rho1)*G.y;

  if (phi.i)
    position (f, phi, G1, Z, true);
  else {
    phi = new_scalar("phi");
    position (f, phi, G1, Z, false);
    _attribute[f.i].phi = phi;
  }
}{end_tracing("acceleration_2","/home/Min/basilisk/src/reduced.h",50);return 0;}end_tracing("acceleration_2","/home/Min/basilisk/src/reduced.h",50);}
#line 11 "reactor.c"
#line 1 "tracer.h"
#line 1 "/home/Min/basilisk/src/tracer.h"
#line 15 "/home/Min/basilisk/src/tracer.h"
extern scalar * tracers;
extern vector uf;
extern double dt;
#line 43
static int tracer_advection_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 41 "/home/Min/basilisk/src/tracer.h"
#line 1 "bcg.h"
#line 42 "/home/Min/basilisk/src/tracer.h"

      static int tracer_advection_1(const int i,const double t,Event *_ev){tracing("tracer_advection_1","/home/Min/basilisk/src/tracer.h",43); {
  advection (tracers, uf, dt
#line 67 "/home/Min/basilisk/src/bcg.h"
, 
NULL
#line 44 "/home/Min/basilisk/src/tracer.h"
);
}{end_tracing("tracer_advection_1","/home/Min/basilisk/src/tracer.h",45);return 0;}end_tracing("tracer_advection_1","/home/Min/basilisk/src/tracer.h",45);}




static int tracer_diffusion_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}





#line 50
static int tracer_diffusion_0(const int i,const double t,Event *_ev){;return 0;}
#line 12 "reactor.c"

#line 1 "/home/Min/basilisk/src/ast/std/math.h"
#include <math.h>
#line 14 "reactor.c"
#line 30 "reactor.c"
double We;


double Re;


double F;


scalar  T={15};
scalar * tracers =((scalar[]) {{15},{-1}});


double sigma_sq_init=0;


double ka = 0;
double sum_sq = 0;
double sum_mean = 0;
double sigma_sq = 0;



vector  hei={{16},{17}};
double interfacial_height(double xp) {
  double dh[30];
  double yh[30];
  double yp;
  double y0 = -0.4, y1 = 0.4;


  for (int i = 0; i < 30; i++) {
    yp = y0 + i*(y1-y0)/(30 -1);
    Point point = locate(xp,yp
#line 806 "/home/Min/basilisk/src/grid/multigrid.h"
, 0
#line 63 "reactor.c"
);int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

    if (val(hei.y,0,0,0) != 1e30) {
      yh[i] = y + height(val(hei.y,0,0,0))*Delta;
      dh[i] = fabs(y-yh[i]);
      return yh[i];
    } else {
      yh[i] = -1000;
      dh[i] = 1000;
    }
  }


  int j = 0;
  for (int i = 1; i < 30; i++) {
    if (dh[i] < dh[j]) { j = i; }
  }

  return yh[j];
}

int main()
{
#line 103
_init_solver();
  
#line 86
size (1);
  DT = 1.;
  origin (-0.5,-0.5
#line 930 "/home/Min/basilisk/src/common.h"
, 0.
#line 88 "reactor.c"
);
  TOLERANCE=1e-5;

  Re = 1000 * 0.5 * 0.05 * 0.05 / 0.003;
  We = 1000 * 0.05 * 0.5 * 0.05 * 0.5 / 0.07;
  F = 0.05 * 0.5 * 0.5 / 9.81;
  N=256;


  mu1 = 1./Re;
  mu2 = 0.0033*mu1;
  _attribute[f.i].sigma = 1./We;

  G.y = 1./F;
  run();
free_solver();

#line 103
}

scalar  un={18};
#line 124
static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary7_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
#line 108
static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t=0)!=0;*ip=i;*tp=t;return ret;}



#line 108
      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","reactor.c",108); {


  do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 111 \"reactor.c\"\nval_out_(phi,0,0,0) = -y + 0.;"}) {_stencil_val_a(phi,0,0,0);    }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = -y + 0.;end_foreach_vertex();} fractions (phi, f
#line 121 "/home/Min/basilisk/src/fractions.h"
,
(
  
#line 122
vector) {0}, 0.
#line 111 "reactor.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);
#line 120 "reactor.c"
  do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 120 \"reactor.c\"\nval_out_(phi,0,0,0) = (pow(x-0.19,2)+pow((y+0.25),2)) < pow(0.1,2);"}) {_stencil_val_a(phi,0,0,0);    }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = (pow(x-0.19,2)+pow((y+0.25),2)) < pow(0.1,2);end_foreach_vertex();} fractions (phi, T
#line 121 "/home/Min/basilisk/src/fractions.h"
,
(
  
#line 122
vector) {0}, 0.
#line 120 "reactor.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);
  do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 121 \"reactor.c\"\nval_out_(phi,0,0,0) = (pow(x+0.19,2)+pow((y+0.25),2)) < pow(0.1,2);"}) {_stencil_val_a(phi,0,0,0);    }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = (pow(x+0.19,2)+pow((y+0.25),2)) < pow(0.1,2);end_foreach_vertex();} fractions (phi, T
#line 121 "/home/Min/basilisk/src/fractions.h"
,
(
  
#line 122
vector) {0}, 0.
#line 121 "reactor.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);


_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[top]=_boundary4,_attribute[u.y.i].boundary_homogeneous[top]=_boundary4_homogeneous;
_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[top]=_boundary5,_attribute[u.x.i].boundary_homogeneous[top]=_boundary5_homogeneous;
_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[embed]=_boundary6,_attribute[u.x.i].boundary_homogeneous[embed]=_boundary6_homogeneous;
_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[embed]=_boundary7,_attribute[u.y.i].boundary_homogeneous[embed]=_boundary7_homogeneous;
  foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"un","scalar",(void *)&un,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 129 \"reactor.c\"\nval_out_(un,0,0,0)=val(u.x,0,0,0);"})
    {_stencil_val(u.x,0,0,0);_stencil_val_a(un,0,0,0);}end_foreach_stencil();
  {
#line 128
foreach()
    val(un,0,0,0)=val(u.x,0,0,0);end_foreach();}


  rho1 = 1.;
  rho2 = rho1*0.001;

}{end_tracing("init_0","reactor.c",135);return 0;}end_tracing("init_0","reactor.c",135);}







static int update_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( t<=200.)!=0;*ip=i;*tp=t;return ret;}static int update_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 143
      static int update(const int i,const double t,Event *_ev){tracing("update","reactor.c",143); {

  double theta_max = (1./3.)*3.14159/4.0;

  double theta = 0;

  double t_prime = 0;
  t_prime = fmod(t,1.0);

  double t_max_prime = 1.;


  double t_1_prime = t_max_prime/4.;


  double t_2_prime = 3.*t_max_prime/4.;


  if(t_prime<=t_1_prime){
    theta = t_prime * theta_max / t_1_prime;
    } else if(t_prime<=t_2_prime){
      theta = theta_max - (4.*(t_prime - t_1_prime)/t_max_prime)*theta_max;
      } else {
        theta = -theta_max + (4.*(t_prime - t_2_prime)/t_max_prime)*theta_max;
      }
  G.x = (1./F) * sin(theta);
  G.y = -(1./F) * cos(theta);


}{end_tracing("update","reactor.c",172);return 0;}end_tracing("update","reactor.c",172);}







static int intermediate_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t+=0.25)!=0;*ip=i;*tp=t;return ret;}








#line 180
      static int intermediate(const int i,const double t,Event *_ev){tracing("intermediate","reactor.c",180); {
  heights(f,hei);
  int j, a, b;
  double u_x = 0;
  double y_step = 0.375;
  double x_step = 0;
  double x_pos = 0;
  double T_out = 0;
  double f_out = 0;
  char str[80];





  sum_sq = 0;
  sum_mean = 0;


  sprintf(str, "mix%.3f.dat", t);

  FILE * fp_3 = fopen(str,"w");

  for (b = 0; b<101; b++){
    for (a=0; a<101; a++) {
      T_out = interpolate(T,-0.5+a*0.01,-0.5+b*0.01
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0., 
true
#line 205 "reactor.c"
);
      f_out = interpolate(f,-0.5+a*0.01,-0.5+b*0.01
#line 837 "/home/Min/basilisk/src/grid/cartesian-common.h"
, 0., 
true
#line 206 "reactor.c"
);
      fprintf(fp_3, "%lf %lf %lf %lf \n", -0.5+a*0.01, -0.5+b*0.01, T_out, f_out);
    }
  }
  fclose(fp_3);
#line 224 "reactor.c"
}{end_tracing("intermediate","reactor.c",224);return 0;}end_tracing("intermediate","reactor.c",224);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
datasize=19*sizeof(double);
  
#line 6
init_solver();
  {
#line 24
multigrid_methods();

    

    
#line 12
{
      
      
    
      {  
#line 895 "/home/Min/basilisk/src/embed.h"
event_register((Event){0,1,metric,{metric_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/embed.h",895,"metric"});  
#line 938
event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/embed.h",938,"defaults"});  
#line 42 "/home/Min/basilisk/src/run.h"
event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/run.h",42,"defaults"});  
#line 126 "/home/Min/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,defaults_1,{defaults_1_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",126,"defaults"});  
#line 187
event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",187,"default_display"});  








event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",196,"init"});  
#line 127 "/home/Min/basilisk/src/vof.h"
event_register((Event){0,1,defaults_2,{defaults_2_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/vof.h",127,"defaults"});  
#line 10 "/home/Min/basilisk/src/two-phase-generic.h"
event_register((Event){0,1,defaults_3,{defaults_3_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/two-phase-generic.h",10,"defaults"});  
#line 30 "/home/Min/basilisk/src/iforce.h"
event_register((Event){0,1,defaults_4,{defaults_4_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/iforce.h",30,"defaults"});  
#line 108 "reactor.c"
event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"reactor.c",108,"init"});  
#line 143
event_register((Event){0,2,update,{update_expr0,update_expr1},((int *)0),((double *)0),"reactor.c",143,"update"});  
#line 180
event_register((Event){0,1,intermediate,{intermediate_expr0},((int *)0),((double *)0),"reactor.c",180,"intermediate"});
	
	
	
      
#line 22 "ast/init_solver.h"
}
#line 1257 "/home/Min/basilisk/src/common.h"
init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1}},"zerof",(double[]){0.,0.,0.});
init_const_vector((vector){{_NVARMAX+2},{_NVARMAX+3}},"unityf",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+4},"unity", 1.);
init_const_scalar((scalar){_NVARMAX+5},"zeroc", 0.);



init_const_vector((vector){{_NVARMAX+6},{_NVARMAX+7}},"unityf0",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+8},"unity0", 1.);  init_scalar((scalar){0},"cs");  init_face_vector((vector){{1},{2}},"fs"); 
#line 461 "/home/Min/basilisk/src/embed.h"
embed=new_bid();  init_scalar((scalar){3},"p");  init_vector((vector){{4},{5}},"u");  init_vector((vector){{6},{7}},"g");  init_scalar((scalar){8},"pf");  init_face_vector((vector){{9},{10}},"uf");  init_scalar((scalar){11},"f");  init_face_vector((vector){{12},{13}},"alphav");  init_scalar((scalar){14},"rhov");  init_scalar((scalar){15},"T");  init_vector((vector){{16},{17}},"hei");  init_scalar((scalar){18},"un");
    
#line 23 "ast/init_solver.h"
}_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;  
#line 50 "/home/Min/basilisk/src/run.h"
event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/run.h",50,"cleanup"});  
#line 222 "/home/Min/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",222,"set_dtmax"});  

event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",224,"stability"});  









event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",234,"vof"});  
event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",235,"tracer_advection"});  
event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",236,"tracer_diffusion"});  






event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",243,"properties"});  
#line 316
event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",316,"advection_term"});  
#line 345
event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",345,"viscous_term"});  
#line 381
event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",381,"acceleration"});  
#line 421
event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",421,"projection"});  
#line 436
event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/navier-stokes/centered.h",436,"end_timestep"});  
#line 140 "/home/Min/basilisk/src/vof.h"
event_register((Event){0,1,stability_0,{stability_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/vof.h",140,"stability"});  
#line 380
event_register((Event){0,1,vof_0,{vof_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/vof.h",380,"vof"});  
#line 50 "/home/Min/basilisk/src/two-phase-generic.h"
event_register((Event){0,1,tracer_advection_0,{tracer_advection_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/two-phase-generic.h",50,"tracer_advection"});  
#line 83
event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/two-phase-generic.h",83,"properties"});  
#line 45 "/home/Min/basilisk/src/iforce.h"
event_register((Event){0,1,acceleration_0,{acceleration_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/iforce.h",45,"acceleration"});  
#line 36 "/home/Min/basilisk/src/tension.h"
event_register((Event){0,1,stability_1,{stability_1_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/tension.h",36,"stability"});  
#line 72
event_register((Event){0,1,acceleration_1,{acceleration_1_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/tension.h",72,"acceleration"});  
#line 36 "/home/Min/basilisk/src/reduced.h"
event_register((Event){0,1,acceleration_2,{acceleration_2_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/reduced.h",36,"acceleration"});  
#line 43 "/home/Min/basilisk/src/tracer.h"
event_register((Event){0,1,tracer_advection_1,{tracer_advection_1_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/tracer.h",43,"tracer_advection"});  






event_register((Event){0,1,tracer_diffusion_0,{tracer_diffusion_0_expr0},((int *)0),((double *)0),"/home/Min/basilisk/src/tracer.h",50,"tracer_diffusion"});
  
#line 24 "ast/init_solver.h"
}
  set_fpe();
}
