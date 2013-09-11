#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

/* Numerical values for bases. */
#define MASKED_BASE_BIT 8
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */
/* Some other type synonyms */
#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */

#define boolean int

/* How big is this array? */
#define ArraySize(a) (sizeof(a)/sizeof((a)[0]))

#ifndef min
#define min(a,b) ( (a) < (b) ? (a) : (b) )
/* Return min of a and b. */
#endif

/* inline functions: To declare a function inline, place the entire function
 * in a header file and prefix it with the INLINE macro.  If used with a
 * compiler that doesn't support inline, change the INLINE marco to be simply
 * `static'.
 */
#ifndef INLINE
#define INLINE static inline
#endif

#define NEEDMEM_LIMIT 500000000


void *needMem(size_t size)
/* Need mem calls abort if the memory allocation fails. The memory
 * is initialized to zero. */
{
  void *pt;
  memset(pt, 0, size);
  return pt;
}


#define AllocVar(pt) (pt = needMem(sizeof(*pt)))
/* Shortcut to allocating a single variable on the heap and
 * assigning pointer to it. */

#define AllocArray(pt, size) (pt = needLargeZeroedMem(sizeof(*pt) * (size)))


INLINE void zeroBytes(void *vpt, int count)
/* fill a specified area of memory with zeroes */
{
  memset(vpt, '\0', count);
}

#define ZeroVar(v) zeroBytes(v, sizeof(*v))

typedef char DNA;

int ntVal[256];
int ntValLower[256];	/* NT values only for lower case. */
int ntValUpper[256];	/* NT values only for upper case. */
int ntVal5[256];
int ntValNoN[256]; /* Like ntVal, but with T_BASE_VAL in place of -1 for nonexistent ones. */
DNA valToNt[(N_BASE_VAL|MASKED_BASE_BIT)+1];

/* convert tables for bit-4 indicating masked */
int ntValMasked[256];
DNA valToNtMasked[256];

struct twoBit
/* Two bit representation of DNA. */
{
    struct twoBit *next;	/* Next sequence in list */
    char *name;			/* Name of sequence. */
    UBYTE *data;		/* DNA at two bits per base. */
    bits32 size;		/* Size of this sequence. */
    bits32 nBlockCount;		/* Count of blocks of Ns. */
    bits32 *nStarts;		/* Starts of blocks of Ns. */
    bits32 *nSizes;		/* Sizes of blocks of Ns. */
    bits32 maskBlockCount;	/* Count of masked blocks. */
    bits32 *maskStarts;		/* Starts of masked regions. */
    bits32 *maskSizes;		/* Sizes of masked regions. */
    bits32 reserved;		/* Reserved for future expansion. */
};

struct dnaSeq
/* A dna sequence in one-character per base format. */
{
    struct dnaSeq *next;  /* Next in list. */
    char *name;           /* Name of sequence. */
    DNA *dna;             /* Sequence base by base. */
    int size;             /* Size of sequence. */
    //Bits* mask;           /* Repeat mask (optional) */
};

boolean inittedNtVal = false;
static void initNtVal()
{
        if (!inittedNtVal)
        {
                int i;
                for (i=0; i<ArraySize(ntVal); i++)
                {
                        ntValUpper[i] = ntValLower[i] = ntVal[i] = -1;
                        ntValNoN[i] = T_BASE_VAL;
                        if (isspace(i) || isdigit(i))
                                ntVal5[i] = ntValMasked[i] = -1;
                        else
                        {
                                ntVal5[i] = N_BASE_VAL;
                                ntValMasked[i] = (islower(i) ? (N_BASE_VAL|MASKED_BASE_BIT) : N_BASE_VAL);
                        }
                }
                ntVal5['t'] = ntVal5['T'] = ntValNoN['t'] = ntValNoN['T'] = ntVal['t'] = ntVal['T'] = 
                ntValLower['t'] = ntValUpper['T'] = T_BASE_VAL;
                ntVal5['u'] = ntVal5['U'] = ntValNoN['u'] = ntValNoN['U'] = ntVal['u'] = ntVal['U'] = 
                ntValLower['u'] = ntValUpper['U'] = U_BASE_VAL;
                ntVal5['c'] = ntVal5['C'] = ntValNoN['c'] = ntValNoN['C'] = ntVal['c'] = ntVal['C'] = 
                ntValLower['c'] = ntValUpper['C'] = C_BASE_VAL;
                ntVal5['a'] = ntVal5['A'] = ntValNoN['a'] = ntValNoN['A'] = ntVal['a'] = ntVal['A'] = 
                ntValLower['a'] = ntValUpper['A'] = A_BASE_VAL;
                ntVal5['g'] = ntVal5['G'] = ntValNoN['g'] = ntValNoN['G'] = ntVal['g'] = ntVal['G'] = 
                ntValLower['g'] = ntValUpper['G'] = G_BASE_VAL;

                valToNt[T_BASE_VAL] = valToNt[T_BASE_VAL|MASKED_BASE_BIT] = 't';
                valToNt[C_BASE_VAL] = valToNt[C_BASE_VAL|MASKED_BASE_BIT] = 'c';
                valToNt[A_BASE_VAL] = valToNt[A_BASE_VAL|MASKED_BASE_BIT] = 'a';
                valToNt[G_BASE_VAL] = valToNt[G_BASE_VAL|MASKED_BASE_BIT] = 'g';
                valToNt[N_BASE_VAL] = valToNt[N_BASE_VAL|MASKED_BASE_BIT] = 'n';

                /* masked values */
                ntValMasked['T'] = T_BASE_VAL;
                ntValMasked['U'] = U_BASE_VAL;
                ntValMasked['C'] = C_BASE_VAL;
                ntValMasked['A'] = A_BASE_VAL;
                ntValMasked['G'] = G_BASE_VAL;

                ntValMasked['t'] = T_BASE_VAL|MASKED_BASE_BIT;
                ntValMasked['u'] = U_BASE_VAL|MASKED_BASE_BIT;
                ntValMasked['c'] = C_BASE_VAL|MASKED_BASE_BIT;
                ntValMasked['a'] = A_BASE_VAL|MASKED_BASE_BIT;
                ntValMasked['g'] = G_BASE_VAL|MASKED_BASE_BIT;

                valToNtMasked[T_BASE_VAL] = 'T';
                valToNtMasked[C_BASE_VAL] = 'C';
                valToNtMasked[A_BASE_VAL] = 'A';
                valToNtMasked[G_BASE_VAL] = 'G';
                valToNtMasked[N_BASE_VAL] = 'N';

                valToNtMasked[T_BASE_VAL|MASKED_BASE_BIT] = 't';
                valToNtMasked[C_BASE_VAL|MASKED_BASE_BIT] = 'c';
                valToNtMasked[A_BASE_VAL|MASKED_BASE_BIT] = 'a';
                valToNtMasked[G_BASE_VAL|MASKED_BASE_BIT] = 'g';
                valToNtMasked[N_BASE_VAL|MASKED_BASE_BIT] = 'n';

                inittedNtVal = true;
        }
}


/* 128*8*1024*1024 == 1073741824 == 2^30 on 32 bit machines,size_t == 4 bytes*/
/* on 64 bit machines, size_t = 8 bytes, 2^30 * 2 * 2 * 2 * 2 = 2^34 == 16 Gb */
static size_t maxAlloc = (size_t)128*8*1024*1024*(sizeof(size_t)/4)*(sizeof(size_t)/4)*(sizeof(size_t)/4*(sizeof(size_t)/4));
void *needLargeMem(size_t size)
/* This calls abort if the memory allocation fails. The memory is
 * not initialized to zero. */
{
  void *pt;
  return pt;
}

static char *cloneStringZExt(const char *s, int size, int copySize)
/* Make a zero terminated copy of string in memory */
{
char *d = (char*)needMem(copySize+1);
copySize = min(size,copySize);
memcpy(d, s, copySize);
d[copySize] = 0;
return d;
}

char *cloneStringZ(const char *s, int size)
/* Make a zero terminated copy of string in memory */
{
return cloneStringZExt(s, strlen(s), size);
}

char *cloneString(const char *s)
/* Make copy of string in dynamic memory */
{
        int size = 0;
        if (s == NULL)
                return NULL;
        size = strlen(s);
        return cloneStringZExt(s, size, size);
}


void *needLargeZeroedMem(size_t size)
/* Request a large block of memory and zero it. */
{
        void *v;
        v = needLargeMem(size);
        memset(v, 0, size);
        return v;
}


static int packedSize(int unpackedSize)
/* Return size when packed, rounding up. */
{
return ((unpackedSize + 3) >> 2);
}


UBYTE packDna4(DNA *in)
/* Pack 4 bases into a UBYTE */
{
  UBYTE out = 0;
  int count = 4;
  int bVal;
  while (--count >= 0) {
    bVal = ntValNoN[(int)*in++];
    out <<= 2;
    out += bVal;
  }
  //printf("\n");
  //printf("%d\n",out);
  return out;
}

struct twoBit *twoBitFromDnaSeq(struct dnaSeq *seq, boolean doMask)
/* Convert dnaSeq representation in memory to twoBit representation.
 * If doMask is true interpret lower-case letters as masked. */
{
        int ubyteSize = packedSize(seq->size);
        UBYTE *pt;
        struct twoBit *tB;
        DNA last4[4];	/* Holds few bases. */
        DNA *dna;
        int i, end;

        /* Allocate structure and fill in name. */
        //AllocVar(twoBit);
        tB = (twoBit*)malloc(sizeof(*tB));
        //pt = AllocArray(twoBit->data, ubyteSize);
        tB->data = (UBYTE*)malloc(sizeof(*(tB->data)) * ubyteSize);//needLargeZeroedMem(sizeof(*pt) * (ubyteSize));
        pt = tB->data;
        //twoBit->name = cloneString(seq->name);
        //twoBit->size = seq->size;

        /* Convert to 4-bases per byte representation. */
        dna = seq->dna;
        end = seq->size - 4;
        for (i=0; i<end; i += 4)
        {
                *pt++ = packDna4(dna+i);
                
        }
        
        /* Take care of conversion of last few bases. */
        last4[0] = last4[1] = last4[2] = last4[3] = 'T';
        memcpy(last4, dna+i, seq->size-i);
        *pt = packDna4(last4);
	
        //printf("%d\n",(*pt));
        
        /* Deal with blocks of N. */
        /*twoBit->nBlockCount = countBlocksOfN(dna, seq->size);
        if (twoBit->nBlockCount > 0)
        {
                AllocArray(twoBit->nStarts, twoBit->nBlockCount);
                AllocArray(twoBit->nSizes, twoBit->nBlockCount);
                storeBlocksOfN(dna, seq->size, twoBit->nStarts, twoBit->nSizes);
        }
        */
        return tB;
}


int main(int argc, char *argv[]) {
  struct twoBit *twoBitList = NULL, *twoBit;
  int i;
  struct dnaSeq seq;
  ZeroVar(&seq);
  boolean noMask = true;
  
  seq.dna = "AAAAAAAACCCCTTTTGGGG"; seq.size = 20; seq.name = "Test1";
  initNtVal();
  twoBit = twoBitFromDnaSeq(&seq, !noMask);
  printf("%d\n",*(twoBit->data+2));
}
