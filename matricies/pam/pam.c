/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/** @file pam.c
 *  Compute log-odds matrix for specified integral PAM distance.
 *
 * @author Stephen F. Altshul, E. Michael Gertz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#define ALPHSIZE 20

/** amino acid alphabet */
static const char alphabet[ALPHSIZE] = "ARNDCQEGHILKMFPSTWYV";


/** Sum the elements in a vector. */
static double 
s_SumVector(double v[], int n)
{
    int i;
    double sum = 0.0;
    for (i = 0;  i < n;  i++)
        sum += v[i];
    return sum;
}

/** Scale the elements of a vector. */
static void
s_ScaleVector(double v[], int n, double alpha)
{
    int i;
    for (i = 0;  i < n;  i++)
        v[i] *= alpha;
}

/** Take the dot product of two vectors. */
static double
s_DotProductVectors(double v[], int n, double w[])
{
    int i;
    double dot = 0.0;
    for (i = 0;  i < n;  i++)
        dot += v[i] * w[i];
    return dot;
}


/** Normalize a vector so it sums to 1 (unless it is the zero vector). */
static void
s_NormalizeVector(double v[], int n)
{
    int i;
    double sum = s_SumVector(v, n);
    if (sum != 0.0) {
        for (i = 0;  i < n;  i++) {
            v[i] /= sum;
        }
    }
}


/** Multiply two matrices of size 20: C = A * B. */
static void
s_MultiplyTrueAaMatrices(double C[][ALPHSIZE],
                         double A[][ALPHSIZE],
                         double B[][ALPHSIZE])
{
    int i, j, k;

    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < ALPHSIZE;  j++) {
            C[i][j] = 0.0;
            for (k = 0;  k < ALPHSIZE;  k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


/** Symmetrize the joint probabilities and normalize them so that they
   sum to one. */
static void
s_NormalizeJointProbs(double jnt[][ALPHSIZE])
{
    int i, j;
    double sum = 0.0;

    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < i;  j++) {
            jnt[i][j] = (jnt[i][j] + jnt[j][i])/2.0;
            jnt[j][i] = jnt[i][j];
        }
    }
    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < ALPHSIZE;  j++) {
            sum += jnt[i][j];
        }
    }
    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < ALPHSIZE;  j++) {
            jnt[i][j] /= sum;
        }
    }
}


/** Compute pam = A^n, where A has size 20. */
static void
s_PowerTrueAaMatrix(double pam[][ALPHSIZE],
                    double A[][ALPHSIZE],
                    int n)
{
    int pam_is_identity = 1;
    double work[ALPHSIZE][ALPHSIZE];
    double powA[ALPHSIZE][ALPHSIZE];

    const size_t size_true_aa_matrix =
        sizeof(double[ALPHSIZE][ALPHSIZE]);

    memcpy(powA, A, size_true_aa_matrix);

    /* We don't handle the n == 0 case, since it would require extra
     * code and no one wants PAM0 */
    assert(n > 0);

    while (n > 0) {
        if (n % 2) {
            if (pam_is_identity) {
                memcpy(pam, powA, size_true_aa_matrix);
                pam_is_identity = 0;
            } else {
                s_MultiplyTrueAaMatrices(work, pam, powA);
                memcpy(pam, work, size_true_aa_matrix);
            }
        }
        n /= 2;
        if (n != 0) {
            /* Square powA */
            s_MultiplyTrueAaMatrices(work, powA, powA);
            memcpy(powA, work, size_true_aa_matrix);
        }
    }
}


/*  Dayhoff mutability data         */
static const double
mutability[ALPHSIZE] = {
    100.0, 65.0, 134.0, 106.0, 20.0, 93.0, 102.0, 49.0, 66.0, 96.0,
    40.0,  56.0, 94.0,  41.0,  56.0, 120.0, 97.0, 18.0, 41.0, 74.0
};

/*  Dayhoff mutation data           */
static const int
mutation[ALPHSIZE * (ALPHSIZE - 1) / 2] = {
 30,
109, 17,
154,  0,532,
 33, 10,  0,  0,
 93,120, 50, 76,  0,
266,  0, 94,831,  0,422,
579, 10,156,162, 10, 30,112,
 21,103,226, 43, 10,243, 23, 10,
 66, 30, 36, 13, 17,  8, 35,  0,  3,
 95, 17, 37,  0,  0, 75, 15, 17, 40,253,
 57,477,322, 85,  0,147,104, 60, 23, 43, 39,
 29, 17,  0,  0,  0, 20,  7,  7,  0, 57,207, 90,
 20,  7,  7,  0,  0,  0,  0, 17, 20, 90,167,  0, 17,
345, 67, 27, 10, 10, 93, 40, 49, 50,  7, 43, 43,  4,  7,
772,137,432, 98,117, 47, 86,450, 26, 20, 32,168, 20, 40,269,
590, 20,169, 57, 10, 37, 31, 50, 14,129, 52,200, 28, 10, 73,696,
  0, 27,  3,  0,  0,  0,  0,  0,  3,  0, 13,  0,  0, 10,  0, 17, 0,
 20,  3, 36,  0, 30,  0, 10,  0, 40, 13, 23, 10,  0,260,  0, 22, 23,  6,
 365, 20, 13, 17, 33, 27, 37, 97, 30,661,303, 17, 77, 10, 50, 43, 186, 0, 17};


/** Calculate the transition matrix for PAM. */
static void
s_CalculatePam1(double pam1[][ALPHSIZE], double bgnd_probs[])
{
    int i, j, k;
    double mutab_divisor;
    double mutab[ALPHSIZE];

    /* Fill the off diagonals of the matrix, setting diagonals to zero */
    k = 0;
    for (i = 0; i < ALPHSIZE; ++i) {
        for(j = 0;  j < i;  ++j) {
            pam1[i][j] = pam1[j][i] = mutation[k++];
        }
        pam1[i][i] = 0.0;
    }
    /* Calculate the background frequencies */
    for (i = 0;  i < ALPHSIZE;  ++i) {
        double sum = 0.0;
        for (j = 0;  j < ALPHSIZE;  ++j) {
            sum += pam1[i][j];
        }
        bgnd_probs[i] = sum / mutability[i];
    }
    s_NormalizeVector(bgnd_probs, ALPHSIZE);

    /* Normalize the mutabilities */
    memcpy(mutab, mutability, sizeof(double[ALPHSIZE]));
    mutab_divisor = 100.0 * s_DotProductVectors(mutab, ALPHSIZE, bgnd_probs);
    if (mutab_divisor != 0.0) {
        s_ScaleVector(mutab, ALPHSIZE, 1.0 / mutab_divisor);
    }
    /** Normalize the off-diagonals and set the diagonals */
    for (i = 0;  i < ALPHSIZE;  ++i) {
        double scale_factor =  mutab[i] / s_SumVector(pam1[i], ALPHSIZE);
        s_ScaleVector(pam1[i], ALPHSIZE, scale_factor);
        pam1[i][i] = 1.0 - mutab[i];
    }
}


/** Calculate the joint and background probabilities for PAM */
static void
s_CalculatePamProbs(double jnt_probs[][ALPHSIZE],
                    double bgnd_probs[], int n)
{
    int i, j;
    double jnt_sum = 0.0;

    double pam1[ALPHSIZE][ALPHSIZE];
    double pam[ALPHSIZE][ALPHSIZE];

    s_CalculatePam1(pam1, bgnd_probs);
    s_PowerTrueAaMatrix(pam, pam1, n);

    jnt_sum = 0.0;
    for (i = 0;  i < ALPHSIZE;  ++i) {
        for (j = 0;  j < ALPHSIZE;  ++j) {
            jnt_probs[i][j] = pam[j][i] * bgnd_probs[j];
            jnt_sum += jnt_probs[i][j];
        }
    }
    s_NormalizeJointProbs(jnt_probs);
}


/** Calculate a scoring matrix from a set of joint and background
    probabilities */
static void
s_CalcMatrixFromProbs(double pam[][ALPHSIZE],
                      double joint_probs[][ALPHSIZE], double bkgd[],
                      double lambda)
{
    int i, j;

    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < ALPHSIZE;  j++) {
            pam[i][j] = log(joint_probs[i][j] / (bkgd[i] * bkgd[j])) / lambda;
        }
    }
}


/* Print the joint probabilities in the style of C language matrices */
static void
s_PrintJointProbsAsCarrays(double jnt_probs[][ALPHSIZE], double bkgd[], int n)
{
    const int per_line = 3;  /* items per line */

    int i, j;

    printf("/** Joint probabilities for PAM%d */\n", n);
    printf("double PAM%d_JOINT_PROBS[20][20] =\n", n);
    for(i = 0; i < ALPHSIZE; i++) {
        for(j = 0; j < ALPHSIZE; j++) {
            double elt = jnt_probs[i][j];
            if (j == 0) {
                /* Start of a new row */
                if (i == 0) {
                    /* First row, open the matrix */
                    printf(" {{%.16e, ", elt);
                } else {
                    /* Open the row array */
                    printf("  {%.16e, ", elt);
                }
            } else if (j == ALPHSIZE - 1) {
                /* End of the current row */
                if (i == ALPHSIZE - 1) {
                    /* Last row, close the matrix */
                    printf("%.16e}};\n", elt);
                } else {
                    /* Close the row array */
                    printf("%.16e},\n", elt);
                }
            } else if((j + 1) % per_line == 0) {
                /* Current output line is full, move to the next */
                printf("%.16e,\n   ", elt);
            } else {
                printf("%.16e, ", elt);
            }
        }
    }
    printf("\n\n");
    printf("/** Background frequencies for PAM%d */\n", n);
    printf("double PAM%d_bg[20] =\n", n);
    for(j = 0; j < ALPHSIZE; j++) {
        double p = bkgd[j];

        if (j == 0) {
            /* Open the array */
            printf(" {%.16e, ", p);
        } else if (j == ALPHSIZE - 1) {
            /* Close the array */
            printf("%.16e};\n", p);
        } else if ((j + 1) % per_line == 0) {
            /* Current output line is full, move to the next */
            printf("%.16e,\n  ", p);
        } else {
            printf("%.16e, ", p);
        }
    }
}


/* Print the joint probabilities */
static void
s_PrintJointProbs(double jnt_probs[][ALPHSIZE])
{
    int i, j;

    for (i = 0;  i < ALPHSIZE;  i++) {
        for (j = 0;  j < ALPHSIZE;  j++) {
            printf("%10.6e   ", jnt_probs[i][j]);
        }
        printf("\n");
    }
}


/* Print the PAM matrix */
static void
s_PrintPamMatrix(double pam[][ALPHSIZE],
double Barray[ALPHSIZE+3],
    double Zarray[ALPHSIZE+3],
    double Xarray[ALPHSIZE+3])
{
    int i, j;
    
        
    
    for (i = 0;  i < ALPHSIZE;  i++) {
        printf("\t");
        for (j = 0;  j < ALPHSIZE;  j++) {
            printf("%4.4f\t", pam[i][j]);
        }
        printf("%4.4f\t", Barray[i]);
        printf("%4.4f\t", Zarray[i]);
        printf("%4.4f\t", Xarray[i]);
        printf("\n");
    }
    
    printf("\t");
    for (i = 0;  i < ALPHSIZE+3;  i++) {
        printf("%4.4f\t", Barray[i]);
    }
    printf("\n\t");
    for (i = 0;  i < ALPHSIZE+3;  i++) {
        printf("%4.4f\t", Zarray[i]);
    }
    printf("\n\t");
    for (i = 0;  i < ALPHSIZE+3;  i++) {
        printf("%4.4f\t", Xarray[i]);
    }
    printf("\n");
    
    for(j = 0; j < ALPHSIZE; ++j) {
        printf("   %c", alphabet[j]);
    }
    printf("\n");
}


/** Convert a string to a positive integer. */
static int
as_positive_int(const char * str, int * status)
{
    int n;
    char * end;
    n = strtol(str, &end, 10);
    if (n > 0) {
        /* We have interpreted a prefix of the string as an integer;
         * allow optional spaces following the integer, but nothing
         * else. */
        while (*end && isspace(*end))
           end++; 
        
        if (*end == '\0') {
            *status = 0;
            return n;
        }
    }
    fprintf(stderr,
            "ERROR: can not interpret %s as a positive integer\n", str);
    *status = 1;
    return 0;
}


/** Convert a string to a positive double */
static double
as_positive_double(const char * str, int * status)
{
    double x;
    char * end;
    x = strtod(str, &end);
    if (x > 0.0) {
        /* We have interpreted a prefix of the string as a double;
         * allow optional spaces following the double, but nothing
         * else. */
        while (*end && isspace(*end))
           end++; 
        
        if (*end == '\0') {
            *status = 0;
            return x;
        }
    }
    fprintf(stderr,
            "ERROR: can not interpret %s as a positive floating "
            "point number\n", str);
    *status = 1;
    return 0.0;
}
        
/* Options specifying whether or how to print joint probabilities */
enum {dontPrintJointProbs = 0,  printJPtoLowPrecision = 1,
      printJPasCarrays = 2};

/** Options structure for this program */
typedef struct PamOptions {
    int n;                   /**< Pam matrix number */
    int print_joint_probs;   /**< Print the joint probabilities, instead of
                                  the matrix */
    double lambda;           /**< the scale of the matrix */
} PamOptions;


/* Parse the command-line options; exit the program on error */
static void
s_ParseOptions(PamOptions * options, int argc, char ** argv)
{
    int i;
    int print_help_and_exit = 0, bad_option = 0;

    options->n = 70;
    options->print_joint_probs = dontPrintJointProbs;
    options->lambda = log(2)/3.0;

    /* Parse options; we don't assume getopt is available */
    for (i = 1;  i < argc;  ++i) {
        if (0 == strcmp("-n", argv[i])) {
            /* Next argument must be an integer */
            if (++i >= argc) {
                bad_option = i - 1;
                break;
            }
            options->n = as_positive_int(argv[i], &bad_option);
            if (bad_option) {
                bad_option = i;
                break;
            }
        } else if (0 == strncmp("-n", argv[i], 2)) {
            /* The rest of the argument must be an integer */
            options->n = as_positive_int(argv[i] + 2, &bad_option);
            if (bad_option) {
                bad_option = i;
                break;
            }
        } else if (0 == strcmp("-s", argv[i])) {
            if (++i >= argc) {
                bad_option = i - 1;
                break;
            }
            options->lambda = as_positive_double(argv[i], &bad_option);
            if (bad_option) {
                bad_option = i;
                break;
            }
        } else if (0 == strncmp("-s", argv[i], 2)) {
            /* The rest of the argument must be a float */
            options->lambda = as_positive_double(argv[i] +  2, &bad_option);
            if (bad_option) {
                bad_option = i;
                break;
            }
        } else if (0 == strcmp("-j", argv[i])) {
            options->print_joint_probs = printJPtoLowPrecision;
        } else if (0 == strcmp("-C", argv[i])) {
            options->print_joint_probs = printJPasCarrays;
        } else if (0 == strcmp("-?", argv[i]) ||
                   0 == strcmp("-h", argv[i]) ||
                   0 == strcmp("--help", argv[i])) {
            print_help_and_exit = 1;
            break;
        } else {
            bad_option = i;
            break;
        }
    }
    if (bad_option)
        fprintf(stderr, "Unrecognized option: %s\n", argv[bad_option]);

    if (print_help_and_exit || bad_option) {
        FILE * outfile;
        if (bad_option) {
            outfile = stderr;
        } else {
            outfile = stdout;
        }
        fprintf(outfile, "\nUsage pam [-j] [-n N] [-s scale]\n\n");
        fprintf(outfile, "Print the Dayhoff and Dayhoff PAM N matrix (used to"
                " score alignments\nof protein sequences).\n\n");
        fprintf(outfile, "  -n N            print PAM N (default 70)\n");
        fprintf(outfile, "  -s scale        scale of the PAM matrix\n");
        fprintf(outfile, "  -j              print joint probabilities "
                "instead of matrix scores\n");
        fprintf(outfile, "  -C              print joint and background "
                "probabilities in the syntax of\n"
                "                  the C programming language, "
                "instead of printing matrix\n"
                "                  scores\n");
        fprintf(outfile, "  -?, -h, --help  display this help and exit\n\n");
        if (bad_option) {
            exit(1);
        } else {
            exit(0);
        }
    }
}

static void s_CalcBZX(double pam[][ALPHSIZE],double Barray[],double Zarray[],double Xarray[],double background_freqs[]){
    int i,j;
    double total;
    
    Barray[ALPHSIZE+2] = 0;
    Zarray[ALPHSIZE+2] = 0;
    Xarray[ALPHSIZE+2] = 0;
    
    for (i = 0;  i < ALPHSIZE;  i++) {
       total += background_freqs[i];
    }
    
    //printf("Total: %lf\n",total);
    
    for (i = 0;  i < ALPHSIZE;  i++) {
       Barray[i] = (background_freqs[2] * pam[i][2] + background_freqs[3] * pam[i][3])/(background_freqs[2]+background_freqs[3]);
       Zarray[i] = (background_freqs[5] * pam[i][5] + background_freqs[6] * pam[i][6])/(background_freqs[5]+background_freqs[6]);
       Xarray[i] = 0;
       total = 0;
       for (j = 0;  j < ALPHSIZE;  j++) {
           Xarray[i] += background_freqs[j] * pam[i][j];
           Xarray[ALPHSIZE+2] += background_freqs[j] * background_freqs[i] * pam[i][j];
       }
       Barray[ALPHSIZE+2] += Barray[i]/ALPHSIZE;
       Zarray[ALPHSIZE+2] += Zarray[i]/ALPHSIZE;
    }
    
    Xarray[ALPHSIZE] = Barray[ALPHSIZE+2];
    Xarray[ALPHSIZE+1] = Zarray[ALPHSIZE+2];
    
    Barray[ALPHSIZE] = background_freqs[2] * background_freqs[2] * pam[2][2];
    Barray[ALPHSIZE] += background_freqs[3] * background_freqs[3] * pam[3][3];
    Barray[ALPHSIZE] += background_freqs[2] * background_freqs[3] * pam[2][3] * 2;
    Barray[ALPHSIZE] /= ((background_freqs[2] + background_freqs[3]));
    Barray[ALPHSIZE] *= Barray[ALPHSIZE];
    
    Zarray[ALPHSIZE+1] = background_freqs[5] * background_freqs[5] * pam[5][5];
    Zarray[ALPHSIZE+1] += background_freqs[6] * background_freqs[6] * pam[6][6];
    Zarray[ALPHSIZE+1] += background_freqs[5] * background_freqs[6] * pam[5][6] * 2;
    Zarray[ALPHSIZE+1] /= ((background_freqs[5] + background_freqs[6]));
    Zarray[ALPHSIZE+1] *= Zarray[ALPHSIZE+1];
    
    Barray[ALPHSIZE+1] = background_freqs[2] * background_freqs[5] * pam[2][5];
    Barray[ALPHSIZE+1] += background_freqs[2] * background_freqs[6] * pam[2][6];
    Barray[ALPHSIZE+1] += background_freqs[3] * background_freqs[5] * pam[3][5];
    Barray[ALPHSIZE+1] += background_freqs[3] * background_freqs[6] * pam[3][6];
    Barray[ALPHSIZE+1] /= ((background_freqs[2] + background_freqs[3])*(background_freqs[5] + background_freqs[6]));
    Barray[ALPHSIZE+1] *= Barray[ALPHSIZE+1];
    Zarray[ALPHSIZE] = Barray[ALPHSIZE+1];
    
}

int main(int argc, char **argv)
{
    PamOptions options;

    double jnt_probs[ALPHSIZE][ALPHSIZE];
    double pam[ALPHSIZE][ALPHSIZE];
    double background_freqs[ALPHSIZE];
    double Barray[ALPHSIZE+3];
    double Zarray[ALPHSIZE+3];
    double Xarray[ALPHSIZE+3];


    s_ParseOptions(&options, argc, argv);

    s_CalculatePamProbs(jnt_probs, background_freqs, options.n);
    switch(options.print_joint_probs) {
    case dontPrintJointProbs:
        s_CalcMatrixFromProbs(pam, jnt_probs, background_freqs,
                              options.lambda);
        s_CalcBZX(pam,Barray,Zarray,Xarray,background_freqs);
        s_PrintPamMatrix(pam,Barray,Zarray,Xarray);
        break;
    case printJPtoLowPrecision:
        s_PrintJointProbs(jnt_probs);
        break;
    case printJPasCarrays:
        s_PrintJointProbsAsCarrays(jnt_probs, background_freqs, options.n);
        break;
    }
    return 0;
}
