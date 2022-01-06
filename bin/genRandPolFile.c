#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "base/base.h"
#include "polynomials/app_rat_poly.h"
#include "flint/fmpq_poly.h"
#include "flint/flint.h"

/* intermediate functions for random generation */
int  isIntinListInt (int el, int * lis, int length);

/* functions for polynomials generation */
void randomDense_polynomial( realRat_poly_t dest, int degree, int bitsize, flint_rand_t state);
void randomSparse_polynomial( realRat_poly_t dest, int degree, int bitsize, int nbterms, flint_rand_t state);


int main(int argc, char **argv){
    
    srand(time(NULL));
    
    if (argc<4){
        printf("usage: %s randomDense  degree [OPTIONS: format bitsize         nbpols location] \n", argv[0]);
        printf("       %s randomSparse degree [OPTIONS: format bitsize nbterms nbpols location] \n", argv[0]);
        printf("                                                                    \n");
        printf("       -f , --format: the format of the output file                 \n");
        printf("                      1: [default] Ccluster (.ccl)                  \n");
        printf("                      2:           MPsolve  (.mpl)                  \n");
        printf("                      3:           ANewDsc  (.dsc)                  \n");
        printf("       -b , --bitsize: the bitsize of the coeffs \n");
        printf("                      8: [default] or a positive integer            \n");
        printf("       -n , --nbterms: the number of non-zero coeffs (for randomSparse)\n");
        printf("                      10: [default] or a positive integer            \n");
        printf("       -p , --nbpols:  the number of polynomials generated\n");
        printf("                      1: [default] or a positive integer            \n");
        printf("       -l , --location: where the files are generated\n");
        printf("                      .: [default] or string            \n");
        return -1;
    }
    
    char poly[100];
    char filename[100];
    int firstArg  = 0;
    char location[100];
    sprintf(location, "." );
    
    /* optional arguments */
    int format    = 1;
    int bitsize   = 8;
    int nbterms   = 10;
    int nbpols    = 1;
    
    char randomDense[] = "randomDense\0";
    char randomSparse[] = "randomSparse\0";
    
    int parse=1;
    
    parse = parse*sscanf(argv[1], "%s", poly);
    parse = parse*sscanf(argv[2], "%d", &firstArg);
    parse = parse*(firstArg>=0);
    if (firstArg<0) {
        printf("%s ERROR: NON-VALID DEGREE/ITERATION (should be positive) \n", argv[0]);
    }
    parse = parse*sscanf(argv[3], "%s", filename);
    
    /* loop on arguments to figure out options */
    for (int arg = 3; arg< argc; arg++) {
        
        if ( (strcmp( argv[arg], "-f" ) == 0) || (strcmp( argv[arg], "--format" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &format);
                if ((format<=0)||(format>3)) {
                    printf("%s ERROR: NON-VALID FORMAT (should be 1, 2 or 3) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-b" ) == 0) || (strcmp( argv[arg], "--bitsize" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &bitsize);
                if (bitsize<=0){
                    printf("%s ERROR: NON-VALID BITSIZE (should be >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-n" ) == 0) || (strcmp( argv[arg], "--nbterms" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &nbterms);
                if (nbterms<=0){
                    printf("%s ERROR: NON-VALID nb of terms (should be >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-p" ) == 0) || (strcmp( argv[arg], "--nbpols" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &nbpols);
                if (nbterms<=0){
                    printf("%s ERROR: NON-VALID NUMBER OF POLS (should be >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-l" ) == 0) || (strcmp( argv[arg], "--location" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%s", location);
//                 if (nbterms<=0){
//                     printf("%s ERROR: NON-VALID PATH \n", argv[0]);
//                     parse = 0;
//                 }
                arg++;
            }
        }
        
    }
    
    if (!parse) {
        printf("%s PARSING ERROR\n", argv[0] );
        return -1;
    }
    
    if ((strcmp(poly, randomDense)!=0)&&(strcmp(poly, randomSparse)!=0)) {
        printf ("%s PARSING ERROR; INVALID POLYNOMIAL: %s\n", argv[0], poly);
        parse = 0;
        return -1;
    }
    
//     printf ("%s PARSING OK, degree: %d, format: %d\n", argv[0], firstArg, format);
    
    realRat_poly_t p;
//     realRat_poly_init(p);
    
    flint_rand_t state;
    flint_randinit(state);
    
    FILE * curFile;
    
    for (int nbp = 1; nbp<=nbpols; nbp++){
        
        realRat_poly_init(p);
        
        if (strcmp(poly, randomDense)==0) {
            if (format==1) {
                sprintf(filename, "%s/%s_%d_%d_%d.ccl", location, poly, firstArg, bitsize, nbp );
            } else if (format==2) {
                sprintf(filename, "%s/%s_%d_%d_%d.mpl", location, poly, firstArg, bitsize, nbp );
            } else if (format==3) {
                sprintf(filename, "%s/%s_%d_%d_%d.dsc", location, poly, firstArg, bitsize, nbp );
            }
        } else {
            if (format==1) {
                sprintf(filename, "%s/%s_%d_%d_%d_%d.ccl", location, poly, firstArg, bitsize, nbterms, nbp );
            } else if (format==2) {
                sprintf(filename, "%s/%s_%d_%d_%d_%d.mpl", location, poly, firstArg, bitsize, nbterms, nbp );
            } else if (format==3) {
                sprintf(filename, "%s/%s_%d_%d_%d_%d.dsc", location, poly, firstArg, bitsize, nbterms, nbp );
            }
        }
        
//         printf("filename: %s \n", filename);
        
        curFile = fopen (filename,"w");
        
        
        if (strcmp(poly, randomDense)==0) {
            randomDense_polynomial(p, firstArg, bitsize, state);
        } else if (strcmp(poly, randomSparse)==0) {
            nbterms = CCLUSTER_MIN( firstArg, nbterms );
            randomSparse_polynomial(p, firstArg, bitsize, nbterms, state);
        }
        
        if (format==1) {
            realRat_poly_fprint(curFile, p);
        }
        else if (format==2) {
            
            realRat_t coeff;
            realRat_init(coeff);
            
            fprintf(curFile, "Sparse;\n");
            fprintf(curFile, "Monomial;\n");
            fprintf(curFile, "Real;\n");
            fprintf(curFile, "Rational;\n");
            fprintf(curFile, "Degree = %ld;\n", p->length -1);
            fprintf(curFile, "\n");
            
            for(slong i = p->length -1; i>=0; i--){
                realRat_poly_get_coeff_realRat(coeff, p, i);
                if (!realRat_is_zero(coeff)) {
                    fprintf(curFile, "%ld ", i);
                    realRat_fprint(curFile, coeff);
                    fprintf(curFile, "\n");
                }
            }
            
            realRat_clear(coeff);
            
        }
        else if (format==3) {
            fmpq_poly_canonicalise(p);
            fprintf(curFile, "%ld\n", p->length -1);
            for (slong i=0; i<p->length; i++){
                fmpz_fprint(curFile, p->coeffs + i);
                fprintf(curFile, "\n");
            }
        }
        else 
            printf("format %d not implemented\n", format);
        
        fclose (curFile);
        realRat_poly_clear(p);    
    }
    
//     realRat_poly_clear(p);
}

void randomDense_polynomial( realRat_poly_t dest, int degree, int bitsize, flint_rand_t state) {
    
    realRat_poly_fit_length(dest, (slong) degree+1);
    realRat_poly_zero(dest);
    
    fmpz_randtest_not_zero(dest->coeffs + degree, state, bitsize);
    fmpz_randtest_not_zero(dest->coeffs + 0,      state, bitsize);
    for (int i=1;i<degree; i++){
        fmpz_randtest_not_zero(dest->coeffs + i,           state, bitsize);
//         fmpz_randtest(dest->coeffs + i,           state, bitsize);
    }
    dest->length = degree +1;
}

int  isIntinListInt (int el, int * lis, int length){
    for (int index=0;index<length;index++)
        if (el==lis[index])
            return 1;
    return 0;
}

void randomSparse_polynomial( realRat_poly_t dest, int degree, int bitsize, int nbterms, flint_rand_t state){
    
//     printf("degree: %d, bitsize: %d, nbterms: %d\n", degree, bitsize, nbterms);
    realRat_poly_fit_length(dest, (slong) degree+1);
    realRat_poly_zero(dest);
    
    fmpz_randtest_not_zero(dest->coeffs + degree, state, bitsize);
    fmpz_randtest_not_zero(dest->coeffs + 0,      state, bitsize);
    
    int * list_coeffs;
    int length = 2;
    int coeff;
    list_coeffs = (int *) malloc (nbterms*sizeof(int));
    list_coeffs[0]=degree;
    list_coeffs[1]=0;
    
    while (length < nbterms) {
        coeff = (rand() % (degree));
        if (!isIntinListInt (coeff, list_coeffs, length)){
            
//             fmpz_randtest(dest->coeffs + coeff,           state, bitsize);
            fmpz_randtest_not_zero(dest->coeffs + coeff,           state, bitsize);
            list_coeffs[length]=coeff;
            length++;
        }
    }
    
    dest->length = degree +1;
}

