#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "base/base.h"
#include "polynomials/app_rat_poly.h"
#include "flint/fmpq_poly.h"

/* intermediate functions for random generation */
int  isIntinListInt (int el, int * lis, int length);

/* functions for polynomials generation */
void regularGrid_polynomial(realRat_poly_t dest, int res);
void Chebyshev1_polynomial(realRat_poly_t dest, int degree);
void Chebyshev2_polynomial(realRat_poly_t dest, int degree);
void Legendre_polynomial(realRat_poly_t dest, int degree);
void randomDense_polynomial( realRat_poly_t dest, int degree, int bitsize);
void randomSparse_polynomial( realRat_poly_t dest, int degree, int bitsize, int nbterms);

void Mandelbrot_polynomial( realRat_poly_t dest, int iterations);
void Runnels_polynomial( realRat_poly_t dest, int iterations);
void Wilkinson_polynomial(realRat_poly_t dest, int degree);

int main(int argc, char **argv){
    
    srand(time(NULL));
    
    if (argc<4){
        printf("usage: %s Bernoulli    degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s RegularGrid  degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Mignotte     degree    filename [OPTIONS: format bitsize] \n", argv[0]);
        printf("       %s Chebyshev1   degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Chebyshev2   degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Legendre     degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s randomDense  degree    filename [OPTIONS: format bitsize] \n", argv[0]);
        printf("       %s randomSparse degree    filename [OPTIONS: format bitsize nbterms] \n", argv[0]);
        printf("       %s Wilkinson    degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Mandelbrot   iteration filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Runnels      iteration filename [OPTIONS: format] \n", argv[0]);
        printf("                                                                    \n");
        printf("       -f , --format: the format of the output file                 \n");
        printf("                      1: [default] Ccluster (.ccl)                  \n");
        printf("                      2:           MPsolve  (.mpl)                  \n");
        printf("                      3:           ANewDsc  (.dsc)                  \n");
        printf("       -b , --bitsize: the bitsize of the coeffs (for Mignotte, randomDense, randomSparse\n");
        printf("                      8: [default] or a positive integer            \n");
        printf("       -n , --nbterms: the bnumber of non-zero coeffs (for randomSparse\n");
        printf("                      10: [default] or a positive integer            \n");
        return -1;
    }
    
    char poly[100];
    char filename[100];
    int firstArg  = 0;
    
    /* optional arguments */
    int format    = 1;
    int bitsize   = 8;
    int nbterms   = 10;
    
    char bernoulli[] = "Bernoulli\0";
    char mignotte[] = "Mignotte\0";
    char regularGrid[] = "RegularGrid\0";
    char Chebyshev1[] = "Chebyshev1\0";
    char Chebyshev2[] = "Chebyshev2\0";
    char Legendre[] = "Legendre\0";
    char randomDense[] = "randomDense\0";
    char randomSparse[] = "randomSparse\0";
    char Wilkinson[] = "Wilkinson\0";
    char Mandelbrot[] = "Mandelbrot\0";
    char Runnels[] = "Runnels\0";
    
    int parse=1;
    
    parse = parse*sscanf(argv[1], "%s", poly);
    parse = parse*sscanf(argv[2], "%d", &firstArg);
    parse = parse*(firstArg>=0);
    if (firstArg<0) {
        printf("%s ERROR: NON-VALID DEGREE/ITERATION (should be positive) \n", argv[0]);
    }
    parse = parse*sscanf(argv[3], "%s", filename);
    
    /* loop on arguments to figure out options */
    for (int arg = 4; arg< argc; arg++) {
        
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
                    printf("%s ERROR: NON-VALID BITSIZE (should >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-n" ) == 0) || (strcmp( argv[arg], "--nbterms" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &nbterms);
                if (nbterms<=0){
                    printf("%s ERROR: NON-VALID BITSIZE (should >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
    }
    
    if (!parse) {
        printf("%s PARSING ERROR\n", argv[0] );
        return -1;
    }
    
    realRat_poly_t p;
    realRat_poly_init(p);
    
    if (strcmp(poly, bernoulli)==0) {
        bernoulli_polynomial( p, firstArg);
    } else
    
    if (strcmp(poly, regularGrid)==0) {
        regularGrid_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, Chebyshev1)==0) {
        Chebyshev1_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, Chebyshev2)==0) {
        Chebyshev1_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, Legendre)==0) {
        Legendre_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, mignotte)==0) {
        mignotte_polynomial(p, firstArg, bitsize);
    } else
    
    if (strcmp(poly, randomDense)==0) {
        randomDense_polynomial(p, firstArg, bitsize);
    } else
    
    if (strcmp(poly, randomSparse)==0) {
        nbterms = CCLUSTER_MIN( firstArg, nbterms );
        randomSparse_polynomial(p, firstArg, bitsize, nbterms);
    } else
    
    if (strcmp(poly, Mandelbrot)==0) {
        Mandelbrot_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, Runnels)==0) {
        Runnels_polynomial(p, firstArg);
    } else
    
    if (strcmp(poly, Wilkinson)==0) {
        Wilkinson_polynomial( p, firstArg);
    } else
    {
        printf ("%s PARSING ERROR; INVALID POLYNOMIAL: %s\n", argv[0], poly);
        parse = 0;
        return -1;
    }

    FILE * curFile;
    printf ("%s PARSING OK; output file: %s\n", argv[0], filename);
    curFile = fopen (filename,"w");
    
    if (format==1) {
//         fmpq_poly_canonicalise(p);
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
//          fmpz_t den, coeff;
//          fmpz_init(den);
//          fmpz_init(coeff);
//          fmpz_set(den, fmpq_poly_denref(p));
         
//          fmpz_print(den);
         fprintf(curFile, "%ld\n", p->length -1);
         for (slong i=0; i<p->length; i++){
//              fmpz_mul(coeff, p->coeffs + i, den);
//              fmpz_fprint(curFile, p->coeffs + i);
//              fprintf(curFile, "     ");
//              fmpz_fprint(curFile, den);
//              fprintf(curFile, "     ");
//              fmpz_fprint(curFile, coeff);
//              fprintf(curFile, "\n");
             fmpz_fprint(curFile, p->coeffs + i);
             fprintf(curFile, "\n");
         }
         
         
//          fmpz_clear(den);
//          fmpz_clear(coeff);
    }
    else 
        printf("format %d not implemented\n", format);
    
    fclose (curFile);
    
    realRat_poly_clear(p);
}

void regularGrid_polynomial(realRat_poly_t preg, int res){
    
    realRat_poly_t ptemp, ptemp2;
    realRat_poly_init2(ptemp,2);
    realRat_poly_init2(ptemp2,2);
    realRat_poly_one(preg);
    realRat_poly_zero(ptemp);
    realRat_poly_zero(ptemp2);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    realRat_poly_set_coeff_si_ui(ptemp2, 2, 1, 1);
    
    for (int i=0; i<=res; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(preg, preg, ptemp);
        
        for (int j=1; j<=res; j++){
            realRat_poly_set_coeff_si_ui(ptemp2, 1, 2*i, 1);
            realRat_poly_set_coeff_si_ui(ptemp2, 0, i*i+j*j, 1);
            realRat_poly_mul(preg, preg, ptemp2);
        }
        
        if (i>0) {
            realRat_poly_set_coeff_si_ui(ptemp, 0, i, 1);
            realRat_poly_mul(preg, preg, ptemp);
            
            for (int j=1; j<=res; j++){
                realRat_poly_set_coeff_si_ui(ptemp2, 1, -2*i, 1);
                realRat_poly_set_coeff_si_ui(ptemp2, 0, i*i+j*j, 1);
                realRat_poly_mul(preg, preg, ptemp2);
            }
        }
            
    }
    
    realRat_poly_clear(ptemp);
    realRat_poly_clear(ptemp2);
}

void Chebyshev1_polynomial(realRat_poly_t pdest, int degree){
    
    realRat_poly_t PNM1, PNM2;
    realRat_poly_init2(PNM1,2);
    realRat_poly_init2(PNM2,2);
    
    if (degree == 0) {
        realRat_poly_one(pdest);
    }
    else if (degree == 1) {
        realRat_poly_zero(pdest);
        realRat_poly_set_coeff_si_ui(pdest, 1, 1, 1);
    }
    else {
        realRat_poly_one(PNM2);
        realRat_poly_zero(PNM1);
        realRat_poly_set_coeff_si_ui(PNM1, 1, 1, 1);
        
        for (int i=2; i<=degree; i++) {
            realRat_poly_shift_left(pdest, PNM1, 1);
            realRat_poly_scalar_mul_si(pdest, pdest, 2);
            realRat_poly_neg(PNM2, PNM2);
            realRat_poly_add(pdest, pdest, PNM2);
            realRat_poly_set(PNM2, PNM1);
            realRat_poly_set(PNM1, pdest);
        }
    }
    
    realRat_poly_clear(PNM1);
    realRat_poly_clear(PNM2);
}

void Chebyshev2_polynomial(realRat_poly_t pdest, int degree){
    
    realRat_poly_t PNM1, PNM2;
    realRat_poly_init2(PNM1,2);
    realRat_poly_init2(PNM2,2);
    
    if (degree == 0) {
        realRat_poly_one(pdest);
    }
    else if (degree == 1) {
        realRat_poly_zero(pdest);
        realRat_poly_set_coeff_si_ui(pdest, 1, 2, 1);
    }
    else {
        realRat_poly_one(PNM2);
        realRat_poly_zero(PNM1);
        realRat_poly_set_coeff_si_ui(PNM1, 1, 2, 1);
        
        for (int i=2; i<=degree; i++) {
            realRat_poly_shift_left(pdest, PNM1, 1);
            realRat_poly_scalar_mul_si(pdest, pdest, 2);
            realRat_poly_neg(PNM2, PNM2);
            realRat_poly_add(pdest, pdest, PNM2);
            realRat_poly_set(PNM2, PNM1);
            realRat_poly_set(PNM1, pdest);
        }
    }
    
    realRat_poly_clear(PNM1);
    realRat_poly_clear(PNM2);
}

void Legendre_polynomial(realRat_poly_t pdest, int degree){
    
    realRat_poly_t PNM1, PNM2;
    realRat_t q;
    realRat_poly_init2(PNM1,2);
    realRat_poly_init2(PNM2,2);
    realRat_init(q);
    
    if (degree == 0) {
        realRat_poly_one(pdest);
    }
    else if (degree == 1) {
        realRat_poly_zero(pdest);
        realRat_poly_set_coeff_si_ui(pdest, 1, 1, 1);
    }
    else {
        realRat_poly_one(PNM2);
        realRat_poly_zero(PNM1);
        realRat_poly_set_coeff_si_ui(PNM1, 1, 1, 1);
        
        for (int i=2; i<=degree; i++) {
            realRat_poly_shift_left(pdest, PNM1, 1);
            realRat_poly_scalar_mul_si(pdest, pdest, (2*(i-1)+1));
            realRat_poly_scalar_mul_si(PNM2, PNM2, i-1);
            realRat_poly_neg(PNM2, PNM2);
            realRat_poly_add(pdest, pdest, PNM2);
            realRat_set_si(q, 1, i);
            realRat_poly_scalar_mul_realRat(pdest, pdest, q);
            realRat_poly_set(PNM2, PNM1);
            realRat_poly_set(PNM1, pdest);
        }
    }
    
    realRat_poly_clear(PNM1);
    realRat_poly_clear(PNM2);
    realRat_clear(q);
}

void randomDense_polynomial( realRat_poly_t dest, int degree, int bitsize) {
    
    realRat_poly_fit_length(dest, (slong) degree+1);
    realRat_poly_zero(dest);
    realRat_poly_set_coeff_si_ui(dest, (slong) degree, (slong) 1, (ulong) 1);
    ulong den = 0x1 << bitsize;
    slong num = 0;
    for (int i=0;i<degree; i++){
        num = (rand() % (2*den)) - den;
        realRat_poly_set_coeff_si_ui(dest, (slong) i, num, den);
    }
    
}

int  isIntinListInt (int el, int * lis, int length){
    for (int index=0;index<length;index++)
        if (el==lis[index])
            return 1;
    return 0;
}

void randomSparse_polynomial( realRat_poly_t dest, int degree, int bitsize, int nbterms){
    
    realRat_poly_fit_length(dest, (slong) degree+1);
    realRat_poly_zero(dest);
    realRat_poly_set_coeff_si_ui(dest, (slong) degree, (slong) 1, (ulong) 1);
    ulong den = 0x1 << bitsize;
    slong num = (rand() % (2*den)) - den;
    realRat_poly_set_coeff_si_ui(dest, 0, num, den);
    
    int * list_coeffs;
    int length = 2;
    int coeff;
    list_coeffs = (int *) malloc (nbterms*sizeof(int));
    list_coeffs[0]=degree;
    list_coeffs[1]=0;
    
    while (length < nbterms) {
        coeff = (rand() % (degree));
        if (!isIntinListInt (coeff, list_coeffs, length)){
            slong num = (rand() % (2*den)) - den;
            realRat_poly_set_coeff_si_ui(dest, (slong) coeff, num, den);
            list_coeffs[length]=coeff;
            length++;
        }
    }
}

void Mandelbrot_polynomial( realRat_poly_t pmand, int iterations){
    
    realRat_poly_t pone, px;
    realRat_poly_init(pone);
    realRat_poly_init(px);
    
    realRat_poly_one(pmand);
    realRat_poly_one(pone);
    realRat_poly_zero(px);
    realRat_poly_set_coeff_si_ui(px, 1, 1, 1);
    
    for (int i = 1; i<=iterations; i++) {
        realRat_poly_pow(pmand, pmand, 2);
        realRat_poly_mul(pmand, pmand, px);
        realRat_poly_add(pmand, pmand, pone);
    }
    
    realRat_poly_clear(pone);
    realRat_poly_clear(px);    
}

void Runnels_polynomial( realRat_poly_t prun, int iterations){
    
    realRat_poly_t prunm1, prunm2, pone, px;
    realRat_poly_init(prunm1);
    realRat_poly_init(prunm2);
    realRat_poly_init(pone);
    realRat_poly_init(px);
    
    realRat_poly_zero(px);
    realRat_poly_set_coeff_si_ui(px, 1, 1, 1);
    
    realRat_poly_one(prunm2);
    realRat_poly_set(prunm1, px);
    realRat_poly_one(prun);
    
    for (int i = 2; i<=iterations; i++) {
        
        realRat_poly_pow(prunm2, prunm2, 4);
        realRat_poly_mul(prunm2, prunm2, px);
        
        realRat_poly_pow(prun, prunm1, 2);
        realRat_poly_add(prun, prun, prunm2);
        
        realRat_poly_set(prunm2, prunm1);
        realRat_poly_set(prunm1, prun);
    }
     
    realRat_poly_clear(prunm1);
    realRat_poly_clear(prunm2);
    realRat_poly_clear(pone);
    realRat_poly_clear(px);
}

void Wilkinson_polynomial(realRat_poly_t pdest, int degree){
    
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    
    realRat_poly_one(pdest);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_mul(pdest, pdest, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}
