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
void regularGrid_polynomial(realRat_poly_t dest, int res);
void Chebyshev1_polynomial(realRat_poly_t dest, int degree);
void Chebyshev2_polynomial(realRat_poly_t dest, int degree);
void Legendre_polynomial(realRat_poly_t dest, int degree);
void randomDense_polynomial( realRat_poly_t dest, int degree, int bitsize);
void randomSparse_polynomial( realRat_poly_t dest, int degree, int bitsize, int nbterms);

void Mandelbrot_polynomial( realRat_poly_t dest, int iterations);
void Runnels_polynomial( realRat_poly_t dest, int iterations);
void Wilkinson_polynomial(realRat_poly_t dest, int degree);

void MignotteGen_polynomial( realRat_poly_t dest, int degree, int bitsize, int power);
void MignotteMul_polynomial( realRat_poly_t dest, int degree, int bitsize, int power);
void WilkRat_polynomial(  realRat_poly_t dest, int degree);
void WilkMul_polynomial(  realRat_poly_t dest, int degree);
void Laguerre_polynomial(  realRat_poly_t dest, int degree);

void genSpiralPolFile( FILE * file, int degree, int prec);
void genClusterPolFile( FILE * file, int iteration, int prec);


int main(int argc, char **argv){
    
    srand(time(NULL));
    
    if (argc<4){
        printf("usage: %s Bernoulli      degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s RegularGrid    degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Mignotte       degree    filename [OPTIONS: format bitsize] \n", argv[0]);
        printf("       %s MignotteGen    degree    filename [OPTIONS: format bitsize power]\n", argv[0]);
        printf("       %s MignotteMul    degree    filename [OPTIONS: format bitsize power]\n", argv[0]);
        printf("       %s Chebyshev1     degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Chebyshev2     degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Legendre       degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s randomDense    degree    filename [OPTIONS: format bitsize] \n", argv[0]);
        printf("       %s randomSparse   degree    filename [OPTIONS: format bitsize nbterms] \n", argv[0]);
        printf("       %s Wilkinson      degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s WilkRat        degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s WilkMul        nbOfRoots filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Mandelbrot     iteration filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Runnels        iteration filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Laguerre       degree    filename [OPTIONS: format] \n", argv[0]);
        printf("       %s Spiral         degree    filename [OPTIONS:        precision]; only for MPsolve\n", argv[0]);
        printf("       %s nestedClusters iteration filename [OPTIONS:        precision]; only for MPsolve\n", argv[0]);
        printf("                                                                    \n");
        printf("       -f , --format: the format of the output file                 \n");
        printf("                      1: [default] Ccluster (.ccl)                  \n");
        printf("                      2:           MPsolve  (.mpl)                  \n");
        printf("                      3:           ANewDsc  (.dsc)                  \n");
        printf("       -b , --bitsize: the bitsize of the coeffs (for Mignotte, randomDense, randomSparse\n");
        printf("                      8: [default] or a positive integer            \n");
        printf("       -n , --nbterms: the number of non-zero coeffs (for randomSparse)\n");
        printf("                      10: [default] or a positive integer            \n");
        printf("       -p , --power: (for MignotteGen: nb of roots in the cluster, \n");
        printf("                     (for MignotteMul: multiplicity of the roots, \n");
        printf("                      2: [default] or a positive integer            \n");
        printf("       -L , --precision: (for NestedClusters and Spiral)            \n"); 
        printf("                      53: [default] or a positive integer            \n");
        return -1;
    }
    
    char poly[100];
    char filename[100];
    int firstArg  = 0;
    
    /* optional arguments */
    int format    = 1;
    int bitsize   = 8;
    int nbterms   = 10;
    int power     = 2;
    int precision = 53;
    
    char bernoulli[] = "Bernoulli\0";
    char mignotte[] = "Mignotte\0";
    char mignotteMul[] = "MignotteMul\0";
    char mignotteGen[] = "MignotteGen\0";
    char regularGrid[] = "RegularGrid\0";
    char Chebyshev1[] = "Chebyshev1\0";
    char Chebyshev2[] = "Chebyshev2\0";
    char Legendre[] = "Legendre\0";
    char randomDense[] = "randomDense\0";
    char randomSparse[] = "randomSparse\0";
    char Wilkinson[] = "Wilkinson\0";
    char WilkMul[] = "WilkMul\0";
    char WilkRat[] = "WilkRat\0";
    char Mandelbrot[] = "Mandelbrot\0";
    char Runnels[] = "Runnels\0";
    char Laguerre[] = "Laguerre\0";
    char Spiral[] = "Spiral\0";
    char nestedClusters[] = "nestedClusters\0";
    
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
        
        if ( (strcmp( argv[arg], "-p" ) == 0) || (strcmp( argv[arg], "--power" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &power);
                if (nbterms<=0){
                    printf("%s ERROR: NON-VALID POWER (should >0) \n", argv[0]);
                    parse = 0;
                }
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-L" ) == 0) || (strcmp( argv[arg], "--precision" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*sscanf(argv[arg+1], "%d", &precision);
                if (nbterms<=0){
                    printf("%s ERROR: NON-VALID PRECISION (should >0) \n", argv[0]);
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
        
    if (strcmp(poly, mignotteMul)==0) {
        MignotteMul_polynomial(p, firstArg, bitsize, power);
    } else
        
    if (strcmp(poly, mignotteGen)==0) {
        MignotteGen_polynomial(p, firstArg, bitsize, power);
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
        
    if (strcmp(poly, WilkMul)==0) {
        WilkMul_polynomial( p, firstArg);
    } else
        
    if (strcmp(poly, WilkRat)==0) {
        WilkRat_polynomial( p, firstArg);
    } else
        
    if (strcmp(poly, Laguerre)==0) {
        Laguerre_polynomial( p, firstArg);
    } else
        
    if ((strcmp(poly, Spiral)==0)&&(parse==1)) {
        FILE * curFile;
        printf ("%s PARSING OK; output file: %s\n", argv[0], filename);
        curFile = fopen (filename,"w");
        genSpiralPolFile( curFile, firstArg, precision );
        fclose (curFile);
        realRat_poly_clear(p);
        return 1;
    } else
        
    if ((strcmp(poly, nestedClusters)==0)&&(parse==1)) {
        FILE * curFile;
        printf ("%s PARSING OK; output file: %s\n", argv[0], filename);
        curFile = fopen (filename,"w");
        genClusterPolFile( curFile, firstArg, precision );
        fclose (curFile);
        realRat_poly_clear(p);
        return 1;
    } else
        
    {
        printf ("%s PARSING ERROR; INVALID POLYNOMIAL: %s\n", argv[0], poly);
        parse = 0;
        realRat_poly_clear(p);
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

void MignotteGen_polynomial( realRat_poly_t dest, int degree, int bitsize, int power){
    mignotte_generalized(dest, (slong) degree, (ulong) power, (slong) bitsize);
}

void MignotteMul_polynomial( realRat_poly_t dest, int degree, int bitsize, int power){
    mignotte_polynomial(dest, (slong) degree, (ulong) bitsize);
    realRat_poly_pow(dest, dest, (ulong) power);
}

void WilkRat_polynomial(  realRat_poly_t dest, int degree){
    
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,2);
    realRat_poly_one(dest);
    realRat_poly_zero(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, ((ulong) degree)+1);
        realRat_poly_mul(dest, dest, ptemp);
    }
    
    realRat_poly_clear(ptemp);
    
}

void WilkMul_polynomial(  realRat_poly_t dest, int degree){
    
    realRat_poly_t ptemp;
    realRat_poly_init2(ptemp,degree+1);
    realRat_poly_one(dest);
    
    for (int i=1; i<=degree; i++){
        realRat_poly_zero(ptemp);
        realRat_poly_set_coeff_si_ui(ptemp, 1, 1, 1);
        realRat_poly_set_coeff_si_ui(ptemp, 0, -i, 1);
        realRat_poly_pow(ptemp, ptemp, (ulong) i);
//         realRat_poly_pow(wilkptemp, ptemp, (ulong) (degree-i)+1);
        realRat_poly_mul(dest, dest, ptemp);
    }
    
    realRat_poly_clear(ptemp);
}

void Laguerre_polynomial(  realRat_poly_t dest, int degree){
    
    realRat_poly_t pone, pzero, ptemp;
    realRat_poly_init(pone);
    realRat_poly_init(pzero);
    realRat_poly_init(ptemp);
    realRat_t coeff;
    realRat_init(coeff);
    
    realRat_poly_one(pzero);
    realRat_poly_one(pone);
    realRat_poly_set_coeff_si_ui(pone, 1, -1, 1);
    realRat_poly_one(ptemp);
    realRat_poly_set_coeff_si_ui(ptemp, 1, -1, 1);
    
    for (int i = 1; i<degree; i++) {
        
        realRat_poly_set_coeff_si_ui(ptemp, 0, 2*i+1, 1);
        
        realRat_poly_mul(dest, ptemp, pone);
        realRat_set_si(coeff, (slong) i, 1);
        realRat_mul_si(coeff, coeff, (slong) -i);
        realRat_poly_scalar_mul_realRat(pzero, pzero, coeff);
        realRat_poly_add(dest, dest, pzero);
        
        realRat_poly_set(pzero, pone);
        realRat_poly_set(pone, dest);
    }
    
    realRat_poly_clear(pone);
    realRat_poly_clear(pzero);
    realRat_poly_clear(ptemp);
}

void genSpiralPolFile( FILE * file, int degree, int prec){
    
    compApp_poly_t dest;
    compApp_poly_init(dest);
    
    realRat_t modu;
    realRat_t argu;
    compApp_t a_modu;
    compApp_t a_argu;
    compApp_t coeff;
    
    realRat_init(modu);
    realRat_init(argu);
    compApp_init(a_modu);
    compApp_init(a_argu);
    compApp_init(coeff);
    
    compApp_poly_t temp;
    compApp_poly_init2(temp,2);
    compApp_poly_set_coeff_si(temp, 1, 1);
    
    compApp_poly_one(dest);
    slong prectemp = degree*prec;
    
    for(int i=1; i<=degree; i++){
        realRat_set_si(modu, -i, (ulong) degree);
        realRat_set_si(argu, 4*i, (ulong) degree);
        compApp_set_realRat( a_modu, modu, prectemp);
        compApp_set_realRat( a_argu, argu, prectemp);
        compApp_exp_pi_i( coeff, a_argu, prectemp);
        compApp_mul( coeff, coeff, a_modu, prectemp);
        compApp_poly_set_coeff_compApp(temp, 0, coeff);
        compApp_poly_mul(dest, dest, temp, prectemp);
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
    }
    
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_poly_clear(temp);
    
    
//     fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Complex;\n");
    fprintf(file, "FloatingPoint;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "Precision = %d;\n", (int) prec);
    fprintf(file, "\n");
    
    char tempstr[100*prec];
    char * temp2;
    
    for(int i = 0; i<=degree; i++){
        compApp_set(coeff, compApp_poly_getCoeff(dest, i));
        
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
        
        temp2 = arb_get_str(compApp_realref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s ", tempstr);
        temp2 = arb_get_str(compApp_imagref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s\n", tempstr);
        
    }

    compApp_clear(coeff);    
    compApp_poly_clear(dest);
}

void clustersIterate( compApp_poly_ptr tabres, compApp_poly_ptr tabprec, int i, slong prec){
    // tabres is a table of 3^i compApp_poly
    // tabres is a table of 3^(i-1) compApp_poly
    realRat_t modu;
    realRat_t argu;
    compApp_t a_modu;
    compApp_t a_argu;
    compApp_t coeff;
    
    realRat_init(modu);
    realRat_init(argu);
    compApp_init(a_modu);
    compApp_init(a_argu);
    compApp_init(coeff);
    int indexInTabRes = 0;
//     printf("pow(3,i-1): %d\n", (int) pow(3,i-1));
    for (int j = 0; j<((int) pow(3,i-1)); j++){
//         printf("%d\n", j);
//         realRat_set_si(modu, -1, (ulong) 0x1<<(4*(i-1)));
        realRat_set_si(modu, -1, (ulong) pow(4, 2*(i-1)));
        compApp_set_realRat( a_modu, modu, prec);
        
        realRat_set_si(argu, 2, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 4, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
        realRat_set_si(argu, 6, 3);
        compApp_set_realRat( a_argu, argu, prec);
        compApp_exp_pi_i( coeff, a_argu, prec);
        compApp_mul( coeff, coeff, a_modu, prec);
        compApp_add( coeff, coeff, compApp_poly_getCoeff(tabprec + j, 0), prec);
        compApp_poly_set( tabres + indexInTabRes, tabprec + j);
        compApp_poly_set_coeff_compApp(tabres + indexInTabRes, 0, coeff);
        indexInTabRes +=1;
        
    }
    realRat_clear(modu);
    realRat_clear(argu);
    compApp_clear(a_modu);
    compApp_clear(a_argu);
    compApp_clear(coeff);    
}

void genClusterPolFile( FILE * file, int iterations, int prec){
    
    compApp_poly_ptr tabprec;
    tabprec = (compApp_poly_ptr) malloc (1*sizeof(compApp_poly));
    compApp_poly_init2( tabprec, 2);
    compApp_poly_zero(tabprec);
    compApp_poly_set_coeff_si(tabprec, 1, 1);
    int degree = 3;
    
    slong prectemp = iterations*10*prec;
    
    for (int i=1; i<=iterations; i++) {
        compApp_poly_ptr tabres;
        tabres = (compApp_poly_ptr) malloc (degree*sizeof(compApp_poly));
        for(int j = 0; j<degree; j++) compApp_poly_init2( tabres + j, 2);
        
        clustersIterate( tabres, tabprec, i, prectemp);
        for(int j = 0; j<((int) (degree/3)); j++) compApp_poly_clear( tabprec + j);
        free(tabprec);
        tabprec = tabres;
        degree = degree*3;
    }
    
    degree = (int) degree/3;
    compApp_poly_t dest;
    compApp_poly_init2(dest, degree+1);
    compApp_poly_one(dest);
    for(int j = 0; j<degree; j++) compApp_poly_mul(dest, dest, tabprec +j, prectemp);
    for(int j = 0; j<degree; j++) compApp_poly_clear( tabprec + j);
    free(tabprec);
    
    compApp_t coeff;
    compApp_init(coeff);
    
//     fprintf(file, "Sparse;\n");
    fprintf(file, "Monomial;\n");
    fprintf(file, "Complex;\n");
    fprintf(file, "FloatingPoint;\n");
    fprintf(file, "Degree = %d;\n", (int) degree);
    fprintf(file, "Precision = %d;\n", (int) prec);
    fprintf(file, "\n");
    
    char tempstr[100*prec];
    char * temp2;
    
    for(int i = 0; i<=degree; i++){
        compApp_set(coeff, compApp_poly_getCoeff(dest, i));
        
//         printf("%s\n", arb_get_str(compApp_realref(coeff), prec, 0));
        temp2 = arb_get_str(compApp_realref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s ", tempstr);
        temp2 = arb_get_str(compApp_imagref(coeff), prec, ARB_STR_NO_RADIUS);
        sprintf(tempstr, "%s", temp2);
        free(temp2);
        if (tempstr[0]=='[') sprintf(tempstr, "0.0");
        fprintf(file, "%s\n", tempstr);
        
    }

    compApp_clear(coeff);    
    compApp_poly_clear(dest);
}
