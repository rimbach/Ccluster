#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
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

int main(int argc, char **argv){
    
    srand(time(NULL));
    
    if (argc<4){
        printf("usage: %s Bernoulli    format degree filename\n", argv[0]);
        printf("       %s Regular Grid format degree filename\n", argv[0]);
        printf("       %s Mignotte     format degree bitsize filename\n", argv[0]);
        printf("       %s Chebyshev1   format degree filename\n", argv[0]);
        printf("       %s Chebyshev2   format degree filename\n", argv[0]);
        printf("       %s Legendre     format degree filename\n", argv[0]);
        printf("       %s randomDense  format degree bitsize filename\n", argv[0]);
        printf("       %s randomSparse format degree bitsize nbterms filename\n", argv[0]);
//         printf("usage: %s MignotteGen format degree bitsize power filename\n", argv[0]);
//         printf("usage: %s MignotteMul format degree bitsize power filename\n", argv[0]);
//         printf("usage: %s Wilkinson   format degree filename\n", argv[0]);
//         printf("usage: %s WilkRat     format degree filename\n", argv[0]);
//         printf("usage: %s WilkMul     format degree filename\n", argv[0]);
//         printf("usage: %s Spiral      format degree prec filename\n", argv[0]);
//         printf("usage: %s NestedClusters format     iterations prec filename\n", argv[0]);
//         printf("usage: %s Mandelbrot  format iterations filename \n", argv[0]);
//         printf("usage: %s Runnels     format iterations filename \n", argv[0]);
//         printf("usage: %s Laguerre    format degree filename \n", argv[0]);
        printf("where format is 1 for ccluster, 2 for mpsolve and 3 for anewdsc\n");
        return -1;
    }
    
    char poly[100];
    char filename[100];
    char bernoulli[] = "Bernoulli\0";
    char mignotte[] = "Mignotte\0";
    char regularGrid[] = "RegularGrid\0";
    char Chebyshev1[] = "Chebyshev1\0";
    char Chebyshev2[] = "Chebyshev2\0";
    char Legendre[] = "Legendre\0";
    char randomDense[] = "randomDense\0";
    char randomSparse[] = "randomSparse\0";
//     char mignotteGen[] = "MignotteGen\0";
//     char mignotteMul[] = "MignotteMul\0";
//     char wilkinson[] = "Wilkinson\0";
//     char wilkRat[] = "WilkRat\0";
//     char wilkMul[] = "WilkMul\0";
//     char spiral[] = "Spiral\0";
//     char cluster[] = "NestedClusters\0";
//     char mandelbrot[] = "Mandelbrot\0";
//     char runnels[] = "Runnels\0";
//     char laguerre[] = "Laguerre\0";
    int format = 0;
    int degree = 0;
    int thirdArg = 0;
    int fourthArg = 0;
//     int fifthArg = 0;
//     FILE * pFile;
    
    sscanf(argv[1], "%s", poly);
    sscanf(argv[2], "%d", &format);
    sscanf(argv[3], "%d", &degree);
    sscanf(argv[argc-1], "%s", filename);
    if (argc>=6)
        sscanf(argv[4], "%d", &thirdArg);
    
    if (argc>=7)
        sscanf(argv[5], "%d", &fourthArg);
    
//     if (argc>=7)
//         sscanf(argv[5], "%d", &fifthArg);
    
    realRat_poly_t p;
    realRat_poly_init(p);
    
    if (strcmp(poly, bernoulli)==0) {
        bernoulli_polynomial( p, degree);
    }
    
    if (strcmp(poly, regularGrid)==0) {
        regularGrid_polynomial(p, degree);
    }
    
    if (strcmp(poly, Chebyshev1)==0) {
        Chebyshev1_polynomial(p, degree);
    }
    
    if (strcmp(poly, Chebyshev2)==0) {
        Chebyshev1_polynomial(p, degree);
    }
    
    if (strcmp(poly, Legendre)==0) {
        Legendre_polynomial(p, degree);
    }
    
    if (strcmp(poly, mignotte)==0) {
        if (argc<5)
            printf("usage: %s Mignotte degree bitsize filename\n", argv[0]);
        else{
            mignotte_polynomial(p, degree, thirdArg);
        }
    }
    
    if (strcmp(poly, randomDense)==0) {
        if (argc<5)
            printf("usage: %s randomDense degree bitsize filename\n", argv[0]);
        else{
            randomDense_polynomial(p, degree, thirdArg);
        }
    }
    
    if (strcmp(poly, randomSparse)==0) {
        if (argc<6)
            printf("usage: %s randomSparse degree bitsize nbterms filename\n", argv[0]);
        else{
            randomSparse_polynomial(p, degree, thirdArg, fourthArg);
        }
    }

    FILE * curFile;
    printf ("output file: %s\n", filename);
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
