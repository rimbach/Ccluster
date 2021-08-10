#include <string.h>
#include <stdio.h>
#include "ccluster/ccluster.h"

#include "parseArgs.h"

compRat_poly_t p_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s <filename> [OPTIONS]", argv[0]);
        printf("                                 \n");
        printf("      -d , --domain: the initial region of interest\n");
        printf("                     global [default] finds all the roots\n");
        printf("                     a box, for instance 0,1/2,100 i.e. the square centered in 0 +i*(1/2) of width 100\n");
        printf("                     if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("      -e , --epsilon: the size of output clusters\n");
        printf("                     +inf [default] output natural clusters with less roots than degree of input polynomial\n");
        printf("                                                                 or size less than 2^(-53)\n");
        printf("                     a positive rational as 1 or 1/100 or a negative power of 2 as -53 for 2^(-53)\n");
        printf("      -o , --output: the way cluster are output; default is NO OUTPUT\n");
        printf("                     0: [default] NO OUTPUT\n");
        printf("                     d>0: d digit precision floating point numbers\n");
        printf("                     -1: rational numbers\n");
        printf("                     -2 or g or G: gnuplot output: can be piped to gnuplot \n");
        printf("                     -3 or gs or GS: gnuplot output with subdivision tree \n");
        printf("      -m, --mode: the version of the algorithm\n");
        printf("                     default [default]: uses improved version of root radii and subdivision; no multi-threading  \n");
        printf("                     CASC2021: uses root radii algo described in CASC2021 and subdivision; no multi-threading  \n");
        printf("                     onlySubd: uses only subdivision\n");
        printf("      -v, --verbose: an integer for verbosity\n");
        printf("                     0: nothing\n");
        printf("                     1 [default]: abstract of input and output\n");
        printf("                     2: detailed reports concerning algorithm\n");
        printf("                     >=3: debugging mode\n");
        printf("      -j, --nbThreads: an positive integer for the number of threads\n");
        printf("                       1 [default]: one thread is used\n");
        printf("                       >1: no compatibility with -o -3 option\n");
        printf("                       >1: no compatibility with -m default\n");
        return -1;
    }
    
    if (argc<=2){ /* display usage */
        printf("usage: %s <filename> [OPTIONS]\n", argv[0]);
        printf("   or: %s to see options\n", argv[0]);
    }
    
    int parse = 1;
//     int degree;
    char * filename;
    char * st;
    int verbosity=1;
    int nbthreads = 1;
    int global = 2; /* by default, search all the roots */
    int infinity = 2;
    int output = 0;
    
    compBox_t bInit;
    realRat_t eps;
    
    char stDefault[] = "default"; 
    st = stDefault;
    
    compBox_init(bInit);
    realRat_init(eps);
    scan_epsilon( "+inf", eps );
    global = scan_initialBox( "global", bInit );
    
    filename = argv[1];
    
    /* loop on arguments to figure out options */
    
    for (int arg = 2; arg< argc; arg++) {
        
        if ( (strcmp( argv[arg], "-v" ) == 0) || (strcmp( argv[arg], "--verbose" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_verbosity(argv[arg+1], &verbosity );
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-d" ) == 0) || (strcmp( argv[arg], "--domain" ) == 0) ) {
            if (argc>arg+1) {
                global = scan_initialBox( argv[arg+1], bInit );
                parse = parse*global;
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-e" ) == 0) || (strcmp( argv[arg], "--epsilon" ) == 0) ) {
            if (argc>arg+1) {
                infinity = scan_epsilon( argv[arg+1], eps );
                parse = parse*infinity;
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-o" ) == 0) || (strcmp( argv[arg], "--output" ) == 0) ) {
            if (argc>arg+1) {
                parse = parse*scan_output(argv[arg+1], &output);
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-j" ) == 0)  ) {
            if (argc>arg+1) {
                parse = parse*scan_nbthreads(argv[arg+1], &nbthreads);
                arg++;
            }
        }
        
        if ( (strcmp( argv[arg], "-m" ) == 0) || (strcmp( argv[arg], "--mode" ) == 0) ) {
            if (argc>arg+1) {
                st = argv[arg+1];
                arg++;
            }
        }
        
    }
    
    realRat_poly_t p;
    realRat_poly_init(p);
    compRat_poly_init(p_global);
    FILE * curFile;
        
    if (parse) {
        
        /* open filename*/
        curFile = fopen (filename,"r");
        if (curFile!=NULL) {
            realRat_poly_fread(curFile, p);
            
            compRat_poly_set_realRat_poly(p_global,p);
            
            if (global==2)
                ccluster_global_interface_realRat_poly( p, eps, st, nbthreads, output, verbosity);
            else
                ccluster_interface_realRat_poly( p, bInit, eps, st, nbthreads, output, verbosity);
            
            fclose (curFile);
        }
        
    }
    
    realRat_poly_clear(p);
    compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
