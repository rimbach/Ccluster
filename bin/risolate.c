#include <string.h>
#include <stdio.h>
#include "ccluster/ccluster.h"
#include "risolate/risolate.h"

#include "parseArgs.h"

int main(int argc, char **argv){
    
    if (argc<2){
        printf("usage: %s <filename> [OPTIONS]", argv[0]);
        printf("                                 \n");
        printf("      -d , --domain: the initial region of interest\n");
        printf("                     global [default] finds all the real roots\n");
        printf("                     a box, for instance 0,1/2,100 i.e. the square centered in 0 +i*(1/2) of width 100\n");
        printf("                     if a bounded box B is given, risolate finds all real roots in B, and possibly some in (5/4)B \n");
        printf("      -e , --epsilon: the size of output isolating intervals\n");
        printf("                     +inf [default] output isolating intervals for real roots\n");
        printf("                     a positive rational as 1 or 1/100 or a negative power of 2 as -53 for 2^(-53)\n");
        printf("      -o , --output: the way roots are output; default is NO OUTPUT\n");
        printf("                     0: [default] NO OUTPUT\n");
        printf("                     d>0: d digit precision floating point numbers\n");
        printf("                     -1: rational numbers\n");
        printf("                     -2 or g or G: gnuplot output: can be piped to gnuplot \n");
        printf("                     -3 or gs or GS: gnuplot output with subdivision tree \n");
        printf("      -m, --mode: the version of the algorithm\n");
        printf("                     default [default]: uses root radii and subdivision\n");
        printf("                     onlySubd: uses only subdivision\n");
        printf("      -v, --verbose: an integer for verbosity\n");
        printf("                     0: nothing\n");
        printf("                     1 [default]: abstract of input and output\n");
        printf("                     2: detailed reports concerning algorithm\n");
        printf("                     3: same as 2 + prints the roots to stdout\n");
        return -1;
    }
    
    if (argc<=2){ /* display usage */
        printf("usage: %s <filename> [OPTIONS]\n", argv[0]);
        printf("   or: %s to see options\n", argv[0]);
    }
    
    int parse = 1;
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
        
        if ( (strcmp( argv[arg], "-m" ) == 0) || (strcmp( argv[arg], "--mode" ) == 0) ) {
            if (argc>arg+1) {
                st = argv[arg+1];
                arg++;
            }
        }
        
    }
    
    realRat_poly_t p;
    realRat_poly_init(p);
    FILE * curFile;
        
    if (parse) {
        
        /* open filename*/
        curFile = fopen (filename,"r");
        if (curFile!=NULL) {
            realRat_poly_fread(curFile, p);
            
            if (global==2)
                risolate_global_interface_poly( p, eps, st, nbthreads, output, verbosity);
            else
                risolate_interface_poly( p, bInit, eps, st, nbthreads, output, verbosity);
            
            fclose (curFile);
        }
        
    }
    
    realRat_poly_clear(p);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
