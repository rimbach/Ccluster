#include <string.h>
#include <stdio.h>
#include "ccluster/ccluster.h"
#include "risolate/risolate.h"

#include "parseArgs.h"

// compRat_poly_t p_global;
// 
// void getApprox(compApp_poly_t dest, slong prec){
//     compApp_poly_set_compRat_poly(dest, p_global, prec);
// }

int main(int argc, char **argv){
    
    if (argc<=2){
        printf("usage: %s [OPTIONS] <filename> ", argv[0]);
        printf("                                 \n");
        printf("      -d , --domain: the initial region of interest\n");
        printf("                     global [default] finds all the roots\n");
        printf("                     a box, for instance 0,1,1,2,100,1 i.e. the square centered in 0/1 +i*(1/2) of width 100/1\n");
        printf("                     if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("      -e , --epsilon: the size of output clusters\n");
        printf("                     +inf [default] output natural clusters wits less roots than degree of input polynomial\n");
        printf("                     a positive number as 1,100 (1/100) or -53 (1/2^(-53))\n");
        printf("      -m, --mode: the version of the algorithm\n");
        printf("                     default value is \"default\"  \n");
        printf("      -v, --verbose: an integer for verbosity\n");
        printf("                     0: nothing\n");
        printf("                     1 [default]: abstract of input and output\n");
        printf("                     2: detailed reports concerning algorithm\n");
        printf("                     3: same as 2 + prints the clusters to stdout\n");
        printf("TODO: nbthreads, display precision of outout\n");
        if (argc<2)
            return -1;
    }
    
    int parse = 1;
//     int degree;
    char * filename;
    char * st;
    int verbosity=1;
    int nbthreads = 1;
    int global = 2; /* by default, search all the roots */
    int infinity = 2;
    
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
        
        if ( (strcmp( argv[arg], "-m" ) == 0) || (strcmp( argv[arg], "--mode" ) == 0) ) {
            if (argc>arg+1) {
//                 parse = parse*scan_strategy( argv[arg+1], st );
                st = argv[arg+1];
                arg++;
            }
        }
        
    }
//     if (argc>=7) {
//         parse = parse*scan_nbthreads(argv[6], &nbthreads );
//     }
    
    realRat_poly_t p;
    realRat_poly_init(p);
//     compRat_poly_init(p_global);
    FILE * curFile;
        
    if (parse) {
        
        /* open filename*/
        curFile = fopen (filename,"r");
        if (curFile!=NULL) {
            realRat_poly_fread(curFile, p);
//             compRat_poly_set_realRat_poly(p_global,p);
            
            if (global==2)
                risolate_global_interface_poly( p, eps, st, nbthreads, verbosity);
            else
                risolate_interface_poly( p, bInit, eps, st, nbthreads, verbosity);
            
            fclose (curFile);
        }
        
    }
    
    realRat_poly_clear(p);
//     compRat_poly_clear(p_global);
    realRat_clear(eps);
    compBox_clear(bInit);
    
    flint_cleanup();
    
    return 0;
}
