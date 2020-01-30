#include <string.h>
#include <stdio.h>
#include "ccluster/ccluster.h"

#include "parseArgs.h"

compRat_poly_t p_global;

void getApprox(compApp_poly_t dest, slong prec){
    compApp_poly_set_compRat_poly(dest, p_global, prec);
}

int main(int argc, char **argv){
    
    if (argc<6){
        printf("usage: %s file ", argv[0]);
        printf("domain epsilon strategy verbosity\n");
        printf("domain:      global (finds all roots) or a box given as:\n");
        printf("             for instance 0,1,1,2,100,1 i.e. the square centered in 0/1 +i*(1/2) of width 100/1\n");
        printf("             if a bounded box B is given, ccluster finds all natural clusters in B, and possibly some in (5/4)B \n");
        printf("eps:         infinity or a positive rational number, as:\n");
        printf("             for instance 1,100 (1/100) or -53 (1/2^(-53))\n");  
        printf("             if infinity, computes a set of natural clusters containing all roots in initial box\n");
        printf("strategy:    default for default strategy\n");
        printf("             test for testing mode\n");
        printf("verbosity:   0: nothing\n");
        printf("             1: abstract of input and output\n");
        printf("             2: detailed reports concerning algorithm\n");
        printf("             3: same as 2 + prints the clusters to stdout\n");
        return -1;
    }
    
    int parse = 1;
//     int degree;
    char * filename;
    char * st;
    int verbosity;
    int nbthreads = 1;
    int global = 0;
    int infinity = 0;
    
    compBox_t bInit;
    realRat_t eps;
    
    compBox_init(bInit);
    realRat_init(eps);
    
    filename = argv[1];
//     parse = parse*scan_degree( argv[1], &degree);
    global = scan_initialBox( argv[2], bInit );
    parse = parse*global;
        
    infinity = scan_epsilon( argv[3], eps );
    parse = parse*infinity;
    
//     parse = parse*scan_strategy(argv[4], &st );
    st = argv[4];
    parse = parse*scan_verbosity(argv[5], &verbosity );
    
    if (argc>=7) {
        parse = parse*scan_nbthreads(argv[6], &nbthreads );
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
                ccluster_global_interface_func( getApprox, eps, st, nbthreads, verbosity);
            else
                ccluster_interface_func( getApprox, bInit, eps, st, nbthreads, verbosity);
            
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
