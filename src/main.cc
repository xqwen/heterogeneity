#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "controller.h"

void show_banner(){

}



int main(int argc, char **argv){

    // creating the grid

    //olist.push_back(0.1);
    //phlist.push_back(0.05);

    char data_file[256];
    double EM_thresh = 0.1;
    
    memset(data_file,0,256);

    int use_zval = 0;

    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
            strcpy(data_file,argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "--zval")==0 ){
            use_zval = 1;
            continue;
        }



        fprintf(stderr, "Error: undefined option %s\n", argv[i]);
        show_banner();
        exit(0);

    }    


    // checking mandatory arguments
    if(strlen(data_file)==0){
        fprintf(stderr,"Error: data file unspecified\n");
        show_banner();
        exit(0);
    }


    // a global variable 
    controller con;
    con.load_data(data_file, use_zval);


}
