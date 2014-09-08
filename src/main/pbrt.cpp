// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"

// main program
int main(int argc, char *argv[]) {
    Options options;
    vector<string> filenames;
    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
                   "[--verbose] [--help] <filename.pbrt> ...\n");
            return 0;
        }
        else filenames.push_back(argv[i]);
    }

    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2014 Matt Pharr and Greg Humphreys.\n");
        printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n");
        printf("See the file LICENSE.txt for the conditions of the license.\n");
        fflush(stdout);
    }

    //systemwide initialization
    pbrtInit(options);
    
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {
        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++)
            if (!ParseFile(filenames[i]))
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
    }

    //final cleanup before system exits
    pbrtCleanup();
    return 0;
}


