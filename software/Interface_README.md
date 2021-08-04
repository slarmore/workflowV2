#Interface_README

##This doc outlines the general format for software interfaces

##Contents:
    ###imports
    ###calculator creation function
        check for valid input
        parse calculation keywords
        construct input file(s)
        define run command
        create calculator object to actually run the calculation
    ###software class
        define attributes (name, input/output file extensions)
        define keywords to search for in output file
        read_output method
        fix_errors method
        resubmit method
    ###Functions for parsing output


    ##A template for building new calculators can be found at interface_template.py
    
