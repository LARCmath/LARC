#!/usr/bin/perl -w   # -w generates warnings
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#   json_scalars_to_strings.pl                                             #
#   Jenny Zito 8/2019                                                      #
#       Prepare to run this program:                                       #
#          * create a directory called "original" containing all the       #
#            json files in the pre-fall2017 format when scalars were       #
#            not saved as strings                                          #
#          * create an empty directory called "modified" where this        #
#            program will put the json files with the string format        #
#            for scalars.                                                  #
#          * create a file called "list_files" which contains a list       #
#            of all the json files in "original" (one per line)            #
#          * either change the paths in this code to point to the          #
#            larc/utils routines where they are, or make a local copy      #
#            of the routines in your directory with a internal correction  #
#            of the path to pylarc so that you can run                     #
#                 canonical_format_json.py                                 #
#                                                                          #
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# NOTE: Remember that the special characters \ | ( ) [ { ^ $ * + ? .       #
#       should be proceeded with backslash in regular expressions/printing #
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
{
    my $list = "list_files";
    my $start_dir = "original";
    my $end_dir = "modified";

    open(LIST, $list) or die "Could not open list of input files $list: $!\n";
    while (my $file_name = <LIST>) {
	chomp($file_name);
	print "The current file is $file_name\n";
	# my $file_name = "test.c";
	my $input_file = sprintf("$start_dir/%s",$file_name);
	my $output_file = sprintf("$end_dir/%s",$file_name);
	print "The path to the input file is $input_file\n";
	print "The path to the output file is $output_file\n";

	# larc/utils/canonical_format_json.py --legacy 
        #     reads in legacy format json files, it then writes
	#     the files out with the string formated scalars
	#     and it also may reduce the matrixIDs since it
	#     used a minimally sized matrix store for this operation.
        system("python canonical_format_json.py --legacy $input_file $output_file");
	    
    }
    close LIST

}

# NOTE: See the program
#          for_each_file_do_something.pl 
#       if you want to also modify the file content using regular expressions.
