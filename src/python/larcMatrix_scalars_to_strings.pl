#!/usr/bin/perl -w   # -w generates warnings

 #*################################################################
 #                                                                #
 # Copyright (C) 2014-2024, Institute for Defense Analyses        #
 # 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           #
 # This material may be reproduced by or for the US Government    #
 # pursuant to the copyright license under the clauses at DFARS   #
 # 252.227-7013 and 252.227-7014.                                 #
 #                                                                #
 # LARC : Linear Algebra via Recursive Compression                #
 # Authors:                                                       #
 #   - Steve Cuccaro (IDA-CCS)                                    #
 #   - John Daly (LPS)                                            #
 #   - John Gilbert (UCSB, IDA adjunct)                           #
 #   - Mark Pleszkoch (IDA-CCS)                                   #
 #   - Jenny Zito (IDA-CCS)                                       #
 #                                                                #
 # Additional contributors are listed in "LARCcontributors".      #
 #                                                                #
 # Questions: larc@super.org                                      #
 #                                                                #
 # All rights reserved.                                           #
 #                                                                #
 # Redistribution and use in source and binary forms, with or     #
 # without modification, are permitted provided that the          #
 # following conditions are met:                                  #
 #   - Redistribution of source code must retain the above        #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer.                                      #
 #   - Redistribution in binary form must reproduce the above     #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer in the documentation and/or other     #
 #     materials provided with the distribution.                  #
 #   - Neither the name of the copyright holder nor the names of  #
 #     its contributors may be used to endorse or promote         #
 #     products derived from this software without specific prior #
 #     written permission.                                        #
 #                                                                #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         #
 # CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    #
 # INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       #
 # MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
 # DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        #
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   #
 # NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   #
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       #
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
 # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, #
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             #
 #                                                                #
 #*################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#   larcMatrix_scalars_to_strings.pl                                       #
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
