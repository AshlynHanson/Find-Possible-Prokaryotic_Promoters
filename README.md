# Find Possible Prokaryotic Promoters

# Main purpose of the program:
This program finds possible matches to prokaryotic promoter sequences (consists of two short sequences of -35 
and -10 positions downstream from the transcription start site) in a single DNA, and return a list of all such 
matches that exceeds a user-defined threshold. The user can select one of the possible matches and have it print 
out the original data with the promoter match highlighted in pink.

# Instructions:
When this program is run, it first displays a general overview of the program and instruction on how to use it.
The user can then enter either a file path to a text(txt) file that contains a string of accurate DNA input. 
This DNA sequence should only contain ACGT, any other nucleotide values will not be accepted and the user will 
be prompted to enter a new input. If the user would rather input the DNA input namually into the command 
window, there is an option to do that as well. 

The user will be asked to enter a threshold for the match score of the nucleotides. This would be measured by
the number of matching nucleotides to the two consensus sequences, -35 and -10. Because a consensus sequence 
can only differ by 1-2 nucleotides to still be considered a consensus sequence, the threshold should most 
likely be found in the 8-12 range.

The program loops through the DNA input and compares each substring of 6 nucleotides with the two consensus 
sequences and stores any possible promoter matches in a dict with the starting position in the DNA input 
string and the match score. The list of all of the possible promoters are displayed to the user along with 
their position and match score. The user would then have the option to see the entire DNA sequence with
one selected possible promoter highlighted.
 
Users can provide the position from a list of all identified possible promoter locations in a selected 
prokaryotic promoter sequence in a single DNA and display it's place within the provided DNA input. 
This program will print the full DNA sequence, highlighting the user selected possible matching promoter 
sequences in a different color and display the overall match score for the matching nucleotides.

The user can choose to see the entire DNA with their choice of selected possible promoter highlighted
as many times as they would like. To exit the program, the user can enter 'X' and it will stop running.

(C) 2021 Ashlyn Hanson
email: ashlyndhanson@gmail.com
