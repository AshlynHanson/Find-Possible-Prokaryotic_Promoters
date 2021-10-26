# ---------------------------------------------------------------------------------------------------------------
# Main purpose of the program:
# This program finds possible matches to prokaryotic promoter sequences (consists of two short sequences of -35 
# and -10 positions downstream from the transcription start site) in a single DNA, and return a list of all such 
# matches that exceeds a user-defined threshold. The user can select one of the possible matches and have it print 
# out the original data with the promoter match highlighted in pink.
#
# Instructions:
# When this program is run, it first displays a general overview of the program and instruction on how to use it.
# The user can then enter either a file path to a text(txt) file that contains a string of accurate DNA input. 
# This DNA sequence should only contain ACGT, any other nucleotide values will not be accepted and the user will 
# be prompted to enter a new input. If the user would rather input the DNA input namually into the command 
# window, there is an option to do that as well. 
#
# The user will be asked to enter a threshold for the match score of the nucleotides. This would be measured by
# the number of matching nucleotides to the two consensus sequences, -35 and -10. Because a consensus sequence 
# can only differ by 1-2 nucleotides to still be considered a consensus sequence, the threshold should most 
# likely be found in the 8-12 range.
#
# The program loops through the DNA input and compares each substring of 6 nucleotides with the two consensus 
# sequences and stores any possible promoter matches in a dict with the starting position in the DNA input 
# string and the match score. The list of all of the possible promoters are displayed to the user along with 
# their position and match score. The user would then have the option to see the entire DNA sequence with
# one selected possible promoter highlighted.
# 
# Users can provide the position from a list of all identified possible promoter locations in a selected 
# prokaryotic promoter sequence in a single DNA and display it's place within the provided DNA input. 
# This program will print the full DNA sequence, highlighting the user selected possible matching promoter 
# sequences in a different color and display the overall match score for the matching nucleotides.
#
# The user can choose to see the entire DNA with their choice of selected possible promoter highlighted
# as many times as they would like. To exit the program, the user can enter 'X' and it will stop running.
#
# (C) 2021 Ashlyn Hanson, Yujun Zhi
# email: adh5584@truman.edu, yz4586@truman.edu
# ---------------------------------------------------------------------------------------------------------------

import re


SEQUENCE_10 = "TATAAT"  # The most common -10 consensus sequence found in most prokaryotic promoters
SEQUENCE_35 = "TTGACA"  # The most common -35 consensus sequence found in most prokaryotic promoters


def check_is_dna(input_dna_sequence):
    '''
    Checks that the input sequence consists of standard DNA nucleotides.
    The DNA input should only contain the ACTG nucleotides: adenine, thymine, cytosine, and guanine. 
    If it contains any other value, the user is told that the input is not a valid DNA sequence and
    the user is prompted to enter a new text string or file. This function also assumes
    that the DNA input does not use nucleotide ambiguity codes.

    :param str input_dna_sequence: the string of DNA inputted by the user
    '''
    if re.search('[^ACTGacgt]', input_dna_sequence):
        print("Please enter a valid DNA sequence!")
        return False
    return True


def fuzzy_match(string1, string2):
    '''
    This function that compares two strings, one of which may use ambiguity codes, and 
    returns the number of matching nucleotides. It assumes that both strings are the same length 
    an that string2 would be the only string to use ambiguity codes

    :param str string1: the first string to compare 
    :param str string2: the second string to compare that can include nucleotide ambiguity codes
    '''
    matching_nucleotides = 0
    string1 = string1.upper()
    string2 = string2.upper()

    for index in range(0, len(string1)):
        if string1[index] == string2[index]:
            matching_nucleotides += 1
        elif string2[index] == 'R':
            if string1[index] == 'A' or string1[index] == 'Y':
                matching_nucleotides += 1
        elif string2[index] == 'Y':
            if string1[index] == 'C' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'S':
            if string1[index] == 'C' or string1[index] == 'G':
                matching_nucleotides += 1
        elif string2[index] == 'W':
            if string1[index] == 'A' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'K':
            if string1[index] == 'G' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'M':
            if string1[index] == 'A' or string1[index] == 'C':
                matching_nucleotides += 1
        elif string2[index] == 'B':
            if string1[index] == 'C' or string1[index] == 'G' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'V':
            if string1[index] == 'A' or string1[index] == 'C' or string1[index] == 'G':
                matching_nucleotides += 1
        elif string2[index] == 'D':
            if string1[index] == 'A' or string1[index] == 'G' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'H':
            if string1[index] == 'A' or string1[index] == 'C' or string1[index] == 'T':
                matching_nucleotides += 1
        elif string2[index] == 'N':
            matching_nucleotides += 1
    return matching_nucleotides


# Function that steps through the input string, analyzing each substring of the appropriate length to be a possible prokaryotic promoter that includes both -35 and -10 elements and calculating an overall match score for that substring
def find_possible_prokaryotic_sequences(input_dna_sequence, threshold):
    '''
    This function steps throught the input string, analyzing each substring of the appropriate length to be a possiblee 
    prokaryotic promoter that includes both -35 and -10 elements and calculates an overall match score for that substring.
    It loops through sections of 6 nucleotides in the inputted DNA and compares it to the -35 consensus. It then checks the 
    next 6 nucleotides, within 16-19 nucleotides downstream and compares it to the -10 consensus sequence. If the difference 
    in nucleotides is with 1-2 nucleotide difference from the consensus and the count of matching nucleotides is greater than
    or equal to the user defined threshold, that sequence gets added to the dict of possible promoters. The possible promoters
    are stored in a dict with the list of the position where the -35 consensus sequence starts, the actual promoter sequence,
    and the match score. 

    :param str input_dna_sequence: the string of DNA inputted from the user
    :param str threshold: the minimum number of nucleotides that should match the consensus sequences in a range of 8-12
    :return dict possible_promoters: a dict of all the possible promoters in the DNA sequence
    '''
    possible_promoters = {
        "positions": [],
        "sequences": [],
        "scores": []
    }
    for index in range(0, len(input_dna_sequence) - 5):
        num_matching_nucleotides = fuzzy_match(SEQUENCE_35, input_dna_sequence[index: index+6])
        if num_matching_nucleotides >= 4:
            next_index = index + 21 #length of the -35 consensus plus the minimum of 16 base pairs between the consensus sequences
            while next_index < len(input_dna_sequence) - 11 and next_index < index + 31:
                # the length of a promoter is anywhere from 16 - 19 bp. It stops checking if it goes over 19 + consensus sequence
                next_matching_nucleotides = fuzzy_match(SEQUENCE_10, input_dna_sequence[next_index: next_index + 6])
                if next_matching_nucleotides >= 4 and num_matching_nucleotides + next_matching_nucleotides >= threshold:
                    if index not in possible_promoters["positions"]:
                        possible_promoters["positions"].append(index)
                        possible_promoters["sequences"].append(input_dna_sequence[index:next_index + 6])
                        possible_promoters["scores"].append(num_matching_nucleotides + next_matching_nucleotides)
                next_index += 1
    return possible_promoters


def print_promoter_sequences(possible_promoters):
    '''
    Prints a list of all the possible promoter in the user provided DNA sequence. The position
    of the possible promoters location in the DNA sequence, the actual sequence of nucleotides, and the match score
    are provided in the output. The output is stucture like the example below:

    Possible Promoters in the DNA sequence:

    Position      Sequence                                            Score   
    [  3          TTGACATGTAAAAAAAAAAAAAAATATAAT                      12   ]
    [  54         ATGACAGATCATTGACAGATGATGGTGATGTAAA                  9    ]
    [  65         TTGACAGATGATGGTGATGTAAAGTATAAT                      12   ]

    :param dict possible_promoters: a dict of all the possible promoters in the DNA sequence
    '''
    print("\nPossible Promoters in the DNA sequence:\n")
    print('{:13s} {:50s}  {:8s}'.format('Position', "Sequence", "Score"))
    for index in range(0, len(possible_promoters["positions"])):
        print('{:2s} {:10s} {:50s}  {:3s} {:2s}'.format("[ ", str(possible_promoters["positions"][index]), 
            possible_promoters["sequences"][index], str(possible_promoters["scores"][index]), " ]"))


def print_each_sequence(input_dna_sequence, possible_promoters, position):
    '''
    This function prints the full DNA sequence, highlighting the user selected possible promoter
    in a different color, pink. The user can provide the position from a list of all identified possible
    promoter locations and this function will display the original DNA input with the selected possible promoter
    in it's original location, highlighted in pink.

    :param str input_dna_sequence: the string of DNA inputted from the user
    :param dict possible_promoters: a dict of all the possible promoters in the DNA sequence
    :param int position: the position of a possible promoter that the user chose to have dispayed
    '''
    if position in possible_promoters["positions"]:
        print("\nDNA Sequence: \n")

        # Print the DNA before the promoter
        print(input_dna_sequence[0: position], end="")

        # Print the possible promoter in pink
        position_index = possible_promoters["positions"].index(position)
        print('\033[95m' + input_dna_sequence[position: position + 
            len(possible_promoters["sequences"][position_index])] + '\033[0;0m', end="")

        # Print the rest of the DNA
        print(input_dna_sequence[(position + len(possible_promoters["sequences"][position_index])):])


def print_instructions():
    '''
    Displays a program overview and general instructions to the user when the program is first started. 
    It tells the user what the program does, restrictions to the DNA input, and instructions to begin.
    '''
    print("\n=====================================================================================")
    print('{:30s} {:50s}'.format(" ","Find Possible Promoters in DNA"))
    print("=====================================================================================\n")

    print("PROGRAM OVERVIEW:\n")
    print("This program reads the DNA downstream from the -35 (TTGACA) consensus sequence to the -10 (TATAAT) consensus sequence. "
        "The possible consensus sequences are assumed to be 6 nucleotides in length. It is assumed that a promoter sequence would " +
        "start with the -35 consensus, which is 6 nucleotides long and then have 16 - 19 nucleotides downstream until the " +
        "-10 consensus sequence, also 6 nucleotides long.\n")

    print("DNA INPUT STRUCTURE:\n") 
    print("Your DNA sequence should contain only Adenine(A), Thymine(T), Cytosine(C), or Guanine(G). " +
        "Nucelotide ambiguity nucleotides will not be accepted for this program.\n")


    print("INSTRUCTIONS:\n")
    print("In order to find any possible promoters, please enter either a text(txt) " +
        "file containing a DNA sequence or manually type your sequence when prompted.\n")

    print("You will have the option to find a possible promoter sequence that is within a specific threshold. The threshold " +
        "is measured by the count of matching nucleotides within the consensus sequences.\n")
    print("To get started, choose either file(F) or text(T) when prompted below then enter your the path to your file or " +
        "your DNA input.\n")


def main():
    '''
    The main function code that calls all of the other functions. Prints informations about the program and 
    provides inctructions to the user for how to use it. Asks the user whether to use a file or manually enter
    the single DNA input. If the input is incorrect, the user is prompted to enter the input again until it is 
    actually valid DNA input, containing only ACGT, no other values are accepted. The user is asked to enter a 
    threshold of the count of matching nucleotides. It should be in a range of 8-12 since a possible promoter
    consensus sequence typically differs in only 1-2 nucleotides, so the total count can only differ by up to 4 
    nucleotides and still count as a possible consensus sequence. A list of the possible promoters are printed out 
    with their position in the DNA input string and their match score. The user can then enter a position and that
    possible promoter is highlighted in pink in the original output. Each selected promoter sequence is printed out
    separately since it is possible to have overlapping possible promoter sequences adn it is easier to visualize 
    when they are outputted separately. To exit the program, the user can enter 'X' and the program will stop running.
    '''
    print_instructions()

    file_or_input = input("Do you want to enter a file or text[F/T]: ")
    if file_or_input.upper() == 'T':
        input_dna_sequence = input("Please enter a non template DNA string: ")
    elif file_or_input.upper() == 'F':
        file_name = input("Please enter a file path: ")
        file = open(file_name, 'r')
        input_dna_sequence = file.read()

    # Check that the input is valid DNA and if not prompt the user to input again until it is valid
    is_valid = check_is_dna(input_dna_sequence)
    while not is_valid:
        if file_or_input.upper() == 'T':
            input_dna_sequence = input("Please enter a non template DNA string: ")
        elif file_or_input.upper() == 'F':
            file_name = input("Please enter the file path to input with correct DNA data: ")
            file = open(file_name, 'r')
            input_dna_sequence = file.read()
        is_valid = check_is_dna(input_dna_sequence)


    if is_valid:
        threshold = input("Please enter a threshold in the range of 8-12 matching nucleotides: ")
        threshold = int(threshold)
        possible_promoters = find_possible_prokaryotic_sequences(input_dna_sequence, threshold)

        print_promoter_sequences(possible_promoters)

        position = input("\nEnter the position of one of the possible sequences to view it inside the DNA sequence (X to exit program): ")
        # Asks the user to choose a promoter's position to display the promoter in the DNA sequence until they enter 'X' to exit
        while position.upper() != "X":
            position = int(position)
            print_each_sequence(input_dna_sequence, possible_promoters, position)
            print_promoter_sequences(possible_promoters)
            position = input("\nEnter the position of one of the possible sequences to view it inside the DNA sequence (X to exit program): ")


if __name__ == "__main__":
    '''
    Calls the main function to start the program
    '''
    main()
	