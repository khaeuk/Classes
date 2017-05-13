import sys, argparse, numpy as np

#*****************************************************************************
#                        L O C A L _ A L I G N M E N T                       *
#*****************************************************************************
#*****************************************************************************
# Name : BackTrack()
# Function : Starts from the coordinate given by start_pt which holds the
#            maximum score within the matrix. It will look at that position's
#            direction ("", "m", "i", "d") and back track until empty
#            direction is met (Same as saying back track until reached 0).
# Parameters :
#       seq1 : First sequence.
#       seq2 : Second sequence.
#       start_pt : Maximum score's coordinate in matrix.
#       matrix : Matrix with calculated scores and directions.
# Return Value : List that holds the aligned sequences.
#*****************************************************************************
def BackTrack ( seq1, seq2, start_pt, matrix ) :
    #To store alignment sequences.
    seq1_aln = ""
    seq2_aln = ""

    curr_node = start_pt

    print("Now performing BackTrack...")
    while True :
        # Terminate alignment when score reached 0 (same as direction with "").
        if matrix[ curr_node[0], curr_node[1] ][1] == "" :
            break
        # Match / Mismatch --> Move diagonally (left, up)
        elif matrix[ curr_node[0], curr_node[1] ][1] == "m" :
            seq1_aln = seq1[curr_node[0]-1] + seq1_aln
            seq2_aln = seq2[curr_node[1]-1] + seq2_aln
            curr_node = (curr_node[0]-1, curr_node[1]-1)
        # Insertion --> Move left
        elif matrix[ curr_node[0], curr_node[1] ][1] == "i" :
            seq1_aln = "-" + seq1_aln
            seq2_aln = seq2[curr_node[1]-1] + seq2_aln
            curr_node = (curr_node[0], curr_node[1]-1)
        # Deletion --> Move up
        elif matrix[ curr_node[0], curr_node[1] ][1] == "d" :
            seq1_aln = seq1[curr_node[0]-1] + seq1_aln
            seq2_aln = "-" + seq2_aln
            curr_node = (curr_node[0]-1, curr_node[1])

    return [seq1_aln, seq2_aln]

#*****************************************************************************
# Name : LocalAlignment()
# Function : Takes in two sequences and peforms local alignment with given
#            user defined match(m), mismatch(s), indel(d) scores.
# Parameters :
#       sequences : List of two sequence to be aligned.
#       m_score : Match score.
#       s_score : Mismatch score.
#       d_score : Indel score.
#       a_flag : Output alignment ONLY if specified.
# Return Value :
#       List of three elements. [ [sequence_alignment], best score, length ]
#*****************************************************************************
def LocalAlignment( seq1, seq2, m_score, s_score, d_score, a_flag ) :
    # Create matrix of tuple of zeros.
    matrix = np.zeros( (len(seq1)+1, len(seq2)+1), dtype = object )
    # Dictionary that holds maximum score and it's coordinate.
    #   Key : (i,j) tuple coordinate   Value : maximum score
    nonzero_scores = {}
    # Stores [Diagonal( Match, Mismatch ), Insertion, Deletion]
    possible_dir = list()
    # List of information to return
    alignment_info = list()

    #DEBUGGING
    tot = len(seq1) * len(seq2)
    counter = 0

    # Initialize matrix with tuples
    for i in range( len(seq1)+1 ) :
        for j in range( len(seq2)+1 ) :
            matrix[i, j] = (0, "")

    print("Initialized Matrix...")
    print("Calculating Scores...")
    #DEBUGGING
    tempPercentage = 0
    print("0%")

    # Finding all scores in matrix
    for i in range( 1, len(seq1)+1 ) :
        #DEBUGGING
        percentage = int((counter * 100)/tot)
        if tempPercentage != percentage :
            if (percentage-tempPercentage) >= 10 :
                tempPercentage = percentage
                print(str(tempPercentage) + "%")

        for j in range( 1, len(seq2)+1 ) :
            #DEBUGGING
            counter = counter + 1

            #clear list
            del possible_dir[:]

            # Diagonal edges -- Match / Mismatch
            if seq1[i-1] == seq2[j-1] :
                possible_dir.append( (matrix[i-1, j-1][0] + m_score, "m") )
            else : possible_dir.append( (matrix[i-1, j-1][0] + s_score, "m") )

            # Left edges -- Insertion (indel)
            possible_dir.append( (matrix[i, j-1][0] + d_score, "i") )

            # Down edges -- Deletion (indel)
            possible_dir.append( (matrix[i-1, j][0] + d_score, "d") )

            score = max(possible_dir)

            # Update score ONLY if > 0 ; Store directions too.
            if int(score[0]) > 0 :
                matrix[i, j] = score
                nonzero_scores[ (i,j) ] = matrix[i, j][0]

    print("100%")
    print("Finished Calculating the Score...")
    # Find maximum score from the entire matrix
    score_max = max(nonzero_scores.values())

    # Get the coordinates for the maximum score.
    for coordinate, score in nonzero_scores.items() :
        if score == score_max : coordinate_max = coordinate

    # Backtrack from highest score until reached 0 to get local alignment.
    alignment = BackTrack( seq1, seq2, coordinate_max, matrix )

    print("Finished BackTracking...")

    # Add all alignment info to the list : best score, length
    alignment_info.append(alignment)
    alignment_info.append(score_max)
    alignment_info.append( len(alignment_info[0][0]) )

    # list of [ [alignment sequences], best score, length ]
    return alignment_info

#*****************************************************************************
# Name : parse_args()
# Function : Parses the arguments with given flags when running the program.
# Parameters/Flags :
#       Input file : FASTA formatted text file.
#       -m : User defined 'match' score.        (Default : 1)
#       -s : User defined 'mismatch' score.     (Default : -10)
#       -d : User defined 'indel' score.        (Default : -1)
#       -a : Output only the alignment itself.  (Default : True)
# Return Value :
#       parser.parse_args() - Returning 'ArgumentParser' object with info.
#*****************************************************************************
def parse_args(args) :
    parser = argparse.ArgumentParser(
                usage = '%(prog)s <seq_file> [-o] <output_file> [-m] <match> '
                        + '[-s] <mismatch> [-d] <indel> [-a]',
                description = 'Local Alignment program that outputs local '
                            + ' alignment of given two sequences.',
                epilog = 'Take a look at \"Readme.txt\" for more information.')
    #Optional Arguments (Flags)
    parser.add_argument('--match', '-m', dest = 'm', required = False,
                        nargs = '?', type = int,
                        default = 1, help = 'Match score.')
    parser.add_argument('--mismatch', '-s', dest = 's', required = False,
                        nargs = '?', type = int,
                        default = -10, help = 'Mismatch score.')
    parser.add_argument('--indel', '-d', dest = 'd', required = False,
                        nargs = '?', type = int,
                        default = -1, help = 'Indel score.')
    parser.add_argument('--alignment', '-a', dest = 'a', required = False,
                        action = 'store_true', help = 'Indel score.')
    parser.add_argument('--output', '-o', dest = 'o', required = False,
                        nargs = '?', type = argparse.FileType('w'),
                        default = argparse.SUPPRESS,
                        help = 'Write out to a file.')

    #Positional Arguments (User input parameters)
    parser.add_argument('file', type = argparse.FileType('r'),
                                help = 'FASTA formatted text file.')

    return parser.parse_args()
#*****************************************************************************
# Name : main()
# Function : To start off the program. Once the algorithm finishes, it
#            outputs the result to the terminal.
# Parameters :
#       sequences : list of two sequences to be aligned.
#       -m : User defined 'match' score.        (Default : 1)
#       -s : User defined 'mismatch' score.     (Default : -10)
#       -d : User defined 'indel' score.        (Default : -1)
#       -a : Output only the alignment itself.  (Default : True)
# Return Value : None
#*****************************************************************************
def main() :
    args = parse_args(sys.argv[1:])
    sequences = list()

    # Read file and store only sequences.
    with args.file as f :
        for line in f :
            if line[0] != '>' and line.strip() :
                    sequences.append(line.strip())

    alignment_result = LocalAlignment( sequences[0], sequences[1],
                                       args.m, args.s, args.d, args.a )

    # Printing out result to the console or to a file.
    if 'o' in args :
        with args.o as outfile :
            if args.a :
                outfile.write("Printing out the Alignment : ")
                outfile.write(alignment_result[0][0] + "\n")
                outfile.write(alignment_result[0][1])
            else :
                outfile.write("Best Score of Local Alignment : "
                              + str(alignment_result[1]) + "\n")
                outfile.write("Length of Best Local Alignment : "
                              + str(alignment_result[-1]) + "\n")
    else :
        if args.a :
            print("Printing out the Alignment : ")
            print(alignment_result[0][0] + "\n")
            print(alignment_result[0][1])
        else :
            print("Best Score of Local Alignment : ", alignment_result[1])
            print("Length of Best Local Alignment : ", alignment_result[-1])


    print("Program Finished.")

#*****************************************************************************
# Description :
#   When python is executed with this file name, parse user's information
#   and pass it to main function.
#*****************************************************************************
if __name__ == '__main__' :
    main()
