import sys, argparse
from graph import Automaton, Node
from exceptions import *

#*****************************************************************************
#                           T R I E _ M A T C H I N G                        *
#*****************************************************************************
#*****************************************************************************
# Name : Evalue()
# Function : This functions calculates e-value of a query.
#            E-value tells the expected number of matches you will find
#            by chance
# Parameters :
#           db : Database string
#           word : query string
# Return Value : interger ( e-value )
#*****************************************************************************
def Evalue( db, word ) :
    e_val = ( (1/4)**len(word) ) * len(db)
    return e_val

#*****************************************************************************
# Name : Create_Automaton()
# Function : This function builds a tree.
#            ( In future, hoping to incorporate the failure function. )
# Parameters :
#           keywords : list of queries
# Return Value : class Automaton
#*****************************************************************************
def Create_Automaton( keywords ) :
    # Initialize automaton
    automaton = Automaton()

    for eachKey in keywords :
        automaton.insert( automaton.root, eachKey )

    return automaton

#*****************************************************************************
# Name : Trie_Search()
# Function : It will search all the queries in the database and count them by
#            using the given automaton.
# Parameters :
#           db : Database
#           automaton : Trie
#           keycount : Dictionary with key (query) : value (count)
#           st_index : index of where to start in the database
# Return Value : dictionary 'keycount'
#*****************************************************************************
def Trie_Search( db, automaton, keycount, st_index ) :
    curr_node = automaton.root
    found = False
    word = ""
    frame_length = len( max( keycount.keys() ) )

    # Iterate over db
    for i in range( st_index, len(db) ) :
        for j in range( i, len(db) ) :
            # Search if there's a match
            for child in curr_node.children :
                if child.key == db[j] :
                    curr_node = child
                    found = True

            if found :
                found = False
                word = word + db[j]
                #print(word)
                # Check if it marks the end of the word
                if curr_node.isEnd :
                    print("FOUND : ", word)
                    keycount[word] = keycount[word] + 1
                    curr_node = automaton.root
                    word = ""
            else :
                word = ""
                curr_node = automaton.root
                break

    return keycount

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
    parser = argparse.ArgumentParser(usage = '%(prog)s [-d] [-o] <query file>',
                        description = "Finds all indices that it" +
                                      "matches from multiple patterns.")
    #Optional Arguments
    parser.add_argument('--file', '-d', dest = 'd', required = True,
                        nargs = '?', type = argparse.FileType('r'),
                        default = argparse.SUPPRESS, help = 'Database file.')
    parser.add_argument('--output', '-o', dest = 'o', required = False,
                        nargs = '?', type = argparse.FileType('w'),
                        default = argparse.SUPPRESS,
                        help = 'Write out to a file.')
    parser.add_argument('--index', '-i', dest = 'i', required = False,
                        nargs = '?', type = int, default = 0,
                        help = 'Index of where to start in database.')

    #Positional Arguments (User input parameters)
    parser.add_argument('file', type = argparse.FileType('r'),
                                help = 'File containing query sequences.')

    return parser.parse_args()

#*****************************************************************************
# Name : main()
# Function : To start off the program. Once the algorithm finishes, it
#            outputs the result to the terminal.
# Parameters : None
# Return Value : None
#*****************************************************************************
def main() :
    # Local Variables
    args = parse_args(sys.argv[1:])
    keywords = list()
    keycount = dict()
    db = ""
    index = 0

    # Parsing out all arguments - Optional Arguments
    if 'd' in args:
        with args.d as dbFile :
            for lines in dbFile :
                if lines[0] == ">" : continue
                else : db = db + lines.strip()
    if args.file :
        with args.file as queryFile :
            for word in queryFile : keywords.append(word.strip())

    # Parsing out all arguments - Positional Arguments
    if 'i' in args :
        try :
            if args.i > len(db) or 0 > args.i : raise IndexError(args.i)
            else : index = args.i
        except Exception as e :
            print("Index : ", e.value, " --> ", e.msg)
            sys.exit(1)

    # Setup dictionary
    for eachkey in keywords : keycount[eachkey] = 0

    # Setup Trie
    automaton = Create_Automaton( keywords )

    # Search
    result_dict = Trie_Search( db, automaton, keycount, index )

    # Print out Results to File
    if 'o' in args :
        with args.o as outfile :
            for k, v in result_dict.items() :
                #outfile.write( str(k) + " : " + str(v) + "\ne-value : " +
                #               str(Evalue( db, k )) + "\n" )
                outfile.write( str(k) + " : " + str(v) + "\n")
    else :
        for k, v in result_dict.items() :
            print( k, " : ", v, "\ne-value : ", Evalue( db, k ) )

#*****************************************************************************
# Description :
#   When python is executed with this file name, parse user's information
#   and pass it to main function.
#*****************************************************************************
if __name__ == '__main__' :
    main()
