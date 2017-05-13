
class Automaton :
    # Constructor
    def __init__ ( self ) :
        self.root = Node("Root")

    def insert ( self, n, keyword ) :
        # Local Variables
        pre = n
        found = False

        # Base Case
        if len(keyword) == 0 :
            pre.isEnd = True
            return

        # If node has no children
        if not pre.children :
#            print("Root has no children")
#            print("Creating node ", keyword[0])
            node = Node( keyword[0] )
            node.state = pre.state + 1                 # Set state of node
            pre.appendNode( node )
            pre = node
#            print("Now currently on node ", pre.key)
            self.insert( pre, keyword[1:] )

        else :
            # Check if keyword in node's children
            for child in pre.children :
                # If exists, move pointer
                if child.key == keyword[0] :
#                    print("Moving pointer from %s to %s " %(pre.key, child.key))
                    pre = child
                    self.insert( pre, keyword[1:] )
                    found = True

            # If not in children
            if not found :
#                print("Creating node ", keyword[0])
                node = Node( keyword[0] )
                node.state = pre.state + 1              # Set state of node
#                print("Connecting node %s to %s" %(pre.key, node.key))
                pre.appendNode( node )
                pre = node
#                print("Now currently on node ", pre.key)
                self.insert( pre, keyword[1:] )

class Node :
    # Constructor
    def __init__ ( self, key ) :
        self.key = key
        self.children = list()   # Holds list of Nodes
        self.state = 0
        self.isEnd = False

    def appendNode( self, node ) :
        self.children.append(node)
