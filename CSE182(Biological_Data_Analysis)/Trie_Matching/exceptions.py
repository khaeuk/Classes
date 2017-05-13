import constants

class IndexError(Exception) :
    def __init__ ( self, val ) :
        self.value = val
        self.msg = constants.ERR_IND_MSG

    def __str__ ( self ) :
        return repr(self.value)
