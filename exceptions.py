WARNING="WARNING: "

class Error(Exception):
    """
    General error class to handle generic errors
    """
    def __init__(self, message):
        self.message = message

class FullWindowException(Exception):

    """
    Exception to throw when attempting to
    add a new observation to an already full window
    """
    def __(self, size):
        self.message = "The Window size has already been reached. Window size: " + str(size)

    def __str__(self):
        return self.message


