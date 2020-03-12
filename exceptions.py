
class FullWindowException(Exception):
    """
    Exception to throw when attempting to
    add a new observation to an already full window
    """
    def __(self, size):
        self.message = "The Window size has already been reached. Window size: " + str(size)

    def __str__(self):
        return self.message
