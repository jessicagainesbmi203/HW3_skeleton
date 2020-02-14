
class Seq:
    """
    A class to hold protein sequences
    """
    def __init__(self, name):
        self.name = name
        self.sequence = ''

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name