

class SpatialGP(object):
    """docstring for SpatialGP."""
    # TODO should Y include all genes or just noe at a time
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y

        self.N = self.Y.shape[0]
        self.P = self.Y.shape[1]

        assert self.X.shape[0] == self.N, 'dimension missmatch'
