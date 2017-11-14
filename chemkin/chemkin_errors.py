class ChemKinError(Exception):
    """ Encapsulates errors encountered in this module. """

    def __init__(self, method, info=None):
        msg = 'Error encountered in chemkin.py method: {0}.'.format(method)
        if info is not None:
            msg = msg + ' ' + info
        Exception.__init__(self, msg)

        self.method = method
        self.info = info
