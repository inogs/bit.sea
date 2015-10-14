class Profile(object): pass

class Instrument(object):
    def profiles(self, var):
        '''
        Return a list of all the profiles for a particular variable
        recorded by the instrument
        '''
        raise NotImplementedError
