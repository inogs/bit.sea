# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
class HarvesterInterface(object):
    def harvest(self, db_path, log):
        raise NotImplementedError

    def rebuild(self, db_path, log):
        raise NotImplementedError
