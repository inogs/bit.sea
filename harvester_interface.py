class HarvesterInterface(object):
    def harvest(self, db_path, log):
        raise NotImplementedError

    def rebuild(self, db_path, log):
        raise NotImplementedError
