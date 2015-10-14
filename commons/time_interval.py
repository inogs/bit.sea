import datetime

class TimeInterval():
    def __init__(self, starttime="19500101", endtime="21000101", dateformat='%Y%m%d'):
        self.start_time = datetime.datetime.strptime(starttime, dateformat)
        self.end_time = datetime.datetime.strptime(endtime   ,dateformat)

    def contains(self, specific_time):
        if (specific_time >= self.start_time) and (specific_time < self.end_time):
            return True
        else:
            return False
