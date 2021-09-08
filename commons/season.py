from datetime import datetime
from dateutil.relativedelta import relativedelta
from commons.time_interval import TimeInterval
#author : epascolo

class season:
    """
    Season Object

    This object provide a tools for season definition and management
    The default setting is the classical oceanographic 4 season, the start day is set as:
    - winter : 0101
    - spring : 0401
    - summer : 0701
    - fall   : 1001
    """


    def __init__(self):
        """
        Init object and set astronomical season
        """
        self.numbers_season = 0
        self.SEASON_LIST = []
        self.SEASON_LIST_NAME = []
        self._reference_year=2000
        self.setseasons(["0101","0401","0701","1001"],["winter","spring","summer","fall"])
        #self.setseasons(["1221","0321","0622","0921"],["winter","spring","summer","fall"])
    

    def setseasons(self,startseason,nameseason):
        """

        Given two arrays where is defined the date when season start and its name
        the subroutine generate the season list. In input take two arrays of string:

        - Example of startseason : ["1221","0321","0622","0921"]
        - Example of nameseason  : ["winter","spring","summer","fall"]

        In the previous example we shown that winter start to 21 december,
        spring to 21 march, summer to 22 june and fall on 21 september.

        """
        self.SEASON_LIST = []
        self.SEASON_LIST_NAME = []
        self.numbers_season = len(startseason)

        if (len(startseason) != len(nameseason)):
            print("ERROR : arrays definitions mistmatch! :-( ")
            exit()

        for i in range(0,self.numbers_season):
            if i==0:
                if startseason[i] == '0101':
                    ref_year = self._reference_year
                else:
                    ref_year = self._reference_year-1
            else:
                ref_year=self._reference_year

            self.SEASON_LIST.append(datetime.strptime(str(ref_year)+startseason[i],'%Y%m%d'))
            self.SEASON_LIST_NAME.append(nameseason[i])

    def get_seasons_number(self):
        """
        Return the number of seasons defined in this object
        """
        return self.numbers_season

    def get_season_dates(self,season_num):
        """
        Given season number, return the range of season dates (start and end)
        and the name of season.
        """

        assert season_num  < self.numbers_season
        start_date=self.SEASON_LIST[season_num]
        if (season_num + 1) == self.numbers_season:
            end_date = self.SEASON_LIST[0] + relativedelta(years = 1)
        else:
            end_date= self.SEASON_LIST[season_num+1]
        TI = TimeInterval.fromdatetimes(start_date, end_date)

        season_name = self.SEASON_LIST_NAME[season_num]

        return TI,season_name

    def findseason(self,date):
        """
        Takes a date as input and return the number and name of season where it is in.

        """
        yearly_date = datetime(self._reference_year, date.month, date.day, date.hour, date.minute, date.second)
        for season_num in range(self.numbers_season):
            ti,_ = self.get_season_dates(season_num)
            if ti.contains(yearly_date):
                return season_num
        yearly_date = datetime(self._reference_year-1, date.month, date.day, date.hour, date.minute, date.second)
        for season_num in range(self.numbers_season):
            ti,_ = self.get_season_dates(season_num)
            if ti.contains(yearly_date):
                return season_num
