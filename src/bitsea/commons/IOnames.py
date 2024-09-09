from xml.dom import minidom
import datetime

class filenamer():
    def __init__(self,prefix='',dateformat='',suffix=''):
        self.prefix=prefix
        self.dateformat=dateformat
        self.suffix=suffix
        lendate = len(datetime.datetime(2000,1,1).strftime(dateformat))
        self.date_startpos = len(self.prefix)
        self.date_endpos   = self.date_startpos+lendate
    def __str__(self):
        return "File name type :  " + self.prefix + self.dateformat + self.suffix


class IOnames():
    def __init__(self,xmlfilename='IOnames.xml'):
        xmldoc = minidom.parse(xmlfilename)
        
        NODE=xmldoc.getElementsByTagName("input")[0]
        prefix =  str(NODE.getAttribute('prefix'))
        dateformat =  str(NODE.getAttribute('dateformat'))
        suffix = str(NODE.getAttribute('suffix'))
        
        self.Input = filenamer(prefix,dateformat,suffix)
        
        NODE=xmldoc.getElementsByTagName("output")[0]
        prefix =  str(NODE.getAttribute('prefix'))
        dateformat =  str(NODE.getAttribute('dateformat'))
        suffix = str(NODE.getAttribute('suffix'))
        
        self.Output = filenamer(prefix,dateformat,suffix)