import pickle
import numpy as np

class container():
    def __init__(self, I,J,values, QI_old, QI_new):
        self.I = I
        self.J = J
        self.values = values
        self.QI_old = QI_old
        self.QI_new = QI_new

fid = open('rejected_only_old.pkl','rb')
DAILY_REJECT_ONLY_OLD=pickle.load(fid)
fid.close()


fid = open('rejected_only_new.pkl','rb')
DAILY_REJECT_ONLY_NEW = pickle.load(fid)
fid.close()

iFrame_start=0
iFrame_end=30

I = np.array([],int)

for C in DAILY_REJECT_ONLY_OLD:
    I = np.concatenate((I,C.I))



