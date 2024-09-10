from xml.dom import minidom
import os,sys

class read_descriptor():
    def __init__(self, VarDescriptorFile):
        
        xmldoc = minidom.parse(VarDescriptorFile)
        
            
        self.NATIVE_VARS=set()
        NODE=xmldoc.getElementsByTagName("vars_for_All_Statistics")[0].getElementsByTagName("native")
        MONTHLY_NATIVE_NODES=NODE[0].getElementsByTagName("var")
        for n in MONTHLY_NATIVE_NODES:
            self.NATIVE_VARS.add( str(n.attributes['name'].value) )
            
                
        self.AGGR_VARS=[]
        NODE=xmldoc.getElementsByTagName("vars_for_All_Statistics")[0].getElementsByTagName("aggregate")
        MONTHLY_NATIVE_NODES=NODE[0].getElementsByTagName("var")
        for n in MONTHLY_NATIVE_NODES:
            self.AGGR_VARS.append(str(n.attributes['name'].value))
            
  
        
        self.SOME_VARS=set()    
        NODE=xmldoc.getElementsByTagName("var_for_Some_Statistics")
        MONTHLY_NATIVE_NODES=NODE[0].getElementsByTagName("var")
        for n in MONTHLY_NATIVE_NODES:
            self.SOME_VARS.add( str(n.attributes['name'].value) )
        
        
        self.LIST_To_Online_PostPROC = self.NATIVE_VARS.union(self.SOME_VARS).union(set(self.AGGR_VARS))  # ottengo la longlist, senza ripetizione
        self.ARCHIVE_VARS=[]
        NODE=xmldoc.getElementsByTagName("toArchive")
        MONTHLY_NATIVE_NODES=NODE[0].getElementsByTagName("var")
        for n in MONTHLY_NATIVE_NODES:
            self.ARCHIVE_VARS.append( str(n.attributes['name'].value) )
            
        
        self.vars2D=set()
        NODE=xmldoc.getElementsByTagName("var2D")[0].getElementsByTagName("vars_for_All")[0].getElementsByTagName("native")
        NATIVE_NODES=NODE[0].getElementsByTagName("var")
        for n in NATIVE_NODES:
            self.vars2D.add(str(n.attributes['name'].value))
        
        self.AGGREGATE_VARS=[]
        self.AGGR_FORMULAS=[]        
        AGGVARNODES=xmldoc.getElementsByTagName("vars_for_All_Statistics")[0].getElementsByTagName("aggregate")[0].getElementsByTagName("aggvar")
        for NODE in AGGVARNODES:
            aggvar=str(NODE.attributes['name'].value)
            self.AGGR_FORMULAS.append(str(NODE.attributes['formula'].value))       
#             N=NODE.getElementsByTagName("var")
#             L=[]
#             for n in N:
#                 L.append( str(n.attributes['name'].value))
            self.AGGREGATE_VARS.append(aggvar)

    def printStatus(self):
        print(len(self.NATIVE_VARS), "native vars")
        print(len(self.AGGR_VARS), " vars which we aggregate")
        print(len(self.SOME_VARS), " vars for Some Statistics")
        print(len(self.LIST_To_Online_PostPROC), "variables for online postproc")
        
                
    def printError(self,group, var):
        print( "*** ERROR in ", group ," variables: ", var , "not defined in model. ")
        sys.exit(1) 
    def printFreqError(self,group, var):
        print( "*** ERROR in ", group ," variables: ", var , "not dumped at high frequency. ")
        sys.exit(1)         
                    
    
    def check(self,group):
        if group ==1: self.check_highfreq()
        
    def check_highfreq(self):
        repeatition=self.SOME_VARS.intersection(self.NATIVE_VARS)
        for i in repeatition :
            print("**** removing", i , "from var_for_Some_Statistics")
            self.SOME_VARS.remove(i)
                    
        SET_ARCHIVE_VARS=set(self.ARCHIVE_VARS)
        if len(SET_ARCHIVE_VARS) != len(self.ARCHIVE_VARS):
            print("duplication in toArchive vars")
            for var in SET_ARCHIVE_VARS:
                if self.ARCHIVE_VARS.count(var) > 1 :
                    print(var, "is duplicated. REMOVE it from xml descriptor ***************")
        # link with MODEL directory
        MODELDIR = os.getenv("CINECA_SCRATCH") + "/" + os.getenv("OPA_HOME") + "/wrkdir/MODEL/"
        if os.path.exists(MODELDIR):
        
            F=file(MODELDIR + "namelist.passivetrc")
            MODELVARS=[]
            for line in F:
                if line.find("ctrcnm") != -1:         
                    apex_1 = line.find("'")
                    MODELVARS.append(line[apex_1+1:apex_1+4])
            F.seek(0)
            HIGHFREQ_MODELVARS=[]
            for line in F:
                if line.find("highfreq_tra") != -1:                     
                    pos = line.find("=")
                    if line[pos+2] == "1":
                        pos1 = line.find("(")
                        pos2 = line.find(")")
                        ind  = int( line[pos1+1:pos2] )
                        
                        HIGHFREQ_MODELVARS.append(MODELVARS[ind-1])
                                
            F.close()
            
            DIA_VARS=[]
            F=file(MODELDIR + "namelist.diagnostics")
            for line in F:
                if line.find("dia") != -1:         
                    apex_1 = line.find("'")
                    DIA_VARS.append(line[apex_1+1:apex_1+4])
            F.seek(0)
            HIGHFREQ_DIAVARS=[]
            for line in F:
                if line.find("highfreq_tra_dia") != -1:                     
                    pos = line.find("=")
                    if line[pos+2] == "1":
                        pos1 = line.find("(")
                        pos2 = line.find(")")
                        ind  = int( line[pos1+1:pos2] )
                        HIGHFREQ_DIAVARS.append(DIA_VARS[ind-1])       
            F.close()
            
            
            MODELVARS.extend(DIA_VARS)
            HIGHFREQ_MODELVARS.extend(HIGHFREQ_DIAVARS)
            
            for var in self.NATIVE_VARS:
                if var not in MODELVARS: self.printError("native", var)                        
                if var not in HIGHFREQ_MODELVARS: self.printFreqError("native", var)
                    
            for var in self.AGGR_VARS:
                if var not in MODELVARS: self.printError("aggregate", var)
                if var not in HIGHFREQ_MODELVARS: self.printFreqError("aggregate", var)
            
            for var in self.ARCHIVE_VARS:
                if var not in MODELVARS: self.printError("archived", var) 
                if var not in HIGHFREQ_MODELVARS: self.printFreqError("archived", var)
                    
            for var in self.SOME_VARS:
                if var not in MODELVARS: self.printError("Some statistics", var) 
                if var not in HIGHFREQ_MODELVARS: self.printFreqError("Some statistics", var)
            
        #    for var in MODELVARS:
        #        if var not in LIST_To_Online_PostPROC: print var ,"produced by model but not processed"



                            
    def dumpfiles(self):
        F=file('longlist_1',"w")
        for line in self.LIST_To_Online_PostPROC:
            F.write(line + "\n")
        F.close()
        F=file('varStoredInAve_1',"w")
        for line in self.ARCHIVE_VARS:
            F.write(line + "\n")
        F.close()


