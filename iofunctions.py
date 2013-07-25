#########################
#
# Some I/O functions
#
# Written by Marcel Zemp
#
#########################

import numpy

# General data reading function

def get_data(FileName):

    file = open('%s'%FileName,'r')

    values = []
    for line in file.readlines():
        if (line.startswith('#')): continue
        values.append([float(value) for value in line.split()])
        
    file.close()

    return numpy.array(values)

# Functions to get data from characteristics file

def get_characteristics_via_ID(FileName,ID):

    file = open('%s'%FileName,'r')
    
    values = []
    while (1):
        line = file.readline()
        if not line: break
        if (line.startswith('#')): continue
        if (ID == int(line.split()[0])):
            values.append([float(value) for value in line.split()])
            break

    file.close()

    return numpy.array(values[0])

def get_line_indices(FileName,d,ID):

    file = open('%s'%FileName,'r')
    
    N = 0
    Nold = 0
    while (1):
        line = file.readline()
        if not line: break
        if (line.startswith('#')): continue
        Nold = N
        N += int(line.split()[9])
        if (d == 2):
            Nz = int(line.split()[12])
            assert(Nz%2 == 0)
            Nz /= 2

        if (ID == int(line.split()[0])): break

    file.close()

    Iupper = N+1
    Ilower = Nold+2

    if (d == 2):
        return Ilower, Iupper, Nz
    else:
        return Ilower, Iupper

# Functions to get profile data

def get_profiles_1d_via_ID(FileName,ID,i1l,i1u):

    file = open('%s'%FileName,'r')

    lines = []
    for i in range(i1u):
        line = file.readline()
        if i >= i1l-1: lines.append(line)

    file.close()
    
    values = []
    for line in lines:
        assert (ID == int(line.split()[0]))
        values.append([float(value) for value in line.split()])

    return numpy.array(values)

def get_profiles_2d_via_ID(FileNameBegin,FileNameEnd,ID,i1l,i1u,Nz):

    M = []
    for i in range(Nz,0,-1):
        FileName = '%s.m%03d.%s'%(FileNameBegin,i,FileNameEnd)
        M.append(get_profiles_1d_via_ID(FileName,ID,i1l,i1u))
    for i in range(1,Nz+1):
        FileName = '%s.p%03d.%s'%(FileNameBegin,i,FileNameEnd)
        M.append(get_profiles_1d_via_ID(FileName,ID,i1l,i1u))

    return numpy.array(M)

def get_ID_etc_and_check(Ccyl,Csph):

    ID = int(Ccyl[0])
    NBin = int(Ccyl[9])
    Nz = int(Ccyl[12])
    assert(Nz%2 == 0)
    Nz /= 2
    assert(ID == int(Csph[0]))

    return ID, NBin, Nz

def get_profiles_2d(FileList,ID,NBin,Nz):

    assert(len(FileList) == 2*Nz)
    M = []
    for File in FileList:
        values = []
        for b in range(NBin):
            profileline = File.readline()
            values.append([float(value) for value in profileline.split()])
            assert(ID == int(values[-1][0]))
        M.append(numpy.array(values))
        
    return numpy.array(M)

def get_profiles_2d_combined(FileListList,ID,NBin,Nz):

    for mt,FileList in zip(range(len(FileListList)),FileListList):
        if (mt == 0):
            PcylCombined = copy(get_profiles_2d(FileList,ID,NBin,Nz))
        else:
            TempPcyl = get_profiles_2d(FileList,ID,NBin,Nz)
            for iz in range(2*Nz):
                for ir in range(NBin):
                    assert(PcylCombined[iz,ir,0] == TempPcyl[iz,ir,0]) # ID
                    assert(PcylCombined[iz,ir,1] == TempPcyl[iz,ir,1]) # ri_1
                    assert(PcylCombined[iz,ir,2] == TempPcyl[iz,ir,2]) # rm_1
                    assert(PcylCombined[iz,ir,3] == TempPcyl[iz,ir,3]) # ro_1
                    assert(PcylCombined[iz,ir,4] == TempPcyl[iz,ir,4]) # ri_2
                    assert(PcylCombined[iz,ir,5] == TempPcyl[iz,ir,5]) # rm_2
                    assert(PcylCombined[iz,ir,6] == TempPcyl[iz,ir,6]) # ro_2
                    PcylCombined[iz,ir,7] += TempPcyl[iz,ir,7] # M
                    assert(PcylCombined[iz,ir,8] == TempPcyl[iz,ir,8]) # N

    return PcylCombined
