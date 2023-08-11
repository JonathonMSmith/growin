import numpy as np

def sami3data(sami3file,sz):
    fid= open(sami3file,'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    if len(sz)==4:
       tempdata1 = np.reshape(tempdata,(sz[0]*sz[1]*sz[2]+2,sz[3]),order='F')
       xx = tempdata1[1:-1,:]
    else:
       xx = tempdata[1:-1]
    samidata= np.reshape(xx,sz,order='F')
    return samidata

def sami3data_grid(sami3file,sz):
    fid= open(sami3file,'rb')
    data = fid.read()
    tempdata = np.frombuffer(data, np.float32)
    fid.close()
    xx = tempdata[1:-1]
    samidata= np.reshape(xx,sz,order='F')
    return samidata
