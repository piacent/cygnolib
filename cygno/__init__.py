# Copyright 2021 dciangot
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Image analisys CYGNO Python Library
# G. Mazzitelli 2017 
# rev oct 2018 - swift direct access 
# rev Oct 2021 clean up and packege
# 

__version__ = '1.0.1'
__all__     = ["cmd", "his", "s3"]


import numpy as np
import glob, os
import re
import sys
from cygno import s3
from cygno import cmd

class myError(Exception):
    pass

__path__ = os.path.dirname(os.path.realpath(__file__))+'/'
## IMPORTING THE CORRECTION TABLES FOR PMT WAVEFORM CORRECTION
__table_path__ = __path__ + 'pmt_correction_tables/'
# LNGS
if(os.path.exists(__table_path__+'table_cell_LNGS.npy')):
    __table_cell_LNGS__ = np.load(__table_path__+'table_cell_LNGS.npy')
else: raise myError('table_cell_LNGS.npy not found')
    
if(os.path.exists(__table_path__+'table_nsample_LNGS.npy')):
    __table_nsample_LNGS__ = np.load(__table_path__+'table_nsample_LNGS.npy')
else: raise myError('table_nsample_LNGS.npy not found')
# LNF
if(os.path.exists(__table_path__+'table_cell_LNF.npy')):   ### NOT PRESENT YES
    __table_cell_LNF__ = np.load(__table_path__+'table_cell_LNF.npy')
else: raise myError('table_cell_LNF.npy not found')
if(os.path.exists(__table_path__+'table_nsample_LNF.npy')):
    __table_nsample_LNF__ = np.load(__table_path__+'table_nsample_LNF.npy')
else: raise myError('table_nsample_LNF.npy not found')
# MAN
if(os.path.exists(__table_path__+'table_cell_MAN.npy')):
    __table_cell_MAN__ = np.load(__table_path__+'table_cell_MAN.npy')
else: raise myError('table_cell_MAN.npy not found')
if(os.path.exists(__table_path__+'table_nsample_MAN.npy')):
    __table_nsample_MAN__ = np.load(__table_path__+'table_nsample_MAN.npy')
else: raise myError('table_nsample_MAN.npy not found')


#
# CYGNO py ROOT Tools
#


class cfile:
    def __init__(self, file, pic, wfm, max_pic, max_wfm, x_resolution, y_resolution):
        self.file         = file
        self.pic          = pic 
        self.wfm          = wfm
        self.max_pic      = max_pic
        self.max_wfm      = max_wfm
        self.x_resolution = x_resolution
        self.y_resolution = y_resolution

def open_mid(run, path='/tmp/',  cloud=True, tag='LNGS', verbose=False):
    import midas.file_reader
    fname = s3.mid_file(run, tag=tag, cloud=cloud, verbose=verbose)
    if verbose: print(fname)
    if not cloud:
        if os.path.exists(path+tag+fname):
            f = midas.file_reader.MidasFile(path+tag+fname)
        else:
            raise myError("openFileError: "+path+tag+fname+" do not exist") 
    else:
        filetmp = cmd.cache_file(fname, cachedir=path, verbose=verbose)
        f = midas.file_reader.MidasFile(filetmp)  
    return f

def open_root(run, path='/tmp/',  cloud=True,  Bari=False,  tag='LAB', verbose=False):
    import ROOT
    import root_numpy as rtnp
    fname = s3.root_file(run, tag=tag, cloud=cloud, Bari=Bari, verbose=verbose)
    if not cloud:
        fname=path+fname
    class cfile:
        def __init__(self, file, pic, wfm, max_pic, max_wfm, x_resolution, y_resolution):
            self.file         = file
            self.pic          = pic 
            self.wfm          = wfm
            self.max_pic      = max_pic
            self.max_wfm      = max_wfm
            self.x_resolution = x_resolution
            self.y_resolution = y_resolution
    try:
        f=ROOT.TFile.Open(fname)
        pic, wfm = rootTH2byname(f)
        image = rtnp.hist2array(f.Get(pic[0])).T
        x_resolution = image.shape[1]
        y_resolution = image.shape[0]
        max_pic = len(pic)
        max_wfm = len(wfm)
    except:
        raise myError("openFileError: "+fname)

    if verbose:
        print ('Open file: '+fname)
        print ('Find Keys: '+str(len(f.GetListOfKeys())))
        print ("# of Images (TH2) Files: %d " % (max_pic))
        print ("# of Waveform (TH2) Files: %d " % (max_wfm))
        print ('Camera X, Y pixel: {:d} {:d} '.format(x_resolution, y_resolution))
    return cfile(f, pic, wfm, max_pic, max_wfm, x_resolution, y_resolution)

    
def open_(run, tag='LAB', posix=False, verbose=False):
    BAKET_POSIX_PATH = '/jupyter-workspace/cloud-storage/cygno-data/'
    if posix:
        path = BAKET_POSIX_PATH+tag+'/'
    else:
        path='/tmp/'
    if tag == 'LNGS' or tag == 'LNF' or tag =='TMP':
        f = open_mid(run, path=path,  cloud=posix,  tag=tag, verbose=False)
    else:
        f = open_root(run, path=path,  cloud=posix,  tag=tag, verbose=False)
    return f

def daq_cam2array(bank):
    test_size=bank.size_bytes/10616832      #10616832=2*5308416   5308416=2304*2304 and 8/16=1/2 (banksize in bytes*8 returns the bits and divided by the 16 bit adc gives the total amount of pixels)
    if test_size<1.1:
        #Fusion,Flash
        shapex = shapey = int(np.sqrt(bank.size_bytes/2))
    else:
        #Quest
        shapex=4096
        shapey=2304

    image = np.reshape(bank.data, (shapey, shapex))
    return image, shapex, shapey

def get_bor_odb(mfile): # function to acquire the begin of run ODB entries from the midas file
    try:
        odb = mfile.get_bor_odb_dump()
    except:
        myError("No begin-of-run ODB dump found")
    
    mfile.jump_to_start()
    return odb

def daq_dgz2header(bank, verbose=False):
    nboard = bank.data[0]
    ich = 1
    for iboard in range(nboard):
        name_board = bank.data[ich]
        ich+=1
        number_samples = bank.data[ich]
        ich+=1
        number_channels =  bank.data[ich]
        ich+=1
        number_events = bank.data[ich]
        ich+=1
        vertical_resulution = bank.data[ich]
        ich+=1
        sampling_rate = bank.data[ich]
        if verbose:
            print ("name_board, number_samples, number_events, vertical_resulution, sampling_rate", 
                   name_board, number_samples, number_events, vertical_resulution, sampling_rate)
        cannaels_offset = [None] * number_channels
        for ichannels in range(number_channels):
            ich+=1
            cannaels_offset[ichannels] = bank.data[ich]
        if verbose:
            print ("cannaels_offset: ", cannaels_offset)
        return number_events, number_channels, number_samples

def daq_dgz2array(bank, header, verbose=False):
    waveform = []
    data_offset = 0
    number_events  = header[0]
    number_channels= header[1]
    number_samples = header[2]
    for ievent in range(number_events):       
        for ichannels in range(number_channels):
            if verbose:
                print ("data_offset, data_offset+number_samples",
                       data_offset, data_offset+number_samples)
                print(bank.data[data_offset:data_offset+number_samples])

            waveform.append(bank.data[data_offset:data_offset+number_samples])
            data_offset += number_samples
    if verbose:
        print(waveform, number_events, number_channels)
    return waveform

class dgtz_header:      # very simple class for the dgtz header
    def __init__(self, a):
        self.ntriggers           = a[0]
        self.nchannels           = a[1]
        self.nsamples            = a[2]
        self.vertical_resulution = a[3]
        self.sampling_rate       = a[4]
        self.offsets             = a[5]
        self.TTT                 = a[6]
        self.SIC                 = a[7]
        if len(a)>8:
            self.nBoards         = a[8]
            self.boardNames      = a[9]
            self.DAQversion      = a[10]
            self.pattern         = a[11]
            
        else:
            self.nBoards        = 1
            self.boardNames     = [1742]
            self.DAQversion     = 0
            self.pattern        = []
            print('WARNING: You are using an older version of the data bank. The analysis of the digitizers might be incomplete.')

        self.itemDict = {}
        self.itemDict["0"]  = self.ntriggers
        self.itemDict["1"]  = self.nchannels
        self.itemDict["2"]  = self.nsamples
        self.itemDict["3"]  = self.vertical_resulution
        self.itemDict["4"]  = self.sampling_rate
        self.itemDict["5"]  = self.offsets
        self.itemDict["6"]  = self.TTT
        self.itemDict["7"]  = self.SIC
        self.itemDict["8"]  = self.nBoards
        self.itemDict["9"]  = self.boardNames
        self.itemDict["10"] = self.DAQversion
        self.itemDict["11"] = self.pattern
    
    def __getitem__(self, index):
        return self.itemDict[str(int(index))]
        
def daq_dgz_full2header(bank, verbose=False):
    # v0.1 full PMT recostruction
    import numpy as np

    if bank.data[0] > 1000: # if DAQ version from CYGNO-04
        DAQversion          = bank.data[0]
        nboard              = bank.data[1]

        full_buffer_size    = len(bank.data)

        if verbose: print("Number of board: {:d}".format(nboard))
        board_index         = bank.data[2]
        name_board          = bank.data[3]
        number_samples      = bank.data[4]
        number_channels     = bank.data[5]
        number_events       = bank.data[6]
        vertical_resulution = bank.data[7]
        sampling_rate       = bank.data[8]

        if verbose:
                print ("board: {:d}, name_board: {:d}, number_samples: {:d}, number_channels: {:d}, number_events: {:d}, vertical_resulution: {:d}, sampling_rate: {:d}, DAQversion: {:d}".format(
                    board_index, name_board, number_samples, number_channels,
                    number_events, vertical_resulution, sampling_rate, DAQversion))
            
        ######### Channels offset reading:
        channels_offset  = bank.data[9:(9+number_channels)]

        ######### TTT reading:
        channels_ttt     = bank.data[9 + number_channels : 9 + number_channels + number_events]

        ######### Pattern reading:
        channels_pattern = bank.data[9 + number_channels + number_events : 9 + number_channels + 2 * number_events]
        
        ######### SIC reading:
        channels_SIC = []
        if name_board == 1742:
            channels_SIC = bank.data[9 + number_channels + 2 * number_events : 9 + number_channels + 3 * number_events]
            
    
    else: # if previous versions
        DAQversion          = 0
        nboard              = bank.data[0]
        
        full_buffer_size    = len(bank.data)
        name_board          = np.empty([nboard], dtype=int)
        number_samples      = np.empty([nboard], dtype=int)
        number_channels     = np.empty([nboard], dtype=int)
        number_events       = np.empty([nboard], dtype=int)
        vertical_resulution = np.empty([nboard], dtype=int)
        sampling_rate       = np.empty([nboard], dtype=int)
        channels_offset     = []
        channels_ttt        = []
        channels_pattern    = []
        channels_SIC        = []
    
        if verbose: print("Number of board: {:d}".format(nboard))
        
        ich=0
        for iboard in range(nboard): ######### cicle over the boards
            ich+=1
            name_board[iboard]          = bank.data[ich]
            ich+=1  
            number_samples[iboard]      = bank.data[ich]
            ich+=1  
            number_channels[iboard]     = bank.data[ich]
            ich+=1  
            number_events[iboard]       = bank.data[ich]
            ich+=1  
            vertical_resulution[iboard] = bank.data[ich]
            ich+=1  
            sampling_rate[iboard]       = bank.data[ich]

            if verbose:
                print ("board: {:d}, name_board: {:d}, number_samples: {:d}, number_channels: {:d}, number_events: {:d}, vertical_resulution: {:d}, sampling_rate: {:d}".format( 
                       iboard, name_board[iboard], number_samples[iboard], number_channels[iboard], number_events[iboard], vertical_resulution[iboard], sampling_rate[iboard]))
        
            ######### Channels offset reading:
            channels_offset_tmp = np.empty([number_channels[iboard]], dtype=int)
            for ichannels in range(number_channels[iboard]):
                ich+=1
                channels_offset_tmp[ichannels] = bank.data[ich]
            if verbose:
                print ("channels_offset: ", channels_offset_tmp, flush=True)
            channels_offset.append(channels_offset_tmp)
        
            ######### TTT reading:
            channels_ttt_tmp = np.empty(number_events[iboard], dtype=int)
            for ttt in range(number_events[iboard]):
                ich+=1
                channels_ttt_tmp[ttt] = bank.data[ich]
            if verbose:
                print ("channels_ttt: ", channels_ttt)
            channels_ttt.append(channels_ttt_tmp)
        
            ######### Start Index Cell reading:  
            if name_board[iboard] == 1742:   
                channels_SIC_tmp = np.empty(number_events[iboard], dtype=int)
                for sic in range(number_events[iboard]):
                    ich+=1
                    channels_SIC_tmp[sic] = bank.data[ich]
                channels_SIC.append(channels_SIC_tmp)
                
    full_header = dgtz_header([number_events, number_channels, number_samples, vertical_resulution, 
                              sampling_rate, channels_offset, channels_ttt, channels_SIC, nboard, name_board, DAQversion, channels_pattern])
    return full_header

def daq_dgz_full2array(bank, header, verbose=False, corrected=True, ch_offset=[], tag='LNGS'):

    
    if header.DAQversion > 1000:
        if verbose: print("There is 1 board, and it is {}".format(header.boardNames))
        data_offset = 0

        channels_to_correct = 8 # FOR NOW WE CORRECT ONLY THE FIRST 8 CHANNELS

        waveform = []
        
        if header.boardNames == 1742:

            number_events   = header.ntriggers
            number_channels = header.nchannels
            number_samples  = header.nsamples

            SIC = header.SIC
            to_correct=[]
            if not corrected:
                for ch in range(channels_to_correct):
                    if ch_offset[ch]<-0.25 and ch_offset[ch]>-0.35:
                        to_correct.append(ch)
            
                if number_events!=len(SIC):       ## Check if the start index cell passed are right
                    raise myError("Number of events does not match")
                    
            for ievent in range(number_events):       
                for ichannels in range(number_channels):
                    if verbose:
                        print ("data_offset, data_offset+number_samples",
                            data_offset, data_offset+number_samples)
                        print(bank.data[data_offset:data_offset+number_samples])

                    waveform.append(bank.data[data_offset:data_offset+number_samples])
                    data_offset += number_samples
            if not corrected: ## Correcting the wavefoms (only the ones with offset at -0.3 of first 8 channels)
                waveform = correct_waveforms(waveform, SIC, number_channels, to_correct=to_correct, tag=tag)
            
            if verbose:
                print(number_channels, number_events, number_channels)
                
        elif header.boardNames == 1720:
            
            number_events   = header.ntriggers
            number_channels = header.nchannels
            number_samples  = header.nsamples
            
            waveform = []
            for ievent in range(number_events):       
                for ichannels in range(number_channels):
                    if verbose:
                        print ("data_offset, data_offset+number_samples",
                            data_offset, data_offset+number_samples)
                        print(bank.data[data_offset:data_offset+number_samples])

                    waveform.append(bank.data[data_offset:data_offset+number_samples])
                    data_offset += number_samples
            if verbose:
                print(number_channels, number_events, number_channels)

        else:
            raise myError("You seem to be trying to use a new digitizer model. You need to update the cygno libs for that.")

        return waveform
            
    else:
        
        if verbose: print("There are {} boards, and they are {}".format(header.nBoards, header.boardNames))
    
        data_offset = 0
        channels_to_correct = 8 # FOR NOW WE CORRECT ONLY THE FIRST 8 CHANNELS
    
        waveform_f = []
        waveform_s = []
    
        for idigi,digitizer in enumerate(header.boardNames):
    
            ## Acquiring the "fast digitizer" data
            if str(digitizer) == '1742':  
    
                number_events   = header[0][idigi]
                number_channels = header[1][idigi]
                number_samples  = header[2][idigi]
                SIC = header.SIC
                to_correct=[]
                
                if not corrected:
                    for ch in range(channels_to_correct):
                        if ch_offset[ch]<-0.25 and ch_offset[ch]>-0.35:
                            to_correct.append(ch)
                
                    if number_events!=len(SIC[0]):       ## Check if the start index cell passed are right
                        raise myError("Number of events does not match")
                
                for ievent in range(number_events):       
                    for ichannels in range(number_channels):
                        if verbose:
                            print ("data_offset, data_offset+number_samples",
                                data_offset, data_offset+number_samples)
                            print(bank.data[data_offset:data_offset+number_samples])
    
                        waveform_f.append(bank.data[data_offset:data_offset+number_samples])
                        data_offset += number_samples
    
                if not corrected:              ## Correcting the wavefoms (only the ones with offset at -0.3 of first 8 channels)
                    waveform_f = correct_waveforms(waveform_f, SIC[0], number_channels, to_correct=to_correct, tag=tag)
                
                if verbose:
                    print(number_channels, number_events, number_channels)
            
            ## Acquiring the "slow digitizer" data
            elif str(digitizer) == '1720':  
    
                number_events   = header[0][idigi]
                number_channels = header[1][idigi]
                number_samples  = header[2][idigi]
                waveform_s = []
                for ievent in range(number_events):       
                    for ichannels in range(number_channels):
                        if verbose:
                            print ("data_offset, data_offset+number_samples",
                                data_offset, data_offset+number_samples)
                            print(bank.data[data_offset:data_offset+number_samples])
    
                        waveform_s.append(bank.data[data_offset:data_offset+number_samples])
                        data_offset += number_samples
                if verbose:
                    print(number_channels, number_events, number_channels)
    
            else:
                raise myError("You seem to be trying to use a new digitizer model. You need to update the cygno libs for that.")

        return waveform_f, waveform_s

def daq_slow2array(bank, verbose=False):
    if verbose:
        print(list(bank.data))
    return bank.data

def write2root(fname, img, id=0, option='update', verbose=False):
    import ROOT
    tf = ROOT.TFile.Open(fname+'.root', option)
    img=img.T
    (nx,ny) = img.shape
    h2 = ROOT.TH2D('pic_run',fname+'_'+str(id),nx,0,nx,ny,0,ny)
    h2.GetXaxis().SetTitle('x')
    h2.GetYaxis().SetTitle('y')
    [h2.SetBinContent(bx,by,img[bx,by]) for bx in range(nx) for by in range(ny)]
    h2.Write()
    tf.Close()
    
def TGraph2array(tgraph, verbose=False):
    import ctypes
    xl = []; yl = []
    for i in range(tgraph.GetN()):
        xi = ctypes.c_double(); yi = ctypes.c_double()
        tgraph.GetPoint(i,xi,yi)
        xl.append(xi.value)
        yl.append(yi.value)
    x = np.array(xl)
    y = np.array(yl)
    return x, y

def rootTH2byname(root_file, verbose=False):
    pic = []
    wfm = []
    for i,e in enumerate(root_file.GetListOfKeys()):
        che = e.GetName()
        if ('pic_run' in str(che)):
            pic.append(che)
        elif ('wfm_run' in str(che)):
            wfm.append(che)
    return pic, wfm


def pic_(cfile, iTr=0, verbose=False):
    import ROOT
    import root_numpy as rtnp
    pic, wfm = rootTH2byname(cfile.file)
    image = rtnp.hist2array(cfile.file.Get(pic[iTr])).T
    return image

def wfm_(cfile, iTr=0, iWf=0, verbose=False):
    import ROOT
    import root_numpy as rtnp
    wfm_module=int(cfile.max_wfm/cfile.max_pic)
    if (iTr > cfile.max_pic) or (iWf > wfm_module):
        raise myError("track or wawform out of ragne {:d} {:d}".format(cfile.max_pic, wfm_module))
    i = iTr*wfm_module+iWf
    pic, wfm = rootTH2byname(cfile.file)
    t,a = TGraph2array(cfile.file.Get(wfm[i]))
    return t,a

def read_(f, iTr=0, verbose=False):
    import ROOT
    import root_numpy as rtnp
    pic, wfm = rootTH2byname(f)
    image = rtnp.hist2array(f.Get(pic[iTr])).T
    return image

def ped_(run, path='./ped/', tag = 'LAB', posix=False, min_image_to_read = 0, max_image_to_read = 0, verbose=False):
    #
    # run numero del run
    # path path lettura/scrittura piedistalli
    # tag subdirectory dei dati
    # min_image_to_read , max_image_to_read  range di imagine sul quale fare i piedistalli 
    # max_image_to_read = 0 EQUIVALE A TUTTE LE IMMAGINI
    #
    import ROOT
    import root_numpy as rtnp
    import numpy as np
    import tqdm
    # funzione per fare i piedistalli se gia' non esistino nella diretory

    fileoutm = (path+"mean_Run{:05d}".format(run))
    fileouts = (path+"sigma_Run{:05d}".format(run))

    if os.path.exists(fileoutm+".root") and os.path.exists(fileouts+".root"): 
        # i file gia' esistono
        m_image = read_(ROOT.TFile.Open(fileoutm+".root"))
        s_image = read_(ROOT.TFile.Open(fileouts+".root"))
        print("RELOAD maen file: {:s} sigma file: {:s}".format(fileoutm, fileouts))
        return m_image, s_image
    else:
        # i file non esistono crea il file delle medie e delle sigma per ogni pixel dell'immagine
        if verbose: print (">>> Pedestal Maker! <<<")
        try:
            cfile = open_(run, tag=tag, posix=posix, verbose=verbose)
        except:
            raise myError("openRunError: "+str(run))
        if max_image_to_read == 0:
            max_image_to_read=cfile.max_pic
        print ("WARNING: pdestal from %d to %d" % (min_image_to_read, max_image_to_read))

        m_image = np.zeros((cfile.x_resolution, cfile.y_resolution), dtype=np.float64)
        s_image = np.zeros((cfile.x_resolution, cfile.y_resolution), dtype=np.float64)

        n0 = 0
        for iTr in tqdm.tqdm(range(min_image_to_read, max_image_to_read)):
            image = rtnp.hist2array(cfile.file.Get(cfile.pic[iTr])).T
            image[image<0]=99 #pach per aclune imagini
            m_image += image
            s_image += image**2 
            if verbose and n0 > 0 and n0 % 10==0:  # print progress and debung info for poit 200, 200...
                print ("Debug Image[200,200]: %d => %.2f %.2f %.2f " % (iTr,
                                                image[200,200],
                                                np.sqrt((s_image[200,200] - 
                                                        m_image[200,200]**2 
                                                          / (n0+1)) / n0),
                                                m_image[200,200]/(n0+1),
                                                ))
            n0 += 1
        m_image = m_image/n0
        
        s_image = np.sqrt((s_image - m_image**2 * n0) / (n0 - 1))
        m_image[np.isnan(s_image)==True]=m_image.mean() # pach per i valori insani di sigma e media
        s_image[np.isnan(s_image)==True]=1024
        
       
        ###### print Info and Save OutPut ######################################
        print("WRITING ...")
        write2root(fileoutm, m_image, id=0, option='recreate')
        write2root(fileouts, s_image, id=0, option='recreate')
        print("DONE OUTPUT maen file: {:s} sigma file: {:s}".format(fileoutm, fileouts))
        return m_image, s_image  
    
def ped_mid(run, path_file='/s3/cygno-data/', path_ped='./ped/', tag = 'LNGS', 
            cloud=False, Bari=False, verbose=False):
    #
    # run numero del run
    # path path lettura/scrittura piedistalli
    # tag subdirectory dei dati
    # min_image_to_read , max_image_to_read  range di imagine sul quale fare i piedistalli 
    # max_image_to_read = 0 EQUIVALE A TUTTE LE IMMAGINI
    #
    import ROOT
    import numpy as np
    import tqdm
    import os
    import midas.file_reader
    # funzione per fare i piedistalli se gia' non esistino nella diretory

    fileoutm = (path_ped+"mean_Run{:05d}".format(run))
    fileouts = (path_ped+"sigma_Run{:05d}".format(run))

    if os.path.exists(fileoutm+".root") and os.path.exists(fileouts+".root"): 
        # i file gia' esistono
        m_image = read_(ROOT.TFile.Open(fileoutm+".root"))
        s_image = read_(ROOT.TFile.Open(fileouts+".root"))
        print("RELOAD maen file: {:s} sigma file: {:s}".format(fileoutm, fileouts))
        return m_image, s_image
    else:
        # i file non esistono crea il file delle medie e delle sigma per ogni pixel dell'immagine
        if verbose: print (">>> Pedestal Maker! <<<")
        try:
            mfile = open_mid(run=run, path=path_file, cloud=cloud, Bari=Bari, tag=tag, verbose=verbose)
        except:
            raise myError("openRunError: "+str(run))
            
        init=True
        for event in mfile:
            if event.header.is_midas_internal_event():
                continue
            bank_names = ", ".join(b.name for b in event.banks.values())
            for bank_name, bank in event.banks.items():
                if bank_name=='CAM0': # CAM image
                    image, shape_x_image, shape_y_image = daq_cam2array(bank)
                    if init:
                        m_image = np.zeros((shape_x_image, shape_y_image), dtype=np.float64)
                        s_image = np.zeros((shape_x_image, shape_y_image), dtype=np.float64)

                        n0 = 0
                        init=False
                    #image[image<0]=99 #pach per aclune imagini
                    m_image += image
                    s_image += image**2 
                    n0 += 1

                    if verbose and n0 > 0 and n0 % 10==0:  # print progress and debung info for poit pixel
                        px=1000
                        print ("Debug Image[200,200]: %d => %.2f %.2f %.2f" % (n0,
                                                        image[px,px],
                                                        m_image[px,px]/n0, 
                                                        np.sqrt((s_image[px,px] - (m_image[px,px]**2) / n0) / (n0+1))))
                    
        m_image = m_image/n0    
        s_image = np.sqrt((s_image - m_image**2 * n0) / (n0 - 1))
        m_image[np.isnan(s_image)==True]=m_image.mean() # pach per i valori insani di sigma e media
        s_image[np.isnan(s_image)==True]=1024
        ###### print Info and Save OutPut ######################################
        print("WRITING ...")
        write2root(fileoutm, m_image, id=0, option='recreate')
        write2root(fileouts, s_image, id=0, option='recreate')
        print("DONE OUTPUT maen file: {:s} sigma file: {:s}".format(fileoutm, fileouts))
        return m_image, s_image
    
###
# log book
###
def read_cygno_logbook(sql=True, tag='lngs', start_run=0, end_run=100000000, verbose=False):
    import pandas as pd
    import numpy as np
    import requests
    if sql:
        url = "https://notebook.cygno.cloud.infn.it/cygno_sql_query.php?site="+tag+"&start="+str(start_run)+"&end="+str(end_run)
        if verbose: print ('url: ', url)
        try:
            r = requests.get(url, verify=False)
            df = pd.read_json(url)
            columns = ["varible", "value"]
        except:
            df = pd.DataFrame([])
    else:
        key="1y7KhjmAxXEgcvzMv9v3c0u9ivZVylWp7Z_pY3zyL9F8" # Log Book
        url_csv_file = "https://docs.google.com/spreadsheet/ccc?key="+key+"&output=csv"
        df = pd.read_csv(url_csv_file)
        df = df[df.File_Number.isnull() == False]
        for name in df.columns:
            if name.startswith('Unnamed:'):
                df=df.drop([name], axis=1)
        isacomment = False
        runp = df.File_Number[0]
        for run in df.File_Number:

            if not run.isnumeric():
                if isacomment == False and verbose: print("To Run {}".format(runp)) 
                isacomment = True
                if verbose: print ("--> General comment: {}".format(run))
                index = df[df.File_Number==run].index[0]
                df=df.drop([index], axis=0)
            else:
                if isacomment and verbose: print("From Run {}".format(run)); isacomment = False
            runp = run
        if verbose: print ('Variables: ', df.columns.values)
    return df

def run_info_logbook(run, tag='lngs', sql=True, verbose=False):
    dataInfo=read_cygno_logbook(sql=sql, tag=tag, start_run=run, end_run=run+1, verbose=verbose)
    if sql:
        #out = dataInfo[dataInfo['Run number']==run]
        out = dataInfo
    else:
        out =  dataInfo[dataInfo.File_Number==str(run)]
    if verbose: print(out.values)
    if len(out.values)==0:
        print("NO RUN "+str(run)+" found in history")
    return out

###
# ROOT cygno tool and image tool
###
def cluster_par(xc, yc, image):
    ph = 0.
    dim = xc.shape[0]
    for j in range(0, dim):
        x = int(xc[j])
        y = int(yc[j])
        ph += (image[y,x]) # waring Y prima di X non e' un errore!
    return ph, dim

def n_std_rectangle(x, y, ax, image = np.array([]), n_std=3.0, facecolor='none', **kwargs):
    from matplotlib.patches import Rectangle
    mean_x = x.mean()
    mean_y = y.mean()
    std_x = x.std()
    std_y = y.std()
    half_width = n_std * std_x
    half_height = n_std * std_y
    if image.any():
        rimage = image*0
        xs = int(mean_x - half_width)+1
        xe = int(mean_x + half_width)+1
        ys = int(mean_y - half_height)+1
        ye = int(mean_y + half_height)+1
        # print(ys,ye, xs,xe)
        rimage[ys:ye, xs:xe]=image[ys:ye, xs:xe]
        # print (rimage)
        # print(rimage.sum())
    else:
        rimage = np.array([])
        
    rectangle = Rectangle(
        (mean_x - half_width, mean_y - half_height),
        2 * half_width, 2 * half_height, facecolor=facecolor, **kwargs)
    return ax.add_patch(rectangle), rimage  

def confidence_ellipse(x, y, ax, image = np.array([]), n_std=3.0, facecolor='none', **kwargs):
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    import numpy as np

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    if image.any():
        # ellsisse e' (x-x0)**2/a**2 + (y-y0)**2/b**2 < 1
        # print (mean_x, mean_y, ell_radius_x*scale_x, ell_radius_y*scale_y)
        rimage = image*0
        ar = abs(pearson)
        for x in range(image.shape[1]):
            for y in range(image.shape[0]):
                xr = (y-mean_y)*np.sin(ar)+(x-mean_x)*np.cos(ar)
                yr = (y-mean_y)*np.cos(ar)-(x-mean_x)*np.sin(ar)
                if (xr)**2/(ell_radius_x*scale_x)**2 + (yr)**2/(ell_radius_y*scale_y)**2 < 1:
                    rimage[y,x]=image[y, x]
        # print (rimage)
        # print(rimage.sum())
    else:
        rimage = np.array([])
    
    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ax.add_patch(ellipse), rimage 

def confidence_ellipse_par(x, y, image = np.array([]), n_std=3.0, facecolor='none', **kwargs):
    import numpy as np

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])

    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
                           
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
                           
    width=scale_x*ell_radius_x * 2 
    height=scale_y*ell_radius_y * 2              
    if image.any():
        # ellsisse e' (x-x0)**2/a**2 + (y-y0)**2/b**2 < 1
        # print (mean_x, mean_y, ell_radius_x*scale_x, ell_radius_y*scale_y)
        rimage = image*0
        ar = abs(pearson)
        for x in range(image.shape[1]):
            for y in range(image.shape[0]):
                xr = (y-mean_y)*np.sin(ar)+(x-mean_x)*np.cos(ar)
                yr = (y-mean_y)*np.cos(ar)-(x-mean_x)*np.sin(ar)
                if (xr)**2/(ell_radius_x*scale_x)**2 + (yr)**2/(ell_radius_y*scale_y)**2 < 1:
                    rimage[y,x]=image[y, x]
        # print (rimage)
        # print(rimage.sum())
    else:
        rimage = np.array([])

    
    return width, height, pearson, rimage.sum(), np.size(rimage[rimage>0])

def cluster_elips(points):
    import numpy as np
    x0i= np.argmin(points[:,1])
    a0 = points[x0i][1]
    x1i= np.argmax(points[:,1])
    a1 = points[x1i][1]
    y0i= np.argmin(points[:,0])
    b0 = points[y0i][0]
    y1i= np.argmax(points[:,0])
    b1 = points[y1i][0]
    #print (a0, a1, b0, b1, x0i, points[x0i])
    a  = (a1 - a0)/2.
    b  = (b1 - b0)/2.
    x0 = (a1 + a0)/2.
    y0 = (b1 + b0)/2.
    theta = np.arctan((points[x1i][0]-points[x0i][0])/(points[x1i][1]-points[x0i][1]))
    return x0, y0, a , b, theta

def poit_3d(points, image):
    ########### if 3D #############
    points_3d = []
    for j in range(len(points)):
        y = points[j,0]
        x = points[j,1]
        z = image[int(y),int(x)] # non è un errore Y al posto di X causa asse invertito in python delle imagine
        points_3d.append([y,x,z]) 

    return points_3d

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def smooth(y, box_pts):
    import numpy as np
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def img_proj(img, vmin, vmax, log=False):
    import matplotlib.pyplot as plt
    import numpy as np
    print('INFO: mean: {:.2f}, sigma: {:.2f}, N out of range: {} < vmin={}, {} > vmax={}, insane: {}'.format(
        img.mean(), img.std(), len(img[img<vmin]), vmin, len(img[img>vmax]),vmax, 
        len(img[np.isnan(img)==True])))
    fig, ax = plt.subplots(2,2, figsize=(10,10))
    ax[0,0].imshow(img,  cmap="jet", vmin=vmin,vmax=vmax, aspect="auto")
    x = np.linspace(img.shape[1], 1, img.shape[1])
    #x = np.linspace(1, img.shape[0], img.shape[0])
    ax[0,1].plot(np.sum(img, axis=1),x , 'b-')
    
    x = np.linspace(1, img.shape[0], img.shape[0])
    ax[1,0].plot(x, np.sum(img, axis=0), 'r-')
    ax[1,1].hist(img.ravel(), bins=vmax-vmin, range=(vmin,vmax))
    if log: ax[1,1].set_yscale('log')
    plt.show()

###
# PMT
###
def get_pmt_w_by_triggers(waveform, header, number_of_w_readed, trigger):
    pmt_data=[]
    if trigger <= header[0]:
 
        for t in range(0, header[0]):
            offset = t*header[1]
            if t == trigger:
                for w in range(0, number_of_w_readed):
                    pmt_data.append(waveform[offset])
                    offset+=1
    return np.array(pmt_data)

def correct_waveforms(wfs_in, SIC, nChannels=32, to_correct=list(range(8)), tag='LNGS'):
    nTriggers=0                                   # for now we are correcting only 8 channels

    if len(wfs_in)%nChannels==0:
        nTriggers=int(len(wfs_in)/nChannels)
    else: raise myError("Number of waveforms not understood.")

    if tag=='LNGS':
        table_cell = __table_cell_LNGS__
        table_nsample = __table_nsample_LNGS__
    elif tag=='LNF':
        table_cell = __table_cell_LNF__
        table_nsample = __table_nsample_LNF__
    elif tag=='MAN':
        table_cell = __table_cell_MAN__
        table_nsample = __table_nsample_MAN__
    else: raise myError("Tag not understood.")
        
    wfs = np.copy(wfs_in)
    for trg in range(nTriggers):
        for ch in range(nChannels):
            if ch in to_correct:        # correct only the channels that have an offset ~ -0.3 (determined before)
                indx=trg*nChannels + ch
                wfs[indx] = np.roll(wfs[indx], SIC[trg])
                wfs[indx] = wfs[indx] - table_cell[ch]
                wfs[indx] = np.roll(wfs[indx], -SIC[trg])
                wfs[indx] = wfs[indx] - table_nsample[ch]
        ## peak correction, need the channels of a specific trigger (for now 8)
        schunk=trg*nChannels
        echunk=(trg+1)*nChannels
        tmp_wfs=PeakCorrection(wfs[schunk:echunk]) ## it returns an array with the corrected wf of the channels
        for ch in range(nChannels):
            if ch in to_correct:
                indx=trg*nChannels + ch
                wfs[indx] = tmp_wfs[ch]
        
    return wfs

def PeakCorrection(wfs_in, Nch = 8):
    wfs = np.copy(wfs_in)                        ##list of waveforms
    sample_size = len(wfs[0])                   ## generalize for eventually slow waveforms. 1024 for normal wfs
    avgs = []
    for ch in range(Nch):                       
        avgs.append(np.mean(wfs[ch]))            ## averages of each channel 
    for i in range(1, sample_size):
        offset  = 0
        offset_plus = 0
        for ch in range(Nch):                   #for over the channels
            if i ==1:                           
                if (wfs[ch][2] - wfs[ch][1])>30:
                    offset += 1
                else:
                    if (wfs[ch][3]-wfs[ch][1])>30 and (wfs[ch][3]-wfs[ch][2])>30:
                        offset += 1
            else:
                if i == (sample_size-1) and (wfs[ch][sample_size-2] - wfs[ch][sample_size-1])>30:
                    offset+=1
                else:
                    if (wfs[ch][i-1]-wfs[ch][i])>30:
                        if (wfs[ch][i+1] - wfs[ch][i])>30:
                            offset += 1
                        elif (i+2)<sample_size-2:
                            if (wfs[ch][i+2] - wfs[ch][i])>30 and (wfs[ch][i+1] - wfs[ch][i])<5:
                                offset += 1
                        else:
                            if i == (sample_size-2) or (wfs[ch][i+2]-wfs[ch][i])>30:
                                offset += 1
                                
            if i < (sample_size-6) and (avgs[ch] - wfs[ch][i])<-30 and \
                (avgs[ch] - wfs[ch][i+1])<-30 and \
                (avgs[ch] - wfs[ch][i+2])<-30 and \
                (avgs[ch] - wfs[ch][i+3])<-30 and \
                (avgs[ch] - wfs[ch][i+4])<-30 and \
                (avgs[ch] - wfs[ch][i+5])<-30:
                    offset_plus += 1
        
        if offset == 8:
            for ch in range(Nch):
                if i ==1:
                    if (wfs[ch][2] - wfs[ch][1])>30:
                        wfs[ch][0] = wfs[ch][2]
                        wfs[ch][1] = wfs[ch][2]
                    else:
                        wfs[ch][0] = wfs[ch][3]
                        wfs[ch][1] = wfs[ch][3]
                        wfs[ch][2] = wfs[ch][3]
                else:
                    if i == (sample_size-1):
                        wfs[ch][sample_size-1] = wfs[ch][sample_size-2]
                    else:
                        if (wfs[ch][i+1]-wfs[ch][i])>30:
                            if (wfs[ch][i+1] - wfs[ch][i])>30:
                                wfs[ch][i]   =  int((wfs[ch][i+1]+ wfs[ch][i-1])/2)

                            elif (i+2)<sample_size-2:
                                if (wfs[ch][i+2] - wfs[ch][i])>30 and (wfs[ch][i+1] - wfs[ch][i])<5:
                                    wfs[ch][i]   =  int((wfs[ch][i+2]+ wfs[ch][i-1])/2)
                                    wfs[ch][i+1] =  int((wfs[ch][i+2]+ wfs[ch][i-1])/2)
                        else:
                            if i == (sample_size-2):
                                wfs[ch][sample_size-2] = wfs[ch][sample_size-3]
                                wfs[ch][sample_size-2] = wfs[ch][1023-3]                         
                            else:
                                wfs[ch][i]   = int((wfs[ch][i+2]+wfs[ch][i-1])/2)
                                wfs[ch][i+1] = int((wfs[ch][i+2]+wfs[ch][i-1])/2)

        if offset_plus==8:                                                                      
            for ch in range(Nch):
                for m in range(6):
                    wfs[ch][i+m] = avgs[ch]
    
    return wfs

    
####
# Storage & SQL
###
def daq_sql_cennection(verbose=False):
    import mysql.connector
    import os
    try:
        connection = mysql.connector.connect(
          host=os.environ['MYSQL_IP'],
          user=os.environ['MYSQL_USER'],
          password=os.environ['MYSQL_PASSWORD'],
          database=os.environ['MYSQL_DATABASE'],
          port=int(os.environ['MYSQL_PORT'])
        )
        if verbose: print(connection)
        return connection
    except:
        return False
    
def daq_update_runlog_replica_checksum(connection, run_number, md5sum, verbose=False):
    if verbose: print("md5sum: ", md5sum)
    return cmd.update_sql_value(connection, table_name="Runlog", row_element="run_number", 
                     row_element_condition=run_number, 
                     colum_element="file_checksum", value=md5sum, 
                     verbose=verbose)
def daq_update_runlog_replica_tag(connection, run_number, TAG, verbose=False):
    if verbose: print("TAG: ", TAG)
    return cmd.update_sql_value(connection, table_name="Runlog", row_element="run_number", 
                     row_element_condition=run_number, 
                     colum_element="file_s3_tag", value=TAG, 
                     verbose=verbose)

def daq_update_runlog_replica_size(connection, run_number, size, verbose=False):
    if verbose: print("size: ", size)
    return cmd.update_sql_value(connection, table_name="Runlog", row_element="run_number", 
                     row_element_condition=run_number, 
                     colum_element="file_size", value=size, 
                     verbose=verbose)


def daq_update_runlog_replica_status(connection, run_number, storage, status=-1, verbose=False):
    if storage=="local":
        storage="storage_local_status"
    elif storage=="cloud":
        storage="storage_cloud_status"
    elif storage=="tape":
        storage="storage_tape_status"
    else:
        return 1
    if verbose: print("Storage: "+storage)
    return cmd.update_sql_value(connection, table_name="Runlog", row_element="run_number", 
                     row_element_condition=run_number, 
                     colum_element=storage, value=status, 
                     verbose=verbose)

def daq_read_runlog_replica_status(connection, run_number, storage, verbose=False):
    if storage=="local":
        storage="storage_local_status"
    elif storage=="cloud":
        storage="storage_cloud_status"
    elif storage=="tape":
        storage="storage_tape_status"
    else:
        return -2
    if verbose: print("Storage: "+storage)
    return cmd.read_sql_value(connection, table_name="Runlog", row_element="run_number", 
                     row_element_condition=str(run_number), 
                     colum_element=storage, 
                     verbose=verbose)

def daq_not_on_tape_runs(connection, days=30, verbose=False):
    import numpy as np
    sql = "SELECT * FROM `Runlog` WHERE DATEDIFF(CURRENT_TIMESTAMP, start_time) < "+str(days)+" \
    AND `storage_tape_status` < 1 AND `storage_cloud_status` = 1;"
    mycursor = connection.cursor()
    mycursor.execute(sql)
    value = mycursor.fetchall()
    if verbose: print(mycursor.rowcount)
    mycursor.close()
    try:
        runs = np.array(list(zip(*value))[0])
    except:
        runs = []
    return runs

def daq_run_info(connection, run, verbose=False):
    import numpy as np
    import pandas as pd
    import mysql.connector
    sql = "SELECT * FROM `Runlog` WHERE `run_number` ="+str(run)+";"
    return pd.read_sql(sql, connection)

