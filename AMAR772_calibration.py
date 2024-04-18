#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 12:48:49 2022
using soundtrap or amar recording file  

piston phone wav file

estimate the sensitivity 

@author: xuj
"""

import numpy as np 

# import wave 
# from scipy.io import wavfile
# samplerate, data = wavfile.read('./output/audio.wav')
import soundfile as sf
import matplotlib.pyplot as plt
from scipy.signal import welch
#from scipy.signal import hann as hanning
import scipy.signal as signal

def readWave(fname):
    myfile = sf.SoundFile(fname,mode='r')

#    print myfile.name
#    print myfile.mode 
#    print myfile.samplerate
#    print myfile.channels
#    print myfile.format
#    print myfile.subtype
#    
#    print myfile.endian

    
    data, Fs = sf.read(fname,dtype='int32')
    """
    why divided by 256??
    """
    data=data/256
    myfile.close()
    
    return data,Fs

if __name__=='__main__':
    
    fproj = 'C:/Users/Adamsmi/Documents/PYTHON/GITHUB/HydrophoneCalibration/'
    
    #AMARWav = fproj+'data/ch1/amar772.1.20220707T151736Z.wav'
    AMARWav = fproj+'data/AMAR/AMAR533.20230727T152213Z.wav'
    pistonWav = fproj+'data/MicrophonePiston/20230725-145320(UTC)-Default setup-0000213064.wav'
 
    AMARdata,AMARfs = readWave(AMARWav)
    pistondata,pistonfs = readWave(pistonWav)   
    
    print('piston Fs:',pistonfs)
    
    t1 = np.arange(0,len(AMARdata))/AMARfs
    t2 = np.arange(0,len(pistondata))/pistonfs
    
    print(len(t1),len(AMARdata))
    print(len(t2),len(pistondata))
    
    
    fig,ax1=plt.subplots()
    ax1.plot(t1,AMARdata,'r')
    ax1.plot(t2,pistondata)
    plt.show()
    plt.close()
    
    # exit(1)
    
    pistonStart = 3
    pistonStop = 13
    
    plt.figure()
    
    plt.plot(t2[pistonStart*pistonfs:pistonStop*pistonfs],pistondata[pistonStart*pistonfs:pistonStop*pistonfs])
    plt.show()
    plt.close()
    
    pistonWAV_int = pistondata[pistonStart*pistonfs:pistonStop*pistonfs]
    pistonWAV_v = pistonWAV_int/2**23*10
    
    # plt.figure()
    
    # plt.plot(pistonWAV_v)
    # plt.show()
    # plt.close()
    
    pistonRMS= np.sqrt(np.var(pistonWAV_v))
    print(pistonRMS)
    
    pistonWAV_pascal_dBAir= 20*np.log10(pistonRMS*1000/51.1*1e6/20)

# Microphone 4191-L-001, 12.5 mv/Pa 

    #pistonWAV_pascal_dBAir= 20*np.log10(pistonRMS*1000/12.5*1e6/20)
    # sound pressure wave in refer to 20 microPascal
    
    print('pistonWAV_pascal_dBAir:', pistonWAV_pascal_dBAir)
    
    pistonWAV_pascal_dBwater = 20*np.log10(pistonRMS*1000/51.1*1e6)

    #pistonWAV_pascal_dBwater = 20*np.log10(pistonRMS*1000/12.5*1e6)
    
    print('pistonWAV_pascal_dB water :', pistonWAV_pascal_dBwater)
    

    nblock = pistonfs
    nblock = 2**14

    overlap = nblock/4*3

    # nblock = 1024    
    # overlap = 1024/2
    
    # nblock = pistonfs
    # overlap = pistonfs/2
#    win = signal.flattop(nblock, True)  # old scipy package 
    win = signal.windows.flattop(nblock, True)
    # win = signal.hann(nblock, True)
    # win = signal.blackman(nblock, True)
    
    # print(pistonfs)
    
    
    # plt.figure()
    
    # plt.psd(pistonWAV_v,NFFT=pistonfs,Fs=pistonfs)
    # plt.show()
    # plt.close()
    
    
    # from numpy.fft import fft, ifft
    
    # X = fft(pistonWAV_v)
    # N = len(X)
    # n = np.arange(N)
    # T = N/pistonfs
    # freq = n/T 
    
    # plt.figure()
    # plt.plot(freq,np.abs(X))
    # plt.show()
    
    # plt.close()
    
    f, Pxxf = welch(pistonWAV_v, pistonfs, window=win, noverlap=overlap, nfft=nblock,return_onesided=True,scaling='spectrum')

    # f, Pxxf = signal.welch(pistonWAV_v, pistonfs, 'flattop',  nperseg=pistonfs, nfft=pistonfs, scaling='spectrum')

    # f, Pxxf = signal.periodogram(pistonWAV_v, pistonfs, window=win, nfft=nblock, return_onesided=True,scaling='spectrum')
    
    # f, Pxxf = signal.welch(pistonWAV_v, pistonfs, 'flattop', nperseg=pistonfs,noverlap=pistonfs*3/4,nfft=pistonfs, scaling='spectrum')
    
    plt.figure()
    
    plt.plot(f, Pxxf, '-o')
    plt.xlim(230,260)
    
    plt.grid()
    plt.show()
    plt.close()
    
    # exit(1)    
    max_p = max(Pxxf)
    
    print('max_p (v^2)=', max_p)
    
    max_p_sqrt = np.sqrt(max_p)
    
    print('max_p (v)=', max_p_sqrt)
    
    pistonWAV_pascal_dBAir = 20*np.log10(max_p_sqrt*1000/51.1*1e6/20)
    pistonWAV_pascal_dBwater = 20*np.log10(max_p_sqrt*1000/51.1*1e6)


    #pistonWAV_pascal_dBAir = 20*np.log10(max_p_sqrt*1000/12.5*1e6/20)
    #pistonWAV_pascal_dBwater = 20*np.log10(max_p_sqrt*1000/12.5*1e6)
    
    # Frequency Domain Calculation, should fairly close to time domain estimation. 

    print('pistonWAV_pascal_dBAir:', pistonWAV_pascal_dBAir)
    print('pistonWAV_pascal_dB water :', pistonWAV_pascal_dBwater)
    
    
    # plt.plot(f, Pxxf, '-o')
    # # plt.xlim(230,260)
    
    # plt.grid()
    # plt.show()
    # plt.close()
    wavStart = 88
    wavStop = 100
    
    # exit(1)
    wavSig = AMARdata[wavStart*AMARfs:wavStop*AMARfs]
    sigt1   = t1[wavStart*AMARfs:wavStop*AMARfs]
    plt.figure()
    plt.plot(sigt1, wavSig)
    plt.show()
    plt.close()
    
    
    # exit(1)
    
    wavSigV = wavSig*9/2**24
    plt.figure()
    plt.plot(sigt1,wavSigV)
    plt.show()
    plt.close()    
    
    # wavSigNormMicroPa = wavSigNorm*10**(175.9/20)

    wavSigNormMicroPa = wavSigV*10**((+129.08)/20)  # hydrophone Sensitivity used in the estimation. 

    # wavSigNormMicroPa = wavSigNorm*10**(180/20)
    
    plt.figure()
    plt.plot(wavSigNormMicroPa)
    # plt.show()
    plt.close()            
    
    nblock = 2*17
    # nblock = 1024    
    overlap = nblock/4
    
    nblock = AMARfs
    overlap = AMARfs/2
    
#    win = signal.blackman(nblock, True)
    win = signal.windows.blackman(nblock, True)
    
    f, Pxxf = welch(wavSigNormMicroPa, AMARfs, window=win, noverlap=overlap, nfft=nblock, return_onesided=True,scaling='spectrum')    

    max_p = max(Pxxf)
    
    print('max_p (v^2)=', max_p)    
    max_p_sqrt = np.sqrt(max_p)
    
    print('max_p (v)=', max_p_sqrt)    

    AMARWAV_pascal_dBwater = 20*np.log10(max_p_sqrt)
    print('AMAR WAV_pascal_dB water :', AMARWAV_pascal_dBwater)

    
    plt.plot(f, Pxxf, '-o')
    plt.xlim(230,300)
    plt.show()
    plt.close()
    
