import re
import os
import torch
import wenetruntime as wenet
import wave
from pydub import AudioSegment
import time

decoder=wenet.Decoder(model_dir="wenet\zipfile\chs",lang='chs',enable_timestamp=True) 

def ASR(filepath):
    if not os.path.exists("temp_audio"):
        os.mkdir("temp_audio")
    audio  = wave.open(filepath,"rb")
    mark=[False,False,False] # 声道数，采样频率，采样
    flag = False
    if audio.getnchannels()!=1:
        mark[0]=True
    if audio.getframerate()!=16000:
        mark[1]=True
    if audio.getsampwidth()!=2:
        mark[2]=True
    if sum(mark)>=1:
        path=audio_convert(filepath,mark) # 进行变换并保存音频
        audio=wave.open(path,"rb")# 覆盖原变量
        flag = True
    start= time.time()
    text = inference(audio)
    end = time.time()
    print(end-start)
    audio.close()
    if flag:
        os.remove(path)
    return text

def detect(audio):
    if (audio.getnchannels() !=1) or (audio.getframerate()!=16000):
        return False
    return True

def inference(audio):
    wav = audio.readframes(audio.getnframes()) # 总长度 16000 是 采样率 及每秒多少帧(暂时这么理解)
    # 已知一秒有24帧 wav / interval = 总时长
    # 两个循环嵌套 第一个循环测试(是否可以并行解决)
    slice = int(1*16000)*180 # 一百八十秒的音频
    text=""
    lst=[]
    for i in range(0,len(wav),slice): # 时长180秒的音频
        temp_wav=wav[i:min(i+slice,len(wav))]
        interval = int(0.5*16000) * 60 # 三十秒的长度
        for j in range(0,len(temp_wav),interval):
            last = False if j + interval < len(temp_wav) else True
            chunk_wav=temp_wav[j:min(j+interval,len(temp_wav))]
            ans = decoder.decode(chunk_wav,last)
        text+=eval(ans)["nbest"][0]["sentence"]
        lst.extend(splitByTimeStamp(eval(ans)["nbest"][0]["word_pieces"],400))
        # print(ans)
    return ",".join(lst)

def audio_convert(filepath,mark): #根据标记对音频进行操作
    sound = AudioSegment.from_wav(filepath)
    if mark[0]:
        sound=sound.set_channels(1) # 返回声道数量
    if mark[1]:
        sound=sound.set_frame_rate(16000) # 返回采样频率
    if mark[2]:
        sound=sound.set_sample_width(2)
    path=audio_rename(filepath)
    sound.export(path,format="wav")
    return path

def audio_rename(filepath): # 生成次文件地址
    file_name=re.search(r'[^\\/:*?"<>|\r\n]+$',filepath)
    if filepath == None:
        raise "检查文件名书写是否正确"
    return os.path.join("temp_audio",file_name.group())

def splitByTimeStamp(dct_lst,slice): # 根据时间戳对语音识别结果进行断句
    sentence=""
    sentence_lst=[]
    for dct in dct_lst:
        start=dct["start"]
        if sentence=="":
            sentence+=dct["word"]
        else:
            if start-end>slice:
                sentence_lst.append(sentence)
                sentence=dct["word"]
            else:
                sentence+=dct["word"]
        end=dct["end"]
    sentence_lst.append(sentence)
    return sentence_lst

# sound = AudioSegment.from_wav("黄河大合唱.wav")
# sound=sound.set_channels(1) # 返回声道数量
# sound=sound.set_frame_rate(16000) # 返回采样频率
# sound.export("黄河大合唱.wav",format="wav")

# frames=f.getnframes() # 
# rate = f.getframerate() 
# wav_length= frames / float(rate)
# print("音频长度：",wav_length,"秒") # 单位是秒


if __name__=="__main__":
    c=ASR(r"mrher_short.wav")
    # c=eval(c)
    print(c)
    # print
    # lst =os.listdir("chunks")
    # for i in lst:
    #     c=ASR(f"chunks/{i}")
    #     print(c)