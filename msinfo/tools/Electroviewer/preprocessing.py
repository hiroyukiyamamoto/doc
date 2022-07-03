import os
import glob
import pandas as pd
import csv
import numpy
import sys
import subprocess
from pymsfilereader import MSFileReader

class MyClass:
    def __init__(self, input_folder, input_compoundlist, output_set, output_peak):
        self.input_folder = input_folder
        self.input_compoundlist = input_compoundlist
        self.output_set = output_set
    
    # サンプル名の取得
    def getSampleNames(self):
        """
        各項目のデータを取得
            返り値 : Peak_ID,mz,Cor_MT,MT,Intensity,Area,Rel_Area,SN,Noise,Height,Left_MT,Right_MT,Conc,Tag
        """
        files = glob.glob(input_folder + '*.raw', recursive=True)
        #print(input_folder)
        #print(files)
        samplenames = []
        for file in files:
            rawfile = os.path.basename(file)
            samplename = os.path.splitext(rawfile)[0]
            samplenames.append(samplename)            
        #samplenames.sort(reverse=True) # 要修正
        self.samplenames = samplenames
        self.files = files
    
    # 設定ファイルの出力
    def outSetting(self):   
        """
        各項目のデータを取得
            返り値 : Peak_ID,mz,Cor_MT,MT,Intensity,Area,Rel_Area,SN,Noise,Height,Left_MT,Right_MT,Conc,Tag
        """
        setting = []
        setting.append(['Sample','rawpath'])
        for i in range(len(self.files)):
            setting.append([self.samplenames[i],self.files[i]])
        f = open(output_set, 'w')
        for x in setting:
            f.write('\t'.join([str(i) for i in x]))
            f.write('\n')
        f.close()

    # m/zの範囲(リスト)を取得(csvファイルから)
    def getMzRange(self):
        compoundlist = self.input_compoundlist
        f = open(compoundlist, 'r') # 化合物のリスト
        reader = csv.reader(f)

        next(f)

        deltamz_all = []
        for row in reader:
            errppm = 10 # エラーppmの設定値(10ppm)
            mz = float(row[1])
            mz_high = mz/(1-errppm*10**-6)
            deltamz = mz_high-mz
            deltamz = deltamz
            lowermz = round(mz-deltamz,5) # 桁数(小数点5桁)
            uppermz = round(mz+deltamz,5) # 桁数(小数点5桁)
            deltamz_all.append(str(lowermz)+'-'+str(uppermz)) # m/zの範囲(リスト)

        f.close()

        return deltamz_all  

    # MTのリストを取得(csvファイルから)
    def getMTFromCompoundlist(self):
        compoundlist = self.input_compoundlist
        f = open(compoundlist, 'r') # 化合物のリスト
        reader = csv.reader(f)
        next(f)

        mt_all = []
        for row in reader:
            mt_all.append(float(row[2])) # MT
        f.close()

        return mt_all  
    
    # エレクトロフェログラムを取得
    def getElectropherogram(self, filename_raw, input_mzrange):        
        rawfile = MSFileReader(filename_raw)
        electropherogram = rawfile.GetChroData(startTime=rawfile.StartTime, endTime=rawfile.EndTime,massRange1=input_mzrange)
        mt = electropherogram[0][0]
        intensity = electropherogram[0][1]
        return mt,intensity

        #rawfile.GetChroData(0, 0, 0, null, mzrange, null, 0.0, ref RTLow, ref RTHigh, 0, 0, ref ChroData, ref PeakFlags, ref size);

    # 実サンプルのピークトップのMTを取得
    def getPeaktopMT(self):
        rawfiles = self.files # rawファイル全部
        mzrangelist = self.getMzRange() # 化合物リストから取得したm/zのリスト
        
        topmt_all = []
        # 各サンプル
        for rawfile in rawfiles:
            # 各化合物(化合物リスト中)
            topmt = []
            for mzrange in mzrangelist:
                mt,intensity = self.getElectropherogram(rawfile, mzrange) # エレクトロフェログラムの取得
                max_index = numpy.argmax(intensity)
                topmt.append(mt[max_index])
            topmt_all.append(topmt)

        return topmt_all

    # MTの補正
    def calc_correctedMT(self, output_folder_path, Rexe_path, Rscript_path, topmt_all, mtall):

        #topmt_all = self.topmt_all # ピークトップのMT(化合物×サンプル)
        #mtall = self.getMTFromCompoundlist() # 化合物リストのMTの配列(補正用の説明変数のデータ)

        # 各種設定
        output_mt = output_folder_path + "mt4R.csv"
        cmd = Rexe_path + ' ' + Rscript_path

        # 各サンプル
        corrected_MT = pd.DataFrame([])
        for i in range(0,len(topmt_all)):
            topmt = topmt_all[i]
            mt = pd.DataFrame([topmt, mtall]).T # ピークトップの実測MT(目的変数)、リストのMT(説明変数)
            mt.to_csv(output_mt, index=False, header=False) # 結果の出力                        
            retcode = subprocess.check_call(cmd) # 補正計算

            # Rで計算が出来なかった時
            if retcode==1:
                print("MTの補正値が計算出来ませんでした。")
                sys.exit()

            # 計算できた時
            output_path = output_folder_path + "mt4py.csv"
            csv_input = pd.read_csv(output_path, encoding="ms932", sep=",")

            # 結果の結合(データフレーム)
            corrected_MT = pd.concat([corrected_MT, csv_input], axis=1)

        # 列名の設定
        corrected_MT.columns = self.samplenames

        # 戻り値
        return corrected_MT  
    
    # ピークファイルの出力
    def outPeak(self, corrected_MT):
        """
        各項目のデータを取得
            返り値 : Peak_ID,mz,Cor_MT,MT,Intensity,Area,Rel_Area,SN,Noise,Height,Left_MT,Right_MT,Conc,Tag
        """

        # 化合物リスト
        f = open(self.input_compoundlist, 'r') # 化合物のリスト
        g = open(output_peak, 'w') # 出力ファイル
        g.write(self.output_set + "\n")

        # ヘッダ
        header=["Peak","mz","dmz"]
        for x in range(len(self.files)):
            header.extend(["Sample","LeftMT","MT","RightMT"])
        for x in header:
            g.write(x + "\t")
        g.write("\n")

        next(f) # headerを飛ばす

        mz = []; mt=[]
        reader = csv.reader(f)

        # // 各物質(各行)
        i=0
        for row in reader:
            g.write(row[0]+"\t") # compound name
            g.write(row[1]+"\t") # m/z(理論値)

            errppm = 10 # エラーppmの設定値(10ppm)
            mz = round(float(row[1]),5)
            mz_high = mz/(1-errppm*10**-6)
            deltamz = round(mz_high-mz,10)
           
            g.write(str(deltamz)+"\t") # delta m/z
           
            ## // 各サンプル
            k = 0
            for sample in self.samplenames:
                g.write(sample+"\t") # sample name   
                mt = corrected_MT.iat[i,k] # 補正後のMT                 
                g.write(str(mt-0.3)+"\t") # left MT
                g.write(str(mt)+"\t") # MT
                g.write(str(mt+0.3)+"\t") # right MT
                k = k + 1
            g.write("\n") # 改行
            i = i + 1
        g.close()            

    def setOutput(self):
        print("test")

# Main
if __name__== '__main__':

    # コマンドラインで実行の場合
    args = sys.argv

    # 長さが8の場合
    if (len(args)==8):
        # 入力する
        input_folder = args[1]
        input_compoundlist = args[2]
        output_set = args[3]
        output_peak = args[4]
        output_folder_path = args[5]
        Rexe_path = args[6]
        Rscript_path = args[7]
    else:
        # // 各種準備
        input_folder = "C:/Users/yamamoto/Desktop/ElectroViwer/EV data/FullMS_Cation_plasma/"
        input_compoundlist = "C:/Users/yamamoto/Desktop/ElectroViwer/EV data/FullMS_Cation_plasma/compound_list.csv"
        output_set = "C:/Users/yamamoto/Desktop/ElectroViwer/EV data/FullMS_Cation_plasma/setfile.txt"
        output_peak = "C:/Users/yamamoto/Desktop/ElectroViwer/EV data/FullMS_Cation_plasma/peakfile.txt"

        a = MyClass(input_folder,input_compoundlist,output_set, output_peak)
        a.getSampleNames() # サンプル名の取得
        a.outSetting() # 設定ファイルの出力
    
        mt_all = a.getMTFromCompoundlist() # リストのMTを取得
        topmt_all = a.getPeaktopMT() # 実サンプルのピークトップのMTを取得

        # // MT補正計算
        output_folder_path = "C:/temp/electroviewer/"
        Rexe_path = '"C:/Users/yamamoto/R-3.5.2/bin/x64/Rscript.exe"'
        Rscript_path = '"C:/Users/yamamoto/Documents/R/test/mtcor2.R"'

        corrected_MT = a.calc_correctedMT(output_folder_path, Rexe_path, Rscript_path, topmt_all, mt_all) 

        # // ピークファイルの出力
        a.outPeak(corrected_MT) 



