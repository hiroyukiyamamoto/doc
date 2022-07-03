using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MSFileReaderLib;
using MathNet.Numerics.Statistics;
using System.Windows.Forms;

namespace WindowsFormsApp5
{
    // -------------------
    //  RawEic(rawfile)
    // -------------------
    public class RawEic
    {
        string rawfilename;
        string mz;
        string delta_mz;
        string topmt;
        string leftmt;
        string rightmt;

        public LocalEic EIC;
        public LocalMtIndex InputMtIndex;
        public LocalMtIndex MtIndex = new LocalMtIndex();
        public LocalMt MT;

        public List<LocalEic> EICS;
        public List<LocalMtIndex> InputMtIndexs;

        public LocalPeak Peak;

        // LocalEic
        public struct LocalEic
        {
            public List<double> mt;
            public List<double> intensity;
        }

        // LocalPeak
        public struct LocalPeak
        {
            public double[] rtrange;
            public double[] intrange;
        }

        // LocalMtIndex
        public class LocalMtIndex
        {
            public int top_index;
            public int left_index;
            public int right_index;
        }

        // LocalMt
        public class LocalMt
        {
            public string leftmt;
            public string rightmt;
            public string topmt;
        }

        // LocalMz
        public struct LocalMz
        {
            public string mz;
            public string delta_mz;
        }

        // SetEic
        public void SetEic(string rawfilename, LocalMz inputmz, LocalMt inputmt)
        {
            this.rawfilename = rawfilename;
            this.mz = inputmz.mz;
            this.delta_mz = inputmz.delta_mz;
            this.leftmt = inputmt.leftmt;
            this.rightmt = inputmt.rightmt;
            this.topmt = inputmt.topmt;

            MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
            rawfile.Open(rawfilename);
            rawfile.SetCurrentController(0, 1);

            double mz_low = double.Parse(mz) - double.Parse(delta_mz);
            double mz_high = double.Parse(mz) + double.Parse(delta_mz);

            string mzrange = mz_low.ToString() + '-' + mz_high.ToString();

            object ChroData = null;
            object PeakFlags = null;

            Double RTLow = 0;
            Double RTHigh = 100;
            int size = 100;

            rawfile.GetChroData(0, 0, 0, null, mzrange, null, 0.0, ref RTLow, ref RTHigh, 0, 0, ref ChroData, ref PeakFlags, ref size);

            List<double> MT = new List<double>();
            List<double> Intensity = new List<double>();

            for (int k = 0; k < size; k++)
            {
                double mt = (double)(ChroData as Array).GetValue(0, k);
                double intensity = (double)(ChroData as Array).GetValue(1, k);

                MT.Add(mt);
                Intensity.Add(intensity);
            }

            rawfile.Close();

            EIC.mt = MT;
            EIC.intensity = Intensity;

            double[] rt = EIC.mt.ToArray();
            double[] diff1 = new double[rt.Length];
            double[] diff2 = new double[rt.Length];
            double[] diff3 = new double[rt.Length];

            // MTとLeftMT、RightMTとの差分
            for (int l = 0; l <= rt.Length - 1; l++)
            {
                diff1[l] = Math.Abs(rt[l] - double.Parse(leftmt));
                diff2[l] = Math.Abs(rt[l] - double.Parse(rightmt));
                diff3[l] = Math.Abs(rt[l] - double.Parse(topmt));
            }

            // 最も近いMTを取得
            double min_left = diff1.Minimum();
            double min_right = diff2.Minimum();
            double min_mt = diff3.Minimum();

            // leftMTとrightMTのindexを取得
            int leftindex = Array.IndexOf(diff1.ToArray(), min_left); // leftMTの位置
            int rightindex = Array.IndexOf(diff2.ToArray(), min_right); // rightMTの位置
            int mtindex = Array.IndexOf(diff3.ToArray(), min_mt); // MTの位置

            // -----------------------------------------------------------------------
            // 　MT範囲を超える場合はMTのインデックスを修正する
            // -----------------------------------------------------------------------
            // 基本的には、ピーク取得の時にそうならないように設定しておく(削除する?)    
            if (double.Parse(leftmt) < MT.Minimum() || double.Parse(rightmt) > MT.Maximum() || double.Parse(topmt) < MT.Minimum() || double.Parse(topmt) > MT.Maximum())
            {

                leftindex = 0;
                rightindex = rt.Length - 1;
                // mtindex = 
            }

            MtIndex.top_index = mtindex;
            MtIndex.left_index = leftindex;
            MtIndex.right_index = rightindex;

            InputMtIndex = MtIndex;
        }

        // SetMtIndex
        public void SetMtIndex(LocalMtIndex DeltaMtIndex)
        {
            MtIndex.left_index = InputMtIndex.left_index + DeltaMtIndex.left_index;
            MtIndex.right_index = InputMtIndex.left_index + DeltaMtIndex.right_index;
            MtIndex.top_index = InputMtIndex.top_index;
        }

        // GetEic
        public LocalEic GetEic()
        {
            return EIC;
        }

        // GetPeak
        public LocalPeak GetPeak()
        {
            double[] rtrange = EIC.mt.GetRange(MtIndex.left_index, MtIndex.right_index - MtIndex.left_index).ToArray();
            double[] intrange = EIC.intensity.GetRange(MtIndex.left_index, MtIndex.right_index - MtIndex.left_index).ToArray();

            Peak.rtrange = rtrange;
            Peak.intrange = intrange;

            return Peak;
        }

        // GetMtIndex
        public LocalMtIndex GetMtIndexFromMt()
        {
            double[] rt = EIC.mt.ToArray();
            double[] diff1 = new double[rt.Length];
            double[] diff2 = new double[rt.Length];
            double[] diff3 = new double[rt.Length];

            // MTとLeftMT、RightMTとの差分
            for (int l = 0; l <= rt.Length - 1; l++)
            {
                diff1[l] = Math.Abs(rt[l] - double.Parse(leftmt));
                diff2[l] = Math.Abs(rt[l] - double.Parse(rightmt));
                diff3[l] = Math.Abs(rt[l] - double.Parse(topmt));
            }

            // 最も近いMTを取得
            double min_left = diff1.Minimum();
            double min_right = diff2.Minimum();
            double min_mt = diff3.Minimum();

            // leftMTとrightMTのindexを取得
            int leftindex = Array.IndexOf(diff1.ToArray(), min_left); // leftMTの位置
            int rightindex = Array.IndexOf(diff2.ToArray(), min_right); // rightMTの位置
            int mtindex = Array.IndexOf(diff3.ToArray(), min_mt); // MTの位置

            LocalMtIndex MtIndex = new LocalMtIndex();

            MtIndex.top_index = mtindex;
            MtIndex.left_index = leftindex;
            MtIndex.right_index = rightindex;

            return MtIndex;
        }

        // GetMtIndex
        public LocalMtIndex GetMtIndex()
        {
            return MtIndex;
        }

    }


    // -------------------------
    //  これ使ってない
    // -------------------------
    // RawMultiEic
    public class RawMultiEic
    {
        List<List<double>> MultiIntensity = new List<List<double>>();
        List<List<double>> Area = new List<List<double>>();
        List<double> AlignError = new List<double>();

        double max_delta_left;
        double max_delta_right;
        double median_delta_left;
        double median_delta_right;

        RawEic.LocalMtIndex pindex = new RawEic.LocalMtIndex();

        // SetMultiPeaks
        public void SetMultiPeaks(List<RawEic.LocalMtIndex> mtindex, List<RawEic.LocalEic> eic)
        {

            double[] delta_left = new double[mtindex.Count()];
            double[] delta_right = new double[mtindex.Count()];

            for (int i = 0; i < mtindex.Count(); i++)
            {
                delta_left[i] = mtindex[i].top_index - mtindex[i].left_index;
                delta_right[i] = mtindex[i].right_index - mtindex[i].top_index;
            }

            max_delta_left = delta_left.Max();
            max_delta_right = delta_right.Max();

            median_delta_left = delta_left.Median();
            median_delta_right = delta_right.Median();

            List<double> eeeic;
            for (int i = 0; i < mtindex.Count(); i++)
            {
                int leftmt_index = (int)mtindex[i].top_index - (int)max_delta_left;
                int rightmt_index = (int)mtindex[i].top_index + (int)max_delta_right;

                eeeic = new List<double>();
                eeeic = eic[i].intensity.GetRange(leftmt_index, rightmt_index - leftmt_index);

                MultiIntensity.Add(eeeic);
            }

        }

        // GetMultiPeaks
        public List<List<double>> GetMultiPeaks()
        {
            return MultiIntensity;
        }

        // GetLRTindex
        public RawEic.LocalMtIndex GetPeakIndex()
        {
            int maxlen = MultiIntensity[0].Count();
            double topindex = maxlen - max_delta_right;

            pindex.top_index = maxlen - (int)max_delta_right;
            pindex.left_index = pindex.top_index - (int)median_delta_left;
            pindex.right_index = pindex.top_index + (int)median_delta_right;

            return pindex;
        }

        // GetAreas
        public List<List<double>> GetAreas()
        {
            return Area;
        }

        // GetAlignError
        public List<double> GetAlignError()
        {
            // Normalized errorを計算
            // 特定のピークIDのクロマトを取得して、平均値を計算する
            // 横幅は全部同じでないとまずいが、今はそうなっていない？
            






            return AlignError;
        }

    }



    // LoadRawで読み込んだクロマトを修正する
    // どういう手順で修正するか？
    // LoadRaw.SetEics
    // リストの値を入れ替える
    // 特定のサンプル、特定のピークIDのクロマトを入れかえる
    // SetEicを作ればよい

    // LoadRaw for the first time
    public class LoadRaw
    {

        string rawfilename;
        string mz;
        string delta_mz;
        string topmt;
        string leftmt;
        string rightmt;

        int peakID;

        public LocalEic EIC;
        public LocalMtIndex InputMtIndex;
        public LocalMtIndex MtIndex;
        public LocalMt MT;

        public List<LocalEic> EICS;
        public List<LocalMtIndex> InputMtIndexs;
        List<double> AlignError;

        public LocalPeak Peak;

        // LocalEic
        public struct LocalEic
        {
            public List<double> mt;
            public List<double> intensity;
        }

        // LocalPeak
        public struct LocalPeak
        {
            public double[] rtrange;
            public double[] intrange;
        }

        // LocalMtIndex
        public struct LocalMtIndex
        {
            public int top_index;
            public int left_index;
            public int right_index;
        }

        // LocalMt
        public struct LocalMt
        {
            public string leftmt;
            public string rightmt;
            public string topmt;
        }

        // LocalMz
        public struct LocalMz
        {
            public string mz;
            public string delta_mz;
        }
        
        // SetEics
        // ピークトップから、左右の点数を固定する        
        public void SetEics(string rawfilename, List<RawEic.LocalMz> inputmz, List<RawEic.LocalMt> inputmt)
        {

            double dmt = 0.2; // 設定値(とりあえず)
            // 適当に広めに設定しておけばよさそう
            // 一般的には、設定値として入力できるようにする

            this.rawfilename = rawfilename;

            LocalEic alleic = new LocalEic();

            EICS = new List<LocalEic>();

            MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
            rawfile.Open(rawfilename);
            rawfile.SetCurrentController(0, 1);

            for (int i = 0; i < inputmz.Count(); i++)
            {
                mz = inputmz[i].mz;
                delta_mz = inputmz[i].delta_mz;
                leftmt = (double.Parse(inputmt[i].topmt) - dmt).ToString();
                rightmt = (double.Parse(inputmt[i].topmt) + dmt).ToString();
                topmt = inputmt[i].topmt;

                double mz_low = double.Parse(mz) - double.Parse(delta_mz);
                double mz_high = double.Parse(mz) + double.Parse(delta_mz);

                string mzrange = mz_low.ToString() + '-' + mz_high.ToString();

                object ChroData = null;
                object PeakFlags = null;

                Double RTLow = 0;
                Double RTHigh = 100;
                int size = 100;

                rawfile.GetChroData(0, 0, 0, null, mzrange, null, 0.0, ref RTLow, ref RTHigh, 0, 0, ref ChroData, ref PeakFlags, ref size);
                
                List<double> MT = new List<double>();
                List<double> Intensity = new List<double>();

                // 全領域のEIC
                for (int k = 0; k < size; k++)
                {
                    double mt = (double)(ChroData as Array).GetValue(0, k);
                    double intensity = (double)(ChroData as Array).GetValue(1, k);

                    MT.Add(mt);
                    Intensity.Add(intensity);
                }
                
                EIC.mt = MT;　// 1次元配列
                EIC.intensity = Intensity;

                double[] rt = EIC.mt.ToArray();
                double[] diff1 = new double[rt.Length];
                double[] diff2 = new double[rt.Length];
                double[] diff3 = new double[rt.Length];

                // MTとLeftMT、RightMTとの差分
                for (int l = 0; l <= rt.Length - 1; l++)
                {
                    diff1[l] = Math.Abs(rt[l] - double.Parse(leftmt));
                    diff2[l] = Math.Abs(rt[l] - double.Parse(rightmt));
                    diff3[l] = Math.Abs(rt[l] - double.Parse(topmt));
                }

                // 最も近いMTを取得
                double min_left = diff1.Minimum();
                double min_right = diff2.Minimum();
                double min_mt = diff3.Minimum();

                // leftMTとrightMTのindexを取得
                int leftindex = Array.IndexOf(diff1.ToArray(), min_left); // leftMTの位置
                int rightindex = Array.IndexOf(diff2.ToArray(), min_right); // rightMTの位置
                int mtindex = Array.IndexOf(diff3.ToArray(), min_mt); // MTの位置

                // EICの領域を取得
                List<double> eeeic = new List<double>();
                List<double> mmmt = new List<double>();

                mmmt = EIC.mt.GetRange(leftindex, rightindex - leftindex);
                eeeic = EIC.intensity.GetRange(leftindex, rightindex - leftindex);

                // -----------------------------------------------------------------------
                // 　MT範囲を超える場合はMT全領域にする
                // -----------------------------------------------------------------------
                // 基本的には、ピーク取得の時にそうならないように設定しておく(削除する?)    
                if (double.Parse(leftmt) < MT.Minimum() || double.Parse(rightmt) > MT.Maximum() || double.Parse(topmt) < MT.Minimum() || double.Parse(topmt) > MT.Maximum())
                {

                    leftindex = 0;
                    rightindex = rt.Length - 1;
                    // mtindex = 

                    mmmt = EIC.mt;
                    eeeic = EIC.intensity;
                }

                // 配列に代入
                alleic.intensity = eeeic;
                alleic.mt = mmmt;

                EICS.Add(alleic); // 2次元配列
            }

            rawfile.Close();
        }

        // SetEic
        // ピークトップから、左右の点数を固定する        
        public void SetEic(string rawfilename,　int input_peakID, RawEic.LocalMz inputmz, RawEic.LocalMt inputmt)
        {

            peakID = input_peakID;

            double dmt = 0.2; // 設定値(とりあえず) 、ピーク領域(左下画面)のMT幅を固定している

            this.rawfilename = rawfilename;

            LocalEic alleic = new LocalEic();

            MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
            rawfile.Open(rawfilename);
            rawfile.SetCurrentController(0, 1);

            mz = inputmz.mz; // m/z(入力値)
            delta_mz = inputmz.delta_mz; // m/zの差分
            leftmt = (double.Parse(inputmt.topmt) - dmt).ToString(); // leftmt
            rightmt = (double.Parse(inputmt.topmt) + dmt).ToString(); // rightmt
            topmt = inputmt.topmt; // peaktopのmt

            double mz_low = double.Parse(mz) - double.Parse(delta_mz); // m/zの下側
            double mz_high = double.Parse(mz) + double.Parse(delta_mz); // m/zの上側

            string mzrange = mz_low.ToString() + '-' + mz_high.ToString(); // m\zの範囲

            // エレクトロフェログラム取得
            object ChroData = null;
            object PeakFlags = null;

            Double RTLow = 0;
            Double RTHigh = 100;
            int size = 100;

            rawfile.GetChroData(0, 0, 0, null, mzrange, null, 0.0, ref RTLow, ref RTHigh, 0, 0, ref ChroData, ref PeakFlags, ref size);
            
            List<double> MT = new List<double>();
            List<double> Intensity = new List<double>();

            // 全領域のEICを変数(MT, Intensity)に代入
            for (int k = 0; k < size; k++)
            {
                double mt = (double)(ChroData as Array).GetValue(0, k);
                double intensity = (double)(ChroData as Array).GetValue(1, k);

                MT.Add(mt);
                Intensity.Add(intensity);
            }

            rawfile.Close();


            EIC.mt = MT;　// 1次元配列
            EIC.intensity = Intensity;

            double[] rt = EIC.mt.ToArray();
            double[] diff1 = new double[rt.Length];
            double[] diff2 = new double[rt.Length];
            double[] diff3 = new double[rt.Length];

            // MTとLeftMT、RightMTとの差分
            for (int l = 0; l <= rt.Length - 1; l++)
            {
                diff1[l] = Math.Abs(rt[l] - double.Parse(leftmt));
                diff2[l] = Math.Abs(rt[l] - double.Parse(rightmt));
                diff3[l] = Math.Abs(rt[l] - double.Parse(topmt));
            }

            // 最も近いMTを取得
            double min_left = diff1.Minimum();
            double min_right = diff2.Minimum();
            double min_mt = diff3.Minimum();

            // leftMTとrightMTのindexを取得
            int leftindex = Array.IndexOf(diff1.ToArray(), min_left); // leftMTの位置
            int rightindex = Array.IndexOf(diff2.ToArray(), min_right); // rightMTの位置
            int mtindex = Array.IndexOf(diff3.ToArray(), min_mt); // MTの位置
          
            // EICの領域を取得
            List<double> eeeic = new List<double>();
            List<double> mmmt = new List<double>();

            mmmt = EIC.mt.GetRange(leftindex, rightindex - leftindex);
            eeeic = EIC.intensity.GetRange(leftindex, rightindex - leftindex);

            // -----------------------------------------------------------------------
            // 　MT範囲を超える場合はMT全領域にする
            // -----------------------------------------------------------------------
            // 基本的には、ピーク取得の時にそうならないように設定しておく(削除する?)            
            if (double.Parse(leftmt) < MT.Minimum() || double.Parse(rightmt) > MT.Maximum() || double.Parse(topmt) < MT.Minimum() || double.Parse(topmt) > MT.Maximum())
            {
                leftindex = 0;
                rightindex = rt.Length - 1;

                mmmt = EIC.mt;
                eeeic = EIC.intensity;
            }

            alleic.intensity = eeeic;
            alleic.mt = mmmt;

            EICS[peakID] = alleic; // 1次元配列
        }
    
        // GetAllEic
        public List<LocalEic> GetAllEic()
        {
            return EICS;
        }

    }   

    public class AlignError
    {

        private List<List<LoadRaw.LocalEic>> AllEIC;
        private int leftindex;
        private int rightindex;
        private int peakID;
        private int pointnum;

        // SetAllEIC
        public void SetAllEIC(List<List<LoadRaw.LocalEic>> InputAllEIC)
        {
            this.AllEIC = InputAllEIC;
        }

        // SetLeftIndex
        public void SetLeftIndex(int input_leftindex)
        {
            this.leftindex = input_leftindex;
        }

        // SetRightIndex
        public void SetRightIndex(int input_rightindex)
        {
            this.rightindex = input_rightindex;
        }

        // SetPeakID
        public void SetPeakID(int input_peakID)
        {
            this.peakID = input_peakID;
        }

        // GetAlignError
        public List<double> GetAlignError()
        {
            List<double> a;
            List<double> b;
            List<List<double>> c;

            pointnum = rightindex - leftindex;

            // Normalized peak
            List<List<double>> Y = new List<List<double>>();

            List<double> y;
            List<double> y1 = new List<double>();
            List<double> y2 = new List<double>();

            for (int i = 0; i < AllEIC.Count(); i++)
            {
                y = new List<double>();
                y = AllEIC[i][peakID].intensity.GetRange(leftindex, pointnum);
                // MTがエレクトロフェログラムの範囲を超えてエラーになっている

                y1 = new List<double>();
                for (int l = 0; l < y.Count(); l++)
                {
                    y1.Add(y[l] - y.Min());
                }
                y2 = new List<double>();
                for (int l = 0; l < y1.Count(); l++)
                {
                    y2.Add(y1[l] / y1.Max());
                }
                Y.Add(y2); // [sample][datapoint]
            }

            // Calculate Mean
            //b = new List<double>();
            c = new List<List<double>>();
            for (int j = 0; j < pointnum; j++)
            {
                a = new List<double>();
                for (int i = 0; i < AllEIC.Count(); i++)
                {
                    a.Add(Y[i][j]);
                }
                //b.Add(a.Mean()); // mean
                c.Add(a); // [datapoint][sample]
            }

            // Calculate Error
            List<double> AE = new List<double>(); // AlignmentError
            List<double> z;
            for (int i = 0; i < AllEIC.Count(); i++)
            {
                z = new List<double>();
                for (int j = 0; j < pointnum; j++)
                {
                    z.Add(Math.Pow(c[j][i] - c[j][0], 2)); // 基準となるサンプル : 0
                }
                AE.Add(z.Sum()); // [sample]
            }
            return AE;
        }
    }
}