using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace WindowsFormsApp5
{
    // TxtPeak
    public class TxtPeak
    {

        public LocalPeakList input = new LocalPeakList();
        public List<LocalPeakList> input2 = new List<LocalPeakList>();
        public List<List<LocalPeakList>> input3 = new List<List<LocalPeakList>>();

        public RawEic.LocalMz inputmz = new RawEic.LocalMz();
        public RawEic.LocalMt inputmt = new RawEic.LocalMt();

        public List<RawEic.LocalMz> inputmz2 = new List<RawEic.LocalMz>();
        public List<RawEic.LocalMt> inputmt2 = new List<RawEic.LocalMt>();

        public List<List<RawEic.LocalMz>> inputmz3 = new List<List<RawEic.LocalMz>>();
        public List<List<RawEic.LocalMt>> inputmt3 = new List<List<RawEic.LocalMt>>();
        
        public List<LocalPeakList> output = new List<LocalPeakList>();

        public List<RawEic.LocalMz> outputmz = new List<RawEic.LocalMz>();
        public List<RawEic.LocalMt> outputmt = new List<RawEic.LocalMt>();

        public int n;
        public int targetid;

        public struct LocalPeakList
        {
            public string peak;
            public string sample;
            public string mz;
            public string delta_mz;
            public string left_mt;
            public string mt;
            public string right_mt;
        }

        // SetPeakFromTxt
        public void SetPeakFromTxt(string inputpeakfile)
        {

            // アライメント後ピークリストの読み込み(ピーク毎)
            StreamReader objReader = new StreamReader(inputpeakfile, false);

            objReader.ReadLine(); // 設定ファイルのパスを捨てる
            objReader.ReadLine(); // headerを捨てる

            List<List<string>> minput2 = new List<List<string>>();
            
            string sLine = "";
            while (sLine != null)
            {
                sLine = objReader.ReadLine();
                List<string> minput = new List<string>();

                if (sLine != null)
                {
                    string[] temp_line0 = sLine.Split('\t');
                                        
                    minput.AddRange(temp_line0);
                    minput2.Add(minput);
                }
            }

            int ncol = minput2[0].Count();
            int nrow = minput2.Count();

            List<List<string>> minput4 = new List<List<string>>();

            for (int i = 0; i < ncol; i++)
            {
                List<string> minput3 = new List<string>();
                for (int j = 0; j < nrow; j++)
                {
                    minput3.Add(minput2[j][i]);
                }
                minput4.Add(minput3);
            }

            // サンプル数
            n = (ncol - 3) / 4;
            
            // アライメント後ピークリストの読み込み(サンプル毎)
            input = new LocalPeakList();
            input2 = new List<LocalPeakList>();
        　　input3 = new List<List<LocalPeakList>>();

            // 各サンプル       
            for (int i = 0; i < n; i++)
            {
                input2 = new List<LocalPeakList>(); // 初期化
                inputmz2 = new List<RawEic.LocalMz>(); // 初期化
                inputmt2 = new List<RawEic.LocalMt>(); // 初期化

                // 各ピーク
                for (int j = 0; j < nrow; j++)
                {
                    input = new LocalPeakList(); // 初期化
                    inputmz = new RawEic.LocalMz(); // 初期化
                    inputmt = new RawEic.LocalMt(); // 初期化

                    input.peak = minput4[0][j];
                    input.sample = minput4[4*(i+1)-1][j];
                    input.mz = minput4[1][j];
                    input.delta_mz = minput4[2][j];
                    input.left_mt = minput4[4*(i+1)][j];
                    input.mt = minput4[4*(i+1)+1][j];
                    input.right_mt = minput4[4*(i+1)+2][j];

                    inputmz.mz = minput4[1][j];
                    inputmz.delta_mz = minput4[2][j];
                    inputmt.leftmt = minput4[4*(i+1)][j];
                    inputmt.topmt = minput4[4*(i+1)+1][j];
                    inputmt.rightmt = minput4[4*(i+1)+2][j];

                    input2.Add(input); // 特定のファイルの全ピーク

                    inputmz2.Add(inputmz);
                    inputmt2.Add(inputmt);

                }

                // ピークID=1
                // MessageBox.Show(input2[0].mt.ToString());

                input3.Add(input2); // 全サンプルの全ピーク

                inputmz3.Add(inputmz2); // 全サンプルの全ピーク
                inputmt3.Add(inputmt2); // 全サンプルの全ピーク

            }

        }     

        // SetPeakTopMT
        public void SetPeakTopMT(int SampleNumber, int PeakId, double PeakTopMt)
        {
            inputmt3[SampleNumber][PeakId].topmt = PeakTopMt.ToString();
            inputmt3[SampleNumber][PeakId].leftmt = (PeakTopMt-0.2).ToString();
            inputmt3[SampleNumber][PeakId].rightmt = (PeakTopMt + 0.2).ToString();            
        }

        // GetPeakNum
        public int GetPeakNum()
        {
            return input2.Count;
        }

        // GetPeakId
        public List<string> GetPeakId()
        {
            List<string> id = new List<string>();
            for (int i = 0; i < input2.Count; i++)
            {
                id.Add(input2[i].peak);
            }
            return id;
        }

        // GetPeakFromTxt
        public List<LocalPeakList> GetPeakFromTxt()
        {
            return input2;
        }

        // GetMz
        public List<List<RawEic.LocalMz>> GetMz()
        {
            return inputmz3;
        }

        // GetMt
        public List<List<RawEic.LocalMt>> GetMt()
        {
            return inputmt3;
        }

        // GetAllPeakFromTxt
        public List<List<LocalPeakList>> GetAllPeakFromTxt()
        {
            return input3;
        }

        // SetPeakId
        public void SetPeakId(int inputid)
        {
            targetid = inputid;
        }

        // GetPeakId
        public int GetCPeakId()
        {
            return(targetid);
        }

        // GetPeakFromId
        public List<LocalPeakList> GetPeakFromId()
        {
            output = new List<LocalPeakList>();

            // input3から必要な情報を取得
            for (int i = 0; i < n; i++)
            {
                output.Add(input3[i][targetid]);
            }
            return output;
        }

        // GetPeakmzFromId
        public List<RawEic.LocalMz> GetPeakmzFromId()
        {
            outputmz = new List<RawEic.LocalMz>();

            // input3から必要な情報を取得
            for (int i = 0; i < n; i++)
            {
                outputmz.Add(inputmz3[i][targetid]);
            }
            return outputmz;
        }

        // GetPeakmtFromId
        public List<RawEic.LocalMt> GetPeakmtFromId()
        {
            outputmt = new List<RawEic.LocalMt>();

            // input3から必要な情報を取得
            for (int i = 0; i < n; i++)
            {
                outputmt.Add(inputmt3[i][targetid]);
            }
            return outputmt;
        }

        // GetDeltaMz
        public List<string> GetDeltaMz()
        {
            List<string> output_delta_mz = new List<string>();
            for (int i = 0; i < input2.Count() - 1; i++)
            {
                output_delta_mz.Add(input2[i].delta_mz);
            }
            return (output_delta_mz);
        }

        // サンプル名の取得
        public List<string> GetSampleName()
        {
            List<string> samplename = new List<string>();

            for (int i = 0; i < output.Count(); i++)
            {
                samplename.Add(output[i].sample);
            }
            return (samplename);
        }




    }
}

// MTのmedianを計算して、出力するようにする