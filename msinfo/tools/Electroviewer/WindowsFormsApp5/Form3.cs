using MSFileReaderLib;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.IO;
using MathNet.Numerics.Statistics;

namespace WindowsFormsApp5
{
    public partial class Form3 : Form
    {
        public String rawfilename;
        public string inputmz;

        public Form3(params string[] argumentValues)
        {
            InitializeComponent();
        }       

        private void Form3_Load(object sender, EventArgs e)
        {
            chart1.ChartAreas.Clear();
            chart2.ChartAreas.Clear();
            chart1.Series.Clear();
            chart2.Series.Clear();            
        }

        // load EIC
        // テキストボックスに入力された値から、データを取得する
        private void button1_Click(object sender, EventArgs e)
        {

            // データが正しく取得された時のみ実行する
            if (textBox1.Text != "" & textBox2.Text != "" & rawfilename != null)
            {

                // 初期化
                chart1.ChartAreas.Clear();
                chart2.ChartAreas.Clear();

                // 初期Areaの設定(描画無し)
                string chart_area0 = "Area";
                string chart_area1 = "Area2";
                               
                chart1.ChartAreas.Add(new ChartArea(chart_area0));

                chart2.ChartAreas.Add(new ChartArea(chart_area1));
                chart2.Series.Clear();

                inputmz = textBox1.Text;
                string eppm = textBox2.Text;
                                               
                chart1.Series.Clear();

                string legend1 = "area1";

                chart1.Series.Add(legend1);

                // グラフの種別を指定
                chart1.Series[legend1].ChartType = SeriesChartType.Line;

                // EICを取得
                MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
                rawfile.Open(rawfilename);
                rawfile.SetCurrentController(0, 1);

                // delta_mzの計算
                Double delta_mz = Double.Parse(inputmz) * (Double.Parse(eppm) / Math.Pow(10, 6));
                
                double mz_low = double.Parse(inputmz) - delta_mz;
                double mz_high = double.Parse(inputmz) + delta_mz;

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

                    chart1.Series[legend1].Points.AddXY(mt, intensity);
                }

                rawfile.Close();

                chart1.Series[legend1].BorderWidth = 1;
                chart1.ChartAreas[0].AxisX.MajorGrid.Enabled = false;
                chart1.ChartAreas[0].AxisY.MajorGrid.Enabled = false;

                chart1.ChartAreas[0].CursorX.IsUserEnabled = true;
                chart1.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;
                chart1.ChartAreas[0].CursorX.Interval = 0.01;

                chart1.ChartAreas[0].CursorY.IsUserEnabled = true;
                chart1.ChartAreas[0].CursorY.IsUserSelectionEnabled = true;

                // ここで表示用のピーク領域を設定している
                chart1.ChartAreas[0].AxisX.ScaleView.ZoomReset();

                // 縦軸の範囲                
                chart1.ChartAreas[0].AxisY.ScaleView.ZoomReset();

                chart1.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart1.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;

                chart1.Series[legend1].Color = Color.FromArgb(0, 0, 0);

                chart1.ChartAreas[0].AxisX.LabelStyle.Format = "#.##";
                chart1.ChartAreas[0].AxisY.LabelStyle.Format = "#";

                chart1.ChartAreas[0].AxisX.ScrollBar.BackColor = Color.White;
                chart1.ChartAreas[0].AxisX.ScrollBar.ButtonColor = Color.Silver;
                chart1.ChartAreas[0].AxisX.ScrollBar.ButtonStyle = ScrollBarButtonStyles.ResetZoom;

                chart1.ChartAreas[0].AxisY.ScrollBar.BackColor = Color.White;
                chart1.ChartAreas[0].AxisY.ScrollBar.ButtonColor = Color.Silver;
                chart1.ChartAreas[0].AxisY.ScrollBar.ButtonStyle = ScrollBarButtonStyles.ResetZoom;

                chart1.ChartAreas[0].BackColor = Color.White;

            }
            else
            {
                MessageBox.Show("No Data");
            }
       
        }

        // グラフをクリックしたら、そのピークトップを表示して、Addする
        // chart2にマススペクトルを表示

        // ピークを加える
        private void Form3_MouseClick(object sender, MouseEventArgs e)
        {

            int mb = 0;

            // 左 or 右クリック
            switch (e.Button)
            {
                case MouseButtons.Left:
                    mb = 1;
                    break;
                case MouseButtons.Middle:
                    mb = 2;
                    break;
                case MouseButtons.Right:
                    mb = 3;
                    break;
            }

            var ca = chart2.ChartAreas[0];

            // 右クリックの処理
            if (mb == 3)
            {
                ca.AxisX.ScaleView.Zoomable = false;
                ca.CursorX.IsUserEnabled = true;
                ca.CursorX.IsUserSelectionEnabled = true;
            }
            // 左クリック
            if (mb == 1)
            {
                ca.AxisX.ScaleView.Zoomable = true;
                ca.CursorX.IsUserEnabled = true;
                ca.CursorX.IsUserSelectionEnabled = true;
            }
        }


        private void chart1_MouseClick(object sender, EventArgs e)
        {

            // ChartAreaが無い時にはエラーになるので、修正する必要がある
            double a = chart1.ChartAreas[0].CursorX.Position;

            if (!Double.IsNaN(a))
            {

                chart2.ChartAreas.Clear();

                string chart_area1 = "Area2";                
                chart2.ChartAreas.Add(new ChartArea(chart_area1));
                               
                // label
                label4.Text = "MT : " + a.ToString() + "min";               
                
                MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
                rawfile.Open(rawfilename);
                rawfile.SetCurrentController(0, 1);

                double RT= a;
                double CentroidPeakWidth = 0.0;
                object MassList = null;
                object PeakFlags = null;
                int ArraySize = 0;

                rawfile.GetMassListFromRT(ref RT, null, 1, 0, 0, 0, ref CentroidPeakWidth, ref MassList, ref PeakFlags, ref ArraySize);
                double[,] mslist = (double[,])MassList;

                // マススペクトル

                chart2.Series.Clear();

                string legend1 ="test";

                chart2.Series.Add(legend1);

                for (int k = 0; k < mslist.GetLength(1); k++)
                {
                    double mz = mslist[0, k];
                    double intensity_mz = mslist[1, k]; 

                    chart2.Series[legend1].Points.AddXY(mz, intensity_mz);
                }

                rawfile.Close();

                // chart2の設定
                chart2.ChartAreas[0].CursorX.IsUserEnabled = true;
                chart2.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;
                chart2.ChartAreas[0].CursorX.Interval = 0.01;

                chart2.ChartAreas[0].CursorY.IsUserEnabled = true;
                chart2.ChartAreas[0].CursorY.IsUserSelectionEnabled = true;


                chart2.Series[legend1].Color = Color.FromArgb(0, 0, 0);

                chart2.ChartAreas[0].AxisX.LabelStyle.Format = "#.####";
                chart2.ChartAreas[0].AxisY.LabelStyle.Format = "#";

                chart2.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart2.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;

                // 入力のm/zの+1、+2、+3辺りの同位体のピークを確認するために拡大する
                chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset();
                chart2.ChartAreas[0].AxisX.ScaleView.Zoom(Double.Parse(inputmz) - 0.005, Double.Parse(inputmz) + 0.005);

                chart2.ChartAreas[0].CursorX.Interval = 0.00001;

                //        double min_right = diff2.Minimum();


                //int ScanNumber = Int32.Parse(args[1]); // scan number

                //double CentroidPeakWidth = 0.0;
                //object MassList = null;
                //object PeakFlags = null;
                //int ArraySize = 0;

                //// ---------------
                ////  Mass Spectra
                //// ---------------
                //rawfile.GetMassListFromScanNum(ref ScanNumber, null, 1, 0, 0, 0, ref CentroidPeakWidth, ref MassList, ref PeakFlags, ref ArraySize);
                //double[,] mslist = (double[,])MassList;

                //MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();

                //string input_raw = args[0];
                //rawfile.Open(input_raw);

                //rawfile.SetCurrentController(0, 1); /* Controller type 0 means mass spec device; Controller 1 means first MS device */

                //int num = 0;
                //rawfile.GetNumSpectra(ref num);

                //// ファイル(output)の存在チェック
                //if (System.IO.File.Exists(@"C:\\temp\\spectra\\list.txt"))
                //{
                //    System.IO.File.Delete(@"C:\\temp\\spectra\\list.txt");
                //}

                //// データ収集
                //for (int j = 1; j <= num; ++j)
                //{

                //    int ScanNumber = j;

                //    // 結果の出力
                //    using (StreamWriter w = new StreamWriter(@"C:\\temp\\spectra\\list.txt", true))
                //    {
                //        double CentroidPeakWidth = 0.0;
                //        object MassList = null;
                //        object PeakFlags = null;
                //        int ArraySize = 0;

                //        rawfile.GetMassListFromScanNum(ref ScanNumber, null, 1, 0, 0, 0, ref CentroidPeakWidth, ref MassList, ref PeakFlags, ref ArraySize);
                //        double[,] mslist = (double[,])MassList;

                //        // MT
                //        double dRT = 0;

                //        rawfile.RTFromScanNum(ScanNumber, ref dRT);

                //        w.Write(ScanNumber);
                //        w.Write("\t");
                //        w.Write(dRT);
                //        w.Write("\n");

                //    }

                //}
            }
        }

        private void OpenRawfileToolStripMenuItem_Click(object sender, EventArgs e)
        {
            // rawfileの読み込み

            DialogResult dr = openFileDialog1.ShowDialog();
            if (dr == System.Windows.Forms.DialogResult.OK)
            {
                rawfilename = openFileDialog1.FileName;
            }

        }
    }
}
