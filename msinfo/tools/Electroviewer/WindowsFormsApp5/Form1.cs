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
using MSFileReaderLib;  /* using namespace after add references "XRawfile2_x64.dll" */
using System.IO;
using MathNet.Numerics.Statistics;

namespace WindowsFormsApp5
{

    public partial class Form1 : Form
    {

        private List<RawEic.LocalEic> EIC = new List<RawEic.LocalEic>();
        private List<RawEic.LocalMtIndex> mtindex = new List<RawEic.LocalMtIndex>();
        private List<LoadRaw.LocalEic> EIC0 = new List<LoadRaw.LocalEic>();      
        private List<List<LoadRaw.LocalEic>> AllEIC = new List<List<LoadRaw.LocalEic>>(); // [sample][peakID]
        private List<LoadRaw> test = new List<LoadRaw>();
        private AlignError ae;

        Form3 f = new Form3();

        List<string> rawfiles = new List<string>();
        List<string> samplenames = new List<string>();
        List<string> filenames = new List<string>();
        List<string> peakfiles = new List<string>();

        List<string> filenames_peak = new List<string>();

        List<RawEic.LocalMtIndex> view_range = new List<RawEic.LocalMtIndex>();
            
        // CurrentPeak class
        public class CurrentPeak {

            private int peakID;
            private int samplenum;

            public void SetPeakId(int input_peakID)
            {
                peakID = input_peakID;
            }
            public void SetSampleNum(int input_samplenum)
            {
                samplenum = input_samplenum;
            }
            public int GetPeakID()
            {
                return peakID;
            }
            public int GetSampleNum()
            {
                return samplenum;
            }
        }

        public CurrentPeak cpeak = new CurrentPeak();

        TxtPeak peaks;

        // Form1
        public Form1()
        {
            InitializeComponent();
        }

        // データの読み込み
        private void rAWToolStripMenuItem_Click(object sender, EventArgs e)
        {
            openFileDialog1.Multiselect = true;
            DialogResult dr = openFileDialog1.ShowDialog();
            if (dr == System.Windows.Forms.DialogResult.OK)
            {
                foreach (string file in openFileDialog1.FileNames)
                {
                    rawfiles.Add(file);
                    filenames.Add(System.IO.Path.GetFileName(file));
                }
            }
        }

        // ピークリストファイルを開く
        private void ピークリストToolStripMenuItem_Click(object sender, EventArgs e)
        {
            openFileDialog2.Multiselect = true;
            DialogResult dr = openFileDialog2.ShowDialog();
            if (dr == System.Windows.Forms.DialogResult.OK)
            {
                foreach (string file in openFileDialog2.FileNames)
                {
                    peakfiles.Add(file);
                    filenames_peak.Add(System.IO.Path.GetFileName(file));
                }
            }
        }

        // フォームのロード
        private void Form1_Load(object sender, EventArgs e)
        {
            dataGridView1.RowHeadersVisible = false;
            dataGridView1.ColumnHeadersVisible = false;
            dataGridView1.AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill;
            dataGridView1.ColumnCount = 1;

            chart1.BackColor = Color.Transparent;
            chart2.BackColor = Color.Transparent;
            chart3.BackColor = Color.Transparent;

            chart1.ChartAreas[0].BackColor = Color.Transparent;
            chart2.ChartAreas[0].BackColor = Color.Transparent;
            chart3.ChartAreas[0].BackColor = Color.Transparent;
        }

        // ------------------------------------------------
        //  読み込みボタン : ファイルの読み込み、EIC生成
        // ------------------------------------------------
        private void OpenToolStripMenuItem_Click(object sender, EventArgs e)
        {

            // ピークファイルの読み込み
            string peakfile = "";

            DialogResult dr = openFileDialog1.ShowDialog();
            if (dr == System.Windows.Forms.DialogResult.OK)
            {
                peakfile = openFileDialog1.FileName;
            }

            if (peakfile != "")
            {

                // 設定ファイルのパスの取得
                StreamReader objReader = new StreamReader(peakfile, false);
                string setfilepath = objReader.ReadLine(); // 1行目

                // 設定ファイルを開く
                StreamReader objReader0 = new StreamReader(setfilepath, false);

                //// ファイル読み込み
                Form2 form2 = new Form2();
                form2.Show();

                // headerを捨てる
                objReader0.ReadLine();

                string sLine = "";
                while (sLine != null)
                {
                    sLine = objReader0.ReadLine();
                    if (sLine != null)
                    {
                        string[] temp_line = sLine.Split('\t');
                        rawfiles.Add(temp_line[1]);
                    }
                }

                // -----------------------------
                //      txtファイルを読み込む
                // -----------------------------
                // rawファイルを全部読み込んで、intensityに代入したあと、捨てる
                // chart1をクリックしたときはrawfileを読み込む

                peaks = new TxtPeak();
                peaks.SetPeakFromTxt(peakfile); // オリジナル

                // List<List<TxtPeak.LocalPeakList>> input3 = peaks.GetAllPeakFromTxt();
                int peaknum = peaks.GetPeakNum();  // ピーク数

                List<RawEic.LocalMz> inputmz = new List<RawEic.LocalMz>();
                List<RawEic.LocalMt> inputmt = new List<RawEic.LocalMt>();

                List<List<RawEic.LocalMz>> inputmz3 = peaks.GetMz();
                List<List<RawEic.LocalMt>> inputmt3 = peaks.GetMt();

                form2.progressBar1.Style = ProgressBarStyle.Continuous;
                form2.progressBar1.Maximum = rawfiles.Count() - 1;

                // ----------------------------------
                //      RawファイルからEICの取得
                // ----------------------------------
                List<List<LoadRaw.LocalEic>> AllEIC0 = new List<List<LoadRaw.LocalEic>>();

                test = new List<LoadRaw>();
                LoadRaw test0 = new LoadRaw();

                // 各サンプル
                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    test0 = new LoadRaw();

                    inputmz = new List<RawEic.LocalMz>();
                    inputmt = new List<RawEic.LocalMt>();

                    // 各ピーク
                    for (int j = 0; j < peaknum; j++)
                    {
                        inputmz.Add(inputmz3[i][j]);
                        inputmt.Add(inputmt3[i][j]);
                    }

                    // ピークID1のmt               
                    test0.SetEics(rawfiles[i], inputmz, inputmt);
                    test.Add(test0);

                    EIC0 = new List<LoadRaw.LocalEic>();
                    EIC0 = test[i].GetAllEic();

                    AllEIC.Add(EIC0);

                    form2.progressBar1.Value = i;
                }

                // ピークID
                List<string> id = peaks.GetPeakId();

                //// Table1：サンプル情報
                dataGridView1.AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill;
                dataGridView1.RowHeadersVisible = false;
                dataGridView1.ColumnHeadersVisible = false;

                dataGridView1.ColumnCount = 1; // 表示
                dataGridView1.RowCount = id.Count(); // 表示

                for (int i = 0; i < id.Count(); i++)
                {
                    dataGridView1.Rows[i].Cells[0].Value = id[i].ToString();
                }

                //// EIC表示
                chart1.Series.Clear();
                chart1.ChartAreas.Clear();
                chart2.Series.Clear();
                chart2.ChartAreas.Clear();
                chart3.Series.Clear();
                chart3.ChartAreas.Clear();

                //// 初期Areaの設定(描画無し)
                string chart_area0 = "Area";
                string chart_area1 = "Area2";
                string chart_area2 = "Area3";

                chart1.ChartAreas.Add(new ChartArea(chart_area0));
                chart2.ChartAreas.Add(new ChartArea(chart_area1));
                chart3.ChartAreas.Add(new ChartArea(chart_area2));

                //chart1.ChartAreas[0].;
                chart1.Series.Add("error");

                chart1.ChartAreas[0].CursorX.IsUserEnabled = true;
                chart1.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;

                chart1.ChartAreas[0].BackColor = Color.Transparent;
                chart2.ChartAreas[0].BackColor = Color.Transparent;
                chart3.ChartAreas[0].BackColor = Color.Transparent;

                form2.Close();

                // chart2の表示範囲の設定
                RawEic.LocalMtIndex view_range0;

                for (int i = 0; i < id.Count(); i++)
                {
                    view_range0 = new RawEic.LocalMtIndex();
                    view_range0.left_index = 1;
                    view_range0.top_index = 20;
                    view_range0.right_index = 40;

                    view_range.Add(view_range0);
                }

                button3.Enabled = true; // peak Add button         

            }
        }

        // 
        //private void eICToolStripMenuItem_Click(object sender, EventArgs e)
        //{
            //// EICを保存
            //// ------------------------
            ////// テキストにEICを出力
            //// ------------------------
            ////SaveFileDialogクラスのインスタンスを作成
            //SaveFileDialog sfd = new SaveFileDialog();

            //String outfilename = null;

            //sfd.Filter = "Text Files | *.txt";
            //sfd.CheckPathExists = true;
            //if (sfd.ShowDialog() == DialogResult.OK)
            //{
            //    outfilename = sfd.FileName;
            //}

            //StreamWriter writer = new StreamWriter(outfilename, true); // 出力先固定

            ////// RT書き込み
            //foreach (double str in RT_All[k])
            //{
            //    writer.Write(String.Format(str.ToString()) + "\t");
            //}
            //writer.WriteLine(); // 改行

            //// Intensity書き込み
            //foreach (double str2 in Int_All[k])
            //{
            //    writer.Write(String.Format(str2.ToString()) + "\t");
            //}
            //writer.WriteLine(); // 改行

            //writer.Close();

        //}

        // ---------------------------
        //      chart1をクリック
        // ---------------------------
        private void Form1_MouseDoubleClick(object sender, MouseEventArgs e)
        {
        }

        //// 
        private void Form1_MouseClick(object sender, MouseEventArgs e)
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
            
        // --------------------------------
        //      datagridのセルをクリック
        // --------------------------------
        private void dataGridView1_CellClick(object sender, DataGridViewCellEventArgs e)
        {            

            // chart3の初期化
            chart3.ChartAreas.Clear();
            chart3.Series.Clear();
            string chart_area0 = "Area";
            chart3.ChartAreas.Add(new ChartArea(chart_area0));

            chart3.ChartAreas[0].BackColor = Color.Transparent;

            int rindex = e.RowIndex; // peak index

            List<List<double>> Intensity = new List<List<double>>();

            // データ部分以外をクリックした時
            if (rindex < 0)
            { 
                return; // 処理終了
            }

            // クリア
            chart1.Series.Clear();
            chart2.Series.Clear();
            chart3.Series.Clear();
            
            // 初期Areaの設定(描画無し)
            //string chart_area0 = "Area";

            // ピークファイルの読み込み
            //Peak = new TxtPeak();
            //Peak.SetPeakFromTxt(peakfiles);

            // クリックしたとき
            if (rindex >= 0)
            {

                cpeak.SetPeakId(rindex); // currentpeak
                peaks.SetPeakId(rindex); // peaksにidを設定している

                // -----------------------------------------------------------

                List<TxtPeak.LocalPeakList> peak = peaks.GetPeakFromId(); // ピーク情報の取得           
                List<RawEic.LocalMt> peakmt = peaks.GetPeakmtFromId(); // ピーク情報の取得

                // Averageのm/zのlabel3への表示
                List<double> AllMz = new List<double>();

                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    AllMz.Add(Double.Parse(peak[i].mz));
                }

                label3.Text = "　　　　　m/z : " + AllMz.Mean().ToString(); // 平均にする必要はないが、過去の名残
                label4.Text = "";

                // ---------------------------
                //      要実装
                // ---------------------------
                // chart2の横軸をピークの位置に合わせて拡大しなおす
                // とりあえずzoomをリセット

                chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset();

                // ------------------------------------


                // 領域選択設定
                chart2.ChartAreas[0].CursorX.IsUserEnabled = true;
                // chart2.ChartAreas[0].CursorY.IsUserEnabled = true;
                chart2.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;

                // 特定のピークの全サンプルのintensity
                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    Intensity.Add(AllEIC[i][rindex].intensity);
                }

                // サンプル毎
                for (int i = 0; i < Intensity.Count(); i++)
                {
                    List<double> y = Intensity[i];
                    string legend2 = i.ToString();

                    chart2.Series.Add(legend2);
                    chart2.Series[legend2].ChartType = SeriesChartType.Line;

                    for (int l = 0; l < y.Count() - 1; l++)
                    {
                        chart2.Series[legend2].Points.AddY(y[l]);
                    }

                    chart2.ChartAreas[0].AxisX.Interval = 10;
                    chart2.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                    chart2.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                    chart2.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                    chart2.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;

                    chart2.ChartAreas[0].AxisX.LabelStyle.Enabled = true;
                    chart2.Series[legend2].BorderWidth = 1;

                    chart2.Series[legend2].Color = Color.FromArgb(0, 0, 0);
                }

                // 表示領域
                //chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset(); // これ必ず必要
                //chart2.ChartAreas[0].AxisX.ScaleView.Zoom(Peaks.GetPeakIndex().left_index, Peaks.GetPeakIndex().right_index);

                // 表示領域
                chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset();
                chart2.ChartAreas[0].AxisX.ScaleView.Zoom(view_range[cpeak.GetPeakID()].left_index, view_range[cpeak.GetPeakID()].right_index);
                
                string legend1 = 0.ToString();

                chart1.ChartAreas.Clear();
                chart1.ChartAreas.Add(new ChartArea(chart_area0));

                chart1.Series.Add(legend1);
                chart1.Series[legend1].BorderWidth = 5;
                chart1.Series[legend1].MarkerSize = 10;
                chart1.Series[legend1].ChartType = SeriesChartType.Point;
                chart1.Series[legend1].MarkerStyle = MarkerStyle.Circle;
                chart1.Series[legend1].Color = Color.Black;

                // AlignmentError
                ae = new AlignError();
                ae.SetAllEIC(AllEIC);
                ae.SetLeftIndex(view_range[cpeak.GetPeakID()].left_index);
                ae.SetRightIndex(view_range[cpeak.GetPeakID()].right_index);
                ae.SetPeakID(rindex);

                List<double> align_error = ae.GetAlignError();

                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    chart1.Series[legend1].Points.AddXY(i + 1, align_error[i]);
                }

                chart1.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart1.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart1.ChartAreas[0].AxisX.Interval = 1;

                chart1.ChartAreas[0].CursorX.IsUserEnabled = true;

                chart1.ChartAreas[0].BackColor = Color.White;
                chart2.ChartAreas[0].BackColor = Color.White;

                // Sample name
                label2.Text = "";

            }
        }

        // -------------------
        // 機能していない
        // -------------------
        // キー入力を禁止
        //private void DataGridView1_KeyDown(object sender, KeyEventArgs e)
        //{
        //    if (e.KeyCode == Keys.Up || e.KeyData == Keys.Down)
        //        e.Handled = true;
        //}

        // -----------------------
        //      chart1_Click
        // -----------------------
        // 
        // chart3でピークトップを変更した後、chart1をクリックした時に、chart3のピーク領域が更新されていない
        // chart3を描画するときに、ピーク領域を正しく取得できていない気がする

        private void chart1_Click(object sender, EventArgs e)
        {
            double b = chart1.ChartAreas[0].CursorX.Position;

            if (!Double.IsNaN(b) & b > 0 & b < rawfiles.Count()+1)
            {
                chart2.ChartAreas[0].CursorX.Position = -1;
                chart3.ChartAreas[0].CursorX.Position = -1;

                chart3.Series.Clear();

                int m = (int)b - 1; // 特定のサンプル
                string legend1 = m.ToString();

                cpeak.SetSampleNum(m);

                chart3.Series.Add(legend1);

                // グラフの種別を指定
                chart3.Series[legend1].ChartType = SeriesChartType.Line; // 折れ線グラフを指定してみます

                // ピーク情報
                List<RawEic.LocalMz> inputmz = peaks.GetPeakmzFromId();
                List<RawEic.LocalMt> inputmt = peaks.GetPeakmtFromId();

                //List<RawEic.LocalMz> inputmz = Peak.GetPeakmzFromId();
                //List<RawEic.LocalMt> inputmt = Peak.GetPeakmtFromId();
                
                // EICを再度取得
                RawEic x = new RawEic();

                x.SetEic(rawfiles[m], inputmz[m], inputmt[m]);

                RawEic.LocalEic g = x.GetEic();

                RawEic.LocalMtIndex mtindex0 = new RawEic.LocalMtIndex();

                // MessageBox.Show(x.GetMtIndex().left_index.ToString());

                mtindex0 = x.GetMtIndex();

                // 代入
                double[] x_values = g.mt.ToArray();
                double[] y_values = g.intensity.ToArray();

                // データをシリーズにセットします
                for (int l = 0; l < y_values.Length; l++)
                {
                    chart3.Series[legend1].Points.AddXY(x_values[l], y_values[l]);
                }

                // サンプル数
                for (int i = 0; i < rawfiles.Count; i++)
                {
                    chart2.Series[i].Color = Color.FromArgb(0, 0, 0);
                }

                chart2.Series[m].Color = Color.FromArgb(255, 0, 0);

                chart3.Series[legend1].BorderWidth = 1;
                chart3.ChartAreas[0].AxisX.MajorGrid.Enabled = false;
                chart3.ChartAreas[0].AxisY.MajorGrid.Enabled = false;
                
                chart3.ChartAreas[0].CursorX.IsUserEnabled = true;
                // chart2.ChartAreas[0].CursorY.IsUserEnabled = true;
                chart3.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;
                chart3.ChartAreas[0].CursorX.Interval = 0.001;

                // ここで表示用のピーク領域を設定している
                chart3.ChartAreas[0].AxisX.ScaleView.ZoomReset();
                chart3.ChartAreas[0].AxisX.ScaleView.Zoom(x_values[mtindex0.left_index], x_values[mtindex0.right_index]);           
                
                // 縦軸の範囲                
                chart3.ChartAreas[0].AxisY.ScaleView.ZoomReset();
                double maxint = y_values.ToList().GetRange(mtindex0.left_index, mtindex0.right_index- mtindex0.left_index).ToArray().Max();
                double minint = y_values.ToList().GetRange(mtindex0.left_index, mtindex0.right_index- mtindex0.left_index).ToArray().Min();
                chart3.ChartAreas[0].AxisY.ScaleView.Zoom(minint, maxint);

                chart3.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                chart3.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                chart3.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart3.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;

                chart3.Series[legend1].Color = Color.FromArgb(0, 0, 0);

                chart3.ChartAreas[0].AxisX.LabelStyle.Format = "#.##";
                chart3.ChartAreas[0].AxisY.LabelStyle.Format = "#";

                chart3.ChartAreas[0].AxisX.ScrollBar.BackColor = Color.White;
                chart3.ChartAreas[0].AxisX.ScrollBar.ButtonColor= Color.Silver;
                chart3.ChartAreas[0].AxisX.ScrollBar.ButtonStyle = ScrollBarButtonStyles.ResetZoom;

                // chart3.BackColor = Color.Transparent;
                // chart3.ChartAreas[0].ShadowColor = Color.PaleGreen;

                //chart3.ChartAreas[0].AxisY.ScrollBar.BackColor = Color.White;
                //chart3.ChartAreas[0].AxisY.ScrollBar.ButtonColor = Color.Silver;
                //chart3.ChartAreas[0].AxisY.ScrollBar.ButtonStyle = ScrollBarButtonStyles.ResetZoom;

                chart3.ChartAreas[0].BackColor = Color.White;

                // label4に表示されるピークトップの時間を変更する                
                label4.Text = "　　　PeakTop : " + inputmt[m].topmt.ToString() + "min";

                // samplename
                label2.Text = Path.GetFileName(rawfiles[m]);

            }
        }

        // ------------------------------
        //      ピークの切り直しボタン
        // ------------------------------
        // 上画面の幅修正は、横軸がデータポイント(index)になっているので、MTに直す必要がある
        // leftMTとrightMTを選択し、indexを実際の時間に直して、datagridviewを修正する必要がある

        // 右クリックでピークトップを変更できるようになると良い(実装したい)

        private void button1_Click_1(object sender, EventArgs e)
        {

            //// 上画面
            //double left0 = chart1.ChartAreas[0].CursorX.SelectionStart;
            //double right0 = chart1.ChartAreas[0].CursorX.SelectionEnd;

            //// 下画面
            //double left = chart2.ChartAreas[0].CursorX.SelectionStart;
            //double right = chart2.ChartAreas[0].CursorX.SelectionEnd;

            //// 上画面
            //if (!Double.IsNaN(left0) && !Double.IsNaN(right0) && right0 - left0 > 0 && flag==1){
                
            //    for (int i = 1; i <= rawfiles.Count(); i++)
            //    {
                    
            //        // indexをMTに戻す
            //        // 左端の時間、右端の時間
                    
                    
                    
                    
                    
            //        // Datagridviewの値を修正する
            //        dataGridView2.Rows[i-1].Cells[3].Value = left0.ToString();
            //        dataGridView2.Rows[i-1].Cells[5].Value = right0.ToString();

            //        // right_MTとleft_MTを修正する
                                        
            //        right_MT_All[i-1][rindex+1] = right0.ToString();
            //        left_MT_All[i-1][rindex+1] = left0.ToString();
                    
            //    }
                
            //    // 表示位置を修正する
            //    chart1.ChartAreas[0].AxisX.Minimum = left0;
            //    chart1.ChartAreas[0].AxisX.Maximum = right0;
            //}
            //// 下画面
            //else if (!Double.IsNaN(left) && !Double.IsNaN(right) && right - left > 0 && flag==2)
            //{

            //    // Datagridviewの値を修正する
            //    dataGridView2.Rows[cindex-1].Cells[3].Value = left.ToString();
            //    dataGridView2.Rows[cindex-1].Cells[5].Value = right.ToString();

            //    right_MT_All[cindex-1][rindex+1] = right.ToString(); // ここがミス
            //    left_MT_All[cindex-1][rindex+1] = left.ToString();

            //    // 表示位置を修正する
            //    chart2.ChartAreas[0].AxisX.Minimum = left;
            //    chart2.ChartAreas[0].AxisX.Maximum = right;

            //}
            //else
            //{
            //    // 処理しない
            //}

        }

        // ------------------
        //    Save
        // ------------------
        // ヘッダとサンプル名を出力するように変更
        private void button3_Click(object sender, EventArgs e)
        {

            // Areaの計算
            List<List<double>> area_All = new List<List<double>>();
            List<List<double>> intensity_All = new List<List<double>>();

            // 出力用
            List<double> PID_sv = new List<double>();
            List<List<string>> sample_All = new List<List<string>>();
            List<List<RawEic.LocalMz>> mz_All = peaks.GetMz();
            List<List<RawEic.LocalMt>> mt_All = peaks.GetMt();
                       
            List<double> MT_sv = new List<double>();

            List<string> peakID = peaks.GetPeakId();
            List<string> sample = peaks.GetSampleName();
            List<double> intensity = new List<double>();
            List<double> area = new List<double>();

            // 各ピーク
            for (int i = 1; i <= dataGridView1.RowCount; i++)
            {

                // 初期化
                area = new List<double>();
                intensity = new List<double>();

                // データの取得(各サンプル)
                for (int j = 1; j <= rawfiles.Count; j++)
                {

                    int index1 = view_range[i - 1].left_index;
                    int index2 = view_range[i - 1].right_index;
                                        
                    double[] rtrange = AllEIC[j - 1][i - 1].mt.GetRange(index1, index2 - index1).ToArray(); // ここでエラーになる
                    double[] intrange = AllEIC[j - 1][i - 1].intensity.GetRange(index1, index2 - index1).ToArray();

                    // Areaの計算
                    double area_raw = 0;
                    for (int k = 0; k <= rtrange.Length - 2; k++)
                    {
                        area_raw = area_raw + (intrange[k] + intrange[k + 1]) * (rtrange[k + 1] - rtrange[k]) / 2;
                    }

                    double backarea = (intrange[0] + intrange[rtrange.Length - 1]) * (rtrange[rtrange.Length - 1] - rtrange[0]) / 2;
                    area.Add(area_raw - backarea);

                    // 最大intensity
                    intensity.Add(intrange.Max());
                }
                intensity_All.Add(intensity);
                area_All.Add(area);
                sample_All.Add(sample);
            }

            //// 保存場所(カレントディレクトリ)
            string stCurrentDir = System.IO.Directory.GetCurrentDirectory();

            //// 結果の出力(csvに保存)
            Encoding sjisEnc = Encoding.GetEncoding("Shift_JIS");
            StreamWriter writer = new StreamWriter(stCurrentDir + "\\Test.txt", false, sjisEnc);

            writer.Write("peakID");
            writer.Write("\t");

            for (int j = 1; j <= rawfiles.Count; j++)
            {
                writer.Write("m/z"+" ("+sample[j-1].ToString() + ")");
                writer.Write("\t");
                writer.Write("MT" + " (" + sample[j - 1].ToString() + ")");
                writer.Write("\t");
                writer.Write("Intensity" + " (" + sample[j - 1].ToString() + ")");
                writer.Write("\t");
                writer.Write("Area" + " (" + sample[j - 1].ToString() + ")");
                writer.Write("\t");
            }
            writer.Write("\n");

            for (int i = 0; i <= peakID.Count - 1; i++) {
                writer.Write(peakID[i].ToString());
                writer.Write("\t");

                for (int j = 1; j <= rawfiles.Count; j++) {
                    writer.Write(mz_All[j - 1][i].mz.ToString());
                    writer.Write("\t");
                    writer.Write(mt_All[j - 1][i].topmt.ToString());
                    writer.Write("\t");
                    writer.Write(intensity_All[i][j-1].ToString());
                    writer.Write("\t");
                    writer.Write(area_All[i][j - 1].ToString());
                    writer.Write("\t");
                }
                writer.Write("\n");
            }

                writer.Close();

                MessageBox.Show(stCurrentDir);
                MessageBox.Show("保存しました");
                
            
        }


        private void chart3_Click(object sender, EventArgs e)
        {
            double a = chart3.ChartAreas[0].CursorX.Position;

            if (!Double.IsNaN(a))
            {
                label4.Text = "　　　　　" + a.ToString() + "min";
            }
        }

        private void chart3_DoubleClick(object sender, EventArgs e)
        {
            double a = chart3.ChartAreas[0].CursorX.Position;

            if (!Double.IsNaN(a))
            {
                label4.Text = "　　　PeakTop : " + a.ToString() + "min";

                // ピークトップの変更
                peaks.SetPeakTopMT(cpeak.GetSampleNum(), cpeak.GetPeakID(), a);                         

                // ピーク領域のEICの変更
                RawEic.LocalMt inputmt;
                inputmt = peaks.GetMt()[cpeak.GetSampleNum()][cpeak.GetPeakID()];
                RawEic.LocalMz inputmz;
                inputmz = peaks.GetMz()[cpeak.GetSampleNum()][cpeak.GetPeakID()];
               
                test[cpeak.GetSampleNum()].SetEic(rawfiles[cpeak.GetSampleNum()], cpeak.GetPeakID(), inputmz, inputmt);

                

                // ------------------
                //    chart2を更新
                // ------------------
                chart2.Series.Clear();

                chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset();

                // 領域選択設定
                chart2.ChartAreas[0].CursorX.IsUserEnabled = true;
                // chart2.ChartAreas[0].CursorY.IsUserEnabled = true;
                chart2.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;

                List<List<double>> Intensity = new List<List<double>>();

                // 特定のピークの全サンプルのintensityを取得
                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    List<LoadRaw.LocalEic> ppp = test[i].GetAllEic(); // 特定サンプルの全ピーク
                    List<double> ggg = ppp[cpeak.GetPeakID()].intensity;
                    Intensity.Add(ggg);
                }                            

                // 表示
                for (int i = 0; i < Intensity.Count(); i++)
                {
                    List<double> y = Intensity[i];
                    string legend2 = i.ToString();

                    chart2.Series.Add(legend2);
                    chart2.Series[legend2].ChartType = SeriesChartType.Line;

                    for (int l = 0; l < y.Count() - 1; l++)
                    {
                        chart2.Series[legend2].Points.AddY(y[l]);
                    }

                    chart2.ChartAreas[0].AxisX.Interval = 10;
                    chart2.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                    chart2.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                    chart2.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                    chart2.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;

                    chart2.ChartAreas[0].AxisX.LabelStyle.Enabled = true;
                    chart2.Series[legend2].BorderWidth = 1;

                    chart2.Series[legend2].Color = Color.FromArgb(0, 0, 0);

                    // 赤色
                    if (i == cpeak.GetSampleNum())
                    {
                        chart2.Series[legend2].Color = Color.FromArgb(255, 0, 0);
                    }
                    
                }

                // ------------------
                //    chart1を更新
                // ------------------
                string legend1 = 0.ToString();

                chart1.ChartAreas.Clear();
                chart1.Series.Clear();
                string chart_area0 = "Area";
                chart1.ChartAreas.Add(new ChartArea(chart_area0));

                chart1.Series.Add(legend1);
                chart1.Series[legend1].BorderWidth = 5;
                chart1.Series[legend1].MarkerSize = 10;
                chart1.Series[legend1].ChartType = SeriesChartType.Point;
                chart1.Series[legend1].MarkerStyle = MarkerStyle.Circle;
                chart1.Series[legend1].Color = Color.Black;
                
                AllEIC = new List<List<LoadRaw.LocalEic>>();
                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    AllEIC.Add(test[i].GetAllEic()); // 特定サンプルの全ピーク
                }
                
                // AlignmentError
                ae = new AlignError();
                ae.SetAllEIC(AllEIC);
                ae.SetLeftIndex(view_range[cpeak.GetPeakID()].left_index);
                ae.SetRightIndex(view_range[cpeak.GetPeakID()].right_index);
                ae.SetPeakID(cpeak.GetPeakID());

                List<double> align_error = ae.GetAlignError();

                for (int i = 0; i < rawfiles.Count(); i++)
                {
                    chart1.Series[legend1].Points.AddXY(i + 1, align_error[i]);
                }

                chart1.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
                chart1.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart1.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
                chart1.ChartAreas[0].AxisX.Interval = 1;

                chart1.ChartAreas[0].CursorX.IsUserEnabled = true;

                chart1.ChartAreas[0].BackColor = Color.White;




                // 表示領域
                //chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset(); // これ必ず必要
                //chart2.ChartAreas[0].AxisX.ScaleView.Zoom(Peaks.GetPeakIndex().left_index, Peaks.GetPeakIndex().right_index);

                //MessageBox.Show(peaks.inputmt3[cpeak.GetSampleNum()][cpeak.GetPeakID()].topmt.ToString());
                //MessageBox.Show("Set Peaktop MT as " + a.ToString() + "min.");

                // chart1の再描画
                // AlignmentErrorを再計算
                // chart2の再描画

            }
        }

        private void peaktopToolStripMenuItem_Click(object sender, EventArgs e)
        {

            // ----------------------------------------
            //  幅修正のボタンの実行の処理を記載する
            // ----------------------------------------

            // chart3_Clickのクリック位置をピークトップに変更する
            // chart2のクロマトのピークトップを変更する
            // chart2を再描画
            // chart1を再描画

            // 表示されているmtを取得

            // chart2.ChartAreas[0].CursorX
            double left = chart2.ChartAreas[0].AxisX.ScaleView.ViewMinimum;
            double right = chart2.ChartAreas[0].AxisX.ScaleView.ViewMaximum;

            // chart2の表示範囲を変更する
            view_range[cpeak.GetPeakID()].left_index = (int)left;
            view_range[cpeak.GetPeakID()].right_index = (int)right;

            chart2.ChartAreas[0].AxisX.ScaleView.ZoomReset();
            chart2.ChartAreas[0].AxisX.ScaleView.Zoom(view_range[cpeak.GetPeakID()].left_index, view_range[cpeak.GetPeakID()].right_index);

            // MessageBox.Show("Set PeakPoint from " + left.ToString() + "to " + right.ToString() + ".");
            ae.SetLeftIndex((int)left);
            ae.SetRightIndex((int)right);

            // chart1
            chart1.ChartAreas.Clear();
            chart1.Series.Clear();

            string chart_area0 = "Area";
            chart1.ChartAreas.Add(new ChartArea(chart_area0));

            string legend1 = 0.ToString();
            chart1.Series.Add(legend1);


            List<double> align_error = ae.GetAlignError();

            for (int i = 0; i < rawfiles.Count(); i++)
            {
                chart1.Series[legend1].Points.AddXY(i + 1, align_error[i]);
            }

            // chart1の表示変更
            chart1.Series[legend1].BorderWidth = 5;
            chart1.Series[legend1].MarkerSize = 10;
            chart1.Series[legend1].ChartType = SeriesChartType.Point;
            chart1.Series[legend1].MarkerStyle = MarkerStyle.Circle;
            chart1.Series[legend1].Color = Color.Black;

            chart1.ChartAreas[0].AxisX.MajorGrid.Enabled = true;
            chart1.ChartAreas[0].AxisY.MajorGrid.Enabled = true;
            chart1.ChartAreas[0].AxisX.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
            chart1.ChartAreas[0].AxisY.MajorGrid.LineDashStyle = ChartDashStyle.Dot;
            chart1.ChartAreas[0].AxisX.Interval = 1;

            chart1.ChartAreas[0].CursorX.IsUserEnabled = true;

            chart1.ChartAreas[0].BackColor = Color.White;

        }

        private void chart2_Click(object sender, EventArgs e)
        {

            var ca2 = chart2.ChartAreas[0];

            ////    // chart2のカーソルの選択を可能にする
            ca2.CursorX.Interval = 0.01;

            ca2.CursorX.IsUserEnabled = true;
            ca2.CursorX.IsUserSelectionEnabled = true;

            ////    // カーソルの位置取得
            var w21 = chart2.ChartAreas[0].CursorX.SelectionStart;
            var w22 = chart2.ChartAreas[0].CursorX.SelectionEnd;

            

        }

        private void tableLayoutPanel1_Paint(object sender, PaintEventArgs e)
        {

        }

        private void dataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void MScheckToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            Form3 f = new Form3();
            f.Show();
            // f.rawfiles = rawfiles;
        }

        // DIAテスト用
        private void dIAToolStripMenuItem_Click(object sender, EventArgs e)
        {
            // 
            IXRawfile5 rawtest = (IXRawfile5)new MSFileReader_XRawfile();

            //IXRawfile rawtest = (IXRawfile)new MSFileReader_XRawfile();

            string testfile = "C:\\data\\dia\\Anion_DIA_pro.raw";
            rawtest.Open(testfile);
            rawtest.SetCurrentController(0, 1);

            //double aaa = -1;
            //double bbb = -1;
            //int ccc = -1;
            //rawtest.GetPrecursorRangeForScanNum(1, 2, ref aaa, ref bbb,ref ccc);
            
            //MessageBox.Show(aaa.ToString());
            //MessageBox.Show(bbb.ToString());
            //MessageBox.Show(ccc.ToString());

            //int first_scan_number = -1;
            //rawtest.GetFirstSpectrumNumber(ref first_scan_number);

            //MessageBox.Show(first_scan_number.ToString());

            //double time = -1;
            //rawtest.GetEndTime(ref time);
            //MessageBox.Show(time.ToString());

            //double aaa=-1;
            //rawtest.GetCollisionEnergyForScanNum(100,2,ref aaa);
            //MessageBox.Show(aaa.ToString());

            //double bbb = -1;
            //rawtest.GetIsolationWidthForScanNum(1, 2, ref bbb);

            // MessageBox.Show(bbb.ToString());

            int cc = -1;
            rawtest.GetNumberOfSourceFragmentsFromScanNum(1, ref cc);
            MessageBox.Show(cc.ToString());


            //rawtest.


            //MessageBox.Show(cc.ToString());

            //rawtest.GetIsolationWidthForScanNum(0, 2, ref bbb);
            //MessageBox.Show(bbb.ToString());

        }
    }
}