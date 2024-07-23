using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Collections.Specialized;
using System.Threading.Tasks;
using System.IO;

namespace Preprocessing_MCMC_sampling
{
    public partial class BaseModule : Form
    {
        #region variables
        string currentPath = @"C:\Users\Hp\OneDrive - Imam university\Documents\Az projects\PathTurbEr\GSE59357 analysis\data\";
        string controlDataFile = "";
        string caseDataFile = "";
        string pathwaysFile = "";

        double th_dic_high = 1.5;
        double th_dic_low = -1.5;
        double th_MeanPosteriorNet = 0.5;
        int th_num_Ordinal_OutputNetwork = 1;

        List<string> MethodSettings = new List<string>();
        List<string> OutputSettings = new List<string>();

        double gamma_prior_a = 0.1;
        double gamma_prior_b = 0.1;
        int nburnIn = 10000;
        int nSamplingIter = 20000;
        //int nburnIn = 10;
        //int nSamplingIter = 20;
        int nMaxOutdegree = 4;
        int nMaxIndegree = 4;

        #endregion

        public BaseModule()
        {
            InitializeComponent();

            // variable that are subject to be loaded from users local disc
            controlDataFile = currentPath + @"sensitive_group_matrix.csv";
            caseDataFile = currentPath + @"resistant_group_matrix.csv";
            pathwaysFile = currentPath + @"KGMLs\KEGG_selected_pathways.csv";

            // -------------- Load data --------------
            dataLoad dl = new dataLoad();
            Dictionary<String, double[]> pMat_raw = dl.LoadMatrix(controlDataFile);
            Dictionary<String, double[]> rMat_raw = dl.LoadMatrix(caseDataFile);
            Dictionary<String, List<String>> allPathways = dl.LoadPathways(pathwaysFile);

            // -------------- Discretized case (Resistant) data matrix --------------
            processData pD = new processData(th_dic_high,th_dic_low);
            Dictionary<String, int[]> rMat_discritized = pD.getDiscritizedMat(pMat_raw, rMat_raw);
            //Dictionary<String, double[]> rMat_standardizedOnly = pD.getStandardizedVal(pMat_raw, rMat_raw);
            //printData(rMat_standardizedOnly);
            Dictionary<String, Dictionary<String, int[]>> pathwayMats = pD.getDiscritizedMatForPathways(allPathways, rMat_discritized);

            // -------------- Infer BN strucrure (for each pathway) --------------
            //Parallel.ForEach(allPathways.Keys, pathway =>{
            //foreach (string pathway in allPathways.Keys)
            //{
            //int nNodes = pathwayMats[pathway].Count;
            string pathway = "Sphingolipid_signaling_pathway";
            //Dictionary<String, double[]> pathwayStdData = pD.getStandardizedDataforAPathway(pathway, allPathways, rMat_standardizedOnly);

            int nNodes = pathwayMats[pathway].Count;

            // -- prepare MCMC sampling config settings
            SetMCMCconfig();
            NameValueCollection InputSettings = new NameValueCollection();
                        // ------------ input settings
            InputSettings.Add("nNodes", nNodes.ToString());
            InputSettings.Add("nIteration", nSamplingIter.ToString());                // [Tasaki et al.] sampling iteration: 5000, burn-in (inclusive): 1000 (20%)
            InputSettings.Add("nBurnIn", nburnIn.ToString());                   // http://www.genetics.org/content/early/2015/01/28/genetics.114.172619
            InputSettings.Add("nMaxOutdegree", nMaxOutdegree.ToString());
            InputSettings.Add("nMaxIndegree", nMaxIndegree.ToString());
            InputSettings.Add("pathway", pathway);

            // -- NS algorithm (comment out when need to use)
            Adaptive_NS obj_NS = new Adaptive_NS(InputSettings, MethodSettings, OutputSettings, pathwayMats[pathway]);
            List<String> BNs_NS = obj_NS.StartAnalysis();
            string bn_NS = Find_CombinedBN(BNs_NS);
            List<double[,]> bnMats = getNetworkSfromString(bn_NS, nNodes);
            printNet(bn_NS, obj_NS.nodeNameList, pathway, "NS");
            bugsSampling bugs_NS = new bugsSampling(bnMats[0], 
                bnMats[1], 
                nNodes, 
                pathway, 
                obj_NS.nodeNameList, 
                "NS", 
                nburnIn,
                nSamplingIter,
                gamma_prior_a,
                gamma_prior_b);       // [0] original matrix, [1] complementary matrix
            // -----------

            // -- HAR algorithm (comment out when need to use)
            Adaptive_HAR obj_HAR = new Adaptive_HAR(InputSettings, MethodSettings, OutputSettings, pathwayMats[pathway]);
            List<String> BNs_HAR = obj_HAR.StartAnalysis();
            string bn_HAR = Find_CombinedBN(BNs_HAR);
            List<double[,]> bnMats_HAR = getNetworkSfromString(bn_HAR, nNodes);
            printNet(bn_HAR, obj_HAR.nodeNameList, pathway, "HAR");
            bugsSampling bugs_HAR = new bugsSampling(bnMats_HAR[0], 
                bnMats_HAR[1], 
                nNodes, 
                pathway, 
                obj_HAR.nodeNameList, 
                "HAR",
                nburnIn,
                nSamplingIter,
                gamma_prior_a,
                gamma_prior_b);        // [0] original matrix, [1] complementary matrix
            // -----------


            // -- MH algorithm (comment out when need to use)
            Adaptive_MH obj_MH = new Adaptive_MH(InputSettings, MethodSettings, OutputSettings, pathwayMats[pathway]);
            List<String> BNs_MH = obj_MH.StartAnalysis();
            string bn_MH = Find_CombinedBN(BNs_MH);
            List<double[,]> bnMats_MH = getNetworkSfromString(bn_MH, nNodes);
            printNet(bn_MH, obj_MH.nodeNameList, pathway, "MH");

            bugsSampling bugs_MH = new bugsSampling(bnMats_MH[0], 
                bnMats_MH[1], 
                nNodes, 
                pathway, 
                obj_MH.nodeNameList, 
                "MH", 
                nburnIn, 
                nSamplingIter, 
                gamma_prior_a, 
                gamma_prior_b);        // [0] original matrix, [1] complementary matrix
            // -----------
            //}
            //});
            
            MessageBox.Show("complete");
        }

        private void printData(Dictionary<string, double[]> rMat_standardizedOnly)
        {
            TextWriter twDat = new StreamWriter("GSE59357_standardizedOnly.csv");
            
            foreach(string gene in rMat_standardizedOnly.Keys)
            {
                string line = gene + ",";
                double[] arr = rMat_standardizedOnly[gene];
                foreach(double val in arr)
                    line += (val + ",");
                twDat.WriteLine(line.Substring(0, line.Length - 1));
            }
            twDat.Close();
        }

        private void printNet(string str_BN, NameValueCollection nodeNameList, string pathway, string method)
        {
            string[] SampleNetString = str_BN.Replace("\"", "").Split("|".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            string line = "";
            for (int i = 0; i < SampleNetString.Count(); i++)
            {
                string[] arr = ConvertNetStringToArray(SampleNetString[i], nodeNameList.AllKeys.Length);
                for (int j = 0; j < arr.Count(); j++)
                {
                    if (arr[j].Equals("1"))
                    {
                        line += (nodeNameList[i] + "," + nodeNameList[j] + "\n");
                    }
                }
            }
            TextWriter tw = new StreamWriter(pathway + "_" + method + "_net.csv");
            tw.WriteLine(line);
            tw.Close();
        }

        private List<double[,]> getNetworkSfromString(string str_BN, int _nNodes)
        {
            List<double[,]> retList = new List<double[,]>();
            
            #region Initialize arrays
            double[,] mat = new double[_nNodes, _nNodes];
            double[,] c_mat = new double[_nNodes, _nNodes];
            #endregion

            #region fill two matrices
            string[] SampleNetString = str_BN.Replace("\"", "").Split("|".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            for (int i = 0; i < SampleNetString.Count(); i++)
            {
                string[] arr = ConvertNetStringToArray(SampleNetString[i], _nNodes);
                for (int j = 0; j < arr.Count(); j++)
                {
                    if (arr[j].Equals("1")){
                        mat[i,j] = 1;
                        c_mat[i,j] = 0;
                    }else {
                        mat[i,j] = 0;
                        c_mat[i,j] = 1;
                    }
                }
            }

            retList.Add(mat);               // [0] original matrix
            retList.Add(c_mat);             // [1] complementary matrix
            #endregion

            return retList;
        }
        private string[] ConvertNetStringToArray(string st, int _nNodes)
        {
            string[] Net = new string[_nNodes];

            char[] ns = st.ToCharArray();
            for (int j = 0; j < ns.Count(); j++)
            {
                Net[j] = ns[j].ToString();
            }

            return Net;
        }
        private string Find_CombinedBN(List<String> bNs)
        {
            // -- Here implement the classifier combination algorithm
            return bNs[0]; // default: return the first BN
        }

        private void SetMCMCconfig()
        {

            // ------------ method settings
            string initialNetwork = "Random";
            string SamplingDistribution = "DirMul";
            string AcceptanceRatio = "MH";
            MethodSettings.Add(initialNetwork);
            MethodSettings.Add(SamplingDistribution);
            MethodSettings.Add(AcceptanceRatio);

            // ------------ output network settings
            OutputSettings.Add(th_MeanPosteriorNet.ToString());
            OutputSettings.Add(th_num_Ordinal_OutputNetwork.ToString());

        }
    }
}
