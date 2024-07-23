using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections.Specialized;
using System.IO;
using System.Diagnostics;
using System.Data;

namespace Preprocessing_MCMC_sampling
{
    class Adaptive_MH
    {
        #region Global variable
        int nNodes;
        ulong nMaxIteration=20;
        ulong nBurnInIteration=10;

        Random rand;
       
        NameValueCollection nParent_th = new NameValueCollection();
        NameValueCollection nChild_th = new NameValueCollection();
        int nMaxParent = 12;
        int nMaxChild = 12;



        int nObservation;


        /* this namevalue collectioin saves the sampled network so far */
        NameValueCollection SampledNetwork = new NameValueCollection();
        NameValueCollection SampledNetwork_Fx = new NameValueCollection();

        /* these collections saves Addable-nonAddable, Deletable-nonDeletable, and Reversible-nonReversible edges */
        Dictionary<string, string> A_list = new Dictionary<string, string>();
        Dictionary<string, string> D_list = new Dictionary<string, string>();
        Dictionary<string, string> R_list = new Dictionary<string, string>();

        Dictionary<String, int[]> pathwayDataMat = new Dictionary<string, int[]>();

        //NameValueCollection Z_forgetList = new NameValueCollection();
        //NameValueCollection Z_neighbourhoodList = new NameValueCollection();
        NameValueCollection DataCategories = new NameValueCollection();
        double[] ProbabilityofAllNodes = null;
        NameValueCollection L_conditionalProb = new NameValueCollection(); // probabilities of changed nodes in Z compared to X

        int[][] DataSet = null;
        TextWriter twFx = null;
        string pathway = "";

        #region public variable to use outside of the class
        public List<double> ListLogFx = new List<double>();
        public Dictionary<double, double> ListSSD = new Dictionary<double, double>();
        public NameValueCollection List_TD_vs_ED = new NameValueCollection();

        public NameValueCollection nodeNameList = new NameValueCollection();
        public DataTable MeanPostNetwork = null;
        public DataTable HighestFreqNetwork = null;
        public ulong iter = 0;
        #endregion


        double _global_LogFx = 0.0;
        double _log_epsilon = -34.538776395; // this value comes from log_e(10e-16), which is 16 digit IEEE double regarding the level of precision we want

        string InputFile = "";
        string DataFileName = "";

        string generalFilePath = "";

        //string targetFunction = "Fx_1_";
        string targetFunction = "";
        Stopwatch watch = null;


        /* declare method parameters */
        string ch_InitNet = "Random";
        string ch_SampDist = "DirMul";
        string ch_AcceptRatio = "MH";


        /* declare output parameters */
        double th_MeanPosteriorNet = 0.5;
        int th_num_Ordinal_OutputNetwork = 1;

        #endregion


        //public Adaptive_MH(int numberOfNodes, List<string> methodSettings, List<string> outputSettings, Dictionary<String, int[]> pathwayDataMat)
        public Adaptive_MH(NameValueCollection inputSettings, List<string> methodSettings, List<string> outputSettings, Dictionary<String, int[]> pathwayDataMat)
        {
            #region Set Pathway data matrix
            //// -------- Pathway data matrix
            //this.pathwayDataMat = pathwayDataMat;
            //#endregion

            //#region Set all the Settings parameter values
            //nNodes = numberOfNodes;

            //generalFilePath = "";

            // -------- Pathway data matrix
            this.pathwayDataMat = pathwayDataMat;
            #endregion

            #region Set all the Settings parameter values
            nNodes = Int32.Parse(inputSettings["nNodes"]);
            nMaxIteration = (ulong)Int32.Parse(inputSettings["nIteration"]);
            nBurnInIteration = (ulong)Int32.Parse(inputSettings["nBurnIn"]);
            nMaxChild = Int32.Parse(inputSettings["nMaxOutdegree"]);
            nMaxParent = Int32.Parse(inputSettings["nMaxIndegree"]);

            generalFilePath = "";
            pathway = inputSettings["pathway"];
            Console.WriteLine(inputSettings["pathway"]);
            Directory.CreateDirectory("BNMCMC_output");
            twFx = new StreamWriter(@"BNMCMC_output\\LogFx_MH_" + pathway + ".csv");

            /* set Methods settings */
            ch_InitNet = methodSettings[0];
            ch_SampDist = methodSettings[1];
            ch_AcceptRatio = methodSettings[2];

            // set output settings
            this.th_MeanPosteriorNet = double.Parse(outputSettings[0]);
            this.th_num_Ordinal_OutputNetwork = Int32.Parse(outputSettings[1]);
            #endregion
        }

        public List<string> StartAnalysis()
        {
            #region Output File Names
            //DataFileName = DataFileName.Substring(0, DataFileName.LastIndexOf("."));
            string freqFileName = generalFilePath + "freq.csv";
            string squredDiffFileName = generalFilePath + "ssd.csv";
            string logValueFileName = generalFilePath + "fx.csv";
            string log_of_Execution = generalFilePath + "logEx.csv";

            
            string lineSquaredDiff = "log_Iterations,squaredDiff\n";
            ProbabilityofAllNodes = new double[nNodes];
            #endregion

            #region Read The data File
            DataSet = new int[nNodes][];
            int dat_i_indx = 0;
            foreach (String pathwayGene in pathwayDataMat.Keys)
            {
                DataCategories.Add(pathwayGene, "3");       // we've only 3 states for each node: -1 (under-expression), 0 (neutral), 1 (over-expression)            
                /* inset the node names int the list */
                nodeNameList.Add(dat_i_indx.ToString(), pathwayGene);
                nObservation = pathwayDataMat[pathwayGene].Count();
                DataSet[dat_i_indx] = new int[nObservation];
                pathwayDataMat[pathwayGene].CopyTo(DataSet[dat_i_indx], 0);
                dat_i_indx++;
            }
            #endregion

            #region Initializing Randomization
            string cpu_cycle = ((int)System.Diagnostics.Stopwatch.GetTimestamp()).ToString();
            int seed = Int32.Parse(cpu_cycle.Substring(cpu_cycle.Length - 3, 2));
            rand = new Random(seed);
            //rand = new Random();
            watch = Stopwatch.StartNew();
            #endregion

            #region Step 1. Construct and Initialize the First Network, "N"

            #region Make a connected DAG graph with nNodes
            int[][] net = new int[nNodes][];
            bool[][] flag = new bool[nNodes][];
            int[][] T_x = new int[nNodes][];


            if (ch_InitNet.Equals("Random"))
            {
                initializeParentChildths();
                FindaRandomNetwork(ref net);
            }
            else
                UseFixedInitialNetwork(ref net);

            #endregion

            #endregion


            #region Step 2. Apply Metropolis-Hasting Sampler with Fixed # of iteration
            int[][] X = new int[nNodes][];
            for (int i = 0; i < nNodes; i++)
                X[i] = new int[nNodes];

            int Mu_X = 0;
            int Mu_Y = 0;
            double fx = 0.0F;
            double fy = 0.0F;

            for (; iter < nMaxIteration; iter++)
            {
                //var watch2 = Stopwatch.StartNew();

                /* 0 as the first argument is for fake "percentage" reporting, but iter is for set as "UserState" argument*/
                //backgroundWorker1.ReportProgress(0, iter);

                // ---- initialize the parent_th and child_th collections
                initializeParentChildths();

                if (iter == 0)
                {
                    X = net;
                    /* populate all the lists initially */
                    PopulateAllLists(X);

                    if (ch_SampDist.Equals("DirMul"))
                        BuildandCalculateProbabilitiesforEachNodes(X);
                }

                #region Handle 'X' & Find 'N_x'

                #region Save all the lists ('A_list', 'D_list', 'R_list') associated with the network 'Y'
                Dictionary<string, string> temp_D_list = new Dictionary<string, string>();
                Dictionary<string, string> temp_A_list = new Dictionary<string, string>();
                Dictionary<string, string> temp_R_list = new Dictionary<string, string>();

                Save_All_Lists(ref temp_D_list, ref temp_A_list, ref temp_R_list);
                #endregion

                int X_totalNumofNeighbours = 0;
                Dictionary<string, string> X_neighbourhoodList = new Dictionary<string, string>();

                /* Insert X itself in its neighbourhood list, as it is a neighbour of itself */
                /* Note, first element is this list is the network itself, so the neighbourhoodID '0' (usually chosen Randomely) indicates itself */
                X_neighbourhoodList.Add("Nil", "No_op");
                X_totalNumofNeighbours++;

                /* sequential executions of "Find_OperableEdges" function 3times will populate all the neighbourhoods in "X_neighbourhoodList" list */
                int num_X_DeletableEdges = Find_OperableEdges(D_list, ref X_neighbourhoodList, "Del");
                int num_X_AddableEdges = Find_OperableEdges(A_list, ref X_neighbourhoodList, "Add");
                int num_X_ReversableableEdges = Find_OperableEdges(R_list, ref X_neighbourhoodList, "Rev");

                /*counting number of neighbourhoods */
                X_totalNumofNeighbours += num_X_DeletableEdges + num_X_AddableEdges + num_X_ReversableableEdges;

                #endregion

                #region Calculate Mu_x for the neighbours of X
                Mu_X = X_totalNumofNeighbours;
                #endregion

                #region Calculate f(x)
                /* will be replaced by some function which calculates f(x)*/
                if (ch_SampDist.Equals("Uniform"))
                    fx = 1.0F;
                else
                    fx = Calcuate_Fx(X);
                #endregion

                #region Draw u_x
                //double upBoundX = (double)(fx / ((double)Mu_x)); // 235 and 236 lines for f(x)=1
                //double Ux = GetRandomeNumber(0.0F, upBoundX);


                double Ux = rand.NextDouble(); // this is for Target function
                #endregion

                #region Step 2.2 Sample "Y" randomly as a neighbourhood of X^(t)

                #region Sample Y_ID
                /* Sample Y from N(X): this random number specifies any integer in the range (0, X_totalNumofNeighbours); Note, +1 indicates the rand function here excludes the last number*/
                int Y_neighbourID = rand.Next(0, X_totalNumofNeighbours);
                #endregion

                #region make Y network as a copy of X^t
                /* make Y network as a copy of X^t */
                int[][] Y = new int[nNodes][];
                for (int i = 0; i < nNodes; i++)
                {
                    Y[i] = new int[nNodes];
                    for (int j = 0; j < nNodes; j++)
                    {
                        if (i == j)
                            Y[i][j] = 0;
                        else
                        {
                            Y[i][j] = X[i][j];
                        }
                    }
                }
                #endregion

                #region Generate Y with proper operation AND Update all the lists: A_list, D_list, R_list
                GenerateNeighbour(ref Y, X_neighbourhoodList, Y_neighbourID);
                #endregion


                #region Handle 'Y' & Find 'N_Y'
                int Y_totalNumofNeighbours = 0;
                Dictionary<string, string> Y_neighbourhoodList = new Dictionary<string, string>();

                /* Insert 'Y' itself in its neighbourhood list, as it is a neighbour of itself */
                /* Note, first element is this list is the network itself, so the neighbourhoodID '0' (usually chosen Randomely) indicates itself */
                Y_neighbourhoodList.Add("Nil", "No_op");
                Y_totalNumofNeighbours++;

                /* sequential executions of "Find_OperableEdges" function 3times will populate all the neighbourhoods in "Y_neighbourhoodList" list */
                int num_Y_DeletableEdges = Find_OperableEdges(D_list, ref Y_neighbourhoodList, "Del");
                int num_Y_AddableEdges = Find_OperableEdges(A_list, ref Y_neighbourhoodList, "Add");
                int num_Y_ReversableableEdges = Find_OperableEdges(R_list, ref Y_neighbourhoodList, "Rev");

                /*counting number of neighbourhoods */
                Y_totalNumofNeighbours += num_Y_DeletableEdges + num_Y_AddableEdges + num_Y_ReversableableEdges;

                #endregion

                #region Calculate Mu_Y for the neighbours of Y
                Mu_Y = Y_totalNumofNeighbours;
                #endregion

                #region Calculate f(y)
                if (ch_SampDist.Equals("Uniform"))
                    fy = 1.0F;
                else
                    fy = Calcuate_Fy(X, X_neighbourhoodList, Y_neighbourID, Y);
                #endregion

                #region Accpetance step: Acceptance or Rejection for the new sample
                if ((fy - fx + Math.Log(Mu_X) - Math.Log(Mu_Y)) >= Math.Log(Ux))
                {
                    if (ch_AcceptRatio.Equals("MH"))
                    {
                        #region Accept the sample: copy tempX into X
                        X = Y;
                        #endregion

                        #region Write the final Sample info

                        /*update the Global list of conditionalProbability values of nodes */
                        for (int indx = 0; indx < nNodes; indx++)
                        {
                            if (L_conditionalProb[indx.ToString()] != null)
                                ProbabilityofAllNodes[indx] = double.Parse(L_conditionalProb[indx.ToString()]);
                        }

                        if (iter >= nBurnInIteration)
                        {

                            string sXFinal = MakeNetworkString(X);
                            SampledNetwork.Add(sXFinal, "S");

                            if (SampledNetwork_Fx[sXFinal] == null)
                                SampledNetwork_Fx.Add(sXFinal, _global_LogFx.ToString());


                            //if (ch_LogPosteriori.Equals("true"))
                            //{
                            //    int networkID = FindNetworkID(SampledNetwork, sXFinal);
                            //    twFx.WriteLine(fy.ToString() + "," + networkID.ToString());
                            //    ListLogFx.Add(fy);
                            //}
                            //ListLogFx.Add(fy);
                        }
                        twFx.WriteLine(fx.ToString() + "," + iter.ToString());
                        ListLogFx.Add(fx);

                        #endregion
                    }
                }
                else
                {
                    #region Write the final Sample info

                    if (iter >= nBurnInIteration)
                    {
                        string sXFinal = MakeNetworkString(X);
                        SampledNetwork.Add(sXFinal, "S");

                        if (SampledNetwork_Fx[sXFinal] == null)
                            SampledNetwork_Fx.Add(sXFinal, fx.ToString());

                        //if (ch_LogPosteriori.Equals("true"))
                        //{
                        //    int networkID = FindNetworkID(SampledNetwork, sXFinal);
                        //    twFx.WriteLine(fx.ToString() + "," + networkID.ToString());
                        //    ListLogFx.Add(fx);
                        //}
                        //ListLogFx.Add(fx);
                    }
                    twFx.WriteLine(fx.ToString() + "," + iter.ToString());
                    ListLogFx.Add(fx);

                    #endregion

                    #region Restores all the lists ('A_list', 'D_list', 'R_list') associated with the network 'Y'
                    Restore_All_Lists(temp_D_list, temp_A_list, temp_R_list);
                    #endregion
                }
                #endregion

                #endregion

                //if (iter % stepSize == 0 && iter >= nBurnInIteration)
                //{
                //    if (ch_SSD.Equals("true"))
                //    {
                //        string str_SSD = FindSquaredDiff();
                //        lineSquaredDiff += Math.Log10((double)iter) + "," + str_SSD + "\n";
                //        if (!double.IsInfinity(Math.Abs(Math.Log10((double)iter))))
                //            ListSSD.Add(Math.Log10((double)iter), double.Parse(str_SSD));
                //    }
                //}

                //backgroundWorker1.ReportProgress(iter);

                


            }
            #endregion


            #region Step 3. Output Time
            #region Find the highest frequency network with its frequency
            List<string> retList = Find_Network_with_Highest_Frequency();
            twFx.Close();
            return retList;
            #endregion
            #endregion

        }

        private void initializeParentChildths()
        {
            nParent_th.Clear();
            nChild_th.Clear();

            for(int i=0; i < nNodes; i++)
            {
                nParent_th.Add(i.ToString(), (nMaxParent + 1).ToString());
                nChild_th.Add(i.ToString(), (nMaxChild + 1).ToString());
            }
        }

        private string[] ConvertNetStringToArray(string st)
        {
            string[] Net = new string[nNodes];

            char[] ns = st.ToCharArray();
            for (int j = 0; j < ns.Count(); j++)
            {
                Net[j] = ns[j].ToString();
            }

            return Net;
        }

        private List<string> Find_Network_with_Highest_Frequency()
        {
            List<string> retList = new List<string>();
            Dictionary<string, int> SampledNetwork2 = new Dictionary<string, int>();

            #region find the frequencies of all the network
            foreach (string key in SampledNetwork.AllKeys)
            {
                int cnt = SampledNetwork[key].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).Count();
                SampledNetwork2.Add(key, cnt);
            }
            #endregion

            /* Sort the dictionary based on the key (in descending order) */
            List<KeyValuePair<string, int>> sortedList = SampledNetwork2.ToList();
            sortedList.Sort((pair1, pair2) => pair2.Value.CompareTo(pair1.Value));

            int topValue = sortedList[0].Value;

            for (int i = 0; i < sortedList.Count; i++){
                if(sortedList[i].Value == topValue) // only have networks with highest frequency
                    retList.Add(sortedList[i].Key);
            }
            return retList;
        }

        private NameValueCollection Find_TD_vs_ED_for_all_sampled_Networks()
        {
            double MaxVal = -9999999.0;
            double MinVal = 9999999.0;
            NameValueCollection ProbFromFrequency = new NameValueCollection();
            NameValueCollection ProbFromFx = new NameValueCollection();
            NameValueCollection normProbFromFx = new NameValueCollection();
            NameValueCollection _TD_ED = new NameValueCollection();
            NameValueCollection SampledNetwork2 = new NameValueCollection();
            Dictionary<string, double> tempDic = new Dictionary<string, double>();  // for sorting the network based on frequency probabilities (ascending order)

            #region Make Sampled Network Table_2
            foreach (string key in SampledNetwork.AllKeys)
            {
                int cnt = SampledNetwork[key].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).Count();
                SampledNetwork2.Add(key, cnt.ToString());
            }
            #endregion

            #region Calculate Probabilities for each network
            foreach (string snet in SampledNetwork2.AllKeys)
            {
                #region Calculate Probablilities from Freq
                double cProb = (double.Parse(SampledNetwork2[snet]) / (double)nMaxIteration);
                ProbFromFrequency.Add(snet, cProb.ToString());
                tempDic.Add(snet, cProb);
                #endregion

                #region Calculate Probabilities from Fx
                //double fx = CalcuateTargetFunction(net);
                double fx = double.Parse(SampledNetwork_Fx[snet]);

                //fx = Math.Exp(fx);
                //fx = Math.Abs(fx);
                if (fx < MinVal)
                    MinVal = fx;
                if (fx > MaxVal)
                    MaxVal = fx;
                ProbFromFx.Add(snet, fx.ToString());
                #endregion
            }

            /* find normalize fx */
            normProbFromFx = FindNormalizedLogliklihood_2(ProbFromFx, MaxVal);
            #endregion

            #region make the final list
            /* Sort the dictionary based on the key (in descending order) */
            List<KeyValuePair<string, double>> sortedList = tempDic.ToList();
            sortedList.Sort((pair1, pair2) => pair1.Value.CompareTo(pair2.Value));

            foreach (KeyValuePair<string, double> kv in sortedList)
            {
                string key = kv.Key;
                _TD_ED.Add(key, ProbFromFrequency[key] + "::" + normProbFromFx[key]);
            }
            #endregion

            return _TD_ED;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ProbFromFx"></param>
        /// <param name="MaxVal"></param>
        /// <returns></returns>
        private NameValueCollection FindNormalizedLogliklihood_2(NameValueCollection ProbFromFx, double MaxVal)
        {
            double _log_n = Math.Log(ProbFromFx.AllKeys.Count());
            double _log_epsilon_by_n = _log_epsilon - _log_n;
            double _sum = 0.0;

            /* find Logfx_val - max_Logfx_val */
            NameValueCollection _probFromfx = new NameValueCollection();
            foreach (string key in ProbFromFx.AllKeys)
            {
                double _val = double.Parse(ProbFromFx[key]) - MaxVal;

                if (_val > _log_epsilon_by_n)
                {
                    _val = Math.Exp(_val);
                    _sum += _val;
                    _probFromfx.Add(key, _val.ToString());
                }
                else
                    _probFromfx.Add(key, "0");  // see file from "Cross-validated" forum that was used for this purpose
            }
            foreach (string key in _probFromfx.AllKeys)
            {
                double _val = double.Parse(_probFromfx[key]) / _sum;
                _probFromfx.Remove(key);
                _probFromfx.Add(key, _val.ToString());
            }

            return _probFromfx;
        }

        #region Other Functions
        private double Calcuate_Fy(int[][] _X, Dictionary<string, string> _X_neighbourhoodList, int _Y_neighbourID, int[][] _Y)
        {
            /* find list of nodes for which nParents has changed */
            Dictionary<int, List<int>> L_changedNodes = FindNodeswithChangedParents(_X, _X_neighbourhoodList, _Y_neighbourID, _Y);

            /* find updated conditional probabilities for those changed nodes */
            L_conditionalProb.Clear();
            L_conditionalProb = FindUpdateConditionalProbabilities(L_changedNodes);


            /* calculate the final Fz */
            double log_F1 = 0.0;
            for (int i = 0; i < nNodes; i++)
            {
                if (L_conditionalProb[i.ToString()] != null)
                    log_F1 += double.Parse(L_conditionalProb[i.ToString()]);
                else
                    log_F1 += ProbabilityofAllNodes[i];
            }
            _global_LogFx = log_F1;


            //double F1 = Math.Exp(log_F1);
            //return F1;
            return log_F1;
        }

        private Dictionary<int, List<int>> FindNodeswithChangedParents(int[][] _X, Dictionary<string, string> _X_neighbourhoodList, int _Y_neighbourID, int[][] _Y)
        {
            Dictionary<int, List<int>> retList = new Dictionary<int, List<int>>();

            List<int> L_z_i = new List<int>();
            List<int> L_z_j = new List<int>();

            /* find i,j,p,q */
            int _i = 0;
            int _j = 0;
            string Y_op = "";

            /* find the Pairs that was applied sequentially to get Z network */
            string Y_edge = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];
            if (!Y_edge.Equals("Nil"))
            {
                Y_op = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1];
                _i = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                _j = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);
            }
            else
                return retList;

            /* save all the parents of nodes i, j, p, q that are associated with them, correspondingly */
            for (int indx = 0; indx < nNodes; indx++)
            {
                if (!Y_edge.Equals("Nil"))
                {
                    if (_Y[indx][_i] == 1)
                        L_z_i.Add(indx);
                    if (_Y[indx][_j] == 1)
                        L_z_j.Add(indx);
                }
            }

            if (!retList.ContainsKey(_i))
                retList.Add(_i, L_z_i);
            if (!retList.ContainsKey(_j))
                retList.Add(_j, L_z_j);


            return retList;
        }

        private bool isTwoNetworksSame(int[][] X, int[][] Z)
        {
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = 0; j < nNodes; j++)
                {
                    if (i != j)
                    {
                        if (X[i][j] != Z[i][j])
                            return false;
                    }
                }
            }
            return true;
        }
        private void pringANetwork(int nZ, int[][] Y)
        {
            TextWriter tw = new StreamWriter("Y_network_nNeithbour_ZCount_" + nZ.ToString() + ".txt");
            for (int i = 0; i < nNodes; i++)
            {
                string line = "";
                for (int j = 0; j < nNodes; j++)
                    line += Y[i][j] + ",";
                line = line.Substring(0, line.Length - 1);
                tw.WriteLine(line);
            }
            tw.Close();
        }


        private double Calcuate_Fz(int[][] _X, Dictionary<string, string> _X_neighbourhoodList, int _Y_neighbourID, int[][] _Z, Dictionary<string, string> _Y_neighbourhoodList, int _Z_neighbourID)
        {

            if (isTwoNetworksSame(_X, _Z))
            {
                int comehere = 0;
            }


            /* find list of nodes for which nParents has changed */
            Dictionary<int, List<int>> L_changedNodes = FindNodeswithChangedParents(_X, _X_neighbourhoodList, _Y_neighbourID, _Z, _Y_neighbourhoodList, _Z_neighbourID);

            /* find updated conditional probabilities for those changed nodes */
            L_conditionalProb.Clear();
            L_conditionalProb = FindUpdateConditionalProbabilities(L_changedNodes);


            /* calculate the final Fz */
            double log_F1 = 0.0;
            for (int i = 0; i < nNodes; i++)
            {
                if (L_conditionalProb[i.ToString()] != null)
                    log_F1 += double.Parse(L_conditionalProb[i.ToString()]);
                else
                    log_F1 += ProbabilityofAllNodes[i];
            }
            _global_LogFx = log_F1;


            //double F1 = Math.Exp(log_F1);
            //return F1;
            return log_F1;
        }

        private NameValueCollection FindUpdateConditionalProbabilities(Dictionary<int, List<int>> L_changedNodes)
        {
            NameValueCollection retList = new NameValueCollection();

            Dictionary<int, NameValueCollection> ConditionalProbTable_List = new Dictionary<int, NameValueCollection>();
            Dictionary<int, NameValueCollection> PriorBeliefTable_List = new Dictionary<int, NameValueCollection>();

            #region Calculate updated Conditional Proabilities of all changedNodes
            foreach (int indx in L_changedNodes.Keys)
            {
                #region build Conditional Proability tables
                /* 'indx' indicates node_i */

                /* for each node create a namevaluecollection where the keyName:= "Combination of parents (, seperated) ;; value is the frequency "*/
                NameValueCollection nodeConditionalProbTable = new NameValueCollection();
                List<int> nodeParents = L_changedNodes[indx];

                if (nodeParents.Count() > 0)
                {
                    #region For the nodes which has got at least one parent

                    #region Generate all possible parental conditions
                    List<int[]> input = new List<int[]>();

                    for (int j = 0; j < nodeParents.Count(); j++)
                    {
                        int key = nodeParents[j];
                        int nCat = Int32.Parse(DataCategories[key]);
                        int[] tempArr = new int[nCat];
                        for (int k = 0; k < nCat; k++)
                            tempArr[k] = k;

                        input.Add(tempArr);
                    }

                    int[] size = new int[input.Count()];

                    /* enumerate all possible parent categories */
                    combine(ref nodeConditionalProbTable, input, size, 0);
                    #endregion

                    #region Explore the whole dataset for finding the frequency
                    foreach (string condition in nodeConditionalProbTable.AllKeys)
                    {
                        for (int j = 0; j < nObservation; j++)
                        {
                            string tempS = "";
                            /* this is just the initialization of the node observation value */
                            int nodeObservationValue = 0;

                            for (int k = 0; k < nNodes; k++)
                            {
                                /* if the index 'k' is in the parent list*/
                                if (nodeParents.Contains(k))
                                    tempS += DataSet[k][j] + ":";

                                /* saving the node observation value */
                                if (k == indx)
                                    nodeObservationValue = DataSet[k][j];
                            }
                            tempS = tempS.Substring(0, tempS.Length - 1);

                            /* if the parental condition list matches with the observation list, then the observation into the CT */
                            if (condition.Equals(tempS))
                                nodeConditionalProbTable.Add(condition, nodeObservationValue.ToString());
                        }
                    }
                    #endregion

                    ConditionalProbTable_List.Add(indx, nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    NameValueCollection nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[indx]));
                    PriorBeliefTable_List.Add(indx, nodePriorBeliefTable);
                    #endregion

                    #endregion
                }
                else if (nodeParents.Count() == 0)
                {
                    #region For the nodes which have no parents

                    #region Explore the whole dataset for finding the frequency

                    for (int j = 0; j < nObservation; j++)
                    {
                        /* this is just the initialization of the node observation value */
                        int nodeObservationValue = 0;

                        for (int k = 0; k < nNodes; k++)
                        {
                            /* saving the node observation value */
                            if (k == indx)
                                nodeObservationValue = DataSet[k][j];
                        }

                        /* since the node_i has no parent, save it's id: 'i' as the key */
                        nodeConditionalProbTable.Add(indx.ToString(), nodeObservationValue.ToString());
                    }
                    #endregion

                    ConditionalProbTable_List.Add(indx, nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    NameValueCollection nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[indx]));
                    PriorBeliefTable_List.Add(indx, nodePriorBeliefTable);
                    #endregion
                    #endregion
                }
                /* for each category of node_i: make all*/
                #endregion
            }

            foreach (int indx in L_changedNodes.Keys)
            {
                #region for each node (Random variable), calculate the value of density function
                NameValueCollection nodeConditionalProbTable = ConditionalProbTable_List[indx];
                NameValueCollection nodePriorBeliefTable = PriorBeliefTable_List[indx];

                double sum_conditionVal = 0.0;
                foreach (string condition in nodeConditionalProbTable.AllKeys)
                {
                    string[] arr_alpha_ij = nodePriorBeliefTable[condition].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

                    #region Calculate alpha_ij
                    double alpha_ij = 0.0;
                    foreach (string s in arr_alpha_ij)
                        alpha_ij += double.Parse(s);
                    #endregion

                    int _nCat = Int32.Parse(DataCategories[indx]);
                    double[] arr_N_ij = FindConditionalProbabilityValues(nodeConditionalProbTable[condition], _nCat);

                    #region Calculate N_ij
                    double N_ij = 0.0;
                    foreach (int val in arr_N_ij)
                        N_ij += (double)val;
                    #endregion

                    double sum_categoricalVal = 0.0;
                    for (int k = 0; k < _nCat; k++)
                    {
                        double _categoricalVal = LogGamma(arr_N_ij[k] + double.Parse(arr_alpha_ij[k])) - LogGamma(double.Parse(arr_alpha_ij[k]));
                        sum_categoricalVal += _categoricalVal;
                    }

                    double _conditionVal = sum_categoricalVal + (LogGamma(alpha_ij) - LogGamma(N_ij + alpha_ij));
                    sum_conditionVal += _conditionVal;


                }
                #endregion

                /* save this */
                retList.Add(indx.ToString(), sum_conditionVal.ToString());
            }
            #endregion

            return retList;
        }

        private Dictionary<int, List<int>> FindNodeswithChangedParents(float[][] _X, Dictionary<string, string> _X_neighbourhoodList, int _Y_neighbourID, int[][] _Z, Dictionary<string, string> _Y_neighbourhoodList, int _Z_neighbourID)
        {
            Dictionary<int, List<int>> retList = new Dictionary<int, List<int>>();


            //List<int> L_z_i = new List<int>();
            //List<int> L_z_j = new List<int>();
            //List<int> L_z_p = new List<int>();
            //List<int> L_z_q = new List<int>();

            ///* find i,j,p,q */
            //int _i = 0;
            //int _j = 0;
            //int _p = 0;
            //int _q = 0;
            //string Y_op = "";
            //string Z_op = "";

            ///* find the Pairs that was applied sequentially to get Z network */
            //string Y_edge = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];
            //if (!Y_edge.Equals("Nil"))
            //{
            //    Y_op = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1];
            //    _i = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            //    _j = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);
            //}
            //string Z_edge = _Y_neighbourhoodList.Keys.ToArray()[_Z_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];
            //if (!Z_edge.Equals("Nil"))
            //{
            //    Z_op = _Y_neighbourhoodList.Keys.ToArray()[_Z_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1];
            //    _p = Int32.Parse(Z_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            //    _q = Int32.Parse(Z_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);
            //}

            //if ((Z_op.Equals("Add") && Y_op.Equals("Del") && (_p == _i) && (_q == _j)) ||
            //        (Z_op.Equals("Del") && Y_op.Equals("Add") && (_p == _i) && (_q == _j)) ||
            //        (Z_op.Equals("Rev") && Y_op.Equals("Rev") && (_p == _j) && (_q == _i)))
            //    return retList;

            ///* save all the parents of nodes i, j, p, q that are associated with them, correspondingly */
            //for (int indx = 0; indx < nNodes; indx++)
            //{
            //    if (!Y_edge.Equals("Nil"))
            //    {
            //        if (_Z[indx][_i] == 1)
            //            L_z_i.Add(indx);
            //        if (_Z[indx][_j] == 1)
            //            L_z_j.Add(indx);
            //    }
            //    if (!Z_edge.Equals("Nil"))
            //    {
            //        if (_Z[indx][_p] == 1)
            //            L_z_p.Add(indx);
            //        if (_Z[indx][_q] == 1)
            //            L_z_q.Add(indx);
            //    }
            //}

            ///* */
            //if (!Y_edge.Equals("Nil"))
            //{
            //    if (Y_op.Equals("Rev"))
            //        retList.Add(_i, L_z_i);
            //    else
            //        retList.Add(_j, L_z_j);
            //}

            //if (!Z_edge.Equals("Nil"))
            //{
            //    if (Z_op.Equals("Rev"))
            //    {
            //        if (!retList.ContainsKey(_p))
            //            retList.Add(_p, L_z_p);
            //    }
            //    else
            //    {
            //        if (!retList.ContainsKey(_q))
            //            retList.Add(_q, L_z_q);
            //    }
            //}

            return retList;
        }

        private Dictionary<int, List<int>> FindNodeswithChangedParents(int[][] _X, Dictionary<string, string> _X_neighbourhoodList, int _Y_neighbourID, int[][] _Z, Dictionary<string, string> _Y_neighbourhoodList, int _Z_neighbourID)
        {
            Dictionary<int, List<int>> retList = new Dictionary<int, List<int>>();


            List<int> L_z_i = new List<int>();
            List<int> L_z_j = new List<int>();
            List<int> L_z_p = new List<int>();
            List<int> L_z_q = new List<int>();

            /* find i,j,p,q */
            int _i = 0;
            int _j = 0;
            int _p = 0;
            int _q = 0;
            string Y_op = "";
            string Z_op = "";

            /* find the Pairs that was applied sequentially to get Z network */
            string Y_edge = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];
            if (!Y_edge.Equals("Nil"))
            {
                Y_op = _X_neighbourhoodList.Keys.ToArray()[_Y_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1];
                _i = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                _j = Int32.Parse(Y_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);
            }

            string Z_edge = _Y_neighbourhoodList.Keys.ToArray()[_Z_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];
            if (!Z_edge.Equals("Nil"))
            {
                Z_op = _Y_neighbourhoodList.Keys.ToArray()[_Z_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1];
                _p = Int32.Parse(Z_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                _q = Int32.Parse(Z_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);
            }

            if ((Y_edge.Equals("Nil") && Z_edge.Equals("Nil")) ||
                (Z_op.Equals("Add") && Y_op.Equals("Del") && (_p == _i) && (_q == _j)) ||
                    (Z_op.Equals("Del") && Y_op.Equals("Add") && (_p == _i) && (_q == _j)) ||
                    (Z_op.Equals("Rev") && Y_op.Equals("Rev") && (_p == _j) && (_q == _i)))
                return retList;

            /* save all the parents of nodes i, j, p, q that are associated with them, correspondingly */
            for (int indx = 0; indx < nNodes; indx++)
            {
                if (!Y_edge.Equals("Nil"))
                {
                    if (_Z[indx][_i] == 1)
                        L_z_i.Add(indx);
                    if (_Z[indx][_j] == 1)
                        L_z_j.Add(indx);
                }
                if (!Z_edge.Equals("Nil"))
                {
                    if (_Z[indx][_p] == 1)
                        L_z_p.Add(indx);
                    if (_Z[indx][_q] == 1)
                        L_z_q.Add(indx);
                }
            }

            /* */
            /*if (!Y_edge.Equals("Nil"))
            {
                if (Y_op.Equals("Rev"))
                    retList.Add(_i, L_z_i);
                else
                    retList.Add(_j, L_z_j);
            }

            if (!Z_edge.Equals("Nil"))
            {
                if (Z_op.Equals("Rev"))
                {
                    if (!retList.ContainsKey(_p))
                        retList.Add(_p, L_z_p);
                }
                else
                {
                    if (!retList.ContainsKey(_q))
                        retList.Add(_q, L_z_q);
                }
            }*/

            if (!retList.ContainsKey(_i))
                retList.Add(_i, L_z_i);
            if (!retList.ContainsKey(_j))
                retList.Add(_j, L_z_j);
            if (!retList.ContainsKey(_p))
                retList.Add(_p, L_z_p);
            if (!retList.ContainsKey(_q))
                retList.Add(_q, L_z_q);


            return retList;
        }

        private double Calcuate_Fx(int[][] X)
        {
            double log_F1 = 0.0;
            for (int i = 0; i < nNodes; i++)
                log_F1 += ProbabilityofAllNodes[i];

            _global_LogFx = log_F1;

            //double F1 = Math.Exp(log_F1);
            //return F1;
            return log_F1;
        }


        private void BuildandCalculateProbabilitiesforEachNodes(int[][] T_x)
        {
            #region Make Conditional Distributions for each nodes
            List<NameValueCollection> ConditionalProbTable_List = new List<NameValueCollection>();
            List<NameValueCollection> PriorBeliefTable_List = new List<NameValueCollection>();

            for (int i = 0; i < nNodes; i++)
            {
                #region build Conditional Proability tables
                /* 'i' indicates node_i */

                /* for each node create a namevaluecollection where the keyName:= "Combination of parents (, seperated) ;; value is the frequency "*/
                NameValueCollection nodeConditionalProbTable = new NameValueCollection();
                NameValueCollection nodePriorBeliefTable = new NameValueCollection();

                NameValueCollection nodeParents = Find_nParents_of_node(T_x, i);

                #region Make node Conditional Probability table and Prior table

                if (nodeParents.AllKeys.Count() > 0)
                {
                    #region For the nodes which has got at least one parent

                    #region Generate all possible parental conditions
                    List<int[]> input = new List<int[]>();

                    for (int j = 0; j < nodeParents.AllKeys.Count(); j++)
                    {
                        string key = nodeParents.AllKeys[j];
                        int nCat = Int32.Parse(DataCategories[Int32.Parse(key)]);
                        int[] tempArr = new int[nCat];
                        for (int k = 0; k < nCat; k++)
                            tempArr[k] = k;

                        input.Add(tempArr);
                    }

                    int[] size = new int[input.Count()];

                    /* enumerate all possible parent categories */
                    combine(ref nodeConditionalProbTable, input, size, 0);
                    #endregion

                    #region Explore the whole dataset for finding the frequency
                    foreach (string condition in nodeConditionalProbTable.AllKeys)
                    {
                        for (int j = 0; j < nObservation; j++)
                        {
                            string tempS = "";
                            /* this is just the initialization of the node observation value */
                            int nodeObservationValue = 0;

                            for (int k = 0; k < nNodes; k++)
                            {
                                /* if the index 'k' is in the parent list*/
                                if (nodeParents[k.ToString()] != null)
                                    tempS += DataSet[k][j] + ":";

                                /* saving the node observation value */
                                if (k == i)
                                    nodeObservationValue = DataSet[k][j];
                            }
                            tempS = tempS.Substring(0, tempS.Length - 1);

                            /* if the parental condition list matches with the observation list, then the observation into the CT */
                            if (condition.Equals(tempS))
                                nodeConditionalProbTable.Add(condition, nodeObservationValue.ToString());
                        }
                    }
                    #endregion

                    //ConditionalProbTable_List.Add(nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[i]));
                    //PriorBeliefTable_List.Add(nodePriorBeliefTable);
                    #endregion

                    #endregion
                }
                else if (nodeParents.AllKeys.Count() == 0)
                {
                    #region For the nodes which have no parents

                    #region Explore the whole dataset for finding the frequency

                    for (int j = 0; j < nObservation; j++)
                    {
                        /* this is just the initialization of the node observation value */
                        int nodeObservationValue = 0;

                        for (int k = 0; k < nNodes; k++)
                        {
                            /* saving the node observation value */
                            if (k == i)
                                nodeObservationValue = DataSet[k][j];
                        }

                        /* since the node_i has no parent, save it's id: 'i' as the key */
                        nodeConditionalProbTable.Add(i.ToString(), nodeObservationValue.ToString());
                    }
                    #endregion

                    //ConditionalProbTable_List.Add(nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[i]));
                    //PriorBeliefTable_List.Add(nodePriorBeliefTable);
                    #endregion
                    #endregion
                }
                #endregion

                #region Make Conditional Probability Value
                double sum_conditionVal = 0.0;
                foreach (string condition in nodeConditionalProbTable.AllKeys)
                {
                    string[] arr_alpha_ij = nodePriorBeliefTable[condition].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    #region Calculate alpha_ij
                    double alpha_ij = 0.0;
                    foreach (string s in arr_alpha_ij)
                        alpha_ij += double.Parse(s);
                    #endregion

                    int _nCat = Int32.Parse(DataCategories[i]);
                    double[] arr_N_ij = FindConditionalProbabilityValues(nodeConditionalProbTable[condition], _nCat);
                    #region Calculate N_ij
                    double N_ij = 0.0;
                    foreach (int val in arr_N_ij)
                        N_ij += (double)val;
                    #endregion

                    double sum_categoricalVal = 0.0;
                    for (int k = 0; k < _nCat; k++)
                    {
                        double _categoricalVal = LogGamma(arr_N_ij[k] + double.Parse(arr_alpha_ij[k])) - LogGamma(double.Parse(arr_alpha_ij[k]));
                        sum_categoricalVal += _categoricalVal;
                    }

                    double _conditionVal = sum_categoricalVal + (LogGamma(alpha_ij) - LogGamma(N_ij + alpha_ij));
                    sum_conditionVal += _conditionVal;
                }
                #endregion

                #region save this */
                ProbabilityofAllNodes[i] = sum_conditionVal;
                #endregion

                #endregion
            }
            #endregion



        }


        private void RevertZbacktoinitialY(ref int[][] tempNet, Dictionary<string, string> _neighbourhoodList, int _neighbourID)
        {
            #region Apply the operation on the Network
            /* find the edges to be changed (operated) */
            string _edge = _neighbourhoodList.Keys.ToArray()[_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

            /* if edge = Nil, it means '_neighbourID'='0', means the network itself, so no operation, just return with the current 'TempNet' */
            if (_edge.Equals("Nil"))
                return;

            int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

            /* Find the operation to be done for generating the Neighbour */
            string op = _neighbourhoodList[_neighbourhoodList.Keys.ToArray()[_neighbourID]];

            /* do the operation */
            if (op.Equals("Del"))
                tempNet[_i][_j] = 1;
            else if (op.Equals("Add"))
                tempNet[_i][_j] = 0;
            else if (op.Equals("Rev"))
            {
                tempNet[_j][_i] = 0;
                tempNet[_i][_j] = 1;
            }
            #endregion
        }

        private void UseFixedInitialNetwork(ref int[][] net)
        {
            TextReader tr_net = new StreamReader("InitialNetwork_" + nNodes.ToString() + "_nodes.txt");
            //TextReader tr_net = new StreamReader("Y_network_nNeithbour_ZCount_1003623481.txt");
            int i_indx = 0;
            for (int i = 0; i < nNodes; i++)
            {
                net[i] = new int[nNodes];
                string[] stall = tr_net.ReadLine().Trim().Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                for (int j = 0; j < nNodes; j++)
                {
                    net[i][j] = Int32.Parse(stall[j]);
                }
            }
        }

        private int FindNetworkTotalDegree(int[][] X)
        {
            int count = 0;
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = 0; j < nNodes; j++)
                {
                    if (X[i][j] == 1)
                        count++;
                }
            }
            return count;
        }

        private string FindSquaredDiff()
        {
            NameValueCollection ProbFromFrequency = new NameValueCollection();
            NameValueCollection ProbFromFx = new NameValueCollection();
            NameValueCollection normProbFromFx = new NameValueCollection();
            NameValueCollection SampledNetwork2 = new NameValueCollection();

            //double _log_epsilon = -34.538776395; // this value comes from log_e(10e-16), which is 16 digit IEEE double regarding the level of precision we want
            double MaxVal = -9999999.0;

            #region Make Sampled Network Table_2
            foreach (string key in SampledNetwork.AllKeys)
            {
                int cnt = SampledNetwork[key].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).Count();
                SampledNetwork2.Add(key, cnt.ToString());
            }
            #endregion

            #region Calculate Probabilities for each network
            foreach (string snet in SampledNetwork.AllKeys)
            {
                #region Make the newtork from the string
                //int[][] net = new int[nNodes][];
                //string[] _snet = snet.Replace("\"", "").Split("|".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                //for (int i = 0; i < nNodes; i++)
                //{
                //    net[i] = new int[nNodes];
                //    char[] _cnet = _snet[i].ToCharArray();
                //    for (int j = 0; j < nNodes; j++)
                //    {
                //        string _ts = _cnet[j].ToString();
                //        net[i][j] = Int32.Parse(_ts);
                //    }
                //}
                #endregion

                #region Calculate Probablilities from Freq
                double cProb = (double.Parse(SampledNetwork2[snet]) / (double)nMaxIteration);
                ProbFromFrequency.Add(snet, cProb.ToString());
                #endregion

                #region Calculate Probabilities from Fx
                //double fx = CalcuateTargetFunction(net);
                double fx = double.Parse(SampledNetwork_Fx[snet]);
                //fx = Math.Exp(fx);
                //fx = Math.Abs(fx);
                if (fx > MaxVal)
                    MaxVal = fx;
                ProbFromFx.Add(snet, fx.ToString());
                #endregion
            }

            /* find normalize fx */
            normProbFromFx = FindNormalizedLogliklihood(ProbFromFx, MaxVal);

            /* find the squared-differences */
            double _sqrdDiff = FindSquaredDifference(ProbFromFrequency, normProbFromFx);

            return _sqrdDiff.ToString();
            #endregion

        }
        private NameValueCollection FindNormalizedLogliklihood(NameValueCollection ProbFromFx, double MaxVal)
        {
            double _log_n = Math.Log(ProbFromFx.AllKeys.Count());
            double _log_epsilon_by_n = _log_epsilon - _log_n;
            double _sum = 0.0;

            /* find Logfx_val - max_Logfx_val */
            NameValueCollection _probFromfx = new NameValueCollection();
            foreach (string key in ProbFromFx.AllKeys)
            {
                double _val = double.Parse(ProbFromFx[key]) - MaxVal;

                if (_val > _log_epsilon_by_n)
                {
                    _val = Math.Exp(_val);
                    _sum += _val;
                    _probFromfx.Add(key, _val.ToString());
                }
            }
            foreach (string key in _probFromfx.AllKeys)
            {
                double _val = double.Parse(_probFromfx[key]) / _sum;
                _probFromfx.Remove(key);
                _probFromfx.Add(key, _val.ToString());
            }

            return _probFromfx;
        }
        private double FindSquaredDifference(NameValueCollection ProbFromFrequency, NameValueCollection normProbFromFx)
        {
            double _val = 0.0;
            foreach (string key in normProbFromFx.AllKeys)
            {
                double _probFreq = double.Parse(ProbFromFrequency[key]);
                double _probFx = double.Parse(normProbFromFx[key]);

                _val += Math.Pow((_probFreq - _probFx), 2.0);
            }
            return _val;
        }

        private void Restore_All_Lists(Dictionary<string, string> temp_D_list, Dictionary<string, string> temp_A_list, Dictionary<string, string> temp_R_list)
        {
            D_list.Clear();
            foreach (var key in temp_D_list.Keys)
                D_list.Add(key, temp_D_list[key]);

            A_list.Clear();
            foreach (var key in temp_A_list.Keys)
                A_list.Add(key, temp_A_list[key]);

            R_list.Clear();
            foreach (var key in temp_R_list.Keys)
                R_list.Add(key, temp_R_list[key]);

        }

        private int FindNetworkID(NameValueCollection SampledNetwork, string sX)
        {
            for (int i = 0; i < SampledNetwork.AllKeys.Count(); i++)
            {
                if (SampledNetwork.AllKeys[i].Equals(sX))
                    return i;
            }
            return -1;
        }

        private string MakeNetworkString(int[][] X)
        {
            string s = "\"";
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = 0; j < nNodes; j++)
                {
                    //if (i != j)
                    s += X[i][j].ToString();
                }
                s += "|";
            }
            s += "\"";
            return s;
        }

        private void Save_All_Lists(ref Dictionary<string, string> temp_D_list, ref Dictionary<string, string> temp_A_list, ref Dictionary<string, string> temp_R_list)
        {
            foreach (var key in D_list.Keys)
                temp_D_list.Add(key, D_list[key]);
            foreach (var key in A_list.Keys)
                temp_A_list.Add(key, A_list[key]);
            foreach (var key in R_list.Keys)
                temp_R_list.Add(key, R_list[key]);
        }

        private void GenerateNeighbour(ref int[][] tempNet, Dictionary<string, string> _neighbourhoodList, int _neighbourID)
        {
            #region Apply the operation on the Network
            /* find the edges to be changed (operated) */
            string _edge = _neighbourhoodList.Keys.ToArray()[_neighbourID].Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

            /* if edge = Nil, it means '_neighbourID'='0', means the network itself, so no operation, just return with the current 'TempNet' */
            if (_edge.Equals("Nil"))
                return;

            int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

            /* Find the operation to be done for generating the Neighbour */
            string op = _neighbourhoodList[_neighbourhoodList.Keys.ToArray()[_neighbourID]];

            /* do the operation */
            if (op.Equals("Del"))
            {
                if (tempNet[_i][_j] == 0)
                {
                    int comehere = 0;
                }
                tempNet[_i][_j] = 0;

            }
            else if (op.Equals("Add"))
            {
                if (tempNet[_i][_j] == 1)
                {
                    int comehere = 0;
                }
                tempNet[_i][_j] = 1;

            }
            else if (op.Equals("Rev"))
            {
                if (tempNet[_j][_i] == 1)
                {
                    int comehere = 0;

                }
                tempNet[_j][_i] = 1;
                tempNet[_i][_j] = 0;

            }
            #endregion

            #region Find Bridlist in the New network, so that it can be useful for next steps
            Dictionary<string, string> tempNet_BridgeList = FindBridgesinNetwork(tempNet);
            #endregion

            #region Update all the lists: (A_list, D_list, R_list) after doing the operation

            Update_D_list(ref tempNet, _edge, op, tempNet_BridgeList);
            Update_A_list(ref tempNet, _edge, op);
            Update_R_list(ref tempNet, _edge, op);

            #endregion
        }

        private void Update_R_list(ref int[][] tempNet, string _edge, string op)
        {
            int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

            if (op.Equals("Del"))
            {
                /* step (1): Remove (i,j) from R_list */
                R_list.Remove(_i.ToString() + ":" + _j.ToString());

                /* step (2): check NON-Reversibles only */
                //Check_Selected_Edges_in_R_list_Only(ref tempNet, op);
                Dictionary<string, string> nonreversible_list = (from entry in R_list where entry.Value.StartsWith("N,") select entry).ToDictionary(x => x.Key, x => x.Value);
                foreach (string elem in nonreversible_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_R_list(tempNet, __i, __j);
                }
            }
            else if (op.Equals("Add"))
            {
                /* step (1): Add (i,j) into R_list as NON-Reversible initially, so that it could be checked with others in next steps */
                R_list.Add(_i.ToString() + ":" + _j.ToString(), "Y");

                /* step (2): Check only Reversibles */
                //Check_Selected_Edges_in_R_list_Only(ref tempNet, op);
                Dictionary<string, string> reversible_list = (from entry in R_list where entry.Value == "Y" select entry).ToDictionary(x => x.Key, x => x.Value);
                foreach (string elem in reversible_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_R_list(tempNet, __i, __j);
                }
            }
            else if (op.Equals("Rev"))
            {
                /* step (1): Remove (i,j) from R_list */
                R_list.Remove(_i.ToString() + ":" + _j.ToString());

                /* step (2): Add (j,i) from R_list for checking */
                R_list.Add(_j.ToString() + ":" + _i.ToString(), "N,check");

                /* step (2): Check the whole list */
                //Check_Selected_Edges_in_R_list_Only(ref tempNet, op);
                foreach (string elem in R_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_R_list(tempNet, __i, __j);
                }
            }
        }

        private void Check_Selected_Edges_in_R_list_Only(ref int[][] tempNet, string op)
        {
            foreach (string _edge in R_list.Keys.ToList())
            {
                int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                if (op.Equals("Del"))
                {
                    /* only the NON-Reversible edges in R_list will be checked here */
                    if (R_list[_edge].StartsWith("N,"))
                        Change_R_list(tempNet, _i, _j);
                }
                else if (op.Equals("Add"))
                {
                    /* only the Reversible edges in R_list will be checked here */
                    if (R_list[_edge].Equals("Y"))
                        Change_R_list(tempNet, _i, _j);
                }
                else if (op.Equals("Rev"))
                {
                    /* check the whole list */
                    Change_R_list(tempNet, _i, _j);
                }
            }
        }

        private void Change_R_list(int[][] tempNet, int _i, int _j)
        {
            R_list.Remove(_i.ToString() + ":" + _j.ToString()); /* just to prevent multiple adding; as this pais is going to be added in next steps anyway */

            #region Reverse
            /* firt remove the edge and the see what happens */
            tempNet[_i][_j] = 0;

            /* check adding the reverse edge creates a cycle or not */
            bool _isCycle = _isCycleExist(tempNet, _j, _i);
            if (!_isCycle)
            {
                bool _is_nParent_Restricted = is_nParent_Restricted(tempNet, _j, _i);
                bool _is_nChild_Restricted = is_nChild_Restricted(tempNet, _j, _i);
                if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                    R_list[_i.ToString() + ":" + _j.ToString()] = "Y";
                else
                    R_list[_i.ToString() + ":" + _j.ToString()] = "N,degree_th_Exceed";
            }
            else
                R_list[_i.ToString() + ":" + _j.ToString()] = "N,cycle_detected";

            /* firt remove the edge and the see what happens */
            tempNet[_i][_j] = 1;

            #endregion
        }


        private void Update_D_list(ref int[][] tempNet, string _edge, string op, Dictionary<string, string> tempNet_BridgeList)
        {
            int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

            if (op.Equals("Del"))
            {
                /* step (1): Remove (i,j) from D_list */
                D_list.Remove(_i.ToString() + ":" + _j.ToString());

                /* step (2): Check only the deletable edges for update*/
                //Check_Selected_Edges_in_D_list_Only(tempNet_BridgeList, op);
                Dictionary<string, string> deletable_list = (from entry in D_list where entry.Value == "Y" select entry).ToDictionary(x => x.Key, x => x.Value);

                foreach (string elem in deletable_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_D_list(tempNet_BridgeList, __i, __j);
                }

            }
            else if (op.Equals("Add"))
            {
                /* step (1): Add (i,j) in D_list as deletable */
                D_list.Add(_i.ToString() + ":" + _j.ToString(), "Y");

                /* step (2): check only NON-deletable */
                //Check_Selected_Edges_in_D_list_Only(tempNet_BridgeList, op);
                Dictionary<string, string> nondeletable_list = (from entry in D_list where entry.Value.StartsWith("N,") select entry).ToDictionary(x => x.Key, x => x.Value);
                foreach (string elem in nondeletable_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_D_list(tempNet_BridgeList, __i, __j);
                }
            }
            else if (op.Equals("Rev"))
            {
                /* step (1): Delete (i,j) from D_list */
                D_list.Remove(_i.ToString() + ":" + _j.ToString());

                /* step (2): Add (j,i) in D_list */
                if (!tempNet_BridgeList.ContainsKey(_j.ToString() + "," + _i.ToString()) && !tempNet_BridgeList.ContainsKey(_i.ToString() + "," + _j.ToString()))
                    D_list.Add(_j.ToString() + ":" + _i.ToString(), "Y");
                else
                    D_list.Add(_j.ToString() + ":" + _i.ToString(), "N,bridge_detected");
            }


        }

        private void Check_Selected_Edges_in_D_list_Only(Dictionary<string, string> tempNet_BridgeList, string op)
        {
            foreach (string _e in D_list.Keys.ToList())
            {
                string _edge = _e.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                if (op.Equals("Del"))
                {
                    if (D_list[_edge].Equals("Y"))
                        Change_D_list(tempNet_BridgeList, _i, _j);
                }
                else if (op.Equals("Add"))
                {
                    if (D_list[_edge].StartsWith("N,"))
                        Change_D_list(tempNet_BridgeList, _i, _j);
                }
            }
        }

        private void Change_D_list(Dictionary<string, string> tempNet_BridgeList, int _i, int _j)
        {
            D_list.Remove(_i.ToString() + ":" + _j.ToString()); /* preventing multiple addition for the same pair */
            if (!tempNet_BridgeList.ContainsKey(_i.ToString() + "," + _j.ToString()) && !tempNet_BridgeList.ContainsKey(_j.ToString() + "," + _i.ToString()))
                D_list[_i.ToString() + ":" + _j.ToString()] = "Y";
            else
                D_list[_i.ToString() + ":" + _j.ToString()] = "N,bridge_detected";
        }

        private Dictionary<string, string> FindBridgesinNetwork(int[][] tempNet)
        {
            int cnt = 0;          // counter
            int[] low;        // pre[v] = order in which dfs examines v
            int[] pre;        // low[v] = lowest preorder of any vertex connected to v
            Dictionary<string, string> _BridgeList = new Dictionary<string, string>();
            low = new int[nNodes];
            pre = new int[nNodes];

            for (int u = 0; u < nNodes; u++)
            {
                low[u] = -1;
                pre[u] = -1;
            }

            for (int v = 0; v < nNodes; v++)
                if (pre[v] == -1)
                    dfs(tempNet, v, v, cnt, low, pre, ref _BridgeList);

            return _BridgeList;

        }

        private void Update_A_list(ref int[][] tempNet, string _edge, string op)
        {
            int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
            int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

            if (op.Equals("Del"))
            {
                /* step (1): Add (i,j) in 'A_list' as addable */
                A_list.Add(_i.ToString() + ":" + _j.ToString(), "Y");

                /* step (2): Add (j,i) in 'A_list' as NON-addable [so that we can check it along with other NON-addable edges in the next step] */
                A_list.Add(_j.ToString() + ":" + _i.ToString(), "N,check");

                /* step (3): Check all NON-addable [they may become addable]*/
                //Check_Selected_Edges_in_A_list_Only(ref tempNet, op);

                /* select non-addable entries from A_list */
                Dictionary<string, string> Nonaddable_list = (from entry in A_list where entry.Value.StartsWith("N,") select entry).ToDictionary(x => x.Key, x => x.Value);
                foreach (string elem in Nonaddable_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_A_list(tempNet, __i, __j);
                }

            }
            else if (op.Equals("Add"))
            {
                /* step (1): Remove (i,j) from A_list */
                A_list.Remove(_i.ToString() + ":" + _j.ToString());

                /* step (2): Remove (i,j) from A_list */
                A_list.Remove(_j.ToString() + ":" + _i.ToString());

                /* step (3): Check Addable edges only in 'A_list' for update */
                //Check_Selected_Edges_in_A_list_Only(ref tempNet, op);
                Dictionary<string, string> addable_list = (from entry in A_list where entry.Value == "Y" select entry).ToDictionary(x => x.Key, x => x.Value);
                foreach (string elem in addable_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_A_list(tempNet, __i, __j);
                }
            }
            else if (op.Equals("Rev"))
            {
                ///* step (1): Remove (i,j) from A_list */
                //bool retStatus = A_list.Remove(_j.ToString() + ":" + _i.ToString());
                //if (retStatus)
                //{
                //    int xxxxx = 0;
                //}

                /* now it has to check the whole list */
                //Check_Selected_Edges_in_A_list_Only(ref tempNet, op);
                foreach (string elem in A_list.Keys.ToList())
                {
                    string __edge = elem.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                    int __i = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                    int __j = Int32.Parse(__edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                    Change_A_list(tempNet, __i, __j);
                }
            }
        }


        private void Check_Selected_Edges_in_A_list_Only(ref int[][] tempNet, string op)
        {
            foreach (string _e in A_list.Keys.ToList())
            {
                string _edge = _e.Split("#".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0];

                int _i = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[0]);
                int _j = Int32.Parse(_edge.Split(":".ToCharArray(), StringSplitOptions.RemoveEmptyEntries)[1]);

                if (op.Equals("Add"))
                {
                    /* only the addable edges in A_list will be checked here */
                    if (A_list[_edge].Equals("Y"))
                        Change_A_list(tempNet, _i, _j);
                }
                else if (op.Equals("Del"))
                {
                    /* only the NON-addable edges in A_list will be checked here */
                    if (A_list[_edge].StartsWith("N,"))
                        Change_A_list(tempNet, _i, _j);
                }
                else if (op.Equals("Rev"))
                {
                    /* all the edges in A_list will be checked here */
                    Change_A_list(tempNet, _i, _j);
                }
            }

        }

        private void Change_A_list(int[][] tempNet, int _i, int _j)
        {
            A_list.Remove(_i.ToString() + ":" + _j.ToString()); /* just to prevent multiple adding; as this pais is going to be added in next steps anyway */

            #region Add operations
            bool _isCycle = _isCycleExist(tempNet, _i, _j);
            if (!_isCycle)
            {
                bool _is_nParent_Restricted = is_nParent_Restricted(tempNet, _i, _j);
                bool _is_nChild_Restricted = is_nChild_Restricted(tempNet, _i, _j);
                if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                {
                    A_list[_i.ToString() + ":" + _j.ToString()] = "Y";
                }
                else
                    A_list[_i.ToString() + ":" + _j.ToString()] = "N,degree_th_Exceed";
            }
            else
                A_list[_i.ToString() + ":" + _j.ToString()] = "N,cycle_detected";
            #endregion
        }


        private double CalcuateTargetFunction(int[][] T_x)
        {
            #region Make Conditional Distributions for each nodes
            List<NameValueCollection> ConditionalProbTable_List = new List<NameValueCollection>();
            List<NameValueCollection> PriorBeliefTable_List = new List<NameValueCollection>();

            for (int i = 0; i < nNodes; i++)
            {
                /* 'i' indicates node_i */

                /* for each node create a namevaluecollection where the keyName:= "Combination of parents (, seperated) ;; value is the frequency "*/
                NameValueCollection nodeConditionalProbTable = new NameValueCollection();
                NameValueCollection nodeParents = Find_nParents_of_node(T_x, i);

                if (nodeParents.AllKeys.Count() > 0)
                {
                    #region For the nodes which has got at least one parent

                    #region Generate all possible parental conditions
                    List<int[]> input = new List<int[]>();

                    for (int j = 0; j < nodeParents.AllKeys.Count(); j++)
                    {
                        string key = nodeParents.AllKeys[j];
                        int nCat = Int32.Parse(DataCategories[Int32.Parse(key)]);
                        int[] tempArr = new int[nCat];
                        for (int k = 0; k < nCat; k++)
                            tempArr[k] = k;

                        input.Add(tempArr);
                    }

                    int[] size = new int[input.Count()];

                    /* enumerate all possible parent categories */
                    combine(ref nodeConditionalProbTable, input, size, 0);
                    #endregion

                    #region Explore the whole dataset for finding the frequency
                    foreach (string condition in nodeConditionalProbTable.AllKeys)
                    {
                        for (int j = 0; j < nObservation; j++)
                        {
                            string tempS = "";
                            /* this is just the initialization of the node observation value */
                            int nodeObservationValue = 0;

                            for (int k = 0; k < nNodes; k++)
                            {
                                /* if the index 'k' is in the parent list*/
                                if (nodeParents[k.ToString()] != null)
                                    tempS += DataSet[k][j] + ":";

                                /* saving the node observation value */
                                if (k == i)
                                    nodeObservationValue = DataSet[k][j];
                            }
                            tempS = tempS.Substring(0, tempS.Length - 1);

                            /* if the parental condition list matches with the observation list, then the observation into the CT */
                            if (condition.Equals(tempS))
                                nodeConditionalProbTable.Add(condition, nodeObservationValue.ToString());
                        }
                    }
                    #endregion

                    ConditionalProbTable_List.Add(nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    NameValueCollection nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[i]));
                    PriorBeliefTable_List.Add(nodePriorBeliefTable);
                    #endregion

                    #endregion
                }
                else if (nodeParents.AllKeys.Count() == 0)
                {
                    #region For the nodes which have no parents

                    #region Explore the whole dataset for finding the frequency

                    for (int j = 0; j < nObservation; j++)
                    {
                        /* this is just the initialization of the node observation value */
                        int nodeObservationValue = 0;

                        for (int k = 0; k < nNodes; k++)
                        {
                            /* saving the node observation value */
                            if (k == i)
                                nodeObservationValue = DataSet[k][j];
                        }

                        /* since the node_i has no parent, save it's id: 'i' as the key */
                        nodeConditionalProbTable.Add(i.ToString(), nodeObservationValue.ToString());
                    }
                    #endregion

                    ConditionalProbTable_List.Add(nodeConditionalProbTable);

                    #region Calculate Prior Table for the node_i and Add that into the list
                    NameValueCollection nodePriorBeliefTable = MakeNodePriorBeliefTable(nodeConditionalProbTable, Int32.Parse(DataCategories[i]));
                    PriorBeliefTable_List.Add(nodePriorBeliefTable);
                    #endregion
                    #endregion
                }
                /* for each category of node_i: make all*/
            }
            #endregion


            #region calculate function value
            double log_F1 = 0.0;

            for (int i = 0; i < nNodes; i++)
            {
                #region for each node (Random variable), calculate the value of density function
                NameValueCollection nodeConditionalProbTable = ConditionalProbTable_List[i];
                NameValueCollection nodePriorBeliefTable = PriorBeliefTable_List[i];

                double sum_conditionVal = 0.0;
                foreach (string condition in nodeConditionalProbTable.AllKeys)
                {
                    string[] arr_alpha_ij = nodePriorBeliefTable[condition].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    #region Calculate alpha_ij
                    double alpha_ij = 0.0;
                    foreach (string s in arr_alpha_ij)
                        alpha_ij += double.Parse(s);
                    #endregion

                    int _nCat = Int32.Parse(DataCategories[i]);
                    double[] arr_N_ij = FindConditionalProbabilityValues(nodeConditionalProbTable[condition], _nCat);
                    #region Calculate N_ij
                    double N_ij = 0.0;
                    foreach (int val in arr_N_ij)
                        N_ij += (double)val;
                    #endregion

                    double sum_categoricalVal = 0.0;
                    for (int k = 0; k < _nCat; k++)
                    {
                        double _categoricalVal = LogGamma(arr_N_ij[k] + double.Parse(arr_alpha_ij[k])) - LogGamma(double.Parse(arr_alpha_ij[k]));
                        sum_categoricalVal += _categoricalVal;
                    }

                    double _conditionVal = sum_categoricalVal + (LogGamma(alpha_ij) - LogGamma(N_ij + alpha_ij));
                    sum_conditionVal += _conditionVal;
                }
                #endregion

                log_F1 += sum_conditionVal;
            }

            #region Return final value for the target function
            _global_LogFx = log_F1;

            //double F1 = Math.Exp(log_F1);
            //return F1;
            return log_F1;
            #endregion

            #endregion


        }

        private double[] FindConditionalProbabilityValues(string probLine, int nCat)
        {
            double[] retArr = new double[nCat];

            string[] Line = probLine.Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            for (int i = 0; i < nCat; i++)
            {
                int count = 0;
                foreach (string s in Line)
                {
                    if (Int32.Parse(s) == i)
                        count++;
                }
                retArr[i] = (double)count;
            }

            return retArr;
        }

        public static double Gamma(double x)
        {
            // We require x > 0

            if (x <= 0.0)
            {
                string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
            }

            // Split the function domain into three intervals:
            // (0, 0.001), [0.001, 12), and (12, infinity)

            ///////////////////////////////////////////////////////////////////////////
            // First interval: (0, 0.001)
            //
            // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
            // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
            // The relative error over this interval is less than 6e-7.

            const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

            if (x < 0.001)
                return 1.0 / (x * (1.0 + gamma * x));

            ///////////////////////////////////////////////////////////////////////////
            // Second interval: [0.001, 12)

            if (x < 12.0)
            {
                // The algorithm directly approximates gamma over (1,2) and uses
                // reduction identities to reduce other arguments to this interval.

                double y = x;
                int n = 0;
                bool arg_was_less_than_one = (y < 1.0);

                // Add or subtract integers as necessary to bring y into (1,2)
                // Will correct for this below
                if (arg_was_less_than_one)
                {
                    y += 1.0;
                }
                else
                {
                    n = (int)(Math.Floor(y)) - 1;  // will use n later
                    y -= n;
                }

                // numerator coefficients for approximation over the interval (1,2)
                double[] p =
                {
                    -1.71618513886549492533811E+0,
                     2.47656508055759199108314E+1,
                    -3.79804256470945635097577E+2,
                     6.29331155312818442661052E+2,
                     8.66966202790413211295064E+2,
                    -3.14512729688483675254357E+4,
                    -3.61444134186911729807069E+4,
                     6.64561438202405440627855E+4
                };

                // denominator coefficients for approximation over the interval (1,2)
                double[] q =
                {
                    -3.08402300119738975254353E+1,
                     3.15350626979604161529144E+2,
                    -1.01515636749021914166146E+3,
                    -3.10777167157231109440444E+3,
                     2.25381184209801510330112E+4,
                     4.75584627752788110767815E+3,
                    -1.34659959864969306392456E+5,
                    -1.15132259675553483497211E+5
                };

                double num = 0.0;
                double den = 1.0;
                int i;

                double z = y - 1;
                for (i = 0; i < 8; i++)
                {
                    num = (num + p[i]) * z;
                    den = den * z + q[i];
                }
                double result = num / den + 1.0;

                // Apply correction if argument was not initially in (1,2)
                if (arg_was_less_than_one)
                {
                    // Use identity gamma(z) = gamma(z+1)/z
                    // The variable "result" now holds gamma of the original y + 1
                    // Thus we use y-1 to get back the orginal y.
                    result /= (y - 1.0);
                }
                else
                {
                    // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
                    for (i = 0; i < n; i++)
                        result *= y++;
                }

                return result;
            }

            ///////////////////////////////////////////////////////////////////////////
            // Third interval: [12, infinity)

            if (x > 171.624)
            {
                // Correct answer too large to display. 
                return double.PositiveInfinity;
            }

            return Math.Exp(LogGamma(x));
        }

        public static double LogGamma(double x)
        {
            // x must be positive
            if (x <= 0.0)
            {
                string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
            }

            if (x < 12.0)
            {
                return Math.Log(Math.Abs(Gamma(x)));
            }

            // Abramowitz and Stegun 6.1.41
            // Asymptotic series should be good to at least 11 or 12 figures
            // For error analysis, see Whittiker and Watson
            // A Course in Modern Analysis (1927), page 252

            double[] c =
            {
		         1.0/12.0,
		        -1.0/360.0,
		        1.0/1260.0,
		        -1.0/1680.0,
		        1.0/1188.0,
		        -691.0/360360.0,
		        1.0/156.0,
		        -3617.0/122400.0
            };
            double z = 1.0 / (x * x);
            double sum = c[7];
            for (int i = 6; i >= 0; i--)
            {
                sum *= z;
                sum += c[i];
            }
            double series = sum / x;

            double halfLogTwoPi = 0.91893853320467274178032973640562;
            double logGamma = (x - 0.5) * Math.Log(x) - x + halfLogTwoPi + series;
            return logGamma;
        }

        private NameValueCollection MakeNodePriorBeliefTable(NameValueCollection nodeConditionalProbTable, int r_i)
        {
            NameValueCollection retNm = new NameValueCollection();
            foreach (string key in nodeConditionalProbTable.AllKeys)
            {
                double q_i = (double)nodeConditionalProbTable.AllKeys.Count();
                double alpha = 1.0;
                double val = (double)(alpha / (q_i * (double)r_i));

                for (int i = 0; i < r_i; i++)
                {
                    retNm.Add(key, val.ToString());
                }
            }
            return retNm;
        }

        private void combine(ref NameValueCollection nodeConditionalProbTable, System.Collections.Generic.List<int[]> input, int[] current, int k)
        {
            if (k == input.Count())
            {
                string tempS = "";
                for (int i = 0; i < k; i++)
                {
                    //Console.Write(current[i] + " ");
                    tempS += current[i].ToString() + ":";

                }
                tempS = tempS.Substring(0, tempS.Length - 1);
                nodeConditionalProbTable.Add(tempS, "");
            }
            else
            {
                for (int j = 0; j < input[k].Length; j++)
                {
                    current[k] = input[k][j];
                    combine(ref nodeConditionalProbTable, input, current, k + 1);
                }

            }
        }

        private int Find_OperableEdges(Dictionary<string, string> _list, ref Dictionary<string, string> NeighbourhoodList, string operation)
        {
            int count = 0;
            foreach (var key in _list.Keys)
            {
                if (_list[key].Equals("Y"))
                {
                    count++;
                    NeighbourhoodList.Add(key + "#" + operation, operation);
                }
            }
            return count;
        }

        private NameValueCollection Find_nParents_of_node(int[][] T_x, int i)
        {
            NameValueCollection nm = new NameValueCollection();
            for (int j = 0; j < nNodes; j++)
            {
                if (T_x[j][i] == 1)
                    nm.Add(j.ToString(), "");
            }
            return nm;
        }

        private void PopulateAllLists(int[][] X)
        {
            PopulateDeletaleEdges(X);
            PopulateAddableEdges(X);
            PopulateReversibleEdges(X);
        }

        private void PopulateReversibleEdges(int[][] X)
        {
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = i + 1; j < nNodes; j++)
                {
                    if (X[i][j] == 0 && X[j][i] == 1)
                    {
                        #region Reverse
                        /* firt remove the edge and the see what happens */
                        X[j][i] = 0;

                        /* check adding the reverse edge creates a cycle or not */
                        bool _isCycle = _isCycleExist(X, i, j);
                        if (!_isCycle)
                        {
                            bool _is_nParent_Restricted = is_nParent_Restricted(X, i, j);
                            bool _is_nChild_Restricted = is_nChild_Restricted(X, i, j);
                            if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                                R_list.Add(j.ToString() + ":" + i.ToString(), "Y");
                            else
                                R_list.Add(j.ToString() + ":" + i.ToString(), "N,degree_th_Exceed");
                        }
                        else
                            R_list.Add(j.ToString() + ":" + i.ToString(), "N,cycle_detected");
                        /* firt remove the edge and the see what happens */
                        X[j][i] = 1;

                        #endregion
                    }
                    else if (X[j][i] == 0 && X[i][j] == 1)
                    {
                        #region Reverse
                        /* firt remove the edge and the see what happens */
                        X[i][j] = 0;

                        /* check adding the reverse edge creates a cycle or not */
                        bool _isCycle = _isCycleExist(X, j, i);
                        if (!_isCycle)
                        {
                            bool _is_nParent_Restricted = is_nParent_Restricted(X, j, i);
                            bool _is_nChild_Restricted = is_nChild_Restricted(X, j, i);
                            if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                                R_list.Add(i.ToString() + ":" + j.ToString(), "Y");
                            else
                                R_list.Add(i.ToString() + ":" + j.ToString(), "N,degree_th_Exceed");
                        }
                        else
                            R_list.Add(i.ToString() + ":" + j.ToString(), "N,cycle_detected");
                        /* firt remove the edge and the see what happens */
                        X[i][j] = 1;

                        #endregion
                    }
                }
            }
        }

        private void PopulateDeletaleEdges(int[][] X)
        {
            Dictionary<string, string> BridgeList = FindBridgesinNetwork(X);

            /* checking all the edges and saving only Deletable Edges */
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = i + 1; j < nNodes; j++)
                {
                    #region Try for (i,j), else for (j,i)
                    if (X[i][j] == 1)
                    {
                        if (!BridgeList.ContainsKey(i.ToString() + "," + j.ToString()) && !BridgeList.ContainsKey(j.ToString() + "," + i.ToString()))
                            D_list.Add(i.ToString() + ":" + j.ToString(), "Y");
                        else
                            D_list.Add(i.ToString() + ":" + j.ToString(), "N,bridge_detected");
                    }
                    else if (X[j][i] == 1)
                    {
                        if (!BridgeList.ContainsKey(j.ToString() + "," + i.ToString()) && !BridgeList.ContainsKey(i.ToString() + "," + j.ToString()))
                            D_list.Add(j.ToString() + ":" + i.ToString(), "Y");
                        else
                            D_list.Add(j.ToString() + ":" + i.ToString(), "N,bridge_detected");
                    }
                    #endregion
                }
            }

        }

        private void dfs(int[][] X, int u, int v, int cnt, int[] low, int[] pre, ref Dictionary<string, string> bridgeList)
        {
            pre[v] = cnt++;
            low[v] = pre[v];

            for (int w = 0; w < nNodes; w++)
            {
                if (X[v][w] == 1 || X[w][v] == 1)
                {
                    if (pre[w] == -1)
                    {
                        dfs(X, v, w, cnt, low, pre, ref bridgeList);
                        low[v] = Math.Min(low[v], low[w]);
                        if (low[w] == pre[w])
                        {
                            //MessageBox.Show(v + "-" + w + " is a bridge");
                            //totalCount++;
                            //bridgeList.Add(totalCount.ToString(), v.ToString() + "," + w.ToString());
                            bridgeList.Add(v.ToString() + "," + w.ToString(), "");
                        }
                    }
                    // update low number - ignore reverse of edge leading to v
                    else if (w != u)
                        low[v] = Math.Min(low[v], pre[w]);
                }
            }

        }

        private void PopulateAddableEdges(int[][] X)
        {
            for (int i = 0; i < nNodes; i++)
            {
                for (int j = i + 1; j < nNodes; j++)
                {
                    /* if the edge can be added and there is no chance of creating cycle after adding that edge */
                    if (X[i][j] == 0 && X[j][i] == 0)
                    {
                        #region Try for (i,j)
                        /* try for (i,j) */
                        bool _isCycle = _isCycleExist(X, i, j);
                        if (!_isCycle)
                        {
                            bool _is_nParent_Restricted = is_nParent_Restricted(X, i, j);
                            bool _is_nChild_Restricted = is_nChild_Restricted(X, i, j);
                            if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                            {
                                A_list.Add(i.ToString() + ":" + j.ToString(), "Y");
                            }
                            else
                                A_list.Add(i.ToString() + ":" + j.ToString(), "N,degree_th_Exceed");
                        }
                        else
                            A_list.Add(i.ToString() + ":" + j.ToString(), "N,cycle_detected");
                        #endregion

                        #region Try for (j,i)
                        /* try for (j,i) */
                        _isCycle = _isCycleExist(X, j, i);
                        if (!_isCycle)
                        {
                            bool _is_nParent_Restricted = is_nParent_Restricted(X, j, i);
                            bool _is_nChild_Restricted = is_nChild_Restricted(X, j, i);
                            if (!_is_nParent_Restricted && !_is_nChild_Restricted)
                            {
                                A_list.Add(j.ToString() + ":" + i.ToString(), "Y");
                            }
                            else
                                A_list.Add(j.ToString() + ":" + i.ToString(), "N,degree_th_Exceed");
                        }
                        else
                            A_list.Add(j.ToString() + ":" + i.ToString(), "N,cycle_detected");
                        #endregion
                    }
                }
            }
        }

        private void FindaRandomNetwork(ref int[][] _net)
        {
            int[][] net = new int[nNodes][];
            bool[][] flag = new bool[nNodes][];

            while (true)
            {
                for (int i = 0; i < nNodes; i++)
                {
                    _net[i] = new int[nNodes];
                    flag[i] = new bool[nNodes];
                    for (int j = 0; j < nNodes; j++)
                    {
                        _net[i][j] = 0;
                        if (i == j)
                            flag[i][j] = true;
                    }
                }
                for (int counter = 0; counter < nNodes * nNodes; )
                {
                    /* sample edge uniformly */
                    int _i = rand.Next(0, nNodes);
                    int _j = rand.Next(0, nNodes);

                    //int _i = n.Next();
                    //int _j = n.Next();

                    if (_i != _j)
                    {
                        if (_net[_i][_j] == 0 && _net[_j][_i] == 0)
                        {
                            double dec = rand.NextDouble();
                            //int dec = rand.Next(0,2);
                            counter++;

                            if (dec <= 0.5)
                            //if (dec == 1)
                            {
                                bool _isCycle = _isCycleExist(_net, _i, _j);
                                if (!_isCycle)
                                {
                                    bool _is_nParent_Restricted = is_nParent_Restricted(_net, _i, _j);
                                    bool _is_nChild_Restricted = is_nChild_Restricted(_net, _i, _j);
                                    if ((!_is_nParent_Restricted && !_is_nChild_Restricted))
                                    //if (!_isCycle && (dec == 1))
                                    {
                                        _net[_i][_j] = 1;
                                    }

                                    //_net[_i][_j] = 1;
                                }
                                else
                                {
                                    _isCycle = _isCycleExist(_net, _j, _i);
                                    if (!_isCycle)
                                    {
                                        bool _is_nParent_Restricted = is_nParent_Restricted(_net, _j, _i);
                                        bool _is_nChild_Restricted = is_nChild_Restricted(_net, _j, _i);
                                        if ((!_is_nParent_Restricted && !_is_nChild_Restricted))
                                        //if (!_isCycle && (dec == 1))
                                        {
                                            _net[_j][_i] = 1;
                                        }

                                        //_net[_j][_i] = 1;
                                    }
                                }
                            }
                        }
                    }
                }

                #region check if the random DAG is a connected one
                bool _isConnected = isConnected(_net, 0);

                if (_isConnected)
                    break;
                else
                    continue;
                #endregion

            }

        }

        private bool is_nChild_Restricted(int[][] _X, int _i, int _j)
        {
            //throw new NotImplementedException();
            _X[_i][_j] = 1;

            bool f = nChild_check(_X);

            _X[_i][_j] = 0;

            return f;
        }

        private bool nChild_check(int[][] _X)
        {
            for (int i = 0; i < nNodes; i++)
            {
                int count = 0;
                for (int j = 0; j < nNodes; j++)
                {
                    if (_X[i][j] == 1)
                        count++;
                }
                if (count > int.Parse(nChild_th[i.ToString()]))
                    return true;
            }
            return false;
        }

        private bool is_nParent_Restricted(int[][] _X, int i, int j)
        {
            //throw new NotImplementedException();
            _X[i][j] = 1;

            bool f = nParent_check(_X);

            _X[i][j] = 0;

            return f;
        }

        private bool nParent_check(int[][] _X)
        {
            for (int i = 0; i < nNodes; i++)
            {
                int count = 0;
                for (int j = 0; j < nNodes; j++)
                {
                    if (_X[j][i] == 1)
                        count++;
                }
                if (count > int.Parse(nParent_th[i.ToString()]))
                    return true;
            }
            return false;
        }

        private bool isEdgeDeletable(int[][] X, int i, int j)
        {

            /* check if that removal disconnects the graph */
            bool _isConnected1 = isConnected(X, i);
            //bool _isConnected2 = isConnected(X, j);

            if (!_isConnected1)
            {
                /* if that removal disconnects the graph, then that edge not deletable */
                return false;
            }
            else
                return true;
        }

        private bool isConnected(int[][] net, int row)
        {
            bool dec = true;

            for (int i = 0; i < nNodes; i++)
            {
                if (i != row)
                {
                    dec = isSimplePathExist(net, row, i);
                    if (!dec)
                        return dec;
                    else
                        continue;
                }
            }

            return dec;
        }

        private bool isSimplePathExist(int[][] net, int start, int end)
        {
            bool dec = false;
            int front, rear;
            int[] que = new int[nNodes];
            front = rear = -1;
            int v = start;
            int n = nNodes;

            bool[] visited = new bool[nNodes];
            for (int i = 0; i < n; i++)
                visited[i] = false;

            visited[v] = true;
            rear++;
            front++;
            que[rear] = v;

            while (front <= rear)
            {
                v = que[front];
                front++;

                for (int i = 0; i < n; i++)
                {
                    if ((net[v][i] == 1 || net[i][v] == 1) && visited[i] == false)
                    {
                        if (i == end)
                        {
                            dec = true;
                            break;
                        }
                        visited[i] = true;
                        rear++;
                        que[rear] = i;
                    }
                }
            }

            return dec;
        }

        private int lenSimplePath(int[][] net, int start, int end)
        {
            bool dec = false;
            int front, rear;
            int[] que = new int[nNodes];
            front = rear = -1;
            int v = start;
            int n = nNodes;
            int count = 0;

            bool[] visited = new bool[nNodes];
            for (int i = 0; i < n; i++)
                visited[i] = false;

            visited[v] = true;
            rear++;
            front++;
            que[rear] = v;

            while (front <= rear)
            {
                v = que[front];
                front++;

                for (int i = 0; i < n; i++)
                {
                    if (net[v][i] == 1 && visited[i] == false)
                    {
                        count++;
                        if (i == end)
                        {
                            dec = true;
                            break;
                        }
                        visited[i] = true;
                        rear++;
                        que[rear] = i;
                    }
                }
            }

            return count;
        }

        private bool _isCycleExist(int[][] _X, int i, int j)
        {
            //throw new NotImplementedException();
            _X[i][j] = 1;

            bool f = cycle_check(_X);

            _X[i][j] = 0;

            return f;
        }
        private bool cycle_check(int[][] _X)
        {
            bool[] visited = new bool[nNodes];
            bool[] recStack = new bool[nNodes];

            for (int i = 0; i < nNodes; i++)
            {
                visited[i] = false;
                recStack[i] = false;
            }
            for (int u = 0; u < nNodes; u++)
                if (isCyclicUtil(_X, u, visited, recStack))
                    return true;

            return false;
        }

        private bool isCyclicUtil(int[][] _X, int v, bool[] visited, bool[] recStack)
        {
            if (visited[v] == false)
            {
                // Mark the current node as visited and part of recursion stack
                visited[v] = true;
                recStack[v] = true;

                for (int i = 0; i < nNodes; ++i)
                {
                    if (_X[v][i] == 1)
                    {
                        if (!visited[i] && isCyclicUtil(_X, i, visited, recStack))
                            return true;
                        else if (recStack[i])
                            return true;
                    }
                }
            }
            recStack[v] = false;  // remove the vertex from recursion stack
            return false;
        }
        #endregion


    }
}
