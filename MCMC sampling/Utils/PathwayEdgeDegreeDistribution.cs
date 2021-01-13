using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Collections.Specialized;

namespace Preprocessing_MCMC_sampling.Utils
{
    class PathwayEdgeDegreeDistribution
    {
        NameValueCollection humanNet_parentList = new NameValueCollection();
        NameValueCollection humanNet_childrenList = new NameValueCollection();

        public PathwayEdgeDegreeDistribution() {
            #region read Human Signalling pathway network
            TextReader tr = new StreamReader(@"C:\Users\aazad\Google Drive\Post PHD matters\Projects\Identifying key players in aberrant pathways\Data\HuamnSignalingNet_WangLab.csv");
            string st = tr.ReadLine();

            while (true) {
                st = tr.ReadLine();
                if (st == null)
                    break;

                string[] stall = st.Trim().Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                humanNet_parentList.Add(stall[3], stall[1]);
                humanNet_childrenList.Add(stall[1], stall[3]);
            }
            tr.Close();
            #endregion

            #region read KEGG pathway genes
            
            #endregion

            #region read Reactome pathway genes
            Dictionary<String, List<String>> Reactome = new Dictionary<string, List<string>>();
            tr = new StreamReader(@"C:\Users\aazad\Google Drive\Post PHD matters\Projects\Identifying key players in aberrant pathways\Data\ReactomeSignalingPathways.csv");
            while (true)
            {
                st = tr.ReadLine();
                if (st == null)
                    break;
                string[] stall = st.Split("\t".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                List<String> tempList = stall.ToList();
                tempList.RemoveAt(0);   // remove pathwayName
                tempList.RemoveAt(0);   // remove dummy text
                Reactome.Add(stall[0], tempList);        // key: pathwayName, Value: geneList
            }
            tr.Close();
            #endregion

            #region read WikiPathway pathway genes
            Dictionary<String, List<String>> WikiPathway = new Dictionary<string, List<string>>();
            tr = new StreamReader(@"C:\Users\aazad\Google Drive\Post PHD matters\Projects\Identifying key players in aberrant pathways\Data\ReactomeSignalingPathways.csv");
            while (true)
            {
                st = tr.ReadLine();
                if (st == null)
                    break;
                string[] stall = st.Split("\t".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                List<String> tempList = stall.ToList();
                tempList.RemoveAt(0);   // remove pathwayName
                tempList.RemoveAt(0);   // remove dummy text
                WikiPathway.Add(stall[0], tempList);        // key: pathwayName, Value: geneList
            }
            tr.Close();
            #endregion

            #region read Pathways and analyse
            String pathwayDB = @"C:\Users\aazad\Google Drive\Post PHD matters\Projects\Identifying key players in aberrant pathways\Data\KEGG_45_SIGNALING.csv";
            ReadPathwaysAndAnalyse(pathwayDB, humanNet_parentList, humanNet_childrenList);
            #endregion

        }

        private void ReadPathwaysAndAnalyse(string pathwayDBPath, NameValueCollection humanNet_parentList, NameValueCollection humanNet_childrenList)
        {
            #region Read pathways
            Dictionary<String, List<String>> pathwayDB = new Dictionary<string, List<string>>();
            TextReader tr = new StreamReader(pathwayDBPath);
            while (true)
            {
                string st = tr.ReadLine();
                if (st == null)
                    break;
                string[] stall = st.Split("\t".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                List<String> tempList = stall.ToList();
                tempList.RemoveAt(0);   // remove pathwayName
                tempList.RemoveAt(0);   // remove dummy text
                pathwayDB.Add(stall[0], tempList);        // key: pathwayName, Value: geneList

                //List<int>
                foreach(String pGene in tempList) {
                    List<String> childrenList = humanNet_childrenList[pGene].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).ToList();
                    List<String> parentList = humanNet_parentList[pGene].Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).ToList();

                    int outDegree = childrenList.ToArray().Intersect(tempList.ToArray()).Count();
                    int inDegree = parentList.ToArray().Intersect(tempList.ToArray()).Count();

                }

            }
            tr.Close();
            #endregion

            
        }
    }
}
