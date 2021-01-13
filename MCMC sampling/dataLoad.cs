using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections.Specialized;
using System.IO;

namespace Preprocessing_MCMC_sampling
{
    class dataLoad
    {
        public Dictionary<String,double[]> LoadMatrix(String fileName)
        {
            Dictionary<String, double[]> retMat = new Dictionary<string, double[]>();
            TextReader tr = new StreamReader(fileName);
            while (true){
                string st = tr.ReadLine();
                if (st == null)
                    break;
                string[] stall = st.Split(",".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                if(!retMat.ContainsKey(stall[0].Trim()))    // decide this later, what will happen if there are multiple data vectors for the same gene
                    retMat.Add(stall[0].Trim(), getArr(stall));
            }
            return retMat;
        }
        public Dictionary<String, List<String>> LoadPathways (String fileName) {
            Dictionary<String, List<String>> retList = new Dictionary<string, List<string>>();
            TextReader tr = new StreamReader(fileName);
            while (true)
            {
                string st = tr.ReadLine();
                if (st == null)
                    break;
                string[] stall = st.Split("\t".ToCharArray(),StringSplitOptions.RemoveEmptyEntries);
                List<String> tempList = stall.ToList();
                tempList.RemoveAt(0);   // remove pathwayName
                tempList.RemoveAt(0);   // remove dummy text
                retList.Add(stall[0], tempList);        // key: pathwayName, Value: geneList
            }
            return retList;
        }
        private double[] getArr(string[] stall)
        {
            List<double> retArr = new List<double>();
            for(int i = 1; i < stall.Length; i++)
                retArr.Add(double.Parse(stall[i]));
            return retArr.ToArray();
        }
    }
}
