using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Preprocessing_MCMC_sampling
{
    class processData
    {
        private double th_dic_high;
        private double th_dic_low;

        public processData(double th_dic_high, double th_dic_low)
        {
            this.th_dic_high = th_dic_high;
            this.th_dic_low = th_dic_low;
        }

        public Dictionary<String, int[]> getDiscritizedMat(Dictionary<String, double[]> controlMat, Dictionary<String, double[]> caseMat) {
            Dictionary<String, int[]> retMat = new Dictionary<string, int[]>();
            foreach(string key in caseMat.Keys){
                double ctl_mean = controlMat[key].Average();
                double ctl_std = getStd(controlMat[key]);

                int[] new_val = getStandardizedVal(caseMat[key], ctl_mean,ctl_std);
                retMat.Add(key, new_val);
            }
            return retMat;
        }

        public Dictionary<String, double[]> getStandardizedVal(Dictionary<String, double[]> controlMat, Dictionary<String, double[]> caseMat)
        {
            Dictionary<String, double[]> retMat = new Dictionary<string, double[]>();
            foreach (string key in caseMat.Keys)
            {
                double ctl_mean = controlMat[key].Average();
                double ctl_std = getStd(controlMat[key]);

                double[] new_val = new double[caseMat[key].Length];
                for (int i = 0; i < caseMat[key].Length; i++)
                    new_val[i] = (((double)caseMat[key][i]) - ctl_mean) / ctl_std;
                retMat.Add(key, new_val);
            }
            return retMat;
        }
        private int[] getStandardizedVal(double[] values, double ctl_mean, double ctl_std)
        {
            int[] retArr = new int[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                double tempVal = (((double)values[i]) - ctl_mean) / ctl_std;
                retArr[i] = (tempVal > th_dic_high) ? 1 : (tempVal < th_dic_low) ? -1 : 0; 
            }
            return retArr;
        }

        private double getStd(double[] values)
        {
            double ret = 0;
            int count = values.Count();
            if (count > 1)
            {
                //Compute the Average
                double avg = values.Average();

                //Perform the Sum of (value-avg)^2
                double sum = values.Sum(d => (d - avg) * (d - avg));

                //Put it all together
                ret = Math.Sqrt(sum / count);
            }
            return ret;
        }

        internal Dictionary<string, Dictionary<string, int[]>> getDiscritizedMatForPathways(Dictionary<string, List<string>> allPathways, Dictionary<string, int[]> rMat_discritized)
        {
            Dictionary<string, Dictionary<string, int[]>> retDat = new Dictionary<string, Dictionary<string, int[]>>(); // key: pathwayName; value: dataMatrix for that pathway
            foreach(String pathway in allPathways.Keys) {
                List<String> pGenes = allPathways[pathway];
                Dictionary<string, int[]> pData = getPathwayData(pGenes, rMat_discritized);
                retDat.Add(pathway, pData);
            }
            return retDat;
        }

        private Dictionary<string, int[]> getPathwayData(List<string> pGenes, Dictionary<string, int[]> rMat_discritized)
        {
            Dictionary<string, int[]> retDat = new Dictionary<string, int[]>();
            foreach(string gene in pGenes) {
                if(rMat_discritized.ContainsKey(gene))          // not all the pathwayGene will have data
                    retDat.Add(gene, rMat_discritized[gene]);
            }
            return retDat;
        }

        //public Dictionary<String, double[]> getStandardizedDataforAPathway(string pathway, Dictionary<string, List<string>> allPathways, Dictionary<String, double[]> rMat_standardizedOnly)
        //{
        //    List<String> pGenes = allPathways[pathway];
        //    Dictionary<string, double[]> retDat = new Dictionary<string, double[]>();
        //    foreach (string gene in pGenes)
        //    {
        //        if (rMat_standardizedOnly.ContainsKey(gene))          // not all the pathwayGene will have data
        //            retDat.Add(gene, rMat_standardizedOnly[gene]);
        //    }
        //    return retDat;
        //}
    }
}
