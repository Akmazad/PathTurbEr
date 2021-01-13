using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using SharpJags;
using SharpJags.Jags;
using System.IO;
using SharpJags.Math;
using SharpJags.Parsing;
using System.Collections.Specialized;

namespace Preprocessing_MCMC_sampling
{

    class bugsSampling
    {

        public bugsSampling(double[,] mat, double[,] c_mat, int _nNodes, string pathwayName, NameValueCollection nodeNameList, string method, int nBurnInSteps, int nSamplingIteration, double gamma_prior_a, double gamma_prior_b)
        {
            SharpJags.Jags.JagsConfig.BinPath = @"C:\Program Files\JAGS\JAGS-4.3.0\x64\bin";

            // -- Set up data for MCMC sampling
            var n = _nNodes;  // should be number of genes
            Matrix<double> _y00a = new double[_nNodes, _nNodes];
            _y00a = c_mat;
            Matrix<double> _y11a = new double[_nNodes, _nNodes];
            _y11a = mat;

            //SharpJags.Math.Vector<double> epsilon = SharpJags.Math.Random.DistributionExtensions.Generate(new SharpJags.Math.Random.Gaussian(0, 1), n).ToArray();
            //SharpJags.Math.Vector<double> delta = SharpJags.Math.Random.DistributionExtensions.Generate(new SharpJags.Math.Random.Uniform(-1, 1), n).ToArray();

            // -- Define and instantiate "Data variables" in Model
            var data = new JagsData{
                { "N", n },
                { "y00a", _y00a },
                { "y11a", _y11a },
            };

            // --- Define Model
            var modelDefintion = new SharpJags.Models.ModelDefinition
            {
                Name = "Az_MCMC_" + pathwayName,
                Definition =
                @"
                model {
                        for (i in 1:N) {
                            for (j in 1:N) {

                                y00a[i,j] ~ dbern(p00a[i,j])
                                log(p00a[i,j]) <- lambdaa[i,j]
                                y11a[i,j] ~ dbern(p11a[i,j])
                                log(p11a[i,j]) <- lambdaa[i,j] + thetaa + a[i] + a[j]
                                lambdaa[i,j] <- -log(1 + exp(thetaa + a[i] + a[j]))

                            }
                        }

                    for (i in 1:N) { a[i] ~ dnorm(0,tau.thetaa) }
                    thetaa ~ dnorm(0,tau.thetaa)
                    tau.thetaa ~ dgamma(" + gamma_prior_a.ToString() + @"," + gamma_prior_b.ToString() + @")
                }
            "
            };

            // -- Define monitor variables
            var yMonitor = new JagsMonitor()
            {
                ParameterName = "y11a",
                Thin = 1
            };
            var alphaMonitor = new JagsMonitor()
            {
                ParameterName = "a",
                Thin = 1
            };

            // -- Set up MCMC run parameters
            var run = new JagsRun
            {
                WorkingDirectory = new FileInfo(Directory.GetCurrentDirectory()).ToString(),
                ModelDefinition = modelDefintion,
                ModelData = data,

                Monitors = new List<JagsMonitor>{
                    alphaMonitor
                },

                Parameters = new MCMCParameters
                {
                    BurnIn = nBurnInSteps,
                    SampleCount = nSamplingIteration,
                    Chains = 1
                }
            };

            // -- Setting Sample Reader
            ISampleReader sR = new SharpJags.Parsing.Coda.CodaDataReader();
            ISampleParser pars = new SharpJags.Parsing.Coda.CodaParser(sR);

            // -- Run MCMC sampling 
            var samples = JagsWrapper.Run(run, pars, pathwayName);

            // -- Read monitor variables from samples
            List<SharpJags.Parsing.IModelParameter> parameterSamples = new List<SharpJags.Parsing.IModelParameter>{
                samples.Get<ModelParameterVector>(alphaMonitor)
            };

            ModelParameterVector alphaVector = (ModelParameterVector)parameterSamples[0];

            TextWriter tw = new StreamWriter(@"JAGS_output\\" + pathwayName + "_alpha_" + method + ".csv");
            tw.WriteLine("geneID, alpha_mean, alpha_SD");
            foreach (SharpJags.Parsing.IModelParameter pm in alphaVector.Parameters)
                tw.WriteLine(nodeNameList[alphaVector.Parameters.IndexOf((ModelParameter)pm).ToString()] + "," + pm.Statistics.Mean.ToString() + "," + pm.Statistics.StandardDeviation.ToString());
            tw.Close();
        }
    }
}
