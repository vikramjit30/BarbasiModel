#region MONO/NET System libraries
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Threading.Tasks;
#endregion

#region Much appreciated thirs party libraries
using MathNet.Numerics.Distributions;
#endregion

#region NETGen libraries
using NETGen.Core;
using NETGen.Visualization;
using NETGen.Dynamics.Synchronization;
using System.IO;
using System.Reflection;
using System.Collections;

namespace Launch
{
    static class LaunchBarbasiModel
    {
        static Kuramoto sync;
        static Network pop;
        static Dictionary<int, double> _clusterOrder = new Dictionary<int, double>();
        static Dictionary<int, bool> pacemaker_mode = new Dictionary<int, bool>();
        static double runTime = 10d;
        
        static Random rnd = new Random();
        // To read the parameter.config files for a type of network 
        public static void read_parameters(String configFile, GlobValues glob)
        {
            string[] allLines;
            ArrayList list = new ArrayList();
           

            using (StreamReader sr = File.OpenText(configFile))
            {
                string s = String.Empty;
                while ((s = sr.ReadLine()) != null)
                {
                    allLines = s.Split(new[] { "=" }, StringSplitOptions.None);
                    list.Add(allLines[1]);
                }

                glob.nodes = Int32.Parse(list[0].ToString());
                glob.minPower = Double.Parse(list[1].ToString());
                glob.maxPower = Double.Parse(list[2].ToString());
                glob.numberOfGraphs = Int32.Parse(list[3].ToString());
                glob.couplingStrength = Double.Parse(list[4].ToString());
                glob.couplingProb = Double.Parse(list[5].ToString());
                glob.runningTime = Double.Parse(list[6].ToString());
                runTime = glob.runningTime;
                
            }

        }

         /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main(string[] args)
        {
           
            string configFile, dir = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), netFile, NetworkFile, resFile,destResultFile;
            string[] path = dir.Split(new string[] { "Launch" }, StringSplitOptions.None);

            string srcResultFile = path[0] +"Launch" + Path.DirectorySeparatorChar + "Launch" + Path.DirectorySeparatorChar + "bin" + Path.DirectorySeparatorChar + "Debug" + Path.DirectorySeparatorChar + "result.dat";
            configFile = path[0] + Path.DirectorySeparatorChar + "Launch" + Path.DirectorySeparatorChar + "config.param.txt";
             
            Console.WriteLine("Starting the program to generate the network based on Barbasi Albert Model..");

            GlobValues glob = new GlobValues();
            try
            {
                // Read parameters from param.config file
                read_parameters(configFile, glob);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                return;
            }

         
       for (int i = 1; i<= glob.numberOfGraphs; i++)
            {
                for (double j = glob.minPower; j <= (glob.maxPower + 0.1); j=j+0.1)
                {
                    // Creating network file
                    netFile = i+ "_BarbasiNetwork_N"+glob.nodes+ "_powerLaw"+j+"_K"+glob.couplingStrength+".edges";
                    NetworkFile =  path[0] + "Launch" + Path.DirectorySeparatorChar + "output" + Path.DirectorySeparatorChar+netFile ;
                    resFile = i + "_res_N" + glob.nodes + "_powerLaw" + j + "_K" + glob.couplingStrength+".dat";
                    destResultFile     =  path[0] + "Launch" + Path.DirectorySeparatorChar + "output" + Path.DirectorySeparatorChar + resFile;
                    try
                        {

                                 // upload the network to run the Kuramoto Model
                                 pop = Network.LoadFromEdgeFile(NetworkFile);
                                
                                // Run the Kuramoto model here and store the results in the output directory
                                NetworkColorizer colorizer = new NetworkColorizer();
                                // Distribution of natural frequencies
                                double mean_frequency = 1d;
                                Normal normal = new Normal(mean_frequency, mean_frequency / 5d);
                                                sync = new Kuramoto(pop,
                                                glob.couplingStrength,
                                                glob.couplingProb,
                                                colorizer,
                                                new Func<Vertex, Vertex[]>(v => { return new Vertex[] { v.RandomNeighbor }; })
                                                );

                                foreach (Vertex v in pop.Vertices)
                                    sync.NaturalFrequencies[v] = normal.Sample();

                                //  foreach (int g in network.ClusterIDs)
                                //    pacemaker_mode[g] = false;


                                sync.OnStep += new Kuramoto.StepHandler(recordOrder);

                                Logger.AddMessage(LogEntryType.AppMsg, "Press enter to start synchronization experiment...");
                                Console.ReadLine();

                                // Run the simulation
                                sync.Run();

                                // Write the time series to the resultfile 
                                if (srcResultFile != null)
                                    sync.WriteTimeSeries(srcResultFile);

                                // Moving results of kuramoto model into output directory
                                System.IO.File.Move(srcResultFile, destResultFile);
                                  
                       }
                            catch (Exception e)
                            {
                                Console.WriteLine("Error: " + e);
                            }
                             
               }
         }

   }
   
    private static void recordOrder(double time)
        {
           
            // Compute and record global order parameter
            double globalOrder = sync.GetOrder(pop.Vertices.ToArray());

            foreach (Vertex v in pop.Vertices)
                sync.AddDataPoint(v.Label, sync.CurrentValues[sync._mapping[v]]);

            //  if (time > 30d)
            if (time > runTime)
                sync.Stop();

            //   Logger.AddMessage(LogEntryType.SimMsg, string.Format("Time = {000000}", time)); //Avg. Cluster Order = {1:0.00}, Global Order = {2:0.00}", time, avgLocalOrder, globalOrder));

        }
      }
}           
 
  #endregion