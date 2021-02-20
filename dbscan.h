
#ifndef _DBSCAN_

#define _DBSCAN_



#include "utils.h"
#include "clusters.h"

namespace NWUClustering {

  class ClusteringAlgo:public Clusters {

    public:

      ClusteringAlgo(){ }

      virtual ~ClusteringAlgo();

      // functions for dbscan algorithm
      void set_dbscan_params(double eps, int minPts, int seeds);
      void writeClusters(ostream& o); // regular dbscan algorithm
      void writeClusters_uf(ostream& o); // union find dbscan algorithm

      void getGrowingPoints(vector<int>& growing_points, int sch, int tid, int lower, int upper);

    public:

      // parameters to run dbscan algorithm

      double  m_epsSquare;
      int   m_minPts;
      int m_seeds;

      // int qualitypoints;

      //int     m_parcent_of_data;
      // noise vector

      vector<bool> m_noise;

      // noise vector

      vector<bool> m_visited;
      vector <int> m_parents;
      vector <int> m_corepoint;
      vector <int> m_member;

      // vector <int> growing_points;

  };  

  void run_dbscan_algo_uf(ClusteringAlgo& dbs); // union find dbscan algorithm
  void run_dbscan_algo(ClusteringAlgo& dbs); // regular dbscan algorithm
  void unionize_neighborhood(ClusteringAlgo& dbs, int& root, int root1, int root2, bool locks, omp_lock_t* actualLock, int tid);
};



#endif
