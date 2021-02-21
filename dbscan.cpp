#include "dbscan.h"
#include <ctime>
#include <set>
#include <iterator>
#include <algorithm>
#include <vector>

namespace NWUClustering {

  void ClusteringAlgo::set_dbscan_params(double eps, int minPts, int seeds) {
    m_epsSquare =  eps * eps;
    m_minPts =  minPts;
    m_seeds = seeds;
  }

  ClusteringAlgo::~ClusteringAlgo() {
    m_noise.clear();
    m_visited.clear();
    m_parents.clear();
    m_corepoint.clear();
    m_member.clear();
  }

  void ClusteringAlgo::writeClusters(ostream& o) {
    // writing point id and cluster id pairs per line, noise has cluster id 0 
    int iMaxID = m_clusters.size(), id, i, j;
    int num_points = m_pts->m_i_num_points;
    for(i = 0; i < num_points; i++) {
      //for (j = 0; j < m_pts->m_i_dims; j++)
      //  o << " " << m_pts->m_points[i][j];
      id = m_pid_to_cid[i];
      // o << i << " " << id << endl;
    }
    int sum_points = 0;
    int noise = 0;
    for(i = 0; i < m_clusters.size(); i++) {
      sum_points += m_clusters[i].size();
      //cout << i << "(" << m_clusters[i].size() << ") ";
    }
    for (i = 0; i < num_points; i++) {
      if(m_noise[i])
        noise++;
    } 
    cout << "Points in clusters " << sum_points << " Noise " << noise << " Total points " << noise + sum_points << endl;
    cout << "Total number of clusters " << m_clusters.size() << endl;
  }

  void ClusteringAlgo::writeClusters_uf(ostream& o) {
    //changed return type to int. Was void. Used for testing purposes
    // writing point id and cluster id pairs per line, noise has cluster id 0 
    vector <int> clusters;
    int num_points = m_pts->m_i_num_points;
    clusters.resize(num_points, 0);
    int i, j, sum_points = 0, noise = 0, root, rootcount = 0, tmp;
    for(i = 0; i < num_points; i++) {
      root = m_parents[i];
      // get the number of trees
      if(i == m_parents[i])
        rootcount++;
      // get the root of the tree containing i  
      while(root != m_parents[root])
        root = m_parents[root];
      // compress the tree to reduce the height of the branch of the tree to 1
      j = i;
      while(m_parents[j] != root) {
        tmp  = m_parents[j];
        m_parents[j] = root;
        j = tmp;
      }
      clusters[root]++;
    }
    //cout << "clusters" << endl;
    int count = 0;
    for(i = 0; i < num_points; i++) {
      if(clusters[i] == 1) {
        // vertex i is a noise
        clusters[i] = 0;
        noise++;
      // } else if(clusters[i] >= m_minPts) { TODO
      } else if(clusters[i] > 1) {
        // This conditional statement determines what is and is not considered a cluster
        // If it's greater than this number than it will be counted as a cluster
        count++;
        sum_points += clusters[i];
        clusters[i] = count;
        // cout << clusters[i] << ": " << i << endl;
      }
    }
    // for(i = 0; i < num_points; i++) {
    //   o << i << " " << clusters[m_parents[i]] << endl; //This was the only one not commented out
    // }
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      for (j = 0; j < m_pts->m_i_dims; j++)
        o << m_pts->m_points[i][j] << " ";  
      // o << i << " " << clusters[m_parents[i]] << endl;
      o << clusters[m_parents[i]] << endl;
    }
    
    cout << "Points in clusters " << sum_points << " Noise " << noise << " Total points " << noise + sum_points << endl;
    cout << "Total number of clusters " << count << endl;
    // o << "Points in clusters " << sum_points << " Noise " << noise << " Total points " << noise + sum_points << endl;
    // o << "Total number of clusters " << count << endl;
    clusters.clear();
  }

  /*
    Determines the initial seed points to be used in a pseudo-random fashion. If all of the points in a thread
    are to be used, the points are inspected linearally to speed up the execution time. 

    vector<int>& growing_points - Vector that will hold the seed points for each thread.
    int sch - the total number number of points per thread
    int tid - the thread ID number
    int lower - the lower-bound index to read the points
    int upper - the upper-bound index to read the points
  */
  void ClusteringAlgo::getGrowingPoints(vector<int>& growing_points, int sch, int tid, int lower, int upper) {
    int sid, number_looking_for = m_seeds;
    kdtree2_result_vector ne;
    vector<int>* ind = m_kdtree->getIndex(); // Sets a vector that contains the index of all points
    vector<int> alreadySeen;
    bool allPoints = false;
    srand(time(NULL));

    if(number_looking_for >= (upper - lower)) {
      allPoints = true;
    }
    if(allPoints) {
      // if all points are to be searched, no point in pseudo-randomly picking them
      for(int h=lower; h < upper; h++) {
        sid = (*ind)[h];
        m_kdtree->r_nearest_around_point(sid, 0, m_epsSquare, ne);
        if(ne.size() >= m_minPts) {
          growing_points.push_back(sid); // adds the point to the growing points vector
          // m_member[sid] = 1; // marks the point as a member of a cluster
        }
        alreadySeen.push_back(sid);
        ne.clear();
      }
    } else {

      for(int h=0; h < number_looking_for; h++) {
        //this loop initializes first n growing points randomly
        do {
          sid = (*ind)[(rand() % sch) + (sch * tid)]; // generates random index in the range of each thread's set of data points
          // `sid` should NOT have already been seen.
        } while (find(alreadySeen.begin(), alreadySeen.end(), sid) != alreadySeen.end());

        // Check if the new point is a core point, and if so add it to `growing_points`
        m_kdtree->r_nearest_around_point(sid, 0, m_epsSquare, ne);
        if(ne.size() >= m_minPts) {
          growing_points.push_back(sid); // adds the point to the growing points vector
          // m_member[sid] = 1; // marks the point as a member of a cluster
        }
        // no matter what, add `sid` to `alreadySeen` and clear out `ne` vector
        alreadySeen.push_back(sid);
        ne.clear();
      }
    }
    alreadySeen.clear();
  }

  // A cluster is determined by the root node. However many root nodes there are, that's how many clusters there are

  void run_dbscan_algo_uf(ClusteringAlgo& dbs) {     
    
    int tid, i, pid, j, k, npid, root, root1, root2, sid;
    int num_points = dbs.m_pts->m_i_num_points;
    vector <int> growing_points;
    // srand(time(NULL));

    // initialize some parameters
    dbs.m_clusters.clear();
  
    kdtree2_result_vector ne;
    kdtree2_result_vector ne2;
    // assign parent to itestf
    dbs.m_parents.resize(num_points, -1);

    dbs.m_member.resize(num_points, 0);
    dbs.m_corepoint.resize(num_points, 0);
    int sch, maxthreads = omp_get_max_threads(); //maximum number of threads available to do work 
    
    if(num_points % maxthreads == 0)
      sch = num_points/maxthreads;
    else
      sch = num_points/maxthreads + 1;
    // TODO bottom 5 lines deal with line 9 of pseudocode
    vector < vector <int > > merge; //initializes two dimensional vector. A vector that holds int vectors
    vector <int> init;
    merge.resize(maxthreads, init); // Merge has a size of maxthreads and each element is initialized with the vector init.
    
    for(tid = 0; tid < maxthreads; tid++)
      merge[tid].reserve(num_points); // ensures that each init vector can hold at least m_i_num_points
    
    vector < int > prID;
    prID.resize(num_points, -1); // initializes prID vector with a size of m_i_num_points and intializes each element with -1
    vector<int>* ind = dbs.m_kdtree->getIndex(); // Sets a vector that contains the index of all points
    double start = omp_get_wtime();
    // cout<< endl;

    #pragma omp parallel private(root, root1, root2, tid, ne, ne2, npid, i, j, pid, growing_points, sid) shared(sch, ind) //, prID) // creates threads
    {
      // private means that each thread will have its own private copy of variable in memory
      // shared means that all threads will share same copy of variable in memory
      int lower, upper; 
      tid = omp_get_thread_num(); // gets tid of each thread to identify it

      lower = sch * tid;  //The range of points that each thread has. Sch is number of points per thread
      upper = sch * (tid + 1);
      
      if(upper > num_points)
        upper = num_points;

      for(i = lower; i < upper; i++) {
        pid = (*ind)[i]; // CAN TRY RANDOMLY SELECTING
        dbs.m_parents[pid] = pid;  // Initialize parents to point to themselves  TODO line 5 of pseudocode
        prID[pid] = tid; //sets each element in prID to the thread number. 
      }
      
      //#pragma omp parallel for
      #pragma omp barrier

      dbs.getGrowingPoints(growing_points, sch, tid, lower, upper); // Each thread gets its seed points  TODO line 7 of pseudocode

      // if(tid == 0) cout << "before growing_points: " << growing_points.size() << endl;
      
      //cout << "made it to the barrier" << endl; 
      #pragma omp barrier // all threads will stop here until every thread has reached this point
      
      for(int i = 0; i < growing_points.size(); i++) { // Iterates through every growing point  TODO line 10 of pseudocode
        
        pid = growing_points[i];
        dbs.m_corepoint[pid] = 1;
        // dbs.m_member[pid] = 1; // mark as a member
        ne.clear();
        dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne); // gets nearest neighbors  TODO line 11 of pseudocode

        // get the root containing pid
        root = pid;

        for (j = 0; j < ne.size(); j++) { // TODO line 12 of pseudocode
          //this loop goes through all of nearest neighbors of a point
          npid= ne[j].idx; // gets index of ne[j] TODO line 13 of pseudocode
          if(npid == pid)
            continue;
          //cout << "prID: " << prID[npid] << " tid: " << tid << endl;
          if(prID[npid] != tid) { // this checks to see if the two points are in the same thread. If not, add them to merge
            // TODO line 28 of pseudocode
            merge[tid].push_back(pid); // The two points gets added to the end of the merge vector and will be merged later
            merge[tid].push_back(npid);
            continue; // goes to end of for loop
          }
          // get the root containing npid
          root1 = npid;
          root2 = root;

          if(dbs.m_corepoint[npid] == 1 || dbs.m_member[npid] == 0) { // TODO line 23 of pseudocode
            // mark as clustered and unionize it
            dbs.m_member[npid] = 1; // TODO line 25 of pseudocode
            omp_lock_t* fakeLocks;
            unionize_neighborhood(dbs, root, root1, root2, false, fakeLocks, tid); // TODO line 24 of pseudocode
          }

          ne2.clear();
          dbs.m_kdtree->r_nearest_around_point(npid, 0, dbs.m_epsSquare, ne2);
              
          if(ne2.size() >= dbs.m_minPts) { // TODO line 17 of pseudocode
            // Also check to see if it has already been added to the growing points vector
            if(find(growing_points.begin(), growing_points.end(), npid) == growing_points.end()) {
              //Now mark point as new growing point
              growing_points.push_back(npid); // TODO line 21 of pseudocode
            }
          } 
        } // end of for loop that goes through all nearest neighbors
      } // end of for loop that goes through each point 
      // if(tid == 0) cout << "after growing_points: " << growing_points.size() << endl;
    }


      
    int v1, v2, size;
    kdtree2_result_vector ne3;
    // merge the trees that have not been merged yet
    double stop = omp_get_wtime() ;
    cout  <<  endl;
    // cout << "Quality points: " << dbs.qualitypoints << endl;
    cout << "Local computation took " << stop - start << " seconds." << endl;
    //allocate and initiate locks
    omp_lock_t *nlocks;
    nlocks = (omp_lock_t *) malloc(num_points * sizeof(omp_lock_t)); // Initialize nlocks array for every data point
    //start = stop;
    start = omp_get_wtime();
    #pragma omp parallel for private(i) shared(nlocks)
      for(i = 0; i < num_points; i++) 
        omp_init_lock(&nlocks[i]); // initialize locks

    #pragma omp parallel for shared(maxthreads, merge, nlocks) private(i, v1, v2, root1, root2, size, tid, ne3) // allows for the spawned threads to split up the loop iterations
      for(tid = 0; tid < maxthreads; tid++) { // TODO line 35 of pseudocode
        size = merge[tid].size()/2;
        for(i = 0; i < size; i++) { // TODO line 36 of pseudocode
          v1 = merge[tid][2 * i]; // merge is a 2D vector. v1 = even numbered elements
          v2 = merge[tid][2 * i + 1]; // v2 = odd numbered elements
          int con = 0;

          ne3.clear();
          dbs.m_kdtree->r_nearest_around_point(v2, 0, dbs.m_epsSquare, ne3); // TODO commented this out. It's not in orriginal code
        
          if(ne3.size() >= dbs.m_minPts) // TODO line 37 of pseudocode
          // if(dbs.m_corepoint[v2] == 1) // TODO
            con = 1;
          else if(dbs.m_member[v2] == 0) { // TODO line 39 of pseudocode
            omp_set_lock(&nlocks[v2]);
            if(dbs.m_member[v2] == 0) { // if v2 is not a member yet   makes sure it's still not a member after it has the lock
              con = 1;
              dbs.m_member[v2] = 1; //marks it as being a member of a cluster  TODO line 40 of pseudocode
            }
            omp_unset_lock(&nlocks[v2]);
          }

          if(con == 1) {
            // lock based approach for merging
            root1 = v1;
            root2 = v2;
            // REMS algorithm with splicing compression techniques
            int fakeInt = -42;
            unionize_neighborhood(dbs, fakeInt, root1, root2, true, nlocks, 1); // TODO lines 38 & 41 of pseudocode
            // // REMS algorithm with splicing compression techniques
            // while (dbs.m_parents[root1] != dbs.m_parents[root2]) {
            //   if (dbs.m_parents[root1] < dbs.m_parents[root2]) {
            //     if(dbs.m_parents[root1] == root1) { // root1 is a root
            //       omp_set_lock(&nlocks[root1]);
            //       int p_set = false;
            //       if(dbs.m_parents[root1] == root1) { // if root1 is still a root
            //         dbs.m_parents[root1] = dbs.m_parents[root2];
            //         p_set = true;
            //       }
            //       omp_unset_lock(&nlocks[root1]);
            //       if (p_set) // merge successful
            //         break;
            //     }
            //     // splicing
            //     int z = dbs.m_parents[root1];
            //     dbs.m_parents[root1] = dbs.m_parents[root2];
            //     root1 = z;
            //   } else {
            //     if(dbs.m_parents[root2] == root2) { // root2 is a root
            //       omp_set_lock(&nlocks[root2]);
            //       int p_set = false;
            //       if(dbs.m_parents[root2] == root2) { // check if root2 is a root
            //         dbs.m_parents[root2] = dbs.m_parents[root1];
            //         p_set = true;
            //       }
            //       omp_unset_lock(&nlocks[root2]);
            //       if (p_set) // merge successful
            //         break;
            //     }
            //     //splicing
            //     int z = dbs.m_parents[root2];
            //     dbs.m_parents[root2] = dbs.m_parents[root1];
            //     root2 = z;
            //   } 
            // }
          }
        }
      }
          
    stop = omp_get_wtime();
    free(nlocks);
    
    cout << "Merging took " << stop - start << " seconds."<< endl;

    for(tid = 0; tid < maxthreads; tid++)
      merge[tid].clear();

    merge.clear();
    ne.clear();
    ne2.clear();
    //outputfile.close();
  }

  void run_dbscan_algo(ClusteringAlgo& dbs) {
    // classical DBSCAN algorithm (only sequential)
    int i, pid, j, k, npid;
    int cid = 1; // cluster id
    int num_points = dbs.m_pts->m_i_num_points;
    vector <int> c;
    c.reserve(num_points);
    // initialize some parameters
    dbs.m_noise.resize(num_points, false);
    dbs.m_visited.resize(num_points, false);   
    dbs.m_pid_to_cid.resize(num_points, 0);
    dbs.m_clusters.clear();
    // get the neighbor of the first point and print them
    //cout << "DBSCAN ALGORITHMS============================" << endl;
    kdtree2_result_vector ne;
    kdtree2_result_vector ne2;
    //kdtree2_result_vector ne3;
    ne.reserve(num_points);
    ne2.reserve(num_points);
    vector<int>* ind = dbs.m_kdtree->getIndex();
    double start = omp_get_wtime() ;    
    for(i = 0; i < num_points; i++) {
      pid = (*ind)[i];
      if (!dbs.m_visited[pid]) {
        dbs.m_visited[pid] = true;
        ne.clear();
        dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne);
        if(ne.size() < dbs.m_minPts)
          dbs.m_noise[pid] = true;
        else {
          // start a new cluster
          c.clear();
          c.push_back(pid);
          dbs.m_pid_to_cid[pid] = cid;
          // traverse the neighbors
          for (j = 0; j < ne.size(); j++) {
            npid= ne[j].idx;
            // not already visited
            if(!dbs.m_visited[npid]) { // if point hasn't been visited then
              dbs.m_visited[npid] = true; // label point as visited
              // go to neighbors
              ne2.clear();
              dbs.m_kdtree->r_nearest_around_point(npid, 0, dbs.m_epsSquare, ne2); // gets neighbors of point
              // enough support
              if (ne2.size() >= dbs.m_minPts) { // If it has enough support, then it is a core point and Union (x,x')
                // join
                for(k = 0; k < ne2.size(); k++) //Union function
                  ne.push_back(ne2[k]); // adds element to the end of the vector
              }
            }
            // not already assigned to a cluster
            if (!dbs.m_pid_to_cid[npid]) {
              c.push_back(npid);
              dbs.m_pid_to_cid[npid]=cid;
              dbs.m_noise[npid] = false;
            }
          }
          dbs.m_clusters.push_back(c);
          cid++;
        }
      }
    }
    double stop = omp_get_wtime();
    cout << "Local computation took " << stop - start << " seconds." << endl;
    cout << "No merging stage in classical DBSCAN"<< endl;
    ind = NULL;
    ne.clear();
    ne2.clear();
  }

  /*
    called in `run_dbscan_algo_uf` function.
    Only called if the amount of points found are equal to, or greater than, the minimum points
    needed to make a cluster. 
    The points are joined via a union opperation using the REMS algorithm.

    ClusteringAlgo& dbs - DBScan object. The parent for EVERYTHING.
    int root - initially set to pid
    int root1 - used to find the actual "root" node
    int root2 - used to find the actual "root" node
    bool locks - flag to determine if locks are to be used or not
    omp_lock_t* actualLock - the actual lock to use if `locks` is set to TRUE. It's a useless variable if `locks` is FALSE.
  */
  void unionize_neighborhood(ClusteringAlgo& dbs, int& root, int root1, int root2, bool locks, omp_lock_t* actualLock, int tid) {
    
    if(locks){
      // REMS algorithm with splicing compression techniques
      while (dbs.m_parents[root1] != dbs.m_parents[root2]) {
        if (dbs.m_parents[root1] < dbs.m_parents[root2]) {
          if(dbs.m_parents[root1] == root1) { // root1 is a root
            omp_set_lock(&actualLock[root1]);
            int p_set = false;
            if(dbs.m_parents[root1] == root1) { // if root1 is still a root
              dbs.m_parents[root1] = dbs.m_parents[root2];
              p_set = true;
            }
            omp_unset_lock(&actualLock[root1]);
            if (p_set) // merge successful
              break;
          }
          // splicing
          int z = dbs.m_parents[root1];
          dbs.m_parents[root1] = dbs.m_parents[root2];
          root1 = z;
        } else {
          if(dbs.m_parents[root2] == root2) { // root2 is a root
            omp_set_lock(&actualLock[root2]);
            int p_set = false;
            if(dbs.m_parents[root2] == root2) { // check if root2 is a root
              dbs.m_parents[root2] = dbs.m_parents[root1];
              p_set = true;
            }
            omp_unset_lock(&actualLock[root2]);
            if (p_set) // merge successful
              break;
          }
          //splicing
          int z = dbs.m_parents[root2];
          dbs.m_parents[root2] = dbs.m_parents[root1];
          root2 = z;
        } 
      }
    } else{
      // REMS algorithm to merge the trees
      while(dbs.m_parents[root1] != dbs.m_parents[root2]) {
        if(dbs.m_parents[root1] < dbs.m_parents[root2]) {
          if(dbs.m_parents[root1] == root1) {
            dbs.m_parents[root1] = dbs.m_parents[root2];
            root = dbs.m_parents[root2];
            break;
          }
          // splicing
          int z = dbs.m_parents[root1];
          dbs.m_parents[root1] = dbs.m_parents[root2];
          root1 = z;
        } else {
          if(dbs.m_parents[root2] == root2) {
            dbs.m_parents[root2] = dbs.m_parents[root1];
            root = dbs.m_parents[root1];
            break;
          }
          // splicing
          int z = dbs.m_parents[root2];
          dbs.m_parents[root2] = dbs.m_parents[root1];         
          root2 = z;
        }
      }
    }
  }
};
