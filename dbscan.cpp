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
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      //for (j = 0; j < m_pts->m_i_dims; j++)
      //  o << " " << m_pts->m_points[i][j];
      id = m_pid_to_cid[i];
      o << i << " " << id << endl;
    }
    int sum_points = 0;
    int noise = 0;
    for(i = 0; i < m_clusters.size(); i++) {
      sum_points += m_clusters[i].size();
      //cout << i << "(" << m_clusters[i].size() << ") ";
    }
    for (i = 0; i < m_pts->m_i_num_points; i++) {
      if(m_noise[i])
        noise++;
    } 
    cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
    cout << "Number of clusters: " << m_clusters.size() << endl;
  }

  void ClusteringAlgo::writeClusters_uf(ostream& o) {
    //changed return type to int. Was void. Used for testing purposes
    // writing point id and cluster id pairs per line, noise has cluster id 0 
    vector <int> clusters;
    clusters.resize(m_pts->m_i_num_points, 0);
    int i, j, sum_points = 0, noise = 0, root, rootcount = 0, tmp;
    for(i = 0; i < m_pts->m_i_num_points; i++) {
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
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      if(clusters[i] == 1) {
        // vertex i is a noise
        clusters[i] = 0;
        noise++;
      }
      else if(clusters[i] >= m_minPts) {
        // This conditional statement determines what is and is not considered a cluster
        // If it's greater than this number than it will be counted as a cluster
        count++;
        sum_points += clusters[i];
        clusters[i] = count;

        cout << clusters[i] << ": " << i << endl;
      }
    }
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      o << i << " " << clusters[m_parents[i]] << endl; //This was the only one not commented out
    }
    cout << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
    cout << "Number of clusters: " << count << endl;
    o << "Total points " << noise + sum_points << " pt_in_cls " << sum_points << " noise " << noise << endl;
    o << "Number of clusters: " << count << endl;
    clusters.clear();
  }

  // A cluster is determined by the root node. However many root nodes there are, that's how many clusters there are

  void run_dbscan_algo_uf(ClusteringAlgo& dbs) {     
    int tid, i, pid, j, k, npid, root, root1, root2, sid, h, qualitypoints=0, test=0;
    srand(time(NULL));

        // initialize some parameters
    dbs.m_clusters.clear();
    
    vector <int> growing_points;
    kdtree2_result_vector ne;
    kdtree2_result_vector ne2;
    // assign parent to itestf
    dbs.m_parents.resize(dbs.m_pts->m_i_num_points, -1);

    dbs.m_member.resize(dbs.m_pts->m_i_num_points, 0);
    dbs.m_corepoint.resize(dbs.m_pts->m_i_num_points, 0);
    int sch, maxthreads = omp_get_max_threads(); //maximum number of threads available to do work 
    if(dbs.m_pts->m_i_num_points % maxthreads == 0)
      sch = dbs.m_pts->m_i_num_points/maxthreads;
    else
      sch = dbs.m_pts->m_i_num_points/maxthreads + 1;
    
    vector < vector <int > > merge; //initializes two dimensional vector. A vector that holds int vectors
    vector <int> init;
    merge.resize(maxthreads, init); // Merge has a size of maxthreads and each element is initialized with the vector init.
    for(tid = 0; tid < maxthreads; tid++)
      merge[tid].reserve(dbs.m_pts->m_i_num_points); // ensures that each init vector can hold at least m_i_num_points
    
    vector < int > prID;
    prID.resize(dbs.m_pts->m_i_num_points, -1); // initializes prID vector with a size of m_i_num_points and intializes each element with -1
    vector<int>* ind = dbs.m_kdtree->getIndex(); // Sets a vector that contains the index of all points
    double start = omp_get_wtime();
    cout<< endl;
    #pragma omp parallel private(root, root1, root2, tid, ne, ne2, npid, i, j, pid, growing_points, sid) shared(sch, ind, h, test, qualitypoints) { //, prID) // creates threads
      // private means that each thread will have its own private copy of variable in memory
      // shared means that all threads will share same copy of variable in memory
      int lower, upper;
      tid = omp_get_thread_num(); // gets tid of each thread to identify it
      lower = sch * tid;  //The range of points that each thread has. Sch is number of points per thread
      upper = sch * (tid + 1);
      if(upper > dbs.m_pts->m_i_num_points)
      upper = dbs.m_pts->m_i_num_points;
      for(i = lower; i < upper; i++) {
        pid = (*ind)[i]; // CAN TRY RANDOMLY SELECTING
        //cout << "pid is: " << pid << " it belongs to thread " << tid << endl;
        dbs.m_parents[pid] = pid;  // Initialize parents to point to themselves
        prID[pid] = tid; //sets each element in prID to the thread number. 
      }

      //#pragma omp parallel for
      #pragma omp barrier
      #pragma omp for

      for(h=0; h < (dbs.m_seeds); h++) {
        //this loop initializes first n growing points randomly

        do {
          sid = (*ind)[(rand() % sch) + (sch * tid)]; // generates random index in the range of each thread's set of data points  
        } while(find(growing_points.begin(), growing_points.end(), sid) != growing_points.end());
        //repeats the do while loop if it is not a core point or the point has already been selected as a growing point
        dbs.m_kdtree->r_nearest_around_point(sid, 0, dbs.m_epsSquare, ne);
        if(ne.size() >= dbs.m_minPts) {
          qualitypoints++;
          growing_points.push_back(sid); // adds the point to the growing points vector
          dbs.m_member[sid] = 1; // marks the point as a member of a cluster
        } 
        // tid is the thread number, and sid is the index of the point
        
        //dbs.m_corepoint
      }
      
      //cout << "made it to the barrier" << endl; 
      #pragma omp barrier // all threads will stop here until every thread has reached this point

      for(int i = 0; i < growing_points.size(); i++) { // Iterates through every growing point
        dbs.m_corepoint[i] = 1;
        //test += 1;
        pid = growing_points[i];
        ne.clear();
        dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne); // gets nearest neighbors
        test += 1;
                
        dbs.m_member[pid] = 1; // mark as a member
        // get the root containing pid
        root = pid;
        for (j = 0; j < ne.size(); j++) {
          //this loop goes through all of nearest neighbors of a point
          npid= ne[j].idx; // gets index of ne[j]
          if(npid == pid)
            continue;
          //cout << "prID: " << prID[npid] << " tid: " << tid << endl;
          if(prID[npid] != tid) { // this checks to see if the two points are in the same thread. If not, add them to merge
            merge[tid].push_back(pid); // The two points gets added to the end of the merge vector and will be merged later
            merge[tid].push_back(npid);
            continue; // goes to end of for loop
          }
          // get the root containing npid
          root1 = npid;
          root2 = root;

          dbs.m_kdtree->r_nearest_around_point(npid, 0, dbs.m_epsSquare, ne2);
          test += 1;
              
          if(ne2.size() >= dbs.m_minPts) {
            // REMS algorithm to merge the trees
            while(dbs.m_parents[root1] != dbs.m_parents[root2]) { // while the parents aren't equal
              if(dbs.m_parents[root1] < dbs.m_parents[root2]) { //if root1's value is less than root2's value (smaller nodes point to larger nodes)
                if(dbs.m_parents[root1] == root1) { //if point is the root
                  dbs.m_parents[root1] = dbs.m_parents[root2]; //Sets the parent of root1 to be the parent of root2
                  root = dbs.m_parents[root2]; //sets root to be the parent of root2
                  break; // root has been found. Break from the loop
                }
                // splicing
                // if not at the root, then set root 1 equal to the parent. Advance up the tree
                int z = dbs.m_parents[root1]; //creates temporary variable
                dbs.m_parents[root1] = dbs.m_parents[root2]; //makes root1's subtree a sibling of root2
                root1 = z; //sets root1 to be the parent of root 1. Advance up the tree in order to find the root
              }
              else {
                //root2 < root 1
                if(dbs.m_parents[root2] == root2) { //if point is the root else
                  dbs.m_parents[root2] = dbs.m_parents[root1]; //Sets the parent of root2 to be the parent of root1
                  root = dbs.m_parents[root1]; //sets root to be the parent of root1
                  break; //root has been found. Break from loop
                }
                // splicing
                int z = dbs.m_parents[root2]; //creates temporary variable
                dbs.m_parents[root2] = dbs.m_parents[root1]; // makes root2's subtree a sibling of root1      
                root2 = z;//sets root2 to be the parent of root 2. Advance up the tree in order to find the root
              }
            } // end of while loop that checks to see if the parents are equal

            if(dbs.m_member[npid] == 0) {
              //check to see if the point belongs to a cluster and if not, add to growing_points and mark as clustered
              // Also check to see if it has already been added to the growing points vector
              if(find(growing_points.begin(), growing_points.end(), npid) == growing_points.end()) {
                dbs.m_member[npid] == 1;
                //Now mark point as new growing point
                growing_points.push_back(npid);
              }
            }

          }else if(ne2.size() < dbs.m_minPts) {
            if(dbs.m_member[npid] == 0) {
              //If point is not a core point but it doesn't belong to any cluster yet, mark as clustered and union it
              dbs.m_member[npid] == 1;
              while(dbs.m_parents[root1] != dbs.m_parents[root2]) { // while the parents aren't equal
                if(dbs.m_parents[root1] < dbs.m_parents[root2]) { //if root1's value is less than root2's value (smaller nodes point to larger nodes)
                  if(dbs.m_parents[root1] == root1) { //if point is the root
                    dbs.m_parents[root1] = dbs.m_parents[root2]; //Sets the parent of root1 to be the parent of root2
                    root = dbs.m_parents[root2]; //sets root to be the parent of root2
                    break; // root has been found. Break from the loop
                  }
                  // splicing
                  // if not at the root, then set root 1 equal to the parent. Advance up the tree
                  int z = dbs.m_parents[root1]; //creates temporary variable
                  dbs.m_parents[root1] = dbs.m_parents[root2]; //makes root1's subtree a sibling of root2
                  root1 = z; //sets root1 to be the parent of root 1. Advance up the tree in order to find the root
                }else {
                  //root2 < root 1
                  if(dbs.m_parents[root2] == root2) { //if point is the root
                    dbs.m_parents[root2] = dbs.m_parents[root1]; //Sets the parent of root2 to be the parent of root1
                    root = dbs.m_parents[root1]; //sets root to be the parent of root1
                    break; //root has been found. Break from loop
                  }
                  // splicing
                  int z = dbs.m_parents[root2]; //creates temporary variable
                  dbs.m_parents[root2] = dbs.m_parents[root1]; // makes root2's subtree a sibling of root1      
                  root2 = z;//sets root2 to be the parent of root 2. Advance up the tree in order to find the root
                }
              }
            }
          } // end of else if statement checking to see if it hasn't been clustered yet
        } // end of for loop that goes through all nearest neighbors
        //} // end of if it's a core point loop
      } // end of for loop that goes through each point
    }
      
    int v1, v2, size, test2 = 0;
    kdtree2_result_vector ne3;
    // merge the trees that have not been merged yet
    double stop = omp_get_wtime() ;
    cout  <<  endl;
    cout << "Quality points: " << qualitypoints << endl;
    cout << "Local computation took " << stop - start << " seconds." << endl;
    //allocate and initiate locks
    omp_lock_t *nlocks;
    nlocks = (omp_lock_t *) malloc(dbs.m_pts->m_i_num_points*sizeof(omp_lock_t)); // Initialize nlocks array for every data point
    //start = stop;
    start = omp_get_wtime();
    #pragma omp parallel for private(i) shared(nlocks)
      for(i = 0; i < dbs.m_pts->m_i_num_points; i++) 
        omp_init_lock(&nlocks[i]); // initialize locks
    #pragma omp parallel for shared(maxthreads, merge, nlocks, test2) private(i, v1, v2, root1, root2, size, tid, ne3) // allows for the spawned threads to split up the loop iterations
      for(tid = 0; tid < maxthreads; tid++) {
      size = merge[tid].size()/2;
      for(i = 0; i < size; i++) {
        v1 = merge[tid][2 * i]; // merge is a 2D vector. v1 = even numbered elements
        v2 = merge[tid][2 * i + 1]; // v2 = odd numbered elements
        int con = 0;

        dbs.m_kdtree->r_nearest_around_point(v2, 0, dbs.m_epsSquare, ne3);
        test2 += 1;
        
        if(ne3.size() >= dbs.m_minPts)
          con = 1;
        else if(dbs.m_member[v2] == 0) {
          omp_set_lock(&nlocks[v2]);
          if(dbs.m_member[v2] == 0) { // if v2 is not a member yet   makes sure it's still not a member after it has the lock
            con = 1;
            dbs.m_member[v2] = 1; //marks it as being a member of a cluster
          }
          omp_unset_lock(&nlocks[v2]);
        }
        if(con == 1) {
          // lock based approach for merging
          root1 = v1;
          root2 = v2;
          // REMS algorithm with splicing compression techniques
          while (dbs.m_parents[root1] != dbs.m_parents[root2]) {
            if (dbs.m_parents[root1] < dbs.m_parents[root2]) {
              if(dbs.m_parents[root1] == root1) { // root1 is a root
                omp_set_lock(&nlocks[root1]);
                int p_set = false;
                if(dbs.m_parents[root1] == root1) { // if root1 is still a root after lock has been set
                  dbs.m_parents[root1] = dbs.m_parents[root2];
                  p_set = true;
                }
                omp_unset_lock(&nlocks[root1]);
                if (p_set) // merge successful
                  break;
              }
              // splicing
              int z = dbs.m_parents[root1];
              dbs.m_parents[root1] = dbs.m_parents[root2];
              root1 = z;
              // root1 = dbs.m_parents[root1];
            } else {
              if(dbs.m_parents[root2] == root2) { // root2 is a root
                omp_set_lock(&nlocks[root2]);
                int p_set = false;
                if(dbs.m_parents[root2] == root2) { // check if root2 is a root
                  dbs.m_parents[root2] = dbs.m_parents[root1];
                  p_set = true;
                }
                omp_unset_lock(&nlocks[root2]);
                if (p_set) // merge successful
                  break;
              }
              //splicing
              int z = dbs.m_parents[root2];
              dbs.m_parents[root2] = dbs.m_parents[root1];
              root2 = z;
              //root2 = dbs.m_parents[root2];
            } 
          }
        }
      }
    }
          
    stop = omp_get_wtime();
    free(nlocks);
    //cout << "Merging neighborhood searches: " << test2 << endl;
    cout << "Merging took " << stop - start << " seconds."<< endl;
    //cout << "Number of Neighborhood searches: " << test << endl;
    //  ofstream output;
    //  output.open("bestrun.txt")
    //  current_time = start_total - stop_total;
    //  output << "Number of seeds is " << n << endl;
    //  output << "Total time took " << current_time << " seconds" << endl << endl;
    //  if(current_time < min_time){
    //    min_time = current_time;
    //    best_n = n;
    //  }
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
    vector <int> c;
    c.reserve(dbs.m_pts->m_i_num_points);
    // initialize some parameters
    dbs.m_noise.resize(dbs.m_pts->m_i_num_points, false);
    dbs.m_visited.resize(dbs.m_pts->m_i_num_points, false);   
    dbs.m_pid_to_cid.resize(dbs.m_pts->m_i_num_points, 0);
    dbs.m_clusters.clear();
    // get the neighbor of the first point and print them
    //cout << "DBSCAN ALGORITHMS============================" << endl;
    kdtree2_result_vector ne;
    kdtree2_result_vector ne2;
    //kdtree2_result_vector ne3;
    ne.reserve(dbs.m_pts->m_i_num_points);
    ne2.reserve(dbs.m_pts->m_i_num_points);
    vector<int>* ind = dbs.m_kdtree->getIndex();
    double start = omp_get_wtime() ;    
    for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
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
};
