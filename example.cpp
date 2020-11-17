#include <iostream>
#include "Delaunay.hpp"
#include <vector>
#include <stdio.h>

#include <string>
#include <sstream>
#include <fstream>

#include <mpi.h>
// #include "/usr/local/opt/metis/include/metis.h"
#include <metis.h>
// #include "metis/kmetis.c"

// To compile and run
// mpic++ -std=c++11 example.cpp -lmetis
// mpirun -np 4 ./a.out

using namespace std;

// Test Data
ifstream inFile ("input.txt");

// Vertices XS - x coordinates & YS - y coordinates
vector<int> XS;
vector<int> YS;

int main(int argc, char** argv)
{
    int x;
    inFile >> x;

    for(int i=0;i<x;i++) {
        int x1,x2;
        inFile >> x1 >> x2;
        XS.push_back(x1);
        YS.push_back(x2);
    }

  // Initialize the Delaunay Triangulation Class
  DelaunayTriangulation DT(99, 99);

  // Adding points
  cout << "Adding points:" << endl;
  for (int i = 0; i < XS.size(); i++)
  {
    std::cout << "i: " << i << "/" << XS.size() << ",  ";
    DT.AddPoint(Point(XS[i], YS[i]));
  }

  cout << endl;
  
  DT.print();

  cout << endl;

  // get Triangles
  // DT.getTriangle();

  //hardcodng the vertices in 2d vector format and also the 4 extreme vertices
  vector<vector<int>> vertices
    {
      {7,4,0},
      {8,4,5},
      {10,3,2},
      {10,2,1},
      {13,5,4},
      {13,4,7},
      {13,7,6},
      {13,6,5},
      {14,8,10},
      {15,0,4},
      {15,4,8},
      {15,8,14},
      {15,14,10},
      {15,10,1},
      {15,1,11},
      {15,11,0}
    };

    // checking for the predicates
    // add metis
    // add mpi codes and compilation tags

    // cout << "MPI" << endl;

    ////////////////////////////////////////////////////////////////
    // mpi process
    // message between and parallelization between threads
    ////////////////////////////////////////////////////////////////
    // To run: mpic++ hello_world_mpi.cpp -o hello_world_mpi
    // To execute: mpirun -np 2 ./hello_world_mpi
    ////////////////////////////////////////////////////////////////
    int process_Rank, size_Of_Cluster, message_Item;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    
    // cout << process_Rank << endl;
    if(process_Rank == 0){
        message_Item = 42;
        MPI_Send(&message_Item, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
        printf("Message Sent from Process: %d %d\n", process_Rank, message_Item);
    }

    else if(process_Rank == 1){
        MPI_Recv(&message_Item, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Message Received to Process: %d %d\n", process_Rank, message_Item);
    }

    MPI_Finalize();

    ////////////////////////////////////////////////////////////////
    // metis example
    // To partition graphs or/and data.
    ////////////////////////////////////////////////////////////////
    // To compile: g++ -std=c++11 metis_example.cpp -lmetis
    // To run: ./metis_example
    ///////////////////////////////////////////////////////////////

    idx_t nVertices = 7;
    idx_t nElements = 5;
    idx_t nParts = 2;

    idx_t objval;
    std::vector<idx_t> epart(nElements, 0);
    std::vector<idx_t> npart(nVertices, 0);

    // Indexes of starting points in adjacent array
    std::vector<idx_t> eptr = {5,0,3,6,9,12};

    // Adjacent vertices in consecutive index order
    std::vector<idx_t> eind = {0,1,2,0,6,2,5,4,6,2,4,6,4,2,3};

    // if all weights are equal then can be set to NULL
    // std::vector<idx_t> vwgt(nVertices * nWeights, 0);
    
    // int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj.data(), adjncy.data(),
		// 		       NULL, NULL, NULL, &nParts, NULL,
    //    				  NULL, NULL, &objval, part.data());

    // std::cout << ret << std::endl;
    
    // for(unsigned part_i = 0; part_i < part.size(); part_i++){
    //   cout << process_Rank << ":";
	  //   std::cout << part_i << " " << part[part_i] << std::endl;
    // }

    // partmeshNodal

    int ret2 = METIS_PartMeshNodal( 
      &nElements, &nVertices, eptr.data(), eind.data(), NULL, NULL,
      &nParts, NULL, NULL, &objval, epart.data(), npart.data());

    for(unsigned part_i = 0; part_i < epart.size(); part_i++){
      cout << "EPart" << process_Rank << ":";
	    std::cout << part_i << " " << epart[part_i] << std::endl;
    }

    for(unsigned part_i = 0; part_i < npart.size(); part_i++){
      cout << process_Rank << ":";
	    std::cout << part_i << " " << npart[part_i] << std::endl;
    }

  return 0;
}