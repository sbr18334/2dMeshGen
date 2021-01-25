#include <iostream>
#include "Delaunay.hpp"
#include <vector>
#include <stdio.h>
#include <math.h>

#include <string>
#include <sstream>
#include <fstream>
#include <unordered_set>

#include <mpi.h>
#include <metis.h>

// To compile and run
// mpic++ -std=c++11 example.cpp -lmetis
// mpirun -np 4 ./a.out

using namespace std;

struct Pnt {
  float x, y;
};

// Test Data
ifstream inFile ("input.txt");
vector<vector<float>> points
  {
    {1,2},{2,3},{-4,6},{-2,3},{4,3},{-10,3},{4,8},
    {5,-3},{-1,7},{0,5},{4,4},{5,-1}
  };

// Vertices XS - x coordinates & YS - y coordinates
vector<int> XS;
vector<int> YS;

float threshold1 = 2;
float threshold2 = 5;

float triangleArea(Pnt p1, Pnt p2, Pnt p3) {         //find area of triangle formed by p1, p2 and p3
   return abs((p1.x*(p2.y-p3.y) + p2.x*(p3.y-p1.y)+ p3.x*(p1.y-p2.y))/2.0);
}

bool inCircle(Pnt p1, Pnt p2, Pnt p3, Pnt p) {     //check whether p is inside or outside
   float area = triangleArea (p1, p2, p3);          //area of triangle ABC
   float area1 = triangleArea (p, p2, p3);         //area of PBC
   float area2 = triangleArea (p1, p, p3);         //area of APC
   float area3 = triangleArea (p1, p2, p);        //area of ABP

   return (area == area1 + area2 + area3);        //when three triangles are forming the whole triangle
}

float* circumcenter(float Px, float Py, float Qx, float Qy, float Rx, float Ry
                    ,float x, float y, float z)
{
  float* arr = new float[2];
  cout << x << " " << y << " " << z << endl;
  float angleA = acos((-x*x+y*y+z*z)/(2*y*z));
  float angleB = acos((x*x-y*y+z*z)/(2*x*z));
  float angleC = acos((x*x+y*y-z*z)/(2*x*y));

  cout << angleA << " " << angleB << " " << angleC << endl;

  arr[0] = (Px * sin(2*angleA) + Qx*sin(2*angleB) + Rx*sin(2*angleC))/
           (sin(2*angleA) + sin(2*angleB) + sin(2*angleC));
  arr[1] = (Py * sin(2*angleA) + Qy*sin(2*angleB) + Ry*sin(2*angleC))/
           (sin(2*angleA) + sin(2*angleB) + sin(2*angleC));
           
  cout << arr[0] << " " << arr[1] << endl;
  return arr;
}

float distance(int x1, int y1, int x2, int y2) 
{
    return sqrt(pow(x2 - x1, 2) +  pow(y2 - y1, 2)); 
} 

float triangle_area(float a, float b, float c)
{
  float p = (a+b+c)/2;
  return sqrt(p*(p-a)*(p-b)*(p-c));
}

float midPoint(float a, float b)
{
  return (a+b)/2;
}

int* checkCommonEdge(int a1, int a2, int a3, int b1, int b2, int b3)
{
  int* arr = new int[2];
  if((a1==b1 && a2 == b2) || (a1 == b2 && a2 == b1) || 
     (a1==b2 && a2 == b3) || (a1 == b3 && a2 == b2) ||
     (a1==b1 && a2 == b3) || (a1 == b3 && a2 == b1)) {
      arr[0] = a1; arr[1] = a2;
      return arr;
  }
  if((a3==b1 && a2 == b2) || (a3 == b2 && a2 == b1) || 
     (a3==b2 && a2 == b3) || (a3 == b3 && a2 == b2) ||
     (a3==b1 && a2 == b3) || (a3 == b3 && a2 == b1)) {
      arr[0] = a2; arr[1] = a3;
      return arr;
  }
  if((a1==b1 && a3 == b2) || (a1 == b2 && a3 == b1) || 
     (a1==b2 && a3 == b3) || (a1 == b3 && a3 == b2) ||
     (a1==b1 && a3 == b3) || (a1 == b3 && a3 == b1)) {
      arr[0] = a1; arr[1] = a3;
      return arr;
  }
  else {
    arr[0] = -999;arr[1] = -999;
    return arr;
  }
}

int main(int argc, char** argv)
{

  //Number of processes
  int nProcesses = 4;
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

  DT.print();

  // get Triangles
  // DT.getTriangle();

  //hardcodng the vertices in 2d vector format and also the 4 extreme vertices
  vector<vector<int>> vertices
    {
      {7,4,0}, {8,4,5}, {10,3,2}, {10,2,1}, {13,5,4},
      {13,4,7}, {13,7,6}, {13,6,5}, {14,8,10}, {15,0,4},
      {15,4,8}, {15,8,14}, {15,14,10}, {15,10,1},
      {15,1,11}, {15,11,0}
    };

  // checking for the predicates
  // add metis
  // add mpi codes and compilation tags

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

  ////////////////////////////////////////////////////////////////
  // metis example
  // To partition graphs or/and data.
  ////////////////////////////////////////////////////////////////
  // To compile: g++ -std=c++11 metis_example.cpp -lmetis
  // To run: ./metis_example
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // partmeshNodal
  ///////////////////////////////////////////////////////////////
  idx_t nVertices = 8;
  idx_t nElements = 8;
  idx_t nParts = nProcesses;

  idx_t objval;
  std::vector<idx_t> epart(nElements, 0);
  std::vector<idx_t> npart(nVertices, 0);
  std::vector<idx_t> recpart(nVertices, 0);

  std::vector<idx_t> eptr = {8,0,3,6,9,12,15,18,21};

  std::vector<idx_t> eind = {0,4,5,4,5,6,4,6,7,4,7,3,0,4,3,0,1,3,1,2,3,2,3,7};

  int ret2 = METIS_PartMeshNodal( 
    &nElements, &nVertices, eptr.data(), eind.data(), NULL, NULL,
    &nParts, NULL, NULL, &objval, epart.data(), npart.data());
    
  // if(process_Rank == 0){
  //   cout << endl <<"----------Elements Partition-------" << endl;
  //   for(unsigned part_i = 0; part_i < epart.size(); part_i++){
  //     std::cout << "Element " << part_i << " Allotted to the P" << epart[part_i] << std::endl;
  //   }

  //   cout << endl << "----------Nodes Partition----------" << endl;
  //   for(unsigned part_i = 0; part_i < npart.size(); part_i++){
  //     std::cout << "Node " << part_i << " Allotted to P" << npart[part_i] << std::endl;
  //   }
  // }

  int num_of_DONE = 0;
  vector<int> partNodes;
  vector<int> partElements;

  for(unsigned part_i = 0; part_i < npart.size(); part_i++){
    if(npart[part_i] == process_Rank) {
      partNodes.push_back(part_i);
    }
    if(epart[part_i] == process_Rank) {
      partElements.push_back(part_i);
    }
  }

  cout << partNodes.size();
  for(int i=0;i<partElements.size()/3;i++) {
    //vertices of that triangle
    int Px = points[eind[i*3]][0];
    int Py = points[eind[i*3]][1];
    int Qx = points[eind[i*3+1]][0];
    int Qy = points[eind[i*3+1]][1];
    int Rx = points[eind[i*3+2]][0];
    int Ry = points[eind[i*3+2]][1];

    float x = distance(Qx,Qy,Rx,Ry);
    float y = distance(Rx,Ry,Px,Py);
    float z = distance(Px,Py,Qx,Qy);

    float shortest = x < y ? (x < z ? x : z) : (y < z ? y : z);
    float area = triangle_area(x,y,z);
    float circumRadius = x*y*z/(4*area);

    ///////////////////////////////////////////////////////////////
    // obtaining bad triangles
    ///////////////////////////////////////////////////////////////
    if(circumRadius/shortest > threshold1 || area > threshold2) {
      //bad triangle, find the circumcenter
      cout << "IN" << endl;
      cout << "Shortest" << shortest << endl;
      cout << "Circumradius" << circumRadius << endl;
      cout << "Circum/shortest" << circumRadius/shortest << endl; 
      cout << "Area" << area << endl;

      //circumcenter co-ordinates
      float* ptr = circumcenter(Px,Py,Qx,Qy,Rx,Ry,x,y,z);

      cout << Px << " " << Py << " " << Qx << " " << Qy 
            << " " << Rx << " " << Ry << endl;

      cout << ptr[0] << " " << ptr[1] << endl;
      //got circumcenter

      //commonedges with processid
      cout << "eind.length/3" << eind.size()/3 << endl;
      //total no. of triangles
      int triangleCount = eind.size()/3;
      //complexity verify
      std::vector<idx_t> edgeList;
          //using trashEdgeList add the triangles
      for(int i=0;i<triangleCount;i++) {
        if(epart[i] == process_Rank){
          for(int j=0;j<triangleCount;j++) {
            //check whether the triangles have a common edge
            if(epart[j] != process_Rank){
              cout << eind[3*i] << eind[3*i+1] << eind[3*i+2] <<
                              eind[3*j] << eind[3*j+1] << eind[3*j+2] << endl;
              int* ptr2 = checkCommonEdge(eind[3*i],eind[3*i+1],eind[3*i+2],
                              eind[3*j],eind[3*j+1],eind[3*j+2]);
              cout << "cEdge[0]:" << ptr2[0] << endl;
              cout << "cEdge[1]:" << ptr2[1] << endl;
              if(!(ptr2[0] == -999 && ptr2[1] == -999)){
                edgeList.push_back(ptr2[0]);
                edgeList.push_back(ptr2[1]);
                edgeList.push_back(process_Rank);
              }
            }
          }
        }
      }

      for(int i=0;i<edgeList.size();i++) {
        cout << i << ":" << edgeList[i] << endl;
      }


      // Local cavity
      Pnt pd = {ptr[0], ptr[1]};

      //find shared edges of this process with other processes
      //take common edges and check whether the circumcenter encroaches the segment
      for(int i=0;i<edgeList.size()/3;i++) {
        float centerX = midPoint(points[edgeList[3*i]][0],points[edgeList[3*i+1]][0]);
        float centerY = midPoint(points[edgeList[3*i]][1],points[edgeList[3*i+1]][1]);
        float radius = distance(centerX,centerY,points[edgeList[3*i]][0],
                                points[edgeList[3*i]][1]);
        if(pow((ptr[0]-centerX),2)+pow((ptr[0]-centerX),2)-pow(radius,2) < 0) {
          cout << "Encroaches" << endl;
          // MPI
          // sending mpi message to the respective process
          int data[2];
          data[0] = edgeList[3*i]; data[1] = edgeList[3*i+1];
          // MPI_Send(data,2,MPI_INT,edgeList[3*i+2],NULL,NULL);

          pd = {centerX, centerY};
          break;

        } else { // in bad triangle // common edges
          // for testing purpose code(+9 lines)
          int data[2];
          data[0] = 1; data[1] = 2;
          cout << "Sending MPI Message to designated processes" << endl;
          cout << "-------------------------------------------" << endl;
          MPI_Send(data,2,MPI_INT,1,0,MPI_COMM_WORLD);
          cout << "-------------------------------------------" << endl;
          cout << endl;
          // end of testing purpose code(+9 lines)
          cout << "Doesnot Encroaches" << endl;

        } // end of doesnot encroaches loop

      } // end of each common edge check 

      // Local cavity
      // vector to keep track of all the triangles to delete
      vector<int> trashTriangles;
      vector<int> trashIndices;
      // do not use eind here
      for(int i=0;i<partElements.size()/3;i++) {
        Pnt pa = {points[partElements[3*i]][0],points[partElements[3*i]][1]};
        Pnt pb = {points[partElements[3*i+1]][0],points[partElements[3*i+1]][1]};
        Pnt pc = {points[partElements[3*i+2]][0],points[partElements[3*i+2]][1]};
        cout << "Pa" << pa.x << " " << pa.y;
        cout << "Pb" << pb.x << " " << pb.y;
        cout << "Pc" << pc.x << " " << pc.y;
        cout << "Pd" << pd.x << " " << pd.y;
        if(inCircle(pa,pb,pc,pd)) {
          cout << "Lies inside";
          int a[3] = {partElements[3*i], partElements[3*i+1], partElements[3*i+2]};
          trashTriangles.insert(trashTriangles.end(), a, a+3);
          trashIndices.push_back(i);
        } else {
          cout << "Lies outside";
        }
      }

      vector<int> newEdgelist;
      for(int i=0;i<trashTriangles.size()/3;i++) {
        
        newEdgelist.push_back(trashTriangles[3*i]);
        newEdgelist.push_back(trashTriangles[3*i+1]);
        newEdgelist.push_back(trashTriangles[3*i+2]);
        newEdgelist.push_back(trashTriangles[3*i]);
        newEdgelist.push_back(trashTriangles[3*i+1]);
        newEdgelist.push_back(trashTriangles[3*i+2]);
        // find all the unique edges between these set of triangles
      }

      vector<int> trashEdges;
      for(int i=0;i<newEdgelist.size()/2;i++) {
        for(int j=i;j<newEdgelist.size()/2;j++) {
          if(i==j) {continue;}
          if((newEdgelist[2*i] == newEdgelist[2*j] && newEdgelist[2*i+1] == newEdgelist[2*j+1])
          || (newEdgelist[2*i] == newEdgelist[2*j+1] && newEdgelist[2*i+1] == newEdgelist[2*j])) {
            trashEdges.push_back(i);
            trashEdges.push_back(j);
          }
        }
      }

      // remove duplicates and put it in descending order
      std::unordered_set<int> s(trashEdges.begin(), trashEdges.end());
      trashEdges.assign(s.begin(), s.end());
      sort(trashEdges.begin(), trashEdges.end(), greater<int>());
      // remove from newEdgeList
      for(int i=0;i<trashEdges.size();i++) {
        newEdgelist.erase(newEdgelist.begin() + trashEdges[2*i], newEdgelist.begin() + trashEdges[2*i+1]);
      }

      // adding new vertex
      int newVertexId = nVertices;
      points[newVertexId] = {pd.x,pd.y};
      nVertices++;
      
      // add the remaining edges to cc to form triangles
      for(int i=0;i<newEdgelist.size()/2;i++) {
        // push back all three vertices
        partElements.push_back(newEdgelist[2*i]);
        partElements.push_back(newEdgelist[2*i+1]);
        partElements.push_back(newVertexId);
      }

      // trashtriangles contains the list of triangles that needs to be removed
      for(int i=0;i<trashIndices.size();i++) {
        cout << trashIndices[i] << endl;
        partElements.erase(partElements.begin() + trashIndices[i], partElements.begin() + trashIndices[i]+2);
      }
      
    } // end of each bad triangle loop check

    cout << endl;

  } // end of each triangle loop check


  int val = 55;
  cout << endl << "Im here" << endl;
  for(int i=0;i<nProcesses;i++) {
    if(i == process_Rank)
      continue;
    MPI_Send(&val, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
  }


  while (1) {
    MPI_Status status;
    int val[2];
    MPI_Recv(val, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("root recev %d,%d from %d with tag = %d\n" , val[0], val[1] , status.MPI_SOURCE , status.MPI_TAG );fflush(stdout);
    // split the edge and construct local cavity for the midpoint

    if (status.MPI_TAG == 2)
    num_of_DONE++;
    printf("num_of_DONE=%d\n" , num_of_DONE);fflush(stdout);
    if(num_of_DONE == nProcesses-1)
    break;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
    
  return 0;
}