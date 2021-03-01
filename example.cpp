#include "Delaunay.hpp"
#include <vector>
#include <iostream>
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
vector<vector<float>> points;

// Vertices XS - x coordinates & YS - y coordinates
vector<int> XS;
vector<int> YS;

float threshold1 = 3;
float threshold2 = 0.5;

float triangleArea(Pnt p1, Pnt p2, Pnt p3) {         //find area of triangle formed by p1, p2 and p3
   return abs((p1.x*(p2.y-p3.y) + p2.x*(p3.y-p1.y)+ p3.x*(p1.y-p2.y))/2.0);
}

bool ccw (float ax, float ay, float bx, float by, float cx, float cy) {
    return (bx - ax)*(cy - ay)-(cx - ax)*(by - ay) > 0;
}

bool inCircle(Pnt p1, Pnt p2, Pnt p3, Pnt p) {     //check whether p is inside or outside
  bool ccwV = ccw(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
  bool var;
  
  float ax_ = p1.x-p.x;
  float ay_ = p1.y-p.y;
  float bx_ = p2.x-p.x;
  float by_ = p2.y-p.y;
  float cx_ = p3.x-p.x;
  float cy_ = p3.y-p.y;
  float val = (ax_*ax_ + ay_*ay_) * (bx_*cy_-cx_*by_) -
      (bx_*bx_ + by_*by_) * (ax_*cy_-cx_*ay_) +
      (cx_*cx_ + cy_*cy_) * (ax_*by_-bx_*ay_);
  if(val > 0) {
    if(ccwV) {
      var = 1;
    } else {
      var = 0;
    }
  } else if (val < 0) {
    if(ccwV){
      var = 0;
    } else {
      var = 1;
    }
  } else {
    var = 1;
  }
  return var;
}

float* circumcenter(float Px, float Py, float Qx, float Qy, float Rx, float Ry
                    ,float x, float y, float z)
{
  float* arr = new float[2];
  float angleA = acos((-x*x+y*y+z*z)/(2*y*z));
  float angleB = acos((x*x-y*y+z*z)/(2*x*z));
  float angleC = acos((x*x+y*y-z*z)/(2*x*y));

  arr[0] = (Px * sin(2*angleA) + Qx*sin(2*angleB) + Rx*sin(2*angleC))/
           (sin(2*angleA) + sin(2*angleB) + sin(2*angleC));
  arr[1] = (Py * sin(2*angleA) + Qy*sin(2*angleB) + Ry*sin(2*angleC))/
           (sin(2*angleA) + sin(2*angleB) + sin(2*angleC));
           
  return arr;
}

float distance(float x1, float y1, float x2, float y2) 
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

void localCavityCreation(vector<int> &partElements, vector<int> &eind, int &nVertices, vector<vector<float>> &points, Pnt pd, vector<int> &epart, int process_Rank)
{
  // Local cavity
  // vector to keep track of all the triangles to delete
  vector<int> trashTriangles;
  vector<int> trashIndices;
  // do not use eind here

  cout << "Point to be added: " << pd.x << "," << pd.y << endl;

  cout << "PartElements size:" << partElements.size() << "from:" << process_Rank << endl;

  for(int i=0;i<partElements.size();i++) {
    Pnt pa = {points[eind[partElements[i]*3]][0],points[eind[partElements[i]*3]][1]};
    Pnt pb = {points[eind[partElements[i]*3+1]][0],points[eind[partElements[i]*3+1]][1]};
    Pnt pc = {points[eind[partElements[i]*3+2]][0],points[eind[partElements[i]*3+2]][1]};
    //comment cout << "PE[i]" << partElements[i] << endl;
    if(inCircle(pa,pb,pc,pd)) {
      cout << "Lies inside" << endl;
      //comment cout << pa.x << " " << pa.y << " " << pb.x 
      // << " " << pb.y << " " <<
      //         pc.x <<  " " << pc.y
      //         << " " << pd.x << " "<< pd.y << endl;
      int a[3] = {eind[partElements[i]*3], eind[partElements[i]*3+1], eind[partElements[i]*3+2]};
      trashTriangles.insert(trashTriangles.end(), a, a+3);
      trashIndices.push_back(i);
    } else {
      //comment cout << "Lies outside";
    }
  }

  vector<int> newEdgelist;
  cout << "TTSize:" << trashTriangles.size() << endl;
  for(int i=0;i<trashTriangles.size()/3;i++) {
    newEdgelist.push_back(trashTriangles[3*i]);
    newEdgelist.push_back(trashTriangles[3*i+1]);
    newEdgelist.push_back(trashTriangles[3*i+2]);
    newEdgelist.push_back(trashTriangles[3*i]);
    newEdgelist.push_back(trashTriangles[3*i+1]);
    newEdgelist.push_back(trashTriangles[3*i+2]);
    // find all the unique edges between these set of triangles
  }

  //comment cout << "New Edge List:" << endl;
  // for(int i=0;i<newEdgelist.size();i++) {
  //   cout << newEdgelist[i] << endl;
  // }
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
    newEdgelist.erase(newEdgelist.begin() + trashEdges[i]*2, newEdgelist.begin() + trashEdges[i]*2+1);
  }
  // adding new vertex
  int newVertexId = nVertices;

  // cout << newVertexId << endl;
  // points[nVertices] = {pd.x,pd.y};
  points.push_back({pd.x, pd.y});
  nVertices++;

  cout << "Number of vertices:" << nVertices << endl;
  cout << "New edgelist size:" << newEdgelist.size() << endl;

  // add the remaining edges to cc to form triangles
  for(int i=0;i<newEdgelist.size()/2;i++) {
    // push back all three vertices
    eind.push_back(newEdgelist[2*i]);
    eind.push_back(newEdgelist[2*i+1]);
    eind.push_back(nVertices-1);
    cout << "Added for:" << process_Rank << (eind.size()/3)-1;
    partElements.push_back((eind.size()/3)-1);
    epart[(eind.size()/3)-1] = process_Rank;
  }
  sort(trashIndices.begin(), trashIndices.end(), greater<int>());
  // trashtriangles contains the list of triangles that needs to be removed
  for(int i=0;i<trashIndices.size();i++) {
    // cout << trashIndices[i] << endl;
    partElements.erase(partElements.begin() + trashIndices[i]);
  }

  //comment for(int i=0;i<partElements.size();i++) {
  //   cout << "PESize:" << partElements[i] << endl;
  // }

  newEdgelist.clear();
  trashTriangles.clear();
  trashIndices.clear();
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

  ///////////////////////////////////////////////////////////////
  // partmeshNodal
  ///////////////////////////////////////////////////////////////
  int nVertices = 29;
  idx_t nElements = 29;
  idx_t nParts = nProcesses;

  idx_t objval;
  std::vector<idx_t> epart(nElements, 0);
  std::vector<idx_t> npart(nVertices, 0);
  std::vector<idx_t> recpart(nVertices, 0);

  std::vector<idx_t> eind;
  std::vector<idx_t> eptr;

  std::ifstream infile("A.1.ele");
  std::string line;
  std::getline(infile, line);
  std::istringstream iss(line);
  int a,b,c;
  if (!(iss >> a >> b >> c)) {}
  cout << a;
  eptr.push_back(a);
  int count = 0;
  while (a > 0)
  {
    std::getline(infile, line);
    std::istringstream iss(line);
    int a2, b, c, d;
    if (!(iss >> a2 >> b >> c >> d)) { break; } // error

    // process pair (a,b)
    // cout << a2 << b << c << endl;
    int v[3] = {b-1, c-1, d-1};
    eind.insert(eind.end(), v, v+3);
    eptr.push_back(3*count);
    a--;count++;
  }

  std::ifstream infile2("A.1.node");
  std::string line2;
  std::getline(infile2, line2);
  std::istringstream iss2(line2);
  int a_c;
  if (!(iss2 >> a_c)) {}
  cout << a;
  while (a_c > 0)
  {
    std::getline(infile2, line2);
    std::istringstream iss2(line2);
    float a2, b, c, d;
    if (!(iss2 >> a2 >> b >> c >> d)) { break; } // error

    // process pair (a,b)
    // cout << a2 << b << c << endl;
    vector<float> v = {b, c};
    points.push_back(v);
    a_c--;
  }

  // remove hardcoding
  // std::vector<idx_t> eptr = {8,0,3,6,9,12,15,18,21};

  // cout << "-----POINTS-------" << endl;
  // for(int i=0;i<points.size();i++) {
  //   cout << points[i][0] << points[i][1] << endl;
  // }
  // cout << "---------------------" << endl;

  // capture from A.1.ele
  // std::vector<idx_t> eind = {0,4,5,4,5,6,4,6,7,4,7,3,0,4,3,0,1,3,1,2,3,2,3,7};

  int ret2 = METIS_PartMeshNodal( 
    &nElements, &nVertices, eptr.data(), eind.data(), NULL, NULL,
    &nParts, NULL, NULL, &objval, epart.data(), npart.data());
    
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

  if(process_Rank == 0){
    cout << endl <<"----------Elements Partition-------" << endl;
    for(unsigned part_i = 0; part_i < epart.size(); part_i++){
      std::cout << "Element " << part_i << " Allotted to the P" << epart[part_i] << std::endl;
    }
  }

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

  // cout << "------------------partElements------------" << endl;

  // for(int i=0;i<partElements.size();i++) {
  //   cout << partElements[i] << ": " << process_Rank << endl;
  // }

  for(int i=0;i<partElements.size();i++) {
    if(partElements[i]>30)
    break;
    //vertices of that triangle
    float Px = points[eind[partElements[i]*3]][0];
    float Py = points[eind[partElements[i]*3]][1];
    float Qx = points[eind[partElements[i]*3+1]][0];
    float Qy = points[eind[partElements[i]*3+1]][1];
    float Rx = points[eind[partElements[i]*3+2]][0];
    float Ry = points[eind[partElements[i]*3+2]][1];
    
    // cout << Px << "," << Py << ","<< Qx << ","<< Qy << ","<< Rx << ","<< Ry << endl;

    float x = distance(Qx,Qy,Rx,Ry);
    float y = distance(Rx,Ry,Px,Py);
    float z = distance(Px,Py,Qx,Qy);

    float shortest = x < y ? (x < z ? x : z) : (y < z ? y : z);
    float area = triangle_area(x,y,z);
    float circumRadius = x*y*z/(4*area);

    cout << "circumRadius/shortest" << circumRadius/shortest << endl;
    cout << "area" << area << endl;

    ///////////////////////////////////////////////////////////////
    // obtaining bad triangles
    ///////////////////////////////////////////////////////////////
    if(circumRadius/shortest > threshold1 || area > threshold2) {
      //bad triangle, find the circumcenter
      // cout << "IN" << endl;
      // cout << "Shortest" << shortest << endl;
      // cout << "Circumradius" << circumRadius << endl;
      // cout << "Circum/shortest" << circumRadius/shortest << endl; 
      // cout << "Area" << area << endl;

      //circumcenter co-ordinates
      float* ptr = circumcenter(Px,Py,Qx,Qy,Rx,Ry,x,y,z);

      // cout << "Triangle:" << Px << " " << Py << " " << Qx << " " << Qy 
      //       << " " << Rx << " " << Ry << endl;
      // cout << ptr[0] << ptr[1] << endl;

      // cout << ptr[0] << " " << ptr[1] << endl;
      //got circumcenter

      //commonedges with processid
      // cout << "eind.length/3" << eind.size()/3 << endl;
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
              // cout << eind[3*i] << eind[3*i+1] << eind[3*i+2] <<
              //                 eind[3*j] << eind[3*j+1] << eind[3*j+2] << endl;
              int* ptr2 = checkCommonEdge(eind[3*i],eind[3*i+1],eind[3*i+2],
                              eind[3*j],eind[3*j+1],eind[3*j+2]);
              // cout << "cEdge[0]:" << ptr2[0] << endl;
              // cout << "cEdge[1]:" << ptr2[1] << endl;
              if(!(ptr2[0] == -999 && ptr2[1] == -999)){
                edgeList.push_back(ptr2[0]);
                edgeList.push_back(ptr2[1]);
                // cout << "epart[j]" << j << epart[j] << endl;
                edgeList.push_back(epart[j]);
              }
            }
          }
        }
      }

      // for(int i=0;i<edgeList.size();i++) {
      //   cout << i << ":" << edgeList[i] << endl;
      // }

      // Local cavity
      Pnt pd = {ptr[0], ptr[1]};

      //find shared edges of this process with other processes
      //take common edges and check whether the circumcenter encroaches the segment
      for(int i=0;i<edgeList.size()/3;i++) {
        float centerX = midPoint(points[edgeList[3*i]][0],points[edgeList[3*i+1]][0]);
        float centerY = midPoint(points[edgeList[3*i]][1],points[edgeList[3*i+1]][1]);
        float radius = distance(centerX,centerY,points[edgeList[3*i]][0],
                                points[edgeList[3*i]][1]);
        if(pow((ptr[0]-centerX),2)+pow((ptr[1]-centerY),2)-pow(radius,2) < 0) {
          // cout << "Encroaches" << endl;
          // MPI
          // sending mpi message to the respective process
          int data[3];
          data[0] = edgeList[3*i]; data[1] = edgeList[3*i+1]; data[2] = edgeList[3*i+2];
          // MPI_Send(data,2,MPI_INT,edgeList[3*i+2],NULL,NULL);

          cout << "Sending MPI Message to designated processes" << endl;
          cout << "-------------------------------------------" << endl;
          cout << "msg is being sent to" << data[2] << endl;
          // needed fix here
          // if(data[2] < nProcesses)
          MPI_Send(data,3,MPI_INT,data[2],1,MPI_COMM_WORLD);
          cout << "-------------------------------------------" << endl;
          cout << endl;

          pd.x = centerX;
          pd.y = centerY;
          break;

        } else { // in bad triangle // common edges

          // cout << "Doesnot Encroaches" << endl;

        } // end of doesnot encroaches loop

      } // end of each common edge check

      // local cavity creation
      cout << "Hello all" << endl;

      // cout << "PE: " << endl;
      // for(int i=0;i<partElements.size();i++) {
      //   cout << partElements[i] << endl;
      // }

      // cout << "eind: " << endl;
      // for(int i=0;i<eind.size();i++) {
      //   cout << eind[i] << endl;
      // }

      // cout << "points: " << endl;
      // for(int i=0;i<points.size();i++) {
      //   cout << points[i][0] << points[i][1] << endl;
      // }

      localCavityCreation(partElements, eind, nVertices, points, pd, epart, process_Rank);
      // MPI_Bcast(&points,points.size(),MPI_FLOAT,process_Rank,MPI_COMM_WORLD);
      // cout << endl << "From BC:" << process_Rank << points.size() << endl;
      // MPI_Bcast(&nVertices,1,MPI_INT,process_Rank,MPI_COMM_WORLD);
      // MPI_Bcast(&eind,eind.size(),MPI_INT,process_Rank,MPI_COMM_WORLD);
      // MPI_Bcast(&epart,epart.size(),MPI_INT,process_Rank,MPI_COMM_WORLD);
      cout << endl << "Out of local cavity" << endl;
      
    } // end of each bad triangle loop check

    cout << endl;

  } // end of each triangle loop check


  int data[3];
  data[0] = 10; data[1] = 11; data[2] = 12;
  for(int i=0;i<nProcesses;i++) {
    if(i == process_Rank)
      continue;
    MPI_Send(data, 3, MPI_INT, i, 2, MPI_COMM_WORLD);
  }

  while (1) {
    MPI_Status status;
    int val[3];
    MPI_Recv(val, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    cout << "----------------------------------------" << endl;
    printf("root recev %d,%d from %d with tag = %d\n" , val[0], val[1] , status.MPI_SOURCE , status.MPI_TAG );fflush(stdout);
    // split the edge and construct local cavity for the midpoint
    // local cavity creation
    Pnt pd;
    if(status.MPI_TAG != 2) {
      pd.x = (points[val[0]][0]+points[val[0]][1])/2;
      pd.y = (points[val[1]][0]+points[val[1]][1])/2;
      localCavityCreation(partElements, eind, nVertices, points, pd, epart, process_Rank);
    }

    if (status.MPI_TAG == 2){
      num_of_DONE++;
      printf("num_of_DONE=%d\n" , num_of_DONE);fflush(stdout);
      if(num_of_DONE == nProcesses-1){
        //send current points, partElements and eind to file
        std::ofstream outfile ("test"+std::to_string(process_Rank) +".txt");
        outfile << partElements.size() << endl;
        for(int i=0;i<partElements.size();i++) {
          outfile << partElements[i] << endl;
        }
        outfile << eind.size() << endl;
        for(int i=0;i<eind.size();i++) {
          outfile << eind[i] << endl;
        }
        outfile << points.size() << endl;
        for(int i=0;i<points.size();i++){
          outfile << points[i][0] << " " << points[i][1] << endl;
        }
        break;
      }
    }
}

  cout << "P:" << points.size() << endl;

  if(process_Rank == 0) {
    for(int i=0;i<points.size();i++) {
      cout << points[i][0] << points[i][1] << endl;
    }
  }

  cout << process_Rank << endl;
  for(int i=0;i<partElements.size();i++) {
    cout << partElements[i] << endl;
  }

  // all processes send their extra points and eind to P0 where P0 does all the calculations
  // comment for(int i=0;i<eind.size();i++) {
  // cout << "eind.size()" << eind[i] << endl;
  // }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}