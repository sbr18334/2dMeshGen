#include "../Delaunay.hpp"
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

#include <sys/time.h>

// To compile and run
// mpic++ -std=c++11 example.cpp -lmetis
// mpirun -np 4 ./a.out

using namespace std;

struct Pnt {
  float x, y;
};

vector<vector<float>> points;

// Vertices XS - x coordinates & YS - y coordinates
vector<int> XS;
vector<int> YS;

float threshold1 = 3;
float threshold2 = 2;

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
    if(ccwV) { var = 1; } 
    else { var = 0; }
  } else if (val < 0) {
    if(ccwV){ var = 0; }
    else { var = 1; }
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
  vector<int> a;
  int aTemp[3] = {a1, a2, a3};
  vector<int> b;
  int bTemp[3] = {b1, b2, b3};
  a.insert(a.end(), aTemp, aTemp+3);
  b.insert(b.end(), bTemp, bTemp+3);
  sort(a.begin(), a.end(), greater<int>());
  sort(b.begin(), b.end(), greater<int>());
  int* arr = new int[2];
  if((a[0]==b[0] && a[1] == b[1]) || (a[0] == b[1] && a[1] == b[2]) ||
     (a[0]==b[0] && a[1] == b[2])) {
      arr[0] = a[0]; arr[1] = a[1];
      return arr;
  }
  if((a[0]==b[0] && a[2] == b[1]) || (a[0]==b[1] && a[2] == b[2]) ||
     (a[0]==b[0] && a[2] == b[2])) {
      arr[0] = a[0]; arr[1] = a[2];
      return arr;
  }
  if((a[1]==b[0] && a[2] == b[1]) || (a[1]==b[1] && a[2] == b[2])) {
      arr[0] = a[0]; arr[1] = a[2];
      return arr;
  }
  else {
    arr[0] = -999;arr[1] = -999;
    return arr;
  }
}

void localCavityCreation(vector<int> &partElements, vector<int> &eind, int &nVertices, vector<vector<float>> &points, Pnt pd, vector<int> &epart, int process_Rank)
{
  vector<int> trashTriangles;
  vector<int> trashIndices;
  for(int i=0;i<partElements.size();i++) {
    Pnt pa = {points[eind[partElements[i]*3]][0],points[eind[partElements[i]*3]][1]};
    Pnt pb = {points[eind[partElements[i]*3+1]][0],points[eind[partElements[i]*3+1]][1]};
    Pnt pc = {points[eind[partElements[i]*3+2]][0],points[eind[partElements[i]*3+2]][1]};
    if(inCircle(pa,pb,pc,pd)) {
      int a[3] = {eind[partElements[i]*3], eind[partElements[i]*3+1], eind[partElements[i]*3+2]};
      trashTriangles.insert(trashTriangles.end(), a, a+3);
      trashIndices.push_back(i);
    }
  }

  vector<int> newEdgelist;
  for(int i=0;i<trashTriangles.size()/3;i++) {
    int a[6] = {trashTriangles[3*i], trashTriangles[3*i+1], trashTriangles[3*i+2],
    trashTriangles[3*i], trashTriangles[3*i+1], trashTriangles[3*i+2]};
    newEdgelist.insert(newEdgelist.end(), a, a+6);
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
    newEdgelist.erase(newEdgelist.begin() + trashEdges[i]*2+1);
    newEdgelist.erase(newEdgelist.begin() + trashEdges[i]*2);
  }
  // adding new vertex
  int newVertexId = nVertices;
  points.push_back({pd.x, pd.y});
  nVertices++;

  // add the remaining edges to cc to form triangles
  for(int i=0;i<newEdgelist.size()/2;i++) {
    // push back all three vertices
    if(midPoint(points[newEdgelist[2*i]][0], points[newEdgelist[2*i+1]][0])
    == pd.x && midPoint(points[newEdgelist[2*i]][1], points[newEdgelist[2*i+1]][1])
    == pd.y) {
      continue;
    }
    int a[3] = {newEdgelist[2*i], newEdgelist[2*i+1], newVertexId};
    eind.insert(eind.end(), a, a+3);
    // cout << "Added for:" << process_Rank << (eind.size()/3)-1;
    partElements.push_back((eind.size()/3)-1);
    epart[(eind.size()/3)-1] = process_Rank;
  }
  sort(trashIndices.begin(), trashIndices.end(), greater<int>());
  // trashtriangles contains the list of triangles that needs to be removed
  for(int i=0;i<trashIndices.size();i++) {
    partElements.erase(partElements.begin() + trashIndices[i]);
}

  newEdgelist.clear();
  trashTriangles.clear();
  trashIndices.clear();
  trashEdges.clear();
}

int main(int argc, char** argv)
{
struct timeval  tv1, tv2;
gettimeofday(&tv1, NULL);
  //Number of processes
  int nProcesses = 4;

  ///////////////////////////////////////////////////////////////
  // partmeshNodal
  ///////////////////////////////////////////////////////////////
  int nVertices = 33343;
  idx_t nElements = 64125;
  idx_t nParts = nProcesses;

  idx_t objval;
  std::vector<idx_t> epart(nElements, 0);
  std::vector<idx_t> npart(nVertices, 0);
  std::vector<idx_t> recpart(nVertices, 0);

  std::vector<idx_t> eind;
  idx_t eind1[nElements*3];
  idx_t eptr[nElements+1];
  std::ifstream infile("../input/A.1.ele");
  std::string line;
  std::getline(infile, line);
  std::istringstream iss(line);
  int a,b,c;
  if (!(iss >> a >> b >> c)) {}
  int count = 0;
  int eindIndex = 0;
  while (a > 0)
  {
    std::getline(infile, line);
    std::istringstream issT(line);
    int a2, b2, c2, d;
    if (!(issT >> a2 >> b2 >> c2 >> d)) { break; }
    int v[3] = {b2-1, c2-1, d-1};
    // eind.insert(eind.end(), v, v+3);
    eind1[eindIndex] = b2-1;eindIndex++;
    eind1[eindIndex] = c2-1;eindIndex++;
    eind1[eindIndex] = d-1;eindIndex++;
    a--;count++;
  }
  for(int i=0;i<=count;i++){
    eptr[i] = i*3;
  }
  // eptr[0] = 210;
  std::ifstream infile2("../input/A.1.node");
  std::string line2;
  std::getline(infile2, line2);
  std::istringstream iss2(line2);
  int a_c;
  if (!(iss2 >> a_c)) {}
//   cout << a_c;
  while (a_c > 0)
  {
    std::getline(infile2, line2);
    std::istringstream iss2T(line2);
    float a2, b3, c3, d;
    if (!(iss2T >> a2 >> b3 >> c3 >> d)) { break; }
    vector<float> v = {b3, c3};
    points.push_back(v);
    a_c--;
  }
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  idx_t *npt;
  npt = (idx_t *)malloc(nVertices*10);
  idx_t *xadj, *adjncy;
  idx_t pnumflag = 0;
  idx_t ncommon = 2;
  idx_t ncon=1;
  idx_t nparts = nProcesses;
  idx_t objvalue;
  xadj =  (idx_t *) malloc(nVertices*2);
  adjncy =  (idx_t *) malloc(20*nVertices);

  METIS_MeshToDual(&nElements, &nVertices, eptr, eind1, &ncommon, &pnumflag, &xadj, &adjncy);

  METIS_PartGraphKway(&nElements, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts,
  NULL, NULL, options, &objvalue, npt);

  for(int i=0;i<nElements*3;i++) {
    eind.push_back(eind1[i]);
  }

  delete xadj;
  delete adjncy;

  ////////////////////////////////////////////////////////////////
  // mpi process
  // message between and parallelization between threads
  ////////////////////////////////////////////////////////////////
  // To run: mpic++ hello_world_mpi.cpp -o hello_world_mpi
  // To execute: mpirun -np 2 ./hello_world_mpi
  ////////////////////////////////////////////////////////////////
  int process_Rank, size_Of_Cluster;

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
    for(unsigned part_i = 0; part_i < eind.size()/3; part_i++){
      std::cout << "Element " << part_i << " Allotted to the P" << npt[part_i] << std::endl;
    }
  }

  int num_of_DONE = 0;
  vector<int> partNodes;
  vector<int> partElements;

  for(unsigned part_i = 0; part_i < npart.size(); part_i++){
    if(npart[part_i] == process_Rank) {
      partNodes.push_back(part_i);
    }
  }
  for(unsigned part_i = 0; part_i < eind.size()/3; part_i++){
    if(npt[part_i] == process_Rank) {
      partElements.push_back(part_i);
    }
  }
  
  delete npt;
  npart.clear();
  recpart.clear();

  int triangleCount = eind.size()/3;

  for(int i=0;i<partElements.size();i++) {
    //vertices of that triangle
    float Px = points[eind[partElements[i]*3]][0];
    float Py = points[eind[partElements[i]*3]][1];
    float Qx = points[eind[partElements[i]*3+1]][0];
    float Qy = points[eind[partElements[i]*3+1]][1];
    float Rx = points[eind[partElements[i]*3+2]][0];
    float Ry = points[eind[partElements[i]*3+2]][1];
    
    float x = distance(Qx,Qy,Rx,Ry);
    float y = distance(Rx,Ry,Px,Py);
    float z = distance(Px,Py,Qx,Qy);

    float shortest = x < y ? (x < z ? x : z) : (y < z ? y : z);
    float area = triangle_area(x,y,z);
    float circumRadius = x*y*z/(4*area);

    cout << "circumRadius/shortest" << circumRadius/shortest << endl;
    cout << "area" << area << endl;

    //int triangleCount = eind.size()/3;

    ///////////////////////////////////////////////////////////////
    // obtaining bad triangles
    ///////////////////////////////////////////////////////////////
    if(circumRadius/shortest > threshold1 || area > threshold2) {

      //circumcenter co-ordinates
      float* ptr = circumcenter(Px,Py,Qx,Qy,Rx,Ry,x,y,z);
      
      // Local cavity
      Pnt pd = {ptr[0], ptr[1]};

      std::vector<idx_t> edgeList;
      for(int j=0;j<triangleCount;j++) {
        //check whether the triangles have a common edge
        if(epart[j] != process_Rank){
          int* ptr2 = checkCommonEdge(eind[3*partElements[i]],eind[3*partElements[i]+1],eind[3*partElements[i]+2],
                        eind[3*j],eind[3*j+1],eind[3*j+2]);
          if(!(ptr2[0] == -999 && ptr2[1] == -999)){
            float centerX = midPoint(points[ptr2[0]][0],points[ptr2[1]][0]);
            float centerY = midPoint(points[ptr2[0]][1],points[ptr2[1]][1]);
            float radius = distance(centerX,centerY,points[ptr2[0]][0],
                                    points[ptr2[0]][1]);
            if(pow((ptr[0]-centerX),2)+pow((ptr[1]-centerY),2)-pow(radius,2) < 0) {
              int data[3];
              data[0] = ptr2[0]; data[1] = ptr2[1]; data[2] = epart[j];
              // MPI_Send(data,2,MPI_INT,edgeList[3*i+2],NULL,NULL);

              cout << "Sending MPI Message to designated processes" << endl;
              cout << "-------------------------------------------" << endl;
              cout << "msg is being sent to" << data[2] << endl;
              MPI_Send(data,3,MPI_INT,data[2],1,MPI_COMM_WORLD);
              cout << "-------------------------------------------" << endl;
              cout << endl;

              pd.x = centerX;
              pd.y = centerY;

              cout << "Encroaches" << endl;
              break;

            }
          }
        }
      }
      edgeList.clear();

      // local cavity creation
      cout << pd.x << pd.y << "up:" << process_Rank << endl;
      bool boolean = true;
      for(int m=0;m<points.size();m++) {
        if(pd.x == points[m][0] && pd.y == points[m][1]) {
          boolean = false;
        }
      }
      if(boolean)
      localCavityCreation(partElements, eind, nVertices, points, pd, epart, process_Rank);
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
    Pnt pd;
    if(status.MPI_TAG != 2) {
      pd.x = midPoint(points[val[0]][0],points[val[1]][0]);
      pd.y = midPoint(points[val[0]][1],points[val[1]][1]);
      cout << pd.x << pd.y << "down:" << process_Rank << endl;
      bool boolean = true;
      for(int i=0;i<points.size();i++) {
        if(pd.x == points[i][0] && pd.y == points[i][1]) {
          boolean = false;
        }
      }
      if(boolean)
      localCavityCreation(partElements, eind, nVertices, points, pd, epart, process_Rank);
    }

    if (status.MPI_TAG == 2){
      num_of_DONE++;
      printf("num_of_DONE=%d\n" , num_of_DONE);fflush(stdout);
      if(num_of_DONE == nProcesses-1){
        std::ofstream outfile ("../output/test"+std::to_string(process_Rank) +".txt");
        cout << partElements.size() << endl;
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
  partElements.clear();
  epart.clear();
  eind.clear();
  points.clear();

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

gettimeofday(&tv2, NULL);

printf ("Total time = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
  return 0;
}