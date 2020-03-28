#include <iostream>
#include "argraph.h"
#include "argedit.h"
#include "vf2_sub_state.h"
#include "match.h"

#define MAXNODES 200
using namespace std;

class NodeAttrComparator: public AttrComparator
{ 
  public:
    NodeAttrComparator() { }
    virtual bool compatible(void *pa, void *pb) { 
      int *a = (int *)pa;
      int *b = (int *)pb;
      return (*a==*b);
    }
};

class NodeAttrDestroyer: public AttrDestroyer
{ 
  public:
       virtual void destroy(void *p) { int * intP = (int *) p; delete intP; }
};

bool my_visitor(int n, node_id ni1[], node_id ni2[], void *usr_data) { 

  for(int i=0; i<n; i++) 
    cout << "(" << ni1[i] << " " << ni2[i] << ")";
  cout << endl;
  // Return false to search for the next matching
  return false;
}

int main() {

  ARGEdit ed1, ed2;
  ed1.InsertNode( new int(10));
  ed1.InsertNode( new int(10));
  ed1.InsertNode( new int(20));
  ed1.InsertEdge(0,1,NULL); 
  ed1.InsertEdge(1,2,NULL); 
  ed1.InsertEdge(0,2,NULL); 

  ed2.InsertNode( new int(10));
  ed2.InsertNode( new int(10));
  ed2.InsertEdge( 1,0,NULL);

  Graph large_graph(&ed1);
  Graph small_graph(&ed2);

  // install the attribute destroyer
  large_graph.SetNodeDestroyer(new NodeAttrDestroyer());
  small_graph.SetNodeDestroyer(new NodeAttrDestroyer());
  
  small_graph.SetNodeComparator(new NodeAttrComparator());
  VF2SubState s0(&small_graph, &large_graph);
  
  int n;
  node_id ni1[MAXNODES], ni2[MAXNODES];

/*
  if (!match(&s0, &n, ni1, ni2)) {
    cout << "No match found\n";
  }

  cout << "found a matching with " << n << " nodes:\n"; 
  for (int i=0; i<n; i++) {
    cout << "\tNode " << ni1[i] << " of graph 1 is paired with node " << ni2[i] << " of graph 2\n";
  }
*/
  match(&s0, my_visitor);
  return 0;
}
