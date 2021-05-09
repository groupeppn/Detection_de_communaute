#include <string>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <deque>
#include <climits>
#include <unordered_map>
#include <algorithm>
#include <unistd.h>
#include <ctime>
#include <chrono>
#include <fstream>
#include <omp.h>

using namespace std;

class Vertex
  {
    int vert; 
    int degree; 
    deque<int> *nbrs; 
  public:
    Vertex() { degree = 0; nbrs = new deque<int>(); }
    Vertex(const Vertex& x) { vert = x.vert; degree = x.degree; nbrs = new deque<int>(x.nbrs->begin(),x.nbrs->end()); }
    ~Vertex() { degree = 0; nbrs->clear(); delete nbrs;}
    void clear(void) { degree = 0; nbrs->clear(); }
    int getV(void) { return vert; } 
    void setV(int x) { vert = x; } 
    int getDegree(void) { return degree; }
    void setDegree(void) { degree = nbrs->size(); }
    deque<int>* getNbrs(void) { return nbrs; }
    bool isAdjacent (int x)
      {
        int i;
        for (i=0;i<degree;i++) if ((*nbrs)[i] == x) break;
        return i<degree; // if less, <x> is in the list of nbrs and <vert> and <x> are connected
      }

    void setNbrs(deque<int> x)
      {
         nbrs->clear();
         nbrs->assign(x.begin(),x.end());
         degree = x.size();
      } 
    void addNbr(int y)
      {
         nbrs->push_back(y);
         degree++;
      } 
    void removeNbr(int y)
      {
        // we are supposed to be connected with vertex y and have to remove it
        bool k = false;
        int j = degree;
        
        #pragma omp parallel for
        for (int i=0;i<degree;i++)
          {
            if (k) continue;
            if ((*nbrs)[i]==y)
              {
                j = i;
                k = true;
              }
          }
        if (j<degree)
          {
            if (j<degree-1) (*nbrs)[j] = (*nbrs)[degree-1];
            nbrs->pop_back();
            degree--;
          }
      } 
  };


class Graph
  {
    int nodesCount; 
    deque<Vertex> adj;
    deque<int> graphSorted;
    unordered_map<string,deque<int>>* cliques;
    bool printState = true;
    bool hashSet = false;
    long int rCount = 0;
  public:
    deque<Vertex> degen;

    void setPrint (void) { printState = true;}
    void setNoPrint (void) { printState = false;}

    Graph(int x)
      { 
        nodesCount = x;
        adj.resize(nodesCount);
        degen.resize(nodesCount);
        graphSorted.resize(nodesCount);
        #pragma omp parallel for shared(adj,degen,graphSorted)
        for (int i = 0; i<nodesCount; i++)
          {
            adj[i].setV(i); 
            degen[i].setV(i);
            graphSorted[i] = i;
          }
      }

    void setHashTable (void)
      {
        if (!hashSet)
          {
            cliques = new unordered_map<string,deque<int>>();
            hashSet = true;
          }
        else cliques->clear();
      }

    void readCSV (string fname, int& vStart)
      {
        // It is a plain CSV reader, no prev info about number of vertices or edges
        FILE* infile;
        int vCount=0;
        int u,v;
        char line[100];
        infile = fopen(fname.c_str(),"r");
        if (infile == NULL)
          {
	    std::cout<<"File does not exist"<<fname.c_str()<<  std::endl;

            exit(2);
          }
        // Read the first header line, it has no information
        if (fgets(line,100,infile)==NULL)
          {
	    std::cout<<"File is empty."<<fname.c_str()<<  std::endl;
            exit(3);
          }
        // The most safe solution is to scan the file first (Pass-1), looking for real number of vertices
        // Also it does not include "back" edges, so, during Pass-2 we have to add (u,v) and (v,u) at once
        do {
          line[0]=(char)0;
          fgets(line,100,infile);
          if ((int)line[0] != 0)
            {
              sscanf(line,"%d,%d",&u,&v);
              if (u>vCount) vCount = u;
              if (v>vCount) vCount = v;
            }
        } while ((int)line[0] != 0);
        // prepare data to store info
        nodesCount = vCount+1-vStart;
        adj.resize(nodesCount);
        degen.resize(nodesCount);
        graphSorted.resize(nodesCount);
        for (int i = 0; i<nodesCount; i++)
          {
            adj[i].setV(i); 
            degen[i].setV(i);
            graphSorted[i] = i;
          }
        // Jump to begining and ignore header
        fseek(infile,0,SEEK_SET);
        fgets(line,100,infile);
        do {
          line[0]=(char)0;
          fgets(line,100,infile);
          if ((int)line[0] != 0)
            {
              sscanf(line,"%d,%d",&u,&v);
              u-=vStart;
              v-=vStart;
              // edge from <u> to <v>        
              adj[u].addNbr(v);
              adj[v].addNbr(u);
              degen[u].addNbr(v);
              degen[v].addNbr(u);
            }
        } while ((int)line[0] != 0);
        // prepare data to store info
        fclose(infile);
      }

    Graph(string fname, int& vStart)
      {
        string s; // To read info from the header
        char line[100];
        int u,v,i,eCount,vCount,iu,iv; // To read pairs of numbers
        FILE* infile;
        size_t found;

        // Check input file name extension. If we have .csv file, it has to be processed by diferent procedure
        found = fname.find(".csv");
        if (found != string::npos)
          {
            readCSV(fname, vStart);
            return;
          }

        infile = fopen(fname.c_str(),"r");
        if (infile == NULL)
          {	    
      std::cout<<"File does not exist.\n"<<fname.c_str()<<  std::endl;

            exit(2);
          }
        // Read the first header line
        if (fgets(line,100,infile)==NULL)
          {
	    std::cout<<"File %s is empty.\n"<<fname.c_str()<<  std::endl;

            exit(3);
          }
        i=0;
        while ((line[i]!=(char)0) && (line[i] != ':')) i++;
        i++;
        while ((line[i]!=(char)0) && (line[i] == ' ')) i++;
        // Read number of nodes 
        while ((line[i]!=(char)0) && (line[i] != ' ')) s+=line[i++];
        // we got number of vertices in graph;
        vCount = stoi(s);
        s.clear();
        while ((line[i]!=(char)0) && (line[i] != ':')) i++;
        i++;
        while ((line[i]!=(char)0) && (line[i] == ' ')) i++;
        // Read number of edges 
        while ((line[i]!=(char)0) && (line[i] != ' ')) s+=line[i++];
        // got number of edges, will be used for check only
        eCount = stoi(s);
        // prepare list
        
        nodesCount = vCount;
        adj.resize(nodesCount);
        degen.resize(nodesCount);
        graphSorted.resize(nodesCount);
        for (int i = 0; i<nodesCount; i++)
          {
            adj[i].setV(i); 
            degen[i].setV(i);
            graphSorted[i] = i;
          }
        
        i = 0; // auxiliary to count number of lines with edge definitions
        do
          {
            line[0]=(char)0;
            fgets(line,100,infile);
            if ( (int)line[0] != 0)
              {
                sscanf(line,"%d\t%d",&u,&v);
                // adjest vertex numbers to start from zero
                u-=vStart;
                v-=vStart;
                // edge from <u> to <v>        
                if ((u<nodesCount)&&(v<nodesCount))
                  {
                    // Correct edge definition
                    adj[u].addNbr(v);
                    degen[u].addNbr(v);
                    i++;
                  }
                else
                  {
                    // Bad edge, gonna report and stop
	    std::cout<<"Bad edge definition: : "<<line<<u<<v<<  std::endl;
                    fclose(infile);
                    adj.clear();
                    degen.clear();
                    exit(5);
                  }
              }
          } while ((int)line[0] != 0);
        /* For sake of safety if the number of actual lines
           with edge definitions does not corespond to the header,
           we report and stop; */
        fclose(infile);
        if (i != eCount)
          {
            std::cout<<"Error reading file Wrong number of edges"<<fname.c_str()<<i<<eCount<<  std::endl;
            adj.clear();
            degen.clear();
            exit(10);
          }
      } 

    ~Graph(void)
      {
        adj.clear();
        degen.clear();
        graphSorted.clear();
        nodesCount = 0;
        if (hashSet) delete cliques;
      }

    void addClique (deque<Vertex> * g)
    {
      string s = "";
      int Vd = g->size();
      int i,v;
      deque<int> clique(Vd);
      // complile in <clique> the list of vertices which compose current graph <g>
      #pragma omp parallel for shared(clique,g) private(i)
      for (i=0; i<Vd; i++) clique[i] = (*g)[i].getV();
      // sort <clique>
      sort(clique.begin(),clique.end());
      // create a hash key string from ordered list of vertices
      for (i=0; i<Vd; i++) s+=to_string(clique[i])+"-";
      cliques->insert({s,clique});
    }

    void printClique (deque<Vertex> * g)
      {
        int sz = g->size();
        std::cout<<"Clique maximale"<< " ";
        for (int j=0; j<sz; j++) 
	 std::cout<<(*g)[j].getV()<<" ";
	 std::cout<<"  "<<std::endl;
       
      }

    void printCliques (void)
    {
      int sz;
      for (const auto & [key, value] : *cliques)
        {
          sz = value.size();
	  std::cout<<"Sous graph"<< " " ;

            for (int j=0; j<sz; j++)  
	        std::cout<<value[j] <<"  ";
	      std::cout<<"  "<<std::endl;
        
        }
    }

    // Finds neighbors of vertex x 
    deque<int>* getNbrs(int x) { return adj[x].getNbrs(); } 

    void getVNbrs(int x, deque<Vertex> *sRes)
      {
        //int vv;
        deque<int>* Nb = adj[x].getNbrs();
        int sz = adj[x].getDegree();
        sRes->resize(sz);
        // here we have to replicate the definitions of <Nb> vertices
        #pragma omp parallel for shared(adj,sRes,sz) //private (vv)
        for (int i=0; i<sz; i++)
          {
            int vv= (*Nb)[i];
            (*sRes)[i].setV(vv);
            (*sRes)[i].setNbrs(*(adj[vv].getNbrs()));
          }
      } 

    void get_intersect(deque<Vertex> *sFirst, deque<Vertex> *sSecond, deque<Vertex> *sRes)
      { 
        int j,vv,sz;
        int whois[nodesCount];
        //Initially mark all as absent
        #pragma parallel omp for shared(whois)
          for (int i=0;i<nodesCount;i++) whois[i]=0;

        for (int i = 0; i < sFirst->size(); i++)
          {
            j = 0;
            vv = (*sFirst)[i].getV();
            // looking to the vertex <v> in the second set
            while ((j< sSecond->size()) && (vv != (*sSecond)[j].getV())) j++;
            // if we found it, it is a common element and must be included into the resulting set
            if (j < sSecond->size()) sRes->push_back((*sFirst)[i]);
          }
        // it is the list of vertices in intersection; we must keep only their mutual nbs
        // marking present vertices
        sz = sRes->size();
        for (int i=0;i<sz;i++) whois[(*sRes)[i].getV()]=1;
        #pragma parallel omp for shared(whois,sRes,sz)
        for (int i=0;i<nodesCount;i++)
          if (whois[i]==0)
            {
              //this vertex can't be in the result
              for (int j=0;j<sz;j++)
                (*sRes)[j].removeNbr(i);
            }
      } 

    void get_union(deque<Vertex> *sFirst, deque<Vertex> *sSecond, deque<Vertex> *sRes)
      { 
        int j,vv,sz;
        int whois[nodesCount];
        //Initially mark all as absent
        #pragma parallel omp for shared(whois)
          for (int i=0;i<nodesCount;i++) whois[i]=0;

        // First set enters fully
        sRes->assign(sFirst->begin(),sFirst->end());
        sz = sFirst->size();
        // but now we can add only different vertices to avoid duplicates - it is prohibited
        for (int i=0; i<sSecond->size(); i++)
          {
            vv = (*sSecond)[i].getV();
            for (j=0; j<sz; j++) if ((*sFirst)[j].getV() == vv) break;
            // if we haven't seen it, adding new element
            if (j==sz) sRes->push_back((*sSecond)[i]);
          }
        // it is the list of vertices in union; we must keep only their mutual nbs
        // marking present vertices
        sz = sRes->size();
        for (int i=0;i<sz;i++) whois[(*sRes)[i].getV()]=1;
        #pragma parallel omp for shared(whois,sRes,sz)
        for (int i=0;i<nodesCount;i++)
          if (whois[i]==0)
            {
              //this vertex can't be in the result
              for (int j=0;j<sz;j++)
                (*sRes)[j].removeNbr(i);
            }
      } 

    void removeNbrsIn(deque<Vertex> *x, Vertex *v)
      { 
        deque<int>* sNb = v->getNbrs();
        for (int i = 0; i < sNb->size(); i++) removeV(x,(*sNb)[i]);
      } 

    long int getRcalls (void) { return rCount; }

    void Bron_KerboschWithPivot(deque<Vertex> * R, deque<Vertex> * P, deque<Vertex> * X, deque<Vertex> * temp)
      { 
        int vv;
        Vertex *u;
        deque<Vertex> * P1;
        deque<Vertex>* P_in;
        deque<Vertex>* X_in;
        
        ++rCount;
//        printf("Recursive call %ld\n",rCount);
        if ((P->size() == 0) && (X->size() == 0))
          {
            if (printState) printClique(R);
            else addClique(R);
            return; 
          } 
        P1  = new deque<Vertex>(P->begin(),P->end());
        P_in = new deque<Vertex>();
        X_in = new deque<Vertex>();
        // Find Pivot 
        get_union(P,X,temp);
        u = new Vertex(getMaxDegreeVertex(temp));

        // P = P / Nbrs(u) 
        removeNbrsIn(P, u); 

        for (int i=0; i<P->size(); i++)
          { 
            int vv = (*P)[i].getV();
            temp->clear();
            getVNbrs(vv, temp);
            #pragma omp parallel
            {
              #pragma omp single //1 seul thread qui executra cette partie 
              {
	//deviser en sous taches
                #pragma omp task shared(R,P,i) 
                  R->push_back((*P)[i]);
                #pragma omp task shared(P1,P_in,temp)
                  get_intersect(P1, temp,P_in);
                #pragma omp task shared(X,X_in,temp)
                  get_intersect(X, temp,X_in);
                #pragma omp taskwait
                {
                  temp->clear();
                  Bron_KerboschWithPivot(R, P_in, X_in,temp); 
                  #pragma omp parallel
                  {
                    #pragma omp task shared(R,vv)
                      removeV(R,vv);
                    #pragma omp task shared(P1,vv)
                      removeV(P1,vv);
                    #pragma omp task shared(P_in,X_in,X,P,i)
                    {
                      X->push_back((*P)[i]);
                      P_in->clear();
                      X_in->clear();
                    }
                  }
                }
              }
            }
          }
        delete u;
        delete P1;
        delete P_in;
        delete X_in;
      } 

    void Bron_KerboschWithPivot2(deque<Vertex> * R, deque<Vertex> * P, deque<Vertex> * X, deque<Vertex> * temp)
      { 
        int vv,vvl;
        int i,j,k,l,pos_vv,cnt;
        Vertex *u;
        deque<int>* localNb;
        deque<Vertex>* P_in;
        deque<Vertex>* X_in;
        deque<Vertex> * P1;

        ++rCount;
        //printf("Recursive call %ld\n",rCount);

        if ((P->size() == 0) && (X->size() == 0))
          {
            // Before keeping a cloque we have to check if it complies
            // with conditions of algorithm 2
            // we have to look at all vetices in <R>
            vv = (*R)[0].getV();
            // First vertex in graph represents G_j, so we must use it to check position of
            // other ones in degeneracy order;
            pos_vv = 0;
            while (graphSorted[pos_vv++] != vv);
            for (i=0;i<R->size();i++)
            {
              localNb = (*R)[i].getNbrs();
              for (j=0;j<(*R)[i].getDegree();j++)
                {
                  vvl = (*localNb)[j];
                  // we have to see if <vvl> is located in degereacy list to the left of <vv_pos>
                  k = 0;
                  while ((k<pos_vv) && ((*localNb)[k++] != vvl));
                  // if k == vv_pos, then <vvl> has higher rank in degen order and no need in extra checks
                  if (k<pos_vv)
                    {
                      // neighbor vertex <vvl> has lower rank then main vertex vv
                      // now we have to see if it is adjacent to all vertices in <P>
                      // We can ahieve that calculating virtual intersection
                      // when each vertex of <R> belongs to neighbors of <vvl>
                      cnt = 0;
                      for (l=0;l<R->size();l++)
                        if ((*R)[l].isAdjacent(vvl)) cnt++;
                      if (cnt == R->size()) return; // it is adjacent to all and we have to finish correctly
                    }
                }             
            }
            if (printState) printClique(R);
            else addClique(R);
            return; 
          } 
        P1 = new deque<Vertex>(P->begin(),P->end());
        P_in = new deque<Vertex>();
        X_in = new deque<Vertex>();
        // Find Pivot 
        get_union(P,X,temp);
        u = new Vertex(getMaxDegreeVertex(temp));

        // P = P / Nbrs(u) 
        removeNbrsIn(P, u); 

        for (int i=0; i<P->size(); i++)
          { 
            int vv = (*P)[i].getV();
            temp->clear();
            getVNbrs(vv, temp);
            #pragma omp parallel
            {
              #pragma omp single
              {
                #pragma omp task shared(R,P,i)
                  R->push_back((*P)[i]);
                #pragma omp task shared(P1,P_in,temp)
                  get_intersect(P1, temp,P_in);
                #pragma omp task shared(X,X_in,temp)
                  get_intersect(X, temp,X_in);
                #pragma omp taskwait
                {
                  temp->clear();
                  Bron_KerboschWithPivot2(R, P_in, X_in,temp); 
                  #pragma omp parallel
                  {
                    #pragma omp task shared(R,vv)
                      removeV(R,vv);
                    #pragma omp task shared(P1,vv)
                      removeV(P1,vv);
                    #pragma omp task shared(P_in,X_in,X,P,i)
                    {
                      X->push_back((*P)[i]);
                      P_in->clear();
                      X_in->clear();
                    }
                  }
                }
              }
            }
          }
        delete u;
        delete P1;
        delete P_in;
        delete X_in;
      } 

    void gSort(deque<Vertex> *g)
      {
         int i,j,min,pos,sz,Vd;
         deque<int> t;

         Vd = g->size();
         for (i=0; i<Vd-1; i++)
           {
             min = (*g)[i].getDegree();
             pos = i;
             #pragma omp parallel for private(j,sz) shared(min,pos,g)
             for (j=i+1; j<Vd; j++)
               {
                  sz = (*g)[j].getDegree();
                  if (sz<min)
                    {
                      min = sz;
                      pos = j;
                    }
               }
             t.assign((*g)[i].getNbrs()->begin(),(*g)[i].getNbrs()->end());
             j=(*g)[i].getV();
             (*g)[i].setNbrs(*(*g)[pos].getNbrs());
             (*g)[i].setV((*g)[pos].getV());
             (*g)[pos].setNbrs(t);
             (*g)[pos].setV(j);
             t.clear();
           }
      }

    Vertex getMaxDegreeVertex(deque<Vertex> *g) { 
        gSort(g); 
        return (*g)[g->size() - 1];
    } 

    int removeFirstVertex(deque<Vertex> *g)
      { 
        int i,j,vv;
        Vertex v = (*g)[0];
        vv = v.getV();
        deque<int> *Nb = v.getNbrs();
        for (i=0; i< v.getDegree(); i++)
          {
            #pragma omp parallel for private(j) shared(i,vv,g,Nb)
            for (j=1; j<g->size();j++)
              if ((*g)[j].getV() == (*Nb)[i])
                (*g)[j].removeNbr(vv);
          }
        g->pop_front(); 
        return vv; 
      } 

    void removeV (deque<Vertex> * L, Vertex x)
      {
        int vv=x.getV();
        int sz = (*L).size();
        bool k = false;
        int j = sz;
        
        #pragma omp parallel for
        for (int i=0;i<sz;i++)
          {
            if (k) continue;
            if ((*L)[i].getV()==vv)
              {
                j = i;
                k = true;
              }
          }
        if (j<sz)
          {
            if (j<sz-1)
              {
                (*L)[j].setV((*L)[sz-1].getV());
                (*L)[j].setNbrs(*((*L)[sz-1].getNbrs()));
              }
            L->pop_back();
          }
      }

    void removeV (deque<Vertex> * L, int x)
      {
        int sz = (*L).size();
        bool k = false;
        int j = sz;
       
        #pragma omp parallel for
        for (int i=0;i<sz;i++)
          {
            if (k) continue;
            if ((*L)[i].getV()==x) { j=i; k = true;}
          }
        if (j<sz)
          {
            if (j<sz-1)
              {
                (*L)[j].setV((*L)[sz-1].getV());
                (*L)[j].setNbrs(*((*L)[sz-1].getNbrs()));
              }
            L->pop_back();
          }
      }

    void Alg1Build()
      {
        int vv,vvn;
        int i=0,j,k,l,localDegree;
        deque<Vertex> * X = new deque<Vertex>();
        deque<Vertex> * R = new deque<Vertex>();
        deque<Vertex> * P = new deque<Vertex>();
        deque<Vertex> * P1 = new deque<Vertex>();
        deque<Vertex>* P_in = new deque<Vertex>();
        deque<Vertex>* X_in = new deque<Vertex>();
        deque<Vertex>* temp = new deque<Vertex>();
        Vertex V;
        deque<int> localNb;
        long int sTime, eTime;
        setNoPrint();
        setHashTable();

       	std::cout<<"---ordre de degenerescence---"<< std::endl;
        sTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
	 std::cout<< " < ";
        while (degen.size() > 0)
          { 
            gSort(&degen);
            vv = removeFirstVertex(&degen); 
            std::cout<<"["<<vv<<"] ";
            graphSorted[i++]=vv;
          } 
	      std::cout<< "> "<<std::endl;

        // we will be generating and processing subgraphs in degeneracy order
        for (j =0; j<nodesCount; j++)
          {
            // we have to process all vertices of initial graph and at least this
            // vertex will become a graph element. But we will need only elements 
            // of the same degree and located after that one in degeneracy list
            vv = graphSorted[j];
            // Create subgraph <G_j> induced by vertex <vv> in variable <P>
            V.setV(vv);
            localDegree = adj[vv].getDegree();
            localNb.assign(adj[vv].getNbrs()->begin(),adj[vv].getNbrs()->end());
            // we have to process all nbrs of initial vertex, but according to degeneracy order and
            // only the ones with the same degree
            for (k=0; k<localDegree;k++)
              {
                vvn = localNb[k];
                // now look if this vertex is located afer <vv> in degeneracy order;
                for (l=j+1; l<nodesCount; l++) if (graphSorted[l] == vvn) break;
                // l == nodesCount if not found
                if ((l<nodesCount) && (adj[vvn].getDegree() >= localDegree)) V.addNbr(vvn);
              }
              // degree should be equal to the actual number of nbrs
              V.setDegree();
              P->push_back(V);
              // initial vertex and was added and we reset <localNb> with a list of
              // valid nbrs since they shoutd be added to the sub-graph with mutual nbrs set
              localNb.clear();
              localNb.assign(V.getNbrs()->begin(),V.getNbrs()->end());
              V.clear();
            for (k=0; k<localNb.size();k++)
              {
                vvn = localNb[k];
                V.setV(vvn);
                // Nbr list will start from initial vertex of G_j
                V.addNbr(vv);
                // add all nbr vertices except myself
                for (l=0; l<localNb.size(); l++) if (localNb[l] != vvn) V.addNbr(localNb[l]);
                // new vertex for G_j with valid nbrs only, adding it 
                V.setDegree();
                P->push_back(V);
                V.clear();
              }          
            // Proceed with Bron-Kerbosch
            if (P->size()>1) // only if we created sub-graph with some nbrs, do B_K
              {
                P1->assign(P->begin(),P->end());
                for (int i=0; i<P->size(); i++)
                  { 
                    int vv = (*P)[i].getV();
                    temp->clear();
                    getVNbrs(vv, temp);
                    #pragma omp parallel
                    {
                      #pragma omp single
                      {
                        #pragma omp task shared(R,P,i)
                          R->push_back((*P)[i]);
                        #pragma omp task shared(P1,P_in,temp)
                          get_intersect(P1, temp,P_in);
                        #pragma omp task shared(X,X_in,temp)
                          get_intersect(X, temp,X_in);
                        #pragma omp taskwait
                        {
                          temp->clear();
                          Bron_KerboschWithPivot(R, P_in, X_in,temp); 
                          #pragma omp parallel
                          {
                            #pragma omp task shared(R,vv)
                              removeV(R,vv);
                            #pragma omp task shared(P1,vv)
                              removeV(P1,vv);
                            #pragma omp task shared(X,P,i,P_in,X_in)
                            {
                              X->push_back((*P)[i]);
                              P_in->clear();
                              X_in->clear();
                            }
                          }
                        }
                      }
                    }
                  }
                // we finished with this G_j and need to take care about "temporary" elements
                // preparing them for B_K of the next G_j
                P->clear();
                R->clear();
                X->clear();
                P1->clear();
                P_in->clear();
                X_in->clear();
                temp->clear();
              }
          }
        eTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
        	std::cout<<"Temps d'execution " << eTime-sTime<<std::endl; //printf("Elapsed time %ld, microseconds\n",eTime-sTime);
        delete temp;
        delete P;
        delete P1;
        delete P_in;
        delete X_in;
        delete R;
        delete X;
        setPrint();
      }

    void Alg2Build()
      {
        int vv,vvn;
        int i=0,j,k,l,localDegree;
        deque<Vertex> * X = new deque<Vertex>();
        deque<Vertex> * R = new deque<Vertex>();
        deque<Vertex> * P = new deque<Vertex>();
        deque<Vertex> * P1 = new deque<Vertex>();
        deque<Vertex> * P_in = new deque<Vertex>();
        deque<Vertex> * X_in = new deque<Vertex>();
        deque<Vertex> * temp = new deque<Vertex>();
        Vertex V;
        deque<int> localNb;
        long int sTime, eTime;
         std::cout<<"Ordre de degenerescence"<<std::endl;
        sTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
       	std::cout<< "<";
        while (degen.size() > 0)
          { 
            gSort(&degen);
            vv = removeFirstVertex(&degen); 
        std::cout<<"["<<vv<<"] ";

            graphSorted[i++]=vv;
          } 

       	      std::cout<< "> "<<std::endl;


        // we will be generating and processing subgraphs in degeneracy order
        for (j =0; j<nodesCount; j++)
          {
            // we have to process all vertices of initial graph and at least this
            // vertex will become a graph element. But we will need only elements 
            // of the same degree and located after that one in degeneracy list
            vv = graphSorted[j];
            // Create subgraph <G_j> induced by vertex <vv> in variable <P>
            V.setV(vv);
            localDegree = adj[vv].getDegree();
            localNb.assign(adj[vv].getNbrs()->begin(),adj[vv].getNbrs()->end());
            // we have to process all nbrs of initial vertex, but according to degeneracy order and
            // only the ones with the same degree
            for (k=0; k<localDegree;k++)
              {
                vvn = localNb[k];

                bool flag = false;
                int pos = nodesCount;
       
                #pragma omp parallel for
                for (int l=j+1;l<nodesCount;l++)
                  {
                    if (flag) continue;
                    if (graphSorted[l] == vvn) { pos=l; flag = true;}
                  }
                // l == nodesCount if not found
                if ((pos<nodesCount) && (adj[vvn].getDegree() >= localDegree)) V.addNbr(vvn);
              }
              // degree should be equal to the actual number of nbrs
              V.setDegree();
              P->push_back(V);
              // initial vertex and was added and we reset <localNb> with a list of
              // valid nbrs since they shoutd be added to the sub-graph with mutual nbrs set
              localNb.clear();
              localNb.assign(V.getNbrs()->begin(),V.getNbrs()->end());
              V.clear();
            for (k=0; k<localNb.size();k++)
              {
                vvn = localNb[k];
                V.setV(vvn);
                // Nbr list will start from initial vertex of G_j
                V.addNbr(vv);
                // add all nbr vertices except myself
                #pragma parallel omp for
                for (int l=0; l<localNb.size(); l++) if (localNb[l] != vvn) V.addNbr(localNb[l]);
                // new vertex for G_j with valid nbrs only, adding it 
                V.setDegree();
                P->push_back(V);
                V.clear();
              }          
            // Proceed with Bron-Kerbosch
            if (P->size()>1) // only if we created sub-graph with some nbrs, do B_K
              {
                P1->assign(P->begin(),P->end());
                for (int i=0; i<P->size(); i++)
                  { 
                    int vv = (*P)[i].getV();
                    temp->clear();
                    getVNbrs(vv, temp);
                    #pragma omp parallel
                    {
                      #pragma omp single
                      {
                        #pragma omp task shared(R,P,i)
                          R->push_back((*P)[i]);
                        #pragma omp task shared(P1,P_in,temp)
                          get_intersect(P1, temp,P_in);
                        #pragma omp task shared(X,X_in,temp)
                          get_intersect(X, temp,X_in);
                        #pragma omp taskwait
                        {
                          temp->clear();
                          Bron_KerboschWithPivot2(R, P_in, X_in,temp); 
                          #pragma omp parallel
                          {
                            #pragma omp task shared(R,vv)
                              removeV(R,vv);
                            #pragma omp task shared(P1,vv)
                              removeV(P1,vv);
                            #pragma omp task shared(X,P,i,P_in,X_in)
                            {
                              X->push_back((*P)[i]);
                              P_in->clear();
                              X_in->clear();
                            }
                          }
                        }
                      }
                    }
                  }
                // we finished with this G_j and need to take care about "temporary" elements
                // preparing them for B_K of the next G_j
                P->clear();
                R->clear();
                X->clear();
                P1->clear();
                P_in->clear();
                X_in->clear();
                temp->clear();
              }
          }
        eTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
	std::cout<<"Temps d'execution "<< eTime-sTime<<std::endl;
        delete temp;
        delete P;
        delete P1;
        delete P_in;
        delete X_in;
        delete R;
        delete X;
      }

    void Bron_KerboschDegenBuild()
      { 
        deque<Vertex> * X = new deque<Vertex>();
        deque<Vertex> * R = new deque<Vertex>();
        deque<Vertex> * P = new deque<Vertex>();
        deque<Vertex> * P1 = new deque<Vertex>();
        deque<Vertex> * P_in;
        deque<Vertex> * X_in;
        deque<Vertex> * temp;
        int vv;
        int is = 0;
        long int sTime, eTime;

      std::cout<<"Ordre de degeneresence"<<std::endl;
        sTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
	std::cout<< "<";
        while (degen.size() > 0)
          { 
            gSort(&degen);
            vv = removeFirstVertex(&degen); 
            std::cout<<"["<<vv<<"] ";
            graphSorted[is++]=vv;
          }

         std::cout<< "> "<<std::endl;

 



        // deque<Vertex> P will represent initial graph in physical degeneracy order
        P->resize(nodesCount);
        P1->resize(nodesCount);
        for (int i =0; i<nodesCount; i++)
          {
            vv = graphSorted[i];
            // Build graph P when physical order of vertices corresponds to
            // the degeneracy order
            (*P)[i].setV(vv);
            (*P)[i].setNbrs(*(adj[vv].getNbrs()));
            (*P1)[i].setV(vv);
            (*P1)[i].setNbrs(*(adj[vv].getNbrs()));
          } 
        P_in = new deque<Vertex>();
        X_in = new deque<Vertex>();
        temp = new deque<Vertex>();
        // Proceed with Bron-Kerbosch
        for (int i=0; i<P->size(); i++)
          { 
            int vv = (*P)[i].getV();
            temp->clear();
            getVNbrs(vv, temp);
            #pragma omp parallel
            {
              #pragma omp single
              {
                #pragma omp task shared(R,P,i)
                  R->push_back((*P)[i]);
                #pragma omp task shared(P1,P_in,temp)
                  get_intersect(P1, temp,P_in);
                #pragma omp task shared(X,X_in,temp)
                  get_intersect(X, temp,X_in);
                #pragma omp taskwait
                {
                  temp->clear();

                  Bron_KerboschWithPivot(R, P_in, X_in,temp);
                  #pragma omp parallel
                  {
                    #pragma omp task shared(R,vv)
                      removeV(R,vv);
                    #pragma omp task shared(P1,vv)
                      removeV(P1,vv);
                    #pragma omp task shared(X,P,i,P_in,X_in)
                    {
                      X->push_back((*P)[i]);
                      P_in->clear();
                      X_in->clear();
                    }
                  }
                }
              }
            }
         }
//}
        eTime = static_cast<long int>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
		std::cout<<"Temps d'execution "<< eTime-sTime<<std::endl;

        delete temp;
        delete P;
        delete P1;
        delete R;
        delete X;
        delete P_in;
        delete X_in;
      } 

    void printGraph(void) { printGraph(&adj); }

    void printGraph(deque<Vertex> *G)
      {
	std::cout<<" nodes - Graph:"<<G->size()<<std::endl;
        for (int i = 0; i<G->size(); i++)
         std::cout<<(*G)[i].getV()<<std::endl;

      }

    void printGraphFull(void) { printGraphFull(&adj); } 

    void printGraphFull(deque<Vertex> *G)
      { 
        deque <int> *Nb;
        int i,j;
	std::cout<<"nodes - Graph full:"<<G->size()<<std::endl;;

        for (i = 0; i<G->size(); i++)
          { 
	   std::cout<<"["<<(*G)[i].getV()<<"] > " ;
            Nb = (*G)[i].getNbrs();
            for (j = 0; j<Nb->size(); j++)
		std::cout<<(*Nb)[j] <<"  >  ";
	       std::cout<<" "<<std::endl;
          }
      } 
  };

void usage (void)
{
  std::cout<<"Usage: ./codep -f file_with_graph_edges_TAB_delimited [-a algorithm_number 0..2] [-s vertex__with_smallest_number] [-t number_of_threads]"<<std::endl;
    std::cout<<"By default it uses generic clique search, vertex numeration starting from zero and 2 parallel threads"<<std::endl;
}

int main(int argc, char **argv)
{
    string f_name;
    // by default only 2 threads for safety
    int alg_sel=0,start_vert=0,num_th=2;
    char c;
    /* Parse command line arguments */ 
    while((c = getopt(argc, argv, "a:f:s:t:")) != -1)
    {
        switch(c)
        {
            case 'a':
                alg_sel = atoi(optarg); 	

                break;
            case 's':
                start_vert = atoi(optarg);
                break;
            case 't':
                {
                  // if we want to set a number of parallel threads, default 2 is not interesting;
                  // we can try to use all available threads, but for experiments and if we do not know
                  // the PC specs the desired Th quantity can't be higher than max
                  int val = atoi(optarg);
                  num_th = omp_get_num_procs();
                  if (val < num_th) num_th = val;
                  break;
                }
            case 'f':
                f_name = optarg;
                break;
            default:
              usage();
              exit(-1);
        }
    }
  omp_set_num_threads(num_th);
  Graph gg(f_name,start_vert);
  gg.printGraphFull();
  switch (alg_sel)
    {
      case 0:
	std::cout<<"Generic search for cliques using degeneracy order"<<std::endl;
    
        gg.Bron_KerboschDegenBuild();
        break;
      case 1:
	std::cout<<"Algorithm 1 for cliques using degeneracy order"<<std::endl;
        gg.Alg1Build();
        gg.printCliques();
        break;
      case 2:
	std::cout<<"Algorithm 2 for cliques using degeneracy order"<<std::endl;
        gg.Alg2Build();
        break;
      default:
	std::cout<<"Invalid selection of algorithm."<<std::endl;
        return(1);
    }
	std::cout<<"Overall number of recursive calls:"<<gg.getRcalls()<<std::endl;

  return 0;
}
