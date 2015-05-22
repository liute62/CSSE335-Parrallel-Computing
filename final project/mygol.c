#include "gol_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void distribution_board(gol_board parent, gol_board* sub,int numproc,int rank,int* m_start,int* m_end){
	int m = parent.m;
	int n = parent.n;
	int mul_m = m / numproc; //100 / 7 = 14
	int mul_n = n / numproc;
	int bal_m = m % numproc; //100 % 7 = 2
	int bal_n = m % numproc;
	if (rank == numproc - 1)
	{
		*m_start = mul_m * rank;
		*m_end = mul_m * (rank+1) + bal_m;
	}else{
		*m_start = mul_m * rank;
		*m_end = mul_m * (rank+1);
	}
	//so proecess rank i calculate matrix m_start to m_end, n_start to n_end
	sub->m = *m_end - *m_start;
	sub->n = parent.n;
	sub->data=malloc(sub->m*sub->n*parent.generations*sizeof(enum STATUS));
	int i,j;
	for (i = 0; i < sub->m ; ++i)
	{
	   for (j = 0; j < sub->n ; ++j)
	   {
			set_status(sub,i,j,0,get_status(&parent,i+*m_start,j,0));
	   }
	}
}

void print_matrix(int rank,gol_board gol,int times){
	int i,j;
	for (i = 0; i < gol.m; ++i)
	{
		for (j = 0; j < gol.n; ++j)
		{
			if(get_status(&gol,i,j,times) == OCCUPIED){
				printf("1 ");
			}else{
				printf("0 ");
			}
		}
		printf("rank: %d\n",rank);
	}
}

int* get_a_row(gol_board gol,int generation,int row){
	int n = gol.n;
	int * result = malloc(n * sizeof(int));
	int i;
	for(i = 0; i != n; i++){
	    if(get_status(&gol,row,i,generation) == OCCUPIED){
	    	result[i] = 1;
	    }else{
	    	result[i] = 0;
	    }
	}
	return result;
}

void sendTheData(int* result,int size,int des,int tag){
  MPI_Send(result,size,MPI_INT,des,tag,MPI_COMM_WORLD); 
}

int* receiveTheData(int size,int des,int tag){
  MPI_Status status;
  int* result = malloc(size * sizeof(int));
  MPI_Recv(result,size,MPI_INT,des,tag, MPI_COMM_WORLD, &status);
  return result;
}

void transition(gol_board* g, int t){
  int i; int j;
  int m = g->m; 
  int n = g->n;
  enum STATUS mystatus;
  enum STATUS newstatus;
  int oc_nbrs;
  for (i = 0; i < m; i++){
    for (j = 0; j < n;j++){
      mystatus=get_status(g,i,j,t);
      oc_nbrs=num_occupied_nbrs(g,i,j,t);
      if (mystatus==OCCUPIED){
	if (oc_nbrs==2 || oc_nbrs==3) newstatus=OCCUPIED;
	if (oc_nbrs>=4) newstatus=EMPTY;
	if (oc_nbrs<2) newstatus=EMPTY;
      }
      if (mystatus==EMPTY){
	if (oc_nbrs==3) newstatus=OCCUPIED;
	else newstatus=EMPTY;
      }
      set_status(g,i,j,t+1,newstatus);
    }

  }

}

/**
 1. Split the input matrix to mutiple sub matrixes and a process manage a sub matrix.
 2. Use transition function in each matrix.
 3. The boudary value of matrix should be sent by other process if neccesary.
 4. Send the result to master after a process finish the final transition and output it.
**/
int main(int argc,char** argv){
   
   int numproc;
   int rank;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numproc);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
   
   gol_board g;
   options o; options defaults;
   strcpy(defaults.inputfile,"pulsar.txt");
   defaults.resolution[0]=512; defaults.resolution[1]=512;
   defaults.generations=100;
   o=parse_commandline(argc,argv,defaults);
  // print_options(o);
   g.generations=o.generations;
   initialize_board(o.inputfile,&g);
  // printf("Rank %d , Loaded a %dx%d board.\n",rank,g.m,g.n);

   int m_start;
   int m_end;
   int generation_max = g.generations;
   int generation_index = 0;
   gol_board sub;
   int g_size = g.n;
   distribution_board(g,&sub,numproc,rank,&m_start,&m_end);
   //print_matrix(rank,sub,0);
   if (rank == 0)
   { 
   	  int * tmp;
   	  int * row_1;
   	  int * row_2;
   	  for (generation_index; generation_index != generation_max; ++generation_index)
   	  {	
   	    //send the (m_end - 1) row of values to ranK + 1
   	    tmp = get_a_row(sub,generation_index,sub.m - 1);
   	    sendTheData(tmp,g_size,rank+1,m_end - 1);
   	    //send the m_start row of values to numproc - 1
   	    tmp = get_a_row(sub,generation_index,0);
   	    sendTheData(tmp,g_size,rank+1,m_start);
   	    //get the (g.n - 1) row of values from numproc - 1
   	    row_1  = receiveTheData(g_size,numproc-1,(g.n - 1));
   	    //get the m_end row of data from rank + 1
   	    row_2  = receiveTheData(g_size,rank+1,m_end); 
   	    //calculate [m_start to m_end) & m_start == 0;

   	  }
   	  
   	   //after X times generation, receive the data from all of processes.
   }
   else if(0 < rank && rank < numproc - 1)
   {
   	  int * tmp;
   	  int * row_1;
   	  int * row_2;
   	  for (generation_index; generation_index != generation_max; ++generation_index)
   	  {	
   	    //get the (m_start - 1) row of values from rank - 1
   	    row_1 = receiveTheData(g_size,rank - 1,m_start - 1);
   	    //send the (m_end - 1) row of values to rank + 1
  	    tmp = get_a_row(sub,generation_index,sub.m - 1);
   	    sendTheData(tmp,g_size,rank + 1,m_end - 1);
   	    //get the (m_end) row of values from rank + 1
   	    row_2 = receiveTheData(gol.n,rank + 1,m_end); 
   	    //send the m_start row of values to rank - 1
   	    tmp = get_a_row(sub,generation_index,0);
   	    sendTheData(tmp,g_size,rank - 1,m_start);
   	    //calculate [m_start to m_end)
   	  }
   	    //after X times generation, send the data to master.
   }
   else if(rank == numproc - 1){
   	  int * tmp;
   	  int * row_1;
   	  int * row_2;
   	 for (generation_index; generation_index != generation_max; ++generation_index)
   	 {
   	    //get the (m_start - 1) row of value from rank - 1 
   	 	row_1 = receiveTheData(gol.n,rank - 1,m_start - 1);
   	 	//get the 0 row of value from rank == 0
   	 	row_2 = receiveTheData(gol.n,0,0);
   	 	//send the m_start row of data to rank - 1
   	 	tmp = get_a_row(sub,generation_index,0);
   	 	sendTheData(tmp,gol.n,rank - 1,m_start);
   	 	//sent the m_end - 1 row of data to 0
   	 	tmp = get_a_row(sub,generation_index,sub.m - 1);
   	 	sendTheData(tmp,gol.n,0,m_end - 1);
   	    //calculate [m_start to m_end) & m_end == g.m
   	 }
   	  //after X times generation, send the data to master.
   	 
   }
   if (rank == 0)
   {
   		//print the final gol_board
   }
   MPI_Finalize();
   return 0;
}

