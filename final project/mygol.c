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

void my_find_3_row_nbrs(int i, int j,int n,int* nbr_r, int* nbr_c,int* numnbrs){
  memset(nbr_r,-1,8);
  memset(nbr_c,-1,8);
  *numnbrs=0;
  int pr[8]={-1,-1,-1, 0,0, 1,1,1};
  int pc[8]={-1, 0, 1,-1,1,-1,0,1};
  int l;
  
  for(l=0;l<8;l++){
    nbr_r[*numnbrs]=(i+pr[l]+3)%(3);
    nbr_c[*numnbrs]=(j+pc[l]+n)%(n);
    *numnbrs=*numnbrs+1;      
  }
}

int find_first_row_occupied_nbrs(int *  before_row,
  int * first_row,int* next_row,int j,int n){

  int num_nbrs;
  int nbr_r[8];
  int nbr_c[8];
  //row 0:  before_row (receive from other rank)
  //row 1:  first_row (row 0 in sub gol_board)
  //row 2:  next_row (row 1 in sub gol_board)
  my_find_3_row_nbrs(1,j,n,nbr_r,nbr_c,&num_nbrs);
  int NO=0;
  int l;
  for (l=0;l<num_nbrs;l++){
     if(nbr_r[l] == 0){
        if(before_row[nbr_c[l]] == 1){
            NO++;
        }
     }else if(nbr_r[l] == 1){
         if(first_row[nbr_c[l]] == 1){
            NO++;
         }
     }else if(nbr_r[l] == 2){
        if(next_row[nbr_c[l]] == 1){
           NO++;
        }
     }
  }
  return NO;
}

int find_last_row_occupied_nbrs(int * before_row,
  int * last_row,int* next_row,int j,int n){

  int num_nbrs;
  int nbr_r[8];
  int nbr_c[8];
  //row 0:  before_row (row sub.m - 2 in sub gol_board)
  //row 1:  last_row (row sub.m - 1 in sub gol_board)
  //row 2:  next_row (0 row in next rank sub gol_board)
  my_find_3_row_nbrs(1,j,n,nbr_r,nbr_c,&num_nbrs);
  int NO=0;
  int l;
  for (l=0;l<num_nbrs;l++){
     if(nbr_r[l] == 0){
        if(before_row[nbr_c[l]] == 1){
            NO++;
        }
     }else if(nbr_r[l] == 1){
         if(last_row[nbr_c[l]] == 1){
            NO++;
         }
     }else if(nbr_r[l] == 2){
        if(next_row[nbr_c[l]] == 1){
           NO++;
        }
     }
  }
  return NO;
}

void transition(gol_board* sub_g, int t,int* first_row,int* last_row){
  int i; int j;
  int m = sub_g->m; 
  int n = sub_g->n;
  enum STATUS mystatus;
  enum STATUS newstatus;
  int oc_nbrs;
  int tmp;
  int* row_0 = get_a_row(*sub_g,t,0);
  int* row_m = get_a_row(*sub_g,t,m - 1);
  //calculate i == 0 row
  for (j = 0; j < n;j++)
  {
    tmp = row_0[j];
    if (tmp == 1) mystatus = OCCUPIED;
    else mystatus = EMPTY;
    oc_nbrs=find_first_row_occupied_nbrs(first_row,row_0,last_row,j,n);
    if (mystatus==OCCUPIED)
    {
        if (oc_nbrs==2 || oc_nbrs==3) newstatus=OCCUPIED;
        if (oc_nbrs>=4) newstatus=EMPTY;
        if (oc_nbrs<2) newstatus=EMPTY;
    }
    if (mystatus==EMPTY)
    {
        if (oc_nbrs==3) newstatus=OCCUPIED;
        else newstatus=EMPTY;
    }
    set_status(sub_g,0,j,t+1,newstatus);
  }
  printf("transition part1\n");
   //calculate i == m - 1 row
  for (j = 0; j < n;j++)
  {
    /* code */
    tmp = row_m[j];
    if (tmp == 1) mystatus = OCCUPIED;
    else mystatus = EMPTY;
    oc_nbrs=find_last_row_occupied_nbrs(first_row,row_m,last_row,j,n);
    if (mystatus==OCCUPIED)
    {
        if (oc_nbrs==2 || oc_nbrs==3) newstatus=OCCUPIED;
        if (oc_nbrs>=4) newstatus=EMPTY;
        if (oc_nbrs<2) newstatus=EMPTY;
    }
    if (mystatus==EMPTY)
    {
        if (oc_nbrs==3) newstatus=OCCUPIED;
        else newstatus=EMPTY;
    }
    set_status(sub_g,m-1,j,t+1,newstatus);
  }
  printf("transition part2\n");
  for (i = 1; i < m - 1; i++){
    for (j = 0; j < n;j++){
       mystatus=get_status(sub_g,i,j,t);
       oc_nbrs=num_occupied_nbrs(sub_g,i,j,t);
       if (mystatus==OCCUPIED)
       {
	         if (oc_nbrs==2 || oc_nbrs==3) newstatus=OCCUPIED;
	         if (oc_nbrs>=4) newstatus=EMPTY;
	         if (oc_nbrs<2) newstatus=EMPTY;
       }
      if (mystatus==EMPTY)
      {
	         if (oc_nbrs==3) newstatus=OCCUPIED;
	         else newstatus=EMPTY;
      }

      set_status(sub_g,i,j,t+1,newstatus);
    }
  }
  printf("transition part3\n");
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
   printf("Rank : %d start: %d end: %d\n",rank,m_start,m_end);
   if (rank == 0)
   { 
   	  int * tmp;
   	  int * row_1;
   	  int * row_2;
   	  for (generation_index; generation_index != generation_max; ++generation_index)
   	  {	
        //send the (m_end - 1) row of values to ranK + 1
   	    tmp = get_a_row(sub,generation_index,sub.m - 1);
        printf("Rank: %d send to %d tag: %d\n",rank,rank+1,m_end - 1);
   	    sendTheData(tmp,g_size,rank+1,m_end - 1);
        //send the m_start row of values to numproc - 1
   	    tmp = get_a_row(sub,generation_index,0);
        printf("Rank: %d send to %d tag: %d\n",rank,numproc-1,m_start);
   	    sendTheData(tmp,g_size,numproc-1,m_start);
   	    //get the (g.n - 1) row of values from numproc - 1
        row_1  = receiveTheData(g_size,numproc-1,(g.n - 1));
        printf("Rank: %d receive from %d tag: %d\n",rank,numproc-1,(g.n - 1));
   	    //get the m_end row of data from rank + 1
        row_2  = receiveTheData(g_size,rank+1,m_end); 
        printf("Rank: %d receive from %d tag: %d\n",rank,rank+1,m_end);
   	    //calculate [m_start to m_end) & m_start == 0;
        transition(&sub,generation_index,row_1,row_2);
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
        printf("Rank: %d receive from %d tag: %d\n",rank,rank - 1,m_start - 1);
        //send the (m_end - 1) row of values to rank + 1
  	    tmp = get_a_row(sub,generation_index,sub.m - 1);
        printf("Rank: %d send to %d tag: %d\n",rank,rank + 1,m_end - 1);
   	    sendTheData(tmp,g_size,rank + 1,m_end - 1);
   	    //send the m_start row of values to rank - 1
        tmp = get_a_row(sub,generation_index,0);
        printf("Rank: %d send to %d tag: %d\n",rank,rank - 1,m_start);
        sendTheData(tmp,g_size,rank - 1,m_start);
        //get the (m_end) row of values from rank + 1
        row_2 = receiveTheData(g_size,rank + 1,m_end); 
        printf("Rank: %d receive from %d tag: %d\n",rank,rank + 1,m_end);
        //calculate [m_start to m_end)
        //transition(&sub,generation_index,row_1,row_2);
   	  }
   	    //after X times generation, send the data to master.
   }
   else if(rank == numproc - 1){
   	  int * tmp;
   	  int * row_1;
   	  int * row_2;
   	 for (generation_index; generation_index != generation_max; ++generation_index)
   	 {
       //get the 0 row of value from rank == 0
       row_2 = receiveTheData(g_size,0,0);
       printf("Rank: %d receive from %d tag: %d\n",rank,0,0);
       //get the (m_start - 1) row of value from rank - 1 
   	 	 row_1 = receiveTheData(g_size,rank - 1,m_start - 1);
       printf("Rank: %d receive from %d tag: %d\n",rank,rank - 1,m_start - 1);
       //send the m_start row of data to rank - 1
   	 	 tmp = get_a_row(sub,generation_index,0);
       printf("Rank: %d send to %d tag: %d\n",rank,rank - 1,m_start);
   	 	 sendTheData(tmp,g_size,rank - 1,m_start);
       //sent the m_end - 1 row of data to 0
   	 	 tmp = get_a_row(sub,generation_index,sub.m - 1);
       printf("Rank: %d send to %d tag: %d\n",rank,0,m_end - 1);
   	 	 sendTheData(tmp,g_size,0,m_end - 1);
        //calculate [m_start to m_end) & m_end == g.m
        //transition(&sub,generation_index,row_1,row_2);
   	 }
   	  //after X times generation, send the data to master.
   }
   if (rank == 0)
   {
   		//print the final gol_board
      print_gol_board(&sub);
   }
   MPI_Finalize();
   return 0;
}

