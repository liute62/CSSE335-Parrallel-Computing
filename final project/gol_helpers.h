#ifndef GOL_HELPERS_H
#define GOL_HELPERS_H

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"gifsave89.c"

enum STATUS{OCCUPIED, EMPTY};

typedef struct options{
  char inputfile[100];
  char outputfile[100];
  int resolution[2];
  int numframes;
  int generations;
}options;

typedef struct gol_board{
  enum STATUS* data;
  int generations;
  int m;
  int n;
}gol_board;

enum STATUS get_status(gol_board* g,int i, int j,int t){
  return g->data[t*(g->m*g->n)+i*g->n+j];
}

void set_status(gol_board* g,int i, int j, int t,enum STATUS s){
  g->data[t*(g->m*g->n)+i*g->n+j]=s;
}

void initialize_board(char* filename, gol_board* g){
  FILE* fp=fopen(filename,"r");
  if (fp==NULL){
    fprintf(stderr,"Error opening %s for read.\n",filename);
    exit(1);
  }

  int i,j;
  if (fscanf(fp,"%d,%d",&g->m,&g->n)!=2){
    fprintf(stderr,"Error reading line 1 of %s.  Expected board size.\n",filename);
    exit(1);
  }

  g->data=malloc(g->m*g->n*g->generations*sizeof(enum STATUS));
  int symb;
  int lineno=0;
  for (i=0;i<g->m;i++){
    lineno++;
    for (j=0;j<g->n;j++){
      if (fscanf(fp,"%d",&symb)!=1){
	fprintf(stderr,"Error reading %s, line %d.\n",filename,lineno);
	exit(1);
      }
      if (symb==0){
	g->data[i*g->n+j]=EMPTY;
      }
      else if(symb==1){
	g->data[i*g->n+j]=OCCUPIED;
      }
      else{
	fprintf(stderr,"Error reading %s, line %d.  Got symbol %d and expected 1 or 0.\n",filename,lineno,symb);
	exit(1);
      }
    }
  }
  

}


void find_nbrs(gol_board* g, int i, int j, int* nbr_r, int* nbr_c,int* numnbrs){
  memset(nbr_r,-1,8);
  memset(nbr_c,-1,8);

  *numnbrs=0;
  int pr[8]={-1,-1,-1, 0,0, 1,1,1};
  int pc[8]={-1, 0, 1,-1,1,-1,0,1};
  int l;
  
  for(l=0;l<8;l++){
    nbr_r[*numnbrs]=(i+pr[l]+g->m)%(g->m);
    nbr_c[*numnbrs]=(j+pc[l]+g->n)%(g->n);
    printf("numbers : %d nbr_r: %d nbr_c: %d \n",*numnbrs, nbr_r[*numnbrs],nbr_c[*numnbrs]);
    *numnbrs=*numnbrs+1;      
  }
  //printf("For reference, -1 % 5 is %d\n", (-1+5) % 5);

}

int num_occupied_nbrs(gol_board* g, int i, int j, int t){

  int num_nbrs;
  int nbr_r[8];
  int nbr_c[8];
  find_nbrs(g,i,j,nbr_r,nbr_c,&num_nbrs);
  int NO=0;
  int l;
  for (l=0;l<num_nbrs;l++){
    if (get_status(g,nbr_r[l],nbr_c[l],t)==OCCUPIED){
      NO++;
    }
  }
  return NO;

}

void print_gol_board(gol_board* g){
  int t,i,j;
  char symb;
  for (t=0;t<g->generations;t++){
    printf("State at time: %d\n",t);
    for (i=0;i<g->m;i++){
      for (j=0;j<g->n;j++){
	symb=get_status(g,i,j,t)==OCCUPIED ? 'O':'E';
	printf("%c ",symb);
      }
      printf("\n");
    }
  }

}

void print_gol_board_annotated(gol_board* g){
  int t,i,j;
  for (t=0;t<g->generations;t++){
    printf("State at time: %d\n",t);
    for (i=0;i<g->m;i++){
      for (j=0;j<g->n;j++){
	printf("%c(%d) ",get_status(g,i,j,t)==OCCUPIED ? 'O':'E',num_occupied_nbrs(g,i,j,t));
      }
      printf("\n");
    }
  }

}

void printusage(options defaults){

  fprintf(stderr,"gol options\n");
  fprintf(stderr,"available options are:\n");
  fprintf(stderr,"-h               Print this message\n");
  fprintf(stderr,"-i=inputfile     Use this input file.  Default is %s\n",defaults.inputfile);
  fprintf(stderr,"-o=outputfile    Use this output file. Default is to derive output filename from input filename.\n");
  fprintf(stderr,"-g=g             Use g generations. Default is %d.\n",defaults.generations);
  fprintf(stderr,"--resolution=mxn Resolution of output gif.  Example usage is\ngol --resolution 512x512\nDefault is %dx%d\n",defaults.resolution[0],defaults.resolution[1]);
  fprintf(stderr,"--numframes=n    Use n frames for animated gif.  Default is one frame per timestep, but that may be too much.\n");

}

options parse_commandline(int argc,char** argv,options defaults){
  options o;
  strcpy(o.inputfile,defaults.inputfile);
  int output_override=0;
  o.generations=defaults.generations;
  o.resolution[0]=defaults.resolution[0];
  o.resolution[1]=defaults.resolution[1];
  int numframes_override=0;
  int i;
  int dummyd1,dummyd2;
  char dummys[100];
  for (i=1;i<argc;i++){
    if (strcmp(argv[i],"-h")==0){
      printusage(defaults);
      exit(0);
    }
    else if (sscanf(argv[i],"-i=%s",dummys)==1){
      strcpy(o.inputfile,dummys);
    }
    else if (sscanf(argv[i],"-o=%s",dummys)==1){
      strcpy(o.outputfile,dummys);
      output_override=1;
    }
    else if (sscanf(argv[i],"-g=%d",&dummyd1)==1){
      o.generations=dummyd1;
    }
    else if(sscanf(argv[i],"--resolution=%dx%d",&dummyd1,&dummyd2)==2){
      o.resolution[0]=dummyd1;o.resolution[1]=dummyd2;
    }
    else if (sscanf(argv[i],"--numframes=%d",&dummyd1)==1){
      o.numframes=dummyd1;
      numframes_override=1;
    }
    else{
      fprintf(stderr,"Error parsing option %s.\n\n\n",argv[i]);
      printusage(defaults);
      exit(1);
    }
  }

  if (!numframes_override){
    o.numframes=o.generations;
  }

  if (!output_override){
    char* lastdot=strrchr(o.inputfile,'.');
    if (lastdot==NULL){
      strcpy(o.outputfile,o.inputfile);
      strcat(o.outputfile,".gif");
    }
    else{
      size_t rootlen=lastdot-o.inputfile;
      //printf("rootlen = %d\n",rootlen);
      strncpy(o.outputfile,o.inputfile,rootlen);
      o.outputfile[rootlen]='\0';
      strcat(o.outputfile,".gif");
    }
  }
  return o;
}

void print_options(options o){

  printf("Input File (initial state): %s\n",o.inputfile);
  printf("Output File: %s\n",o.outputfile);
  printf("Generations: %d\n",o.generations);
  printf("Resolution: %dx%d\n",o.resolution[0],o.resolution[1]);
  printf("Number of frames requested: %d\n",o.numframes);
  printf("Number of frames output: %d\n",(o.generations/(o.generations/o.numframes)));

}


int max(int i, int j){
  return i>j?i:j;
}

int min(int i, int j){
  return i<j?i:j;
}

void insert_rect(unsigned char * pixels, int m, int n, int top_left_r, int top_left_c, int bottom_right_r, int bottom_right_c,int colorindex){
  int i;
  for (i=top_left_r; i<=bottom_right_r;i++){

    memset(pixels+i*n+top_left_c, colorindex, bottom_right_c-top_left_c+1);

  }

}


void insert_v_line(unsigned char* pixels, int m, int n, int c, int t, int b, int thickness, int col){
  if (thickness>0) insert_rect(pixels,m,n,t,c,b,c+(thickness-1),col);

}

void insert_h_line(unsigned char* pixels, int m, int n, int r, int left, int right, int thickness, int col){
  //printf("hline in row %d\n",r);
  if (thickness>0)  insert_rect(pixels,m,n,r,left,r+(thickness-1),right,col);
}

void gol_board_to_gif(gol_board* g, options o){

  if (o.resolution[0]<g->m || o.resolution[1]<g->n){
    fprintf(stderr,"There is no way I can reasonably output this to a file.  You don't even have 1 pixel per cell.\n");
    exit(1);
  }
  int colortable[]={255,255,255,0,0,0,128,128,255,0,0,0,-1};
  int bgindex=4;
  unsigned char* pixels=calloc(o.resolution[0]*o.resolution[1],sizeof(unsigned char));
  int frameskip=max(g->generations/o.numframes,1);
  int i,j,t;

  unsigned char* gifimage=NULL;
  void* gsdata=NULL;
  gsdata = newgif((void**) &gifimage,o.resolution[0],o.resolution[1],colortable,bgindex);

  if (gsdata==NULL){
    fprintf(stderr,"Error creating the gif file. Abort\n");
    exit(2);
  }

  animategif(gsdata,1,0,-1,2);
  int* sqr_hieght=malloc(g->m*sizeof(int));
  int* sqr_width=malloc(g->n*sizeof(int));
  for (i=0;i<g->m;i++){
    sqr_hieght[i]=o.resolution[0]/g->m;
    if (i< (o.resolution[0] % g->m)){
      sqr_hieght[i]++;
    }
  }
  for (i=0;i<g->n;i++){
    sqr_width[i]=o.resolution[1]/g->n;
    if (i<o.resolution[1] % g->n){
      sqr_width[i]++;
    }
  }
  int sqr_top_row[g->m*sizeof(int)];
  int sqr_bottom_row[g->m*sizeof(int)];
  int sqr_left_col[g->n*sizeof(int)];
  int sqr_right_col[g->n*sizeof(int)];

  sqr_top_row[0]=0; sqr_left_col[0]=0;
  sqr_bottom_row[g->m-1]=o.resolution[0]-1;
  sqr_right_col[g->n-1]=o.resolution[1]-1;
  
  for (i=0;i<max(g->n,g->m);i++){
    if (i<g->n)   sqr_right_col[i]=sqr_left_col[i]+sqr_width[i]-1;
    if (i<g->n-1) sqr_left_col[i+1]=sqr_right_col[i]+1;
    if (i<g->m)   sqr_bottom_row[i]=sqr_top_row[i]+sqr_hieght[i]-1;
    if (i<g->m-1) sqr_top_row[i+1]=sqr_bottom_row[i]+1;
  }

  int hline_thickness=o.resolution[0]/(2*(g->m-1))/10;
  int vline_thickness=o.resolution[1]/(2*(g->n-1))/10;
  for (t=0;t<g->generations;t+=frameskip){
    memset(pixels,(unsigned char)0,o.resolution[0]*o.resolution[1]);

    for (i=0;i<g->m;i++){
      for (j=0;j<g->n;j++){
	if (i>0) insert_h_line(pixels, o.resolution[0], o.resolution[1], sqr_top_row[i], sqr_left_col[j], sqr_right_col[j], hline_thickness, 1);

	if (j>0) insert_v_line(pixels,o.resolution[0], o.resolution[1], sqr_left_col[j],sqr_top_row[i],sqr_bottom_row[i],vline_thickness,1);
	if (i<g->m-1) insert_h_line(pixels, o.resolution[0], o.resolution[1], sqr_bottom_row[i]-(hline_thickness-1), sqr_left_col[j], sqr_right_col[j], hline_thickness, 1);
	if (j<g->n-1) insert_v_line(pixels,o.resolution[0], o.resolution[1], sqr_right_col[j]-(vline_thickness-1),sqr_top_row[i],sqr_bottom_row[i],vline_thickness,1);
	if (get_status(g,i,j,t)==OCCUPIED) insert_rect(pixels,o.resolution[0],o.resolution[1],sqr_top_row[i]+hline_thickness,sqr_left_col[j]+vline_thickness,sqr_bottom_row[i]-hline_thickness,sqr_right_col[j]-vline_thickness,2);
      }
    }
    
    putgif(gsdata,pixels);
  }
  int nbytes=endgif(gsdata);

  FILE* fp = fopen(o.outputfile,"wb");
  if (fp==NULL){
    fprintf(stderr,"Error opening %s for writing\n",o.outputfile);
    exit(1);
  }
  fwrite(gifimage,sizeof(unsigned char),nbytes,fp);
  fclose(fp);
}

#endif
