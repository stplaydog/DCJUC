#include "dis_master.h"

/**
 * constructor and destructors
 * **/
DisMaster::DisMaster(char *f){
	file = f;
}

DisMaster::~DisMaster(){
	free(v_idx[0]);
	free(v_idx[1]);
	free(e_idx[0]);
	free(e_idx[1]);
}

/**
 * Initiate the dis_master first
 * **/
void
DisMaster::init_graph(){
	//scan the graph to get the basic information
	scan_graph();
	//read the real content of the graph
	read_graph();
}

/**
 * scan the graph to get some basic information
 * **/
void
DisMaster::scan_graph(){
	//declare variables
	FILE *stream;
	char file_buf[100];
	char str[100];
	int v_num;
	int e_num;
	int g_num;
	int v_id=0;
	int v_to=0;
	int color=0;
	int i,j,sum=0;

	//read the head information, and allocate according variables
	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	v_size = v_num;
	v_idx = (int**)malloc(sizeof(int*)*4);
	e_idx = (int**)malloc(sizeof(int*)*4);
	for(i=0;i<2;i++){
		v_idx[i] = (int*)malloc(sizeof(int)*v_size);
		e_idx[i] = (int*)malloc(sizeof(int)*e_num);
		e_size[i] = 0;
	}
	//scan the real content
	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
			printf("error here\n");
		e_size[color-1] += 1;
		v_idx[color-1][v_id] +=1;
	}
	for(i=0;i<2;i++){
		sum=0;
		for(j=0;j<v_num;j++){
			sum += v_idx[i][j];
			v_idx[i][j]=sum;
		}
	}
	fclose(stream);
}


void
DisMaster::read_graph(){
	FILE *stream;
	char file_buf[100];
	char str[100];
	int v_num;
	int e_num;
	int g_num;
	int v_id=0;
	int v_to=0;
	int color=0;
	int i,j;
	int e_num_read = 0;
	//read basic information
	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	//this is for the start of the different postitions
	int **idx = (int**)malloc(sizeof(int*)*3);
	for(i=0;i<3;i++){
		idx[i]=(int*)malloc(sizeof(int)*v_size);
		for(j=0;j<v_size;j++){
			idx[i][j] = j==0?0:v_idx[i][j-1];
		}
	}
	//real read part
	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
			printf("error here\n");
		int pos = idx[color-1][v_id];
		e_idx[color-1][pos] = v_to;
		idx[color-1][v_id]++;
	}
	//frees
	for(i=0; i<3; i++)
		free(idx[i]);
	free(idx);
	fclose(stream);
}
