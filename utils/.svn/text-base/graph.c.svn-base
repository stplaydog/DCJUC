#include "graph.h"

void
from_gr_to_adj(pg g, pga a);
void
from_adj_to_gr(pg g, pga a);
void
free_adj(pga a);
void
resume_zero(pg g);
void scan_csr_idx(pg g, char *file);
void read_graph_csr(pg g, char *file);
void scan_csr_idx(pg g, char *file)
{
	FILE *stream;
	char file_buf[100];
	char str[100];
	int v_num;
	int e_num;
	int count=0;

	int v_id=0;
	int v_to=0;
	int color=0;
	int i,j,sum=0;
	int e_num_read = 0;

	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d\n", &v_num, &e_num)==EOF)
		printf("error here\n");
	g->v_size = v_num;

	g->v_idx = (int**)malloc(sizeof(int*)*4);
	g->v_deg = (int**)malloc(sizeof(int*)*4);
	g->e_idx = (int**)malloc(sizeof(int*)*4);
	g->e_check = (int**)malloc(sizeof(int*)*4);
	g->e_size = (int*)malloc(sizeof(int)*4);
	g->v_check = (int*)malloc(sizeof(int)*g->v_size);
	g->cycle = 0;
	

	for(i=0;i<4;i++){
		g->v_idx[i] = (int*)malloc(sizeof(int)*g->v_size);
		g->v_deg[i] = (int*)malloc(sizeof(int)*g->v_size);
		g->e_idx[i] = (int*)malloc(sizeof(int)*e_num);
		g->e_check[i] = (int*)malloc(sizeof(int)*e_num);
		g->e_size[i] = 0;
		for(j=0;j<v_num;j++)
			g->v_deg[i][j]=0;
	}
	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
			printf("error here\n");
		//printf("%d %d %d\n", v_id, v_to, color);
		g->e_size[color-1] += 1;
		g->v_deg[color-1][v_id] +=1;
	}
	//for(i=0; i<g->v_size; i++)
	//	printf("%d ", g->v_deg[0][i]);
	//printf("\n");
	for(i=0;i<3;i++){
		sum=0;
		for(j=0;j<v_num;j++){
			sum += g->v_deg[i][j];
			g->v_idx[i][j]=sum;
		}
	}
	fclose(stream);
}



void read_graph_csr(pg g, char *file){
	FILE *stream;
	char file_buf[100];
	char str[100];
	int v_num;
	int e_num;
	int count=0;

	int v_id=0;
	int v_to=0;
	int e_weight=0;
	int i,j,sum=0;
	int e_num_read = 0;

	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d\n", &v_num, &e_num)==EOF)
		printf("error here\n");
	//this is for the start of the different postitions
	int **idx = (int**)malloc(sizeof(int*)*3);	
	for(i=0;i<3;i++){
		idx[i]=(int*)malloc(sizeof(int)*g->v_size);
		for(j=0;j<g->v_size;j++){
			idx[i][j] = j==0?0:g->v_idx[i][j-1];
		}
	}

	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &e_weight)==EOF)
			printf("error here\n");
		int pos = idx[e_weight-1][v_id];
		g->e_idx[e_weight-1][pos] = v_to;
		idx[e_weight-1][v_id]++;
	}
	for(i=0; i<g->v_size; i++)
		g->v_check[i] = TRUE;
	for(i=0; i<3; i++)
		free(idx[i]);
	free(idx);
	fclose(stream);
}

void init_gr(pg g, char *file){
	//init global variables
	int i, j, k;
	scan_csr_idx(g, file);
	read_graph_csr(g, file);
}

int check_gap(pg g, int from){
	int start1 = from==0?0:g->v_idx[0][from-1];
	int start2 = from==0?0:g->v_idx[1][from-1];
	int start3 = from==0?0:g->v_idx[2][from-1];
	int end1 = g->v_idx[0][from];
	int end2 = g->v_idx[1][from];
	int end3 = g->v_idx[2][from];
	int gap1 = end1-start1;
	int gap2 = end2-start2;
	int gap3 = end3-start3;
	
	if(gap1==1 &&gap2==1 && gap3==1)
		return 3;
	else if(gap1==1 &&gap2==1 && gap3==0)
		return 2;
	else if(gap1==1 &&gap2==0 && gap3==1)
		return 2;
	else if(gap1==0 &&gap2==1 && gap3==1)
		return 2;
	return FALSE;
}

int check_another(pg g, int another){
	if(g->v_deg[0][another]==1 && g->v_deg[1][another]==1 && g->v_deg[2][another]==0 ||
		g->v_deg[0][another]==1 && g->v_deg[1][another]==0 && g->v_deg[2][another]==1 ||
		g->v_deg[0][another]==0 && g->v_deg[1][another]==1 && g->v_deg[2][another]==1)
		return TRUE;
	else
		return FALSE;
}

void rename_graph_by_car(pg g) {
	int i,j,x,y,z;
	int color;
	for(i=0;i<g->v_size;i++)
		g->v_check[i] = FALSE;
	pga a = (pga)malloc(sizeof(tga));
	from_gr_to_adj(g, a);
	g->cycle = 0;
	g->zero = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0; i<g->v_size; i++)
		g->zero[i]=-1;
	//graph_vis(g, "before.dot");
	//detect ASs
	for(i=0; i<a->v_size; i++){
		if(a->idx[0][i]==0 || a->idx[1][i]==0 || a->idx[2][i]==0)
			continue;
		if(a->idx[0][i]>1 || a->idx[1][i]>1 || a->idx[2][i]>1)
			continue;
		if(a->adj[0][i][0] == a->adj[1][i][0] && a->adj[1][i][0] == a->adj[2][i][0]){
			int to = a->adj[0][i][0];
			if(a->idx[0][to]>1 || a->idx[1][to]>1 || a->idx[2][to]>1)
				continue;
			a->idx[0][i] = a->idx[1][i] = a->idx[2][i] =0;
			a->idx[0][to] = a->idx[1][to] = a->idx[2][to] =0;
			g->cycle += 3;
			g->zero[i] = a->adj[0][i][0];
			g->zero[a->adj[0][i][0]] = i;
		}
	}
	for(color=0; color<3; color++){
		int c1, c2, c3;
		if(color==0){c1=0; c2=1; c3=2;}
		if(color==1){c1=0; c2=2; c3=1;}
		if(color==2){c1=1; c2=2; c3=0;}
		for(i=0; i<a->v_size; i++){
			if(a->idx[c1][i]>1 || a->idx[c2][i]>1 || a->idx[c3][i]>1)
				continue;
			if(a->idx[c1][i]==0 || a->idx[c2][i]==0)
				continue;
			int to1 = a->adj[c1][i][0]; int to2 = a->adj[c2][i][0]; int to3 = a->adj[c3][i][0];
			if(to1 == CAP || to2 == CAP)
				continue;
			if(to1 == to2){
				if(a->idx[c1][to1]>1 || a->idx[c2][to1]>1 || a->idx[c3][to1]>1)
					continue;
				if((a->idx[c3][i]==0 && a->idx[c3][to1]!=0) ||
						(a->idx[c3][i]!=0 && a->idx[c3][to1]==0))
					continue;
				if(a->idx[c3][i]==0 && a->idx[c3][to1]==0){
					a->idx[c1][i] = 0; a->idx[c2][i] = 0; a->idx[c3][i] = 0;
					a->idx[c1][to1] = 0; a->idx[c2][to1] = 0; a->idx[c3][to1] = 0;
					g->cycle+=2;
					g->zero[i] = to1;
					g->zero[to1] = i;
					continue;
				}
				int neighbor = a->adj[c3][i][0];
				int to_neighbor = a->adj[c3][to1][0];
				//printf("v %d neigbor %d to %d to_neighbor\n", i, neighbor, to1, to_neighbor);
				if(neighbor==CAP || to_neighbor==CAP)
					continue;
				for(j=0; j<a->idx[c3][neighbor]; j++)
					if(a->adj[c3][neighbor][j] == i)
						a->adj[c3][neighbor][j] = to_neighbor;
				for(j=0; j<a->idx[c3][to_neighbor]; j++)
					if(a->adj[c3][to_neighbor][j] == to1){
						a->adj[c3][to_neighbor][j] = neighbor;
					}
				a->idx[c1][i] = 0; a->idx[c2][i] = 0; a->idx[c3][i] = 0;
				a->idx[c1][to1] = 0; a->idx[c2][to1] = 0; a->idx[c3][to1] = 0;
				g->cycle+=2;
				g->zero[i] = to1;
				g->zero[to1] = i;
			}
		}
	}
	from_adj_to_gr(g, a);
	//graph_vis(g, "after.dot");
	free_adj(a);
	//for(i=0; i<g->v_size;i++){
	//	int gap_result = check_gap(g, i);
	//	if(g->v_check[i]==1 && gap_result!=FALSE){
	//		int found = FALSE;
	//		int from = i;
	//		int to;
	//		int start1 = from==0?0:g->v_idx[0][from-1];
	//		int start2 = from==0?0:g->v_idx[1][from-1];
	//		int start3 = from==0?0:g->v_idx[2][from-1];
	//		int end1 = g->v_idx[0][from];
	//		int end2 = g->v_idx[1][from];
	//		int end3 = g->v_idx[2][from];
	//		int gap1 = end1-start1;
	//		int gap2 = end2-start2;
	//		int gap3 = end3-start3;
	//		if(g->e_idx[0][start1] == g->e_idx[1][start2] && g->e_idx[0][start1]!=65533){
	//			int another = g->e_idx[0][start1];
	//			int start_another = another==0?0:g->v_idx[2][another-1];
	//			int end_another = g->v_idx[2][another];
	//			int gap_another = check_gap(g, another);
	//			if(gap_another ==3 && gap_result ==3){
	//				//this is the edge shrinking step
	//				int one = g->e_idx[2][start3];
	//				int two = g->e_idx[2][start_another];
	//				int start_one = one==0?0:g->v_idx[2][one-1];
	//				int start_two = two==0?0:g->v_idx[2][two-1];
	//				int end_one = g->v_idx[2][one];
	//				int end_two = g->v_idx[2][two];
	//				for(j=start_one; j<end_one; j++)
	//					if(g->e_idx[2][j] == from)
	//						g->e_idx[2][j] = two;
	//				for(j=start_two; j<end_two; j++)
	//					if(g->e_idx[2][j] == another)
	//						g->e_idx[2][j] = one;

	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//			else if(check_another(g, another)==TRUE && gap3 == 0){
	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//		}
	//		else if(g->e_idx[0][start1] == g->e_idx[2][start3] && g->e_idx[0][start1]!=65533){
	//			int another = g->e_idx[0][start1];
	//			int start_another = another==0?0:g->v_idx[1][another-1];
	//			int end_another = g->v_idx[1][another];
	//			int gap_another = check_gap(g, another);
	//			if(gap_another==3 && gap_result == 3){
	//				//this is the edge shrinking step
	//				int one = g->e_idx[1][start2];
	//				int two = g->e_idx[1][start_another];
	//				int start_one = one==0?0:g->v_idx[1][one-1];
	//				int start_two = two==0?0:g->v_idx[1][two-1];
	//				int end_one = g->v_idx[1][one];
	//				int end_two = g->v_idx[1][two];
	//				for(j=start_one; j<end_one; j++)
	//					if(g->e_idx[1][j] == from)
	//						g->e_idx[1][j] = two;
	//				for(j=start_two; j<end_two; j++)
	//					if(g->e_idx[1][j] == another)
	//						g->e_idx[1][j] = one;

	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//			else if(check_another(g, another)==TRUE && gap2 == 0){
	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//		}
	//		else if(g->e_idx[2][start3] == g->e_idx[1][start2] && g->e_idx[1][start2]!=65533){
	//			int another = g->e_idx[2][start3];
	//			int start_another = another==0?0:g->v_idx[0][another-1];
	//			int end_another = g->v_idx[0][another];
	//			int gap_another = check_gap(g, another);
	//			if(gap_another==3 && gap_result == 3){
	//				//this is the edge shrinking step
	//				int one = g->e_idx[0][start1];
	//				int two = g->e_idx[0][start_another];
	//				int start_one = one==0?0:g->v_idx[0][one-1];
	//				int start_two = two==0?0:g->v_idx[0][two-1];
	//				int end_one = g->v_idx[0][one];
	//				int end_two = g->v_idx[0][two];
	//				for(j=start_one; j<end_one; j++)
	//					if(g->e_idx[0][j] == from)
	//						g->e_idx[0][j] = two;
	//				for(j=start_two; j<end_two; j++)
	//					if(g->e_idx[0][j] == another)
	//						g->e_idx[0][j] = one;

	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//			else if(check_another(g, another)==TRUE && gap1 == 0){
	//				g->v_check[from] = NUL;
	//				g->v_check[another] = NUL;
	//			}
	//		}
	//	}
	//}
}


void from_gr_to_adj(pg g, pga a){
	int i, j, k;
	a->adj = (int***)malloc(sizeof(int**)*3);
	a->idx = (int**)malloc(sizeof(int*)*3);
	a->e_size = (int*)malloc(sizeof(int)*3);
	a->v_size = g->v_size;
	a->c_size = 3;
	for(i=0; i<3; i++){
		a->adj[i] = (int**)malloc(sizeof(int*)*g->v_size);
		a->idx[i] = (int*)malloc(sizeof(int)*g->v_size);
		a->e_size[i] = g->e_size[i];
		for(j=0; j<g->v_size; j++){
			a->adj[i][j] = (int*)malloc(sizeof(int)*g->v_size);
			a->idx[i][j] = 0;
		}
	}
	for(i=0; i<3; i++){
		for(j=0; j<g->v_size; j++){
			int start = j==0?0:g->v_idx[i][j-1];
			int end = g->v_idx[i][j];
			for(k=start; k<end; k++){
				int from = j;
				int to = g->e_idx[i][k];
				if(to != CAP){
					a->adj[i][from][a->idx[i][from]++] = to;
				}
			}
		}
	}
}

void from_adj_to_gr(pg g, pga a){
	int i, j, k;
	//v_idx
	for(i=0; i<a->c_size; i++){
		int sum=0;
		for(j=0; j<a->v_size; j++){
			sum += a->idx[i][j];
			g->v_idx[i][j] = sum;
			g->v_deg[i][j] = a->idx[i][j];
		}
	}
	//e_idx
	for(i=0; i<a->c_size; i++){
		int pos =0;
		for(j=0; j<a->v_size; j++){
			for(k=0; k<a->idx[i][j]; k++){
				g->e_idx[i][pos++] = a->adj[i][j][k];
			}
		}
	}
}

void
free_adj(pga a){
	int i, j, k;
	for(i=0; i<a->c_size; i++){
		for(j=0; j<a->v_size; j++)
			free(a->adj[i][j]);
		free(a->adj[i]);
		free(a->idx[i]);
	}
	free(a->e_size);
	free(a->adj);
	free(a->idx);
	free(a);
}

void
resume_zero(pg g){
	int i, j, k;
	int **adj = (int**)malloc(sizeof(int*)*g->v_size);	
	int *idx = (int*)malloc(sizeof(int)*g->v_size); 
	for(j=0; j<g->v_size; j++){
		idx[j] = 0;
		adj[j] = (int*)malloc(sizeof(int)*g->v_size);
		int start = j==0?0:g->v_idx[3][j-1];
		int end = g->v_idx[3][j-1];
		for(k=start; k<end; k++){
			int from = j;
			int to = g->e_idx[3][k];
			adj[from][idx[from]++] = to;
			adj[to][idx[to]++] = from;
		}
	}
	for(i=0; i<g->v_size; i++){
		if(g->zero[i] != -1){
			int from = i;
			int to = g->zero[i];
			adj[from][idx[from]++] = to;
			adj[to][idx[to]++] = from;
		}
	}
	//v_idx
	int sum=0;
	for(j=0; j<g->v_size; j++){
		sum += idx[j];
		g->v_idx[3][j] = sum;
	}
	//e_idx
	int pos =0;
	for(j=0; j<g->v_size; j++){
		for(k=0; k<idx[j]; k++){
			g->e_idx[3][pos++] = adj[j][k];
		}
	}
	for(j=0; j<g->v_size; j++)
		free(adj[j]);
	free(adj);
	free(idx);
}

void prepare_dis(pg g){
	int color,i;
	for(color=0; color<4; color++){
		for(i=0;i<g->e_size[color]; i++)
			g->e_check[color][i] = FALSE;
	}
}

void graph_vis(pg g, char *file){
	int i=0,j,from;
	FILE *writer = fopen(file, "w");
	fprintf(writer, "graph G{\n");
	int color;
	for(color=0;color<3;color++){
		char cstr[10];
		switch (color){
			case 0:
				sprintf(cstr, "red");
				break;
			case 1:
				sprintf(cstr, "blue");
				break;
			case 2:
				sprintf(cstr, "green");
				break;
			case 3:
				sprintf(cstr, "black");
				break;
		}
		for(from=0;from<g->v_size; from++){
			if(g->v_check[from] != NUL){
				int start = from==0?0:g->v_idx[color][from-1];
				int end = g->v_idx[color][from];
				for(j=start;j<end;j++){
					int to = g->e_idx[color][j];
					//if(from<=to){
						fprintf(writer, "%d -- %d [color=%s];\n", from, to, cstr);
					//}
				}
			}
		}
		//for(i=0;i<g->v_size; i++){
		//	int start = i==0?0:g->v_idx[color][i-1];
		//	int end = g->v_idx[color][i];
		//	for(j=start;j<end;j++)
		//		printf("%d|%d ", i, g->e_idx[color][j]);
		//	printf("%d ", i);
		//}
		//printf("\n");
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void graph_vis_two(pg g, char *file, int c1, int c2){
	int i=0,j,from;
	FILE *writer = fopen(file, "w");
	fprintf(writer, "graph G{\n");
	for(from=0;from<g->v_size; from++){
		//if(g->v_check[from] != NUL){
		int start = from==0?0:g->v_idx[c1][from-1];
		int end = g->v_idx[c1][from];
		for(j=start;j<end;j++){
			int to = g->e_idx[c1][j];
			fprintf(writer, "%d -- %d [color=red];\n", from, to);
		}   
		//}
	}
	for(from=0;from<g->v_size; from++){
		//if(g->v_check[from] != NUL){
		int start = from==0?0:g->v_idx[c2][from-1];
		int end = g->v_idx[c2][from];
		for(j=start;j<end;j++){
			int to = g->e_idx[c2][j];
			fprintf(writer, "%d -- %d [color=blue];\n", from, to);
		}   
		//}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

/*****************************************************************
 * this function is for reducing the search space for exemplar distance
 * computation.
 * ***************************************************************/
void prepare_graph_for_distance(pg g1, pg g2){
	g2->v_size = g1->v_size;
	g2->e_size = (int*)malloc(sizeof(int)*4);
	g2->v_check = (int*)malloc(sizeof(int)*g2->v_size);
	g2->v_idx = (int**)malloc(sizeof(int*)*4);
	g2->v_deg = (int**)malloc(sizeof(int*)*4);
	g2->e_idx = (int**)malloc(sizeof(int*)*4);
	g2->e_check = (int**)malloc(sizeof(int*)*4);
	//get v_deg
	int color, i,j,k;
	for(color=0;color<4;color++){
		g2->v_idx[color] = (int*)malloc(sizeof(int)*g2->v_size);
		g2->v_deg[color] = (int*)malloc(sizeof(int)*g2->v_size);
		for(i=0;i<g1->v_size;i++){
			g2->v_deg[color][i]=0;
			int start = i==0?0:g1->v_idx[color][i-1];
			int end = g1->v_idx[color][i];
			for(j=start;j<end;j++){
				int is_dup = FALSE;
				if(g1->e_idx[color][j]==i)
					is_dup = TRUE;
				for(k=start;k<j;k++)
					if(g1->e_idx[color][j]==g1->e_idx[color][k])
						is_dup = TRUE;
				if(is_dup==FALSE)
					g2->v_deg[color][i]++;		
			}
		}
	}
	//printf("--------------------\n");
	//for(i=0;i<g2->v_size;i++)
	//	printf("%d: %d %d\n", i, g2->v_deg[1][i], g1->v_deg[1][i]);
	//printf("\n");
	//get v_idx
	int **tmp = (int**)malloc(sizeof(int*)*4);
	for(color=0;color<4;color++){
		int sum=0;
		for(i=0;i<g2->v_size;i++){
			sum+=g2->v_deg[color][i];
			g2->v_idx[color][i] = sum;
		}
		g2->e_size[color]=sum;
		g2->e_idx[color] = (int*)malloc(sizeof(int)*sum);
		g2->e_check[color] = (int*)malloc(sizeof(int)*sum);
		tmp[color] = (int*)malloc(sizeof(int)*g2->v_size);
		for(i=0;i<g2->v_size;i++)
			tmp[color][i] = i==0?0:g2->v_idx[color][i-1];
	}
	//copy back
	for(color=0;color<4;color++){
		for(i=0;i<g1->v_size;i++){
			int start = i==0?0:g1->v_idx[color][i-1];
			int end = g1->v_idx[color][i];
			for(j=start;j<end;j++){
				int is_dup = FALSE;
				for(k=start;k<j;k++)
					if(g1->e_idx[color][j]==g1->e_idx[color][k])
						is_dup = TRUE;
				if(is_dup!=TRUE){
					int pos = tmp[color][i];
					g2->e_idx[color][pos] = g1->e_idx[color][j];		
					tmp[color][i]++;
				}
			}
		}
	}
	//copy v_check
	for(i=0;i<g1->v_size;i++)
		g2->v_check[i] = g1->v_check[i];
	free(tmp[0]);
	free(tmp[1]);
	free(tmp[2]);
	free(tmp[3]);
	free(tmp);
}

/*****************************************************************
 * this function is for reducing the search space for exemplar distance
 * computation. this is the version used in the nj dis_mat initialization
 * ***************************************************************/
void prepare_graph_for_distance_two(pg g1, pg g2){
	g2->v_size = g1->v_size;
	g2->e_size = (int*)malloc(sizeof(int)*2);
	g2->v_check = (int*)malloc(sizeof(int)*g2->v_size);
	g2->v_idx = (int**)malloc(sizeof(int*)*2);
	g2->v_deg = (int**)malloc(sizeof(int*)*2);
	g2->e_idx = (int**)malloc(sizeof(int*)*2);
	g2->e_check = (int**)malloc(sizeof(int*)*2);
	//get v_deg
	int color, i,j,k;
	for(color=0;color<2;color++){
		g2->v_idx[color] = (int*)malloc(sizeof(int)*g2->v_size);
		g2->v_deg[color] = (int*)malloc(sizeof(int)*g2->v_size);
		for(i=0;i<g1->v_size;i++){
			g2->v_deg[color][i]=0;
			int start = i==0?0:g1->v_idx[color][i-1];
			int end = g1->v_idx[color][i];
			for(j=start;j<end;j++){
				int is_dup = FALSE;
				if(g1->e_idx[color][j]==i)
					is_dup = TRUE;
				for(k=start;k<j;k++)
					if(g1->e_idx[color][j]==g1->e_idx[color][k])
						is_dup = TRUE;
				if(is_dup==FALSE)
					g2->v_deg[color][i]++;		
			}
		}
	}
	//get v_idx
	int **tmp = (int**)malloc(sizeof(int*)*2);
	for(color=0;color<2;color++){
		int sum=0;
		for(i=0;i<g2->v_size;i++){
			sum+=g2->v_deg[color][i];
			g2->v_idx[color][i] = sum;
		}
		g2->e_size[color]=sum;
		g2->e_idx[color] = (int*)malloc(sizeof(int)*sum);
		g2->e_check[color] = (int*)malloc(sizeof(int)*sum);
		tmp[color] = (int*)malloc(sizeof(int)*g2->v_size);
		for(i=0;i<g2->v_size;i++)
			tmp[color][i] = i==0?0:g2->v_idx[color][i-1];
	}
	//copy back
	for(color=0;color<2;color++){
		for(i=0;i<g1->v_size;i++){
			int start = i==0?0:g1->v_idx[color][i-1];
			int end = g1->v_idx[color][i];
			for(j=start;j<end;j++){
				int is_dup = FALSE;
				for(k=start;k<j;k++)
					if(g1->e_idx[color][j]==g1->e_idx[color][k])
						is_dup = TRUE;
				if(is_dup!=TRUE){
					int pos = tmp[color][i];
					g2->e_idx[color][pos] = g1->e_idx[color][j];		
					tmp[color][i]++;
				}
			}
		}
	}
	//copy v_check
	for(i=0;i<g1->v_size;i++)
		g2->v_check[i] = g1->v_check[i];
	free(tmp[0]);
	free(tmp[1]);
	free(tmp);
}

void 
free_graph(pg g){
	int i=0;
	free(g->v_check);
	free(g->e_size);
	for(i=0; i<4; i++){
		free(g->v_idx[i]);
		free(g->v_deg[i]);
		free(g->e_idx[i]);
		free(g->e_check[i]);
	}
	free(g->v_idx);
	free(g->v_deg);
	free(g->e_idx);
	free(g->e_check);
	free(g);
}

void 
free_graph_two(pg g){
	int i=0;
	free(g->v_check);
	free(g->e_size);
	for(i=0; i<2; i++){
		free(g->v_idx[i]);
		free(g->v_deg[i]);
		free(g->e_idx[i]);
		free(g->e_check[i]);
	}
	free(g->v_idx);
	free(g->v_deg);
	free(g->e_idx);
	free(g->e_check);
	free(g);
}

void 
indentify_comps(int** comp_id, int **comp_type, pg g){
	int color;
	int i, j, k;
	int *ambi_list = (int*)malloc(sizeof(int)*g->v_size);
	int ambi_idx = 0;
	int *visited = (int*)malloc(sizeof(int)*g->v_size);
	for(color=0; color<3; color++){
		for(i=0; i<g->v_size; i++)
			visited[i] = FALSE;
		int comp_num = 0;
		for(i=0; i<g->v_size; i++){
			if(g->v_check[i] == NUL || visited[i] == TRUE){
				visited[i] = TRUE;
				continue;
			}
			comp_id[color][i] = comp_num;
			comp_type[color][comp_num]= TYPE_SING;
			ambi_list[ambi_idx] = i;
			ambi_idx++;
			visited[i] = TRUE;
			//bfs to find connected components
			while(ambi_idx > 0){
				//pop bfs list
				int v = ambi_list[ambi_idx-1];
				ambi_idx--;
				//
				int start_one = v==0?0:g->v_idx[color][v-1];
				int end_one = g->v_idx[color][v];
				int deg_one = end_one - start_one;
				if(deg_one>1)
					comp_type[color][comp_num]= TYPE_MUL;
				for(j=start_one; j<end_one; j++){
					int to = g->e_idx[color][j];
					if(to != NUL && to != CAP && visited[to]==FALSE){
						comp_id[color][to] = comp_num;
						ambi_list[ambi_idx] = to;
						ambi_idx++;
						visited[to] = TRUE;
					}
					else{
						comp_id[color][to] == -1;
					}
				}
				int start_two = v==0?0:g->v_idx[3][v-1];
				int end_two = g->v_idx[3][v];
				int deg_two = end_two - start_two;
				if(deg_two>1)
					comp_type[color][comp_num]= TYPE_MUL;
				for(j=start_two; j<end_two; j++){
					int to = g->e_idx[3][j];
					if(to != NUL && to != CAP && visited[to] == FALSE){
						comp_id[color][to] = comp_num;
						ambi_list[ambi_idx] = to;
						ambi_idx++;
						visited[to] = TRUE;
					}
					else{
						comp_id[color][to] == -1;
					}
				}
			}
			comp_num++;
		}
	}
	free(ambi_list);
	free(visited);
}
