#include "graph.h"
#include "tree.h"
#include "nj.h"
#include <unistd.h>
#include <getopt.h>

int total_gap;
static int verbose_flag;

int main(int argc, char* argv[]){
	////read the file and init the graph	
	int compute_median = FALSE;
	int compute_tree = FALSE;
	int simulate_data = FALSE;
	int is_experiment = FALSE;
	int compute_super_tree = FALSE;
	int is_test_qr = FALSE;
	int num_gene = 0;
	int num_gnm = 0;
	double theta = 0;
	double gamma = 0;
	double phi = 0;
	char *file;
	int c;
	while (1){
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"verbose", no_argument,       &verbose_flag, 1},
			{"brief",   no_argument,       &verbose_flag, 0},
			/* These options don't set a flag.
			   We distinguish them by their indices. */
			{"dist",     no_argument,       0, 'e'},
			{"median",  no_argument,       0, 'm'},
			{"sim",  no_argument,       0, 's'},
			{"tree",  no_argument,       0, 't'},
			{"super",  no_argument,       0, 'b'},
			{"testqr",  no_argument,       0, 'q'},
			{"file",  required_argument, 0, 'f'},
			{"theta",  required_argument, 0, 'h'},
			{"gamma",    required_argument, 0, 'g'},
			{"phi",    required_argument, 0, 'p'},
			{"ngnm",    required_argument, 0, 'n'},
			{"ngene",    required_argument, 0, 'l'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long (argc, argv, "bemtsqf:h:g:p:l:n:", long_options, &option_index);
		if (c == -1)
			break;
		switch (c)
		{
			case 'e':
				is_experiment = TRUE;
				break;
			case 'b':
				compute_super_tree = TRUE;
				break;
			case 'm':
				compute_median = TRUE;
				break;
			case 't':
				compute_tree = TRUE;
				break;
			case 's':
				simulate_data = TRUE;
				break;
			case 'q':
				is_test_qr = TRUE;
				break;
			case 'f':
				file = optarg;
				break;
			case 'l':
				num_gene = atoi(optarg);
				break;
			case 'n':
				num_gnm = atoi(optarg);
				break;
			case 'h':
				theta = atof(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
			case 'p':
				phi = atof(optarg);
				break;
			case '?':
				fprintf (stderr,
						"Unknown option character `\\x%x'.\n",
						optopt);
				return 1;
			default:
				abort ();
		}
	}
	if(compute_median == TRUE){
		//init_for_dis();
		char f_buf[100];
		pg g = (pg)malloc(sizeof(tg));	
		init_gr(g, file);
		rename_graph_by_car(g);
		//get the initial median
		init_median(g);
		//graph_vis_two(g, "before.dot", 0, 3);
		//start local search
		lk_opt(g, 2, 3, FALSE);
		graph_vis_two(g, "after.dot", 0, 3);
		free_graph(g);
	}
	if(compute_tree == TRUE){
		//this is to compute the upper bound for the tree
		//run_example();
		int upper_bound = run_nj(file);
		//int upper_bound = 1;
		printf("the upperbound is %d \n", upper_bound);
		po o = (po)malloc(sizeof(to));	
		init_orders(file, o);
		graph_vis_orders_two(o, "order.dot", 0, 1);
		//int i;
		//for(i=0; i<o->num_genome; i++)
		//	printf("%d ", o->valid[i]);
		//printf("\n");
		int lower_bound = run_lower_bound(o);
		printf("after neighbor joinning the upper_bound is %d the lower_bound is %d\n", upper_bound, lower_bound);
		//this is to find the maximum parsimony tree
		bnb_mp_tree(o, upper_bound, lower_bound);
	}
	if(compute_super_tree == TRUE){
		run_super_tree(file);
	}
	if(is_test_qr == TRUE){
		pdism d_mat = (pdism)malloc(sizeof(tdism));
		init_dis_mat(d_mat, 4);
		d_mat->q_mat[0][0] = 52;
		d_mat->q_mat[0][1] = 30;
		d_mat->q_mat[0][2] = 49;
		d_mat->q_mat[0][3] = 28;
		d_mat->q_mat[1][0] = 30;
		d_mat->q_mat[1][1] = 50;
		d_mat->q_mat[1][2] = 8;
		d_mat->q_mat[1][3] = 44;
		d_mat->q_mat[2][0] = 49;
		d_mat->q_mat[2][1] = 8;
		d_mat->q_mat[2][2] = 46;
		d_mat->q_mat[2][3] = 16;
		d_mat->q_mat[3][0] = 28;
		d_mat->q_mat[3][1] = 44;
		d_mat->q_mat[3][2] = 16;
		d_mat->q_mat[3][3] = 22;
		//householder(d_mat->dis_mat, d_mat->num_spc);
		//init_dis_mat(d_mat, 3);
		//d_mat->q_mat[0][0] = 12;
		//d_mat->q_mat[0][1] = -51;
		//d_mat->q_mat[0][2] = 4;
		//d_mat->q_mat[1][0] = 6;
		//d_mat->q_mat[1][1] = 167;
		//d_mat->q_mat[1][2] = -68;
		//d_mat->q_mat[2][0] = -4;
		//d_mat->q_mat[2][1] = 24;
		//d_mat->q_mat[2][2] = -41;
		double *vector = (double*)malloc(sizeof(double)*4);
		compute_second_eigen_vector(d_mat, vector);	
	}
	if(simulate_data == TRUE){
		char outfile[100];
		sprintf(outfile, "data/sim/%d_%d_%3.2f_%3.2f_%3.2f", num_gnm, num_gene, theta, gamma, phi);
		printf("start generating simulated file to %s\n", outfile);
		run_simulate(outfile, num_gene, num_gnm, theta, gamma, phi);
	}
	if(is_experiment == TRUE){
		total_gap = 0;
		pg g = (pg)malloc(sizeof(tg));	
		printf("start computing distance\n");
		cpec_mem_save(g, 0, 1, 100);
		printf("end computing distance\n");
		free(g);
	}

	return 0;	
}
