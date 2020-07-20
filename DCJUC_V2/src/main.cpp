#include "list.h"
#include "insknap.h"
#include "insmed.h"
#include "insmc.h"
#include "listint.h"
#include "instance.h"
#include "plist.h"
#include "insdis.h"
#include <getopt.h>


static int verbose_flag;
void
print_help();
int 
read_opt(char *opt_file);
int
ede1 ( int invdist, int ngene );

int parallel_mode = -1;

int which_med;


#define DIS_MODE_EXEM 1
#define DIS_MODE_MATC 2
#define ede_a   0.595639
#define ede_c   0.457748

int main(int argc , char *argv[]){

    int seq_len;
    bool is_CC_separation = false;
    bool is_dis           = false;
    bool is_dcj           = false;
    bool is_med           = false;
    bool is_knap          = false;
    bool is_mc            = false;
    bool is_cal_bij       = false;
    int heu_level         = 2;
    int term_move         = 3;
    bool is_opt           = false;
    char *input_file      = NULL;
    char *output_file     = NULL;
    char *output_dir      = NULL;
    char *bij_file        = NULL;
    char *dict_file       = NULL;
    char *tmp_folder      = NULL;
    int dis_mode          = 1;
    int num_t             = 1;
    int c                 = -1;
    char *opt_file        = NULL;
    char *logger          = NULL;

    while (1)
    {
        static struct option long_options[] =
        {   
            /* These options set a flag. */
            {"verbose",     no_argument,       &verbose_flag, 1}, 
            {"brief",       no_argument,       &verbose_flag, 0}, 
            /* These options don't set a flag.
             * We distinguish them by their indices. */
            {"dis",         no_argument,             0, 'D'},
            {"dcj",         no_argument,             0, 'd'},
            {"median",      no_argument,             0, 'm'},
            {"knap",        no_argument,             0, 'k'},
            {"cal_bij",     no_argument,             0, 'j'},
            {"mc",          no_argument,             0, 'C'},
            {"CC",          no_argument,             0, 'c'},
            {"dis_mode",	required_argument, 	     0, 'M'},
            {"num_t",       required_argument,       0, 't'},
            {"input_file",  required_argument,       0, 'f'},
            {"output_file", required_argument,       0, 'F'},
            {"bij_file",    required_argument,       0, 'J'},
            {"dict_file",   required_argument,       0, 'T'},
            {"output_dir",  required_argument,       0, 'r'},
            {"p_mode",      required_argument,       0, 'p'},
            {"heu_level",   required_argument,       0, 'l'},
            {"term_move",   required_argument,       0, 'v'},
            {"is_opt",      required_argument,       0, 'o'},
            {"opt_file",    required_argument,       0, 'O'},
            {"seq_len",     required_argument,       0, 'L'},
            {"tmp_folder",  required_argument,       0, '1'},
            {"logger",      required_argument,       0, '2'},
            {"help",        no_argument,             0, 'h'},
            {0,             0,                       0,  0 }
        }; 
        int option_index = 0;
        c = getopt_long (argc, argv, "df:hjJ:kl:mo:p:t:T:v:1:", 
                long_options, &option_index);
        if (c == -1)
            break;	
        switch (c)
        {
            case 'r':
                {
                    output_dir = optarg;
                    break;
                }
            case 'D':
                {
                    is_dis = true;
                    break;
                }
            case 'L':
                {
                    seq_len = atoi(optarg);
                    break;
                }
            case 'M':
                {
                    dis_mode = atoi(optarg);
                    break;
                }
            case 'd':
                {
                    is_dcj = true;
                    break;
                }
            case 'f':
                {
                    input_file = optarg;
                    break;
                }
            case 'j':
                {
                    is_cal_bij = true;
                    break;
                }
            case 'J':
                {
                    bij_file = optarg;
                    break;
                }
            case 'F':
                {
                    output_file = optarg;
                    break;
                }
            case 'O':
                {
                    opt_file = optarg;
                    break;
                }
            case 'h':
                {
                    print_help();
                    return 1;
                }
            case 'k':
                {
                    is_knap = true;
                    break;
                }
            case 'l':
                {
                    heu_level = atoi(optarg);
                    break;
                }
            case 'm':
                {
                    is_med = true;
                    break;
                }
            case 'C':
                {
                    is_mc = true;
                    break;
                }
            case 'c':
                {
                    is_CC_separation = true;
                    break;
                }
            case 'o':
                {
                    is_opt = atoi(optarg)==0?true:false; 
                    break;
                }
            case 'p':
                {
                    parallel_mode = atoi(optarg);
                    break;
                }
            case 't':
                {
                    num_t = atoi(optarg);
                    break;
                }
            case 'T':
                {
                    dict_file = optarg;
                    break;
                }
            case 'v':
                {
                    term_move = atoi(optarg);
                    break;
                }
            case '1':
                {
                    tmp_folder = optarg;
                    break;
                }
            case '2':
                {
                    logger = optarg;
                    break;
                }
            case '?':
                {
                    print_help();
                    return 1;
                }
            default:
                {
                    abort ();
                }
        }
    }
    if(is_CC_separation == true)
    {
        Instance *ins = new InsMC(input_file, output_file, output_dir);
        printf("finished separation of Connected Components\n");
    }
    if(is_mc == true)
    {
        printf("Start initialization of instances!\n");

        Instance** ins = (Instance**)malloc(sizeof(InsMC*)*2);
        if(logger != NULL)
        {
            char lf_insmc[OPTKIT_FILE_SIZE];
            char input_f[OPTKIT_FILE_SIZE];
            //Utils::get_file_name(logger, input_f, OPTKIT_FILE_SIZE);
            snprintf(lf_insmc, OPTKIT_FILE_SIZE, "%s%s", logger, input_f);
            ins[0] = new InsMC(input_file, NULL, lf_insmc);
        }
        ins[1] = new InsMC(input_file);

        printf("finished initialization of instances!\n");

        int buck_size    = ins[0]->upper_bound - ins[0]->lower_bound + 2;
        int list_size    = 100000;
        int base         = ins[0]->lower_bound-1;
        int num_t        = 1;
        int num_elem     = ins[0]->num_count;
        bool is_ub       = true;
        bool is_enum_all = true;
        IntList *list = new IntList(buck_size, list_size, 
                                    base,      num_t, 
                                    is_ub,     num_elem, 
                                    &is_enum_all);
        exit(1);
        list->bnb(ins, 0);
    }
    if(is_dis == true)
    {
        printf("is_cal_bij %s\n", is_cal_bij?"true":"false");
        Instance** ins = (Instance**)malloc(sizeof(InsDis*)*2);
        if(dis_mode == DIS_MODE_EXEM)
        {
            for(int i=0; i<2; i++)
            {
                if(is_cal_bij)
                {
                    bool is_cal_bij = true;
                    ins[i] = new InsDis(input_file, true, NULL, dict_file, bij_file, &is_cal_bij);
                }
                else
                {
                    ins[i] = new InsDis(input_file, true);
                }
            }
        }
        else if(dis_mode == DIS_MODE_MATC)
        {
            for(int i=0; i<2; i++)
            {
                if(is_cal_bij)
                {
                    ins[i] = new InsDis(input_file, false, NULL, dict_file, bij_file, &is_cal_bij);
                }
                else
                {
                    ins[i] = new InsDis(input_file, false);
                }
            }
        }
        printf("finished initialization of instances!\n");
        if(ins[0]->upper_bound == ins[0]->lower_bound)
        {
            printf("num_extern_comp: %d\n", 
                  (ins[0]->lower_bound+ins[0]->upper_bound)/2);
            exit(1);
        }
        int buck_size = ins[0]->upper_bound - ins[0]->lower_bound + 2;
        int list_size = 100000;
        int base = ins[0]->lower_bound-1;
        int num_t = 1;
        int num_elem = ins[0]->num_count;
        /*TODO check this patch*/
        if(num_elem>1000)
        {
            num_elem = 1000;
        }
        IntList *list = new IntList(buck_size, list_size, 
                        base, num_t, true, num_elem);
        list->bnb(ins, 0);
        FILE *writer = fopen(opt_file, "a");
        fprintf(writer, "num_extern_comp: %d\n", 
                    (list->lower_bound+list->upper_bound)/2);
        printf("num_extern_comp: %d\n", 
                    (list->lower_bound+list->upper_bound)/2);
        fclose(writer);
    }
    if(is_med == true)
    {
        if(parallel_mode == MODE_PAR)
        {
            // TODO change the parallel mode
            //test sequential algorithm for dcj distance
            //test 10 times
            int dis = 10000000;
            int dis0;
            for(int j=0; j<10; j++)
            {
                printf("start parallel bucket processing!\n");
                Instance** ins = (Instance**)malloc(sizeof(InsMed*)*num_t+1);
                for(int i=0; i<num_t+1; i++)
                {
                    ins[i] = new InsMed(input_file, tmp_folder, i);
                }
                printf("initial upper_bound %d lower_bound %d!\n", 
                            ins[0]->upper_bound, ins[0]->lower_bound);
                int buck_size = term_move+2; //the bucket size is smaller
                int list_size = 100000;
                int base = 0; //based will be 0 as there is no need to control position by upper or lower bound
                int num_elem = ins[0]->num_count;
                IntList *list = new IntList(buck_size, list_size, 
                                base, num_t, true, num_elem);
                list->upper_bound = ins[0]->upper_bound;
                if(j==0)
                    ins[0]->print_encode();
                dis0= ins[0]->score;
                list->lk(ins, heu_level, term_move, is_opt);
                if(list->upper_bound< dis)
                    dis = list->upper_bound;
            }
            printf("dis0 %d dis1 %d \n", dis0, dis);
            printf("end parallel bucket processing!\n");
        }
        if(parallel_mode == MODE_SEQ)
        {
            /* test sequential algorithm for dcj distance */
            printf("start sequential bucket processing! num_t %d\n", num_t);
            Instance** ins = (Instance**)malloc(sizeof(InsMed*)*(num_t+1));
            for(int i=0; i<num_t+1; i++)
            {
                printf("init %d\n", i);
                which_med = i;
                bool uh = true;
                int th = 2;
                ins[i] = new InsMed(input_file, tmp_folder, i, NULL, &dis_mode, &uh, &th);
            }
            printf("initial upper_bound %d lower_bound %d!\n", 
                    ins[0]->upper_bound, ins[0]->lower_bound);
            /* the bucket size is smaller */
            int buck_size = term_move+2; 
            int list_size = 100000;
            /* base will be 0 as there is no need to 
               control position by upper or lower bound */
            int base = 0; 
            int num_elem = ins[0]->num_count;
            IntList *list = new IntList(buck_size, list_size, 
                            base, num_t, true, num_elem);
            list->upper_bound = ins[0]->upper_bound;
            /* compute the real score */
            int dis0 = ins[0]->score;
            printf("start computing median using lk method!\n");
            list->lk(ins, heu_level, term_move, is_opt);
            printf("the optimal solution distance is %d\n", list->upper_bound);
        }
    }
    if(is_knap == true)
    {
        if(parallel_mode == MODE_SEQ)
        {
            //test sequential algorithm
            printf("start sequential processing!\n");
            Instance** ins = (Instance**)malloc(sizeof(Instance*)*2);
            for(int i=0; i<2; i++)
                ins[i] = new InsKnapsack(input_file, 0);
            int buck_size = ins[0]->upper_bound - ins[0]->lower_bound + 2;
            int list_size = 100000;
            int base = ins[0]->lower_bound-1;
            int num_elem = ins[0]->num_count; 
            IntList *list = new IntList(buck_size, list_size, 
                                base, num_t, true, num_elem);	
            list->bnb(ins, 0);	
            //printf("time taken %f  %f %f %f %f %f \n", 
            //        list->t.list[0].time, list->t.list[1].time, 
            //        list->t.list[2].time, list->t.list[3].time, 
            //        list->t.list[4].time, list->t.list[5].time);
#ifdef USE_SPC
            printf("search space %d\n", list->search_space);
#endif
            printf("end sequential processing!\n");
        }
        else if(parallel_mode == MODE_PAR)
        {
            //test parallel algorithm
            printf("satrt parallel processing!\n");
            Instance** p_ins = (Instance**)malloc(sizeof(Instance*)*(num_t+1));
            for(int i=0; i<=num_t; i++)
                p_ins[i] = new InsKnapsack(input_file, i);
            printf("finished initialization!\n");
            int buck_size = p_ins[0]->upper_bound - p_ins[0]->lower_bound + 2;
            int list_size = 1000000;
            int base = p_ins[0]->lower_bound-1; 
            int num_elem = p_ins[0]->num_count;
            IntList *p_list = new IntList(buck_size, list_size, 
                                base, num_t, true, num_elem);
            p_list->parallel_base = 100;
            p_list->bnb_parallel_bucket(p_ins);	
            printf("time taken %f  %f %f %f %f %f \n", 
                        p_list->t.list[0].time, p_list->t.list[1].time, 
                        p_list->t.list[2].time, p_list->t.list[3].time, 
                        p_list->t.list[4].time, p_list->t.list[5].time);
#ifdef USE_SPC
            printf("search space %d\n", p_list->search_space);
#endif
#ifdef USE_BW
            long long bytes=0;
            bytes += (2*p_list->read_cnt + p_list->read_num +
                    3*p_list->write_cnt + 2*p_list->write_num);
            printf("bytes %lld\n", bytes*4);
#endif
            printf("end parallel processing!\n");
        }
        else if(parallel_mode == MODE_TPAR)
        {
            printf("start parallel thread processing!\n");
            Instance** p_ins = (Instance**)malloc(sizeof(Instance*)*(num_t+1));
            for(int i=0; i<=num_t; i++) 
                p_ins[i] = new InsKnapsack(input_file, i); 
            int buck_size = p_ins[0]->upper_bound - p_ins[0]->lower_bound + 2;
            int list_size = 1000000;
            int base = p_ins[0]->lower_bound-1;
            int num_elem = p_ins[0]->num_count; 
            PList *pl = new PList(p_ins[0]->lower_bound, 
                        p_ins[0]->upper_bound, num_t);
            pl->base = base;
            pl->num_threads = num_t;
            pl->p_lists = (List**)malloc(sizeof(List*)*num_t);
            for(int i=0; i<num_t; i++)
            {
                pl->p_lists[i] = new IntList(buck_size, 
                                    list_size, base, num_t, true, num_elem);
                pl->p_lists[i]->parallel_base = 100;
            }
            pl->bnb_parallel_threads_ub(p_ins);
#ifdef USE_SPC
            int space = 0;
            for(int i=0; i<num_t; i++)
            {
                space += pl->p_lists[i]->search_space;
            }
            printf("search space %d\n", space);
#endif
#ifdef USE_BW
            long long int bytes=0;
            for(int i=0; i<num_t; i++)
            {
                bytes += (2*pl->p_lists[i]->read_cnt + pl->p_lists[i]->read_num +
                        3*pl->p_lists[i]->write_cnt + 2*pl->p_lists[i]->write_num);
            }
            printf("bytes %lld\n", bytes*4);
#endif
            printf("end parallel thread processing!\n");
        }
    }
}

void
print_help(){
    printf("\n\nOPTKit: Optimization Tool-kit for Prallellizing Discrete Combinatoric Problems in Emerging Platforms (2013-2015)\n\n"); 
    printf("usage:\n"
           "optkit <--dis                                  |\n" 
           "        --dcj                                  |\n" 
           "        --median                               |\n" 
           "        --knap                                 |\n" 
           "        --cal_bij                              |\n" 
           "        --mc --input_file <graph_file>         |\n" 
           "        --CC                                   >\n"
           "       [--logger <log folder>]\n\n\n"
          );
    printf("--dis       is calculate dcj distance with unequal content genomes\n");
    printf("--dcj       is calculate dcj\n");
    printf("--median    is calculate median\n");
    printf("--knap      is calculate knapsack problem\n");
    printf("--cal_bij   is calculate bijection accuracy\n");
    printf("--mc        is calculate maximum clique problem\n");
    printf("--CC        is calculate connected component problem\n");
    printf("--dis_mode <1|2>            0:exemplar distance 1:matching distance\n");
    printf("--num_t <num threads>       assign number of threads\n");
    printf("--input_file <file name>    assign input file\n");
    printf("--output_file <file name>   assign output file\n");
    printf("--bij_file <file_name>      assign bijection output file\n");
    printf("--dict_file <file_name>     assign mapping dictionary file\n");
    printf("--output_dir <dir_name>     assign output directory name\n");
    printf("--logger <logger_folder>    logger folder position\n");
    printf("--p_mode <1|2|3>            1:sequential " 
            "2:parallel_bucket " 
            "3:parallel_threading\n");
    printf("--heu_level <heuristic level> "   
            "uesed in lin-kernighan algorithm\n");
    printf("--term_move <terminal move>"
            "used in lin-kernighan algorithm\n");
    printf("--is_opt <true|false>       used in lin-kernighan algorithm\n");
    printf("--opt_file <file_name>      write optimal result to file\n");
    printf("--seq_len <length of seq>   assign length of the sequence\n");
    printf("--tmp_folder <tmp folder>   folder to keep temporary files\n");
    printf("--help      print this help file\n\n");

    printf("Examples: \n");
    printf("optkit --dis --cal_bij --input_file <input_file> --dis_mode <1|2>\n"
            "       --dict_file <dict_file> --opt_file <opt_file> --p_mode <1|2|3>\n"
            "       --seq_len <seq_len> --bij_file <bij_file>\n\n" );
    printf("optkit --dis --input_file <input_file> --dis_mode <1|2>\n"
            "       --opt_file <opt_file> --p_mode <1|2|3>\n"
            "       --seq_len <seq_len>\n\n" );
    printf("optkit --median --input_file <input_file> --dis_mode <1|2>\n"
            "      --p_mode <1|2|3> --seq_len <seq_len> --tmp_folder <tmp_folder> \n\n" );
}

int 
read_opt(char *opt_file){
    FILE *reader = fopen(opt_file, "r");
    char str[1000];
    float num;
    float result = 0;
    fscanf(reader, "%s %f\n", &str, &num);
    result += num;
    fscanf(reader, "%s %f\n", &str, &num);
    result += num;
    fclose(reader);
    return (int)result;
}

int 
extimate_ede(int dis){
    double d = (double)dis;
    double ede_dis = (pow(d,2)+0.5956*d)/(pow(d,2)+0.4577*d+0.5956);
    return (int)ede_dis;
}

    int
ede1 ( int invdist, int ngene )
{
    double ll, tt, kk, pp, dval;
    int newvalue;

    kk = invdist / ( ngene + 0.0 );

    if ( kk >= 0.999999999999 )
    {                           /* the distance correction has singularity at 1 */
        kk = 0.999999999999;
    }
    if ( kk <= 1 - ede_c )
        return invdist;

    ll = ede_c * kk - ede_a;
    tt = 4 * ede_a * ( 1 - kk ) * kk + ll * ll;
    tt = ll + sqrt ( tt );
    pp = tt / ( 2 * ( 1 - kk ) );
    pp *= ngene;

    dval = pp;
    newvalue = ( int ) ceil ( dval );
    /*if (newvalue-dval > 0) return newvalue-1; */
    return newvalue;
}
