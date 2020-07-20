#pragma once

#ifndef _H_UTILS
#define _H_UTILS
#include <cstdio>
#include <cstdlib>
#include <cstring>


#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

#define red_p(cl, ...)   printf(RED); \
                printf(cl, ##__VA_ARGS__);\
                printf(RESET "\n")
#define blue_p(cl, ...)   printf(BLUE  ##__VA_ARGS__  RESET "\n")
#define ERROR_PRINT() {printf("Error on (%d) line in (%s) file\n", __LINE__, __FILE__); exit(123);}

#define OPTKIT_INVALID -1
/**
 *  @enum       GGPARAM_COMP
 *
 *  @brief      Component type that a parameter is applicable to.
 *
 *  @details    This is a one-2-one mapping of type tcomp in ggparam.xsd.
 *              It enumerates all OGG applications. GGPARAM_COMP and tcomp
 *              must stay in sync with actual applications in the OGG
 *              installation package. Each enum value is defined in a way
 *              so that they can used in a bit map, since parameters can
 *              to more than one component.
**/
typedef enum
{
    OPTKIT_NULL = -1,
    OPTKIT_ZERO = 0, 
    OPTKIT_FILE_SIZE = 512 
}
OPTKIT_DEF;

using namespace std;

class Utils
{
    public:
        static char *readLine(FILE *file);
        static int  str_split(char* a_str, const char a_delim, char **result);
        static void q_sort_two(int *key, int *val, int left, int right);
        static void q_sort(int *number, int left, int right);
        static void get_file_name(const char *full_path_name, char *file_name, int16_t size);
        static void get_file_name(const char *full_path_name, int16_t size1, char *file_name, int16_t size2);
        static void concate_path(const char *path, const char *full_path_name, int16_t size2, const char* concated, int size3);
};

#endif
