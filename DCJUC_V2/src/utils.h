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

using namespace std;

class Utils
{
public:
	static char *readLine(FILE *file);
	static int str_split(char* a_str, const char a_delim, char **result);
	static void q_sort_two(int *key, int *val, int left, int right);
	static void q_sort(int *number, int left, int right);
};

#endif
