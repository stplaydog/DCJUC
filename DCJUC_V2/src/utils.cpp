#include "utils.h"

/**
 * read one line from a file
 *
 * @param[in]       file            the file pointer
 *
 * @return      a line of string
 *
**/
char* Utils::readLine(FILE *file) 
{
	if (file == NULL) 
    {
		printf("Error: file pointer is null.");
		exit(1);
	}

	int maximumLineLength = 128;
	char *lineBuffer = new char[maximumLineLength];
	if (lineBuffer == NULL) 
    {
		printf("Error allocating memory for line buffer.");
		exit(1);
	}

	char ch = getc(file);
	int count = 0;

	while ((ch != '\n') && (ch != EOF)) {
		if (count == maximumLineLength) {
			maximumLineLength += 128;
			lineBuffer = (char*)realloc(lineBuffer, maximumLineLength);
			if (lineBuffer == NULL) {
				printf("Error reallocating space for line buffer.");
				exit(1);
			}
		}
		lineBuffer[count] = ch;
		count++;

		ch = getc(file);
	}

	lineBuffer[count] = '\0';
	char line[count + 1];
	for(int i=0; i<count+1; i++)
		line[i] = lineBuffer[i];
	free(lineBuffer);
	char *constLine = line;
	return constLine;
}

int  
Utils::str_split(char* a_str, const char a_delim, char** result)
{
	//char** result    = 0;
	size_t count     = 0;
	char* tmp        = a_str;
	char* last_comma = 0;
	char delim[2];
	delim[0] = a_delim;
	delim[1] = 0;

	/* Count how many elements will be extracted. */
	while (*tmp)
	{
		if (a_delim == *tmp)
		{
			count++;
			last_comma = tmp;
		}
		tmp++;
	}

	/* Add space for trailing token. */
	count += last_comma < (a_str + strlen(a_str) - 1);

	/* Add space for terminating null string so caller
	 *        knows where the list of returned strings ends. */
	count++;

	result = (char**)malloc(sizeof(char*) * count);

	if (result)
	{
		size_t idx  = 0;
		char* token = strtok(a_str, delim);

		while (token)
		{
			*(result + idx++) = strdup(token);
			token = strtok(0, delim);
		}
		*(result + idx) = 0;
	}

	return count;
}

	void 
Utils::q_sort_two(int *key, int *val, int left, int right)
{
	int pivot, l_hold, r_hold, val_pivot;
	l_hold = left;
	r_hold = right;
	pivot = key[left];
	val_pivot = val[left];
	while (left < right)
	{
		while ((key[right] >= pivot) && (left < right))
			right--;
		if (left != right)
		{
			key[left] = key[right];
			val[left] = val[right];
			left++;
		}
		while ((key[left] <= pivot) && (left < right))
			left++;
		if (left != right)
		{
			key[right] = key[left];
			val[right] = val[left];
			right--;
		}
	}
	key[left] = pivot;
	val[left] = val_pivot;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort_two(key, val, left, pivot-1);
	if (right > pivot)
		q_sort_two(key, val, pivot+1, right);
}

/**
 * sort only one array
**/
void 
Utils::q_sort(int *numbers, int left, int right)
{
	int pivot, l_hold, r_hold;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right))
        {
			right--;
        }

		if (left != right)
		{
			numbers[left] = numbers[right];
			left++;
		}

		while ((numbers[left] <= pivot) && (left < right))
        {
			left++;
        }

		if (left != right)
		{
			numbers[right] = numbers[left];
			right--;
		}
	}

	numbers[left] = pivot;
	pivot = left;
	left = l_hold;
	right = r_hold;

	if (left < pivot)
    {
		q_sort(numbers, left, pivot-1);
    }

	if (right > pivot)
    {
		q_sort(numbers, pivot+1, right);
    }
}

/**
 * @brief   From a file name (not) contains full path,
 *          find the real name.
 *
 * @param[in]       full_path_name  The file with(out) full path.
 * @param[in]       size1           The full_path_name buffer size.
 * @param[out]      file_name       The file deprived of full path.
 * @param[in]       size2           The file_name buffer size.
 *
 * @return      File name deprived of full path.
 *
 * TODO     Treat windows and linux diffrerently.
**/
void Utils::get_file_name(const char *full_path_name, int16_t size1, char *file_name, int16_t size2)
{
    //assert(size1 <= size2);
    //assert(size2 <= OPTKIT_FILE_SIZE);

    //char    ret[OPTKIT_FILE_SIZE];
    //int16_t last_slash_pos = OPTKIT_ZERO;
    //size_t  full_path_size = strlen(full_path_name);

    //for(int16_t i=full_path_size-1; i>=0; i--)
    //{
    //    if(full_path_name[i] == '/')
    //    {
    //        last_slash_pos = i;
    //    }
    //}

    //snprintf(ret, size, "%s", full_path_name+last_slash_pos);
}

/**
 * @brief   Given a path and a full path name, concate this path,
 *          and the real file name, to form a new full path name.
 *
 * @param[in]       full_path_name  The file with(out) full path.
 * @param[in]       size1           The full_path_name buffer size.
 * @param[out]      file_name       The file deprived of full path.
 * @param[in]       size2           The file_name buffer size.
 * @param[out]      file_name       The file deprived of full path.
 * @param[in]       size2           The file_name buffer size.
 *
 * @return      N/A
**/
void Utils::concate_path(const char *path, const char *full_path_name, int16_t size2, const char* concated, int size3)
{
    //char file_name[OPTKIT_FILE_SIZE];
    //Utils::get_file_name(full_path_name, size2, file_name, OPTKIT_FILE_SIZE);
    //
    //assert(strlen(file_name)+strlen(path) < size3);
    
}
