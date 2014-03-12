#ifndef _ERROR_HEADER_H
#define _ERROR_HEADER_H

#include <stdio.h>
#include <stdlib.h>

#define ARGUMENTS_ERROR 1
#define FILE_ERROR 2
#define FILE_SIZE_ERROR 3
#define LEX_1 4
#define LEX_2 5
#define LEX_3 6
#define LEX_4 7
#define LEX_5 8
#define LEX_6 9
#define SYNTAX_1 10
#define SYNTAX_2 11
#define SYNTAX_3 12
#define MISSING_GROUND 13
#define PLOT_1 14

void error( int errcode, int line_num);




#endif
