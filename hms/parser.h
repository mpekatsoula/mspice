#ifndef _PARSER_H
#define _PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdint.h>

void parser(char *filename);
void lexical_analysis( char *line );
int line_num;
unsigned int nodes_counter;
unsigned int voltage_elements_counter;
unsigned int i_elements_counter;
unsigned int resistances;

#define OPTIONS_CASE 0xBEEF
#define DC_CASE 0xBEF0
#define PLOT_CASE 0xBEF1
#define TRAN_CASE 0xBEF2
#define MAX_PLOT_VALUES 256
#define MAX_PWL_PAIRS 8

#define _EXP 0x30
#define _SIN 0x31
#define _PULSE 0x32
#define _PWL 0x33

typedef struct _options {

  int OPTIONS;
  int SPD;
  int DC;
  int PLOT;
  int TRAN;
  int ITER;
  int SPARSE;
  int TR;
  int BE;
  int num_of_dc_plots;
  int num_of_tran_plots;
  char *PLOT_scan[MAX_PLOT_VALUES];
  char *PLOT_tran[MAX_PLOT_VALUES];
  char *DC_source_name;
  char *DC_type;
  float ITOL;
  float DC_start;
  float DC_end;
  float DC_step;
  float TRAN_step;
  float TRAN_fin_time;

} options;

typedef struct _element {

  char *type;
  char *name;
  int node_1;
  int node_2;
  int node_3;
  int node_4;
  double value;
  int area;
  char *model_name;
  double L, W;
  struct _element *next;
  int EXP;
  int SIN;
  int PULSE;
  int PWL;
  double td1, tc1, td2, tc2;
  double ia, fr, td, df, ph, tr, tf, pw, per;
  double t1, i1, i2;
  int n;
  double pwl_t[MAX_PWL_PAIRS], pwl_i[MAX_PWL_PAIRS];
} element;

/* Flag to check if node 0 is present */
int ground_flag;

/* global flag */
int flag;


element **b_info;

/* Global element list that contains circuit elements */
element *Element_list;

/* Last list pointer to optimize searching */
element *last_pointer;

/* Global options struct */
options *NetOptions;

#endif
