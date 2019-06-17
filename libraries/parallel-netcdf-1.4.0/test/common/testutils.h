/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: testutils.h 1468 2013-10-26 16:53:18Z wkliao $ */


#ifndef _UTILS_H
#define _UTILS_H

#include <limits.h>

typedef struct {
	char infname[PATH_MAX];
	char outfname[PATH_MAX];
} params;

void parse_read_args(int argc, char **argv, int rank, params *p);
void parse_write_args(int argc, char **argv, int rank, params *p);
#endif
