/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "ectrans/transi.h"

#include "transi_test.h"


int main ( int arc, char **argv )
{
  int start = allocated();

  print_mem("Initially allocated: ",allocated()-start);

  if( strcmp(ectrans_version(),"transi/contrib") == 0 ) // equals
    trans_set_handles_limit(2);
  else /// Older versions cannot reuse existing handles, so allocate enough.
    trans_set_handles_limit(3);

  trans_use_mpi(false);

  trans_init();

  int iter=0;
  int iter_max=50;
  int mem_leak=1;
  //int start_loop = allocated();
  for( iter=0; iter<iter_max; ++iter )
  {
    printf("iteration %d\n",iter+1);
    int start_iter = allocated();
    struct Trans_t trans;

    trans_new(&trans);
    trans_set_resol_lonlat(&trans,320,161);
    trans_set_trunc(&trans,159);

    trans_setup(&trans);

    print_mem("Allocated in iteration:       ",allocated()-start_iter);

    trans_delete(&trans);

    int allocated_now = allocated();
    print_mem("Possibly leaked in iteration: ",allocated_now-start_iter);
    //print_mem("total leaked in loop:   ",allocated()-start_loop);

    // No memory leaks in subsequent iterations
    if( allocated_now - start_iter == 0 ) {
      mem_leak=0;
      break;
    }
  }
  if( mem_leak != 0 ) {
        printf("ERROR: Memory leaks present between trans_setup() and trans_delete(), tested after %d attempts\n",iter);
  }
  printf( "trans_finalize()\n" );
  trans_finalize();

  print_mem("Total leaked: ",allocated()-start);


  display_mallinfo();

  return mem_leak;
}

