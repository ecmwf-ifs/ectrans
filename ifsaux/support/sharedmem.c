
#include <stdlib.h>

void sharedmem_malloc_bytes (void** ptr, size_t bytes)
{ 
  *ptr = malloc(bytes);
}

void sharedmem_free(void** ptr)
{ 
  free(*ptr);
}

void sharedmem_advance_bytes  (void** ptr, size_t bytes)
{ 
  *ptr += bytes; 
}
