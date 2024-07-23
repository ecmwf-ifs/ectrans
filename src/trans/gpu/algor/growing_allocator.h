#pragma once

extern "C" void growing_allocator_register_free_c(void *,
                                                  void (&)(float *, size_t));
