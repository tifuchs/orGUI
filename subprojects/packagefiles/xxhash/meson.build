project(
  'xxHash',
  'c',
  version: '0.8.2',
  license: 'BSD-2-Clause',
)

xxhash = static_library(
  'xxhash',
  'xxhash.c',
  install: false,
)

xxhash_dep = declare_dependency(
  link_with: xxhash,
  include_directories: include_directories('.'),
)

meson.override_dependency('libxxhash', xxhash_dep)
