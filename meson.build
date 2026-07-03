project(
  'orGUI',
  'cpp',
  version: run_command(
    find_program('python'),
    '-m',
    'setuptools_scm',
    check: true
  ).stdout().strip(),
  default_options: ['cpp_std=c++17']
)

py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11', required: true)

if get_option('native_optimization')
  add_project_arguments('-march=native', language: 'cpp')
endif

# CUDA support is intentionally disabled for now. If we add GPU kernels later,
# re-enable the CUDA language/dependency here, following the debye build.
# cuda_dep = dependency('cuda', required: true)
# add_project_arguments(
#   ['-arch=compute_75', '-code=sm_75', '-lineinfo'],
#   language: 'cuda'
# )

py.extension_module(
  '_CTRcalc_cpp',
  files('orgui/datautils/xrayutils/cpp/CTRcalc_cpp.cpp'),
  dependencies: [pybind11_dep],
  install: true,
  subdir: 'orgui/datautils/xrayutils'
)

py.extension_module(
  '_roi_sum_cpp',
  files('orgui/app/cpp/roi_sum_cpp.cpp'),
  dependencies: [pybind11_dep],
  install: true,
  subdir: 'orgui/app'
)

install_subdir(
  'orgui',
  install_dir: py.get_install_dir(),
  exclude_directories: ['__pycache__']
)
