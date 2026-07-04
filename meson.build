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
cpp = meson.get_compiler('cpp')
native_args = []

if get_option('native_optimization')
  native_args = cpp.get_supported_arguments('-march=native')

  if native_args.length() == 0 and cpp.get_id() == 'msvc'
    if meson.can_run_host_binaries()
      arch_probe = cpp.run(
        files('buildtools/detect_msvc_arch.cpp')[0],
        name: 'MSVC native /arch detection'
      )
      if arch_probe.compiled() and arch_probe.returncode() == 0
        foreach detected_arch : arch_probe.stdout().strip().split()
          supported_arch_args = cpp.get_supported_arguments('/arch:' + detected_arch)
          if supported_arch_args.length() > 0
            native_args = supported_arch_args
            break
          endif
        endforeach
      else
        warning('MSVC native CPU probe failed; building without native CPU flags')
      endif
    else
      warning('native_optimization requested for MSVC, but host binaries cannot be run during setup')
    endif
  endif

  if native_args.length() == 0
    warning('native_optimization requested, but this C++ compiler does not support a known native CPU flag')
  else
    add_project_arguments(native_args, language: 'cpp')
  endif
endif

build_config = configuration_data()
target_cpu = 'generic'
allowed_instruction_set = []

gcc_instruction_options = [
  '-mabm',
  '-madx',
  '-maes',
  '-mavx',
  '-mavx2',
  '-mavx5124fmaps',
  '-mavx5124vnniw',
  '-mavx512bf16',
  '-mavx512bitalg',
  '-mavx512bw',
  '-mavx512cd',
  '-mavx512dq',
  '-mavx512er',
  '-mavx512f',
  '-mavx512ifma',
  '-mavx512pf',
  '-mavx512vbmi',
  '-mavx512vbmi2',
  '-mavx512vl',
  '-mavx512vnni',
  '-mavx512vp2intersect',
  '-mavx512vpopcntdq',
  '-mavxvnni',
  '-mbmi',
  '-mbmi2',
  '-mcrc32',
  '-mf16c',
  '-mfma',
  '-mfma4',
  '-mfxsr',
  '-mmmx',
  '-mpclmul',
  '-mpopcnt',
  '-mprfchw',
  '-mrdpid',
  '-mrdrnd',
  '-mrdseed',
  '-msha',
  '-msse',
  '-msse2',
  '-msse2avx',
  '-msse3',
  '-msse4',
  '-msse4.1',
  '-msse4.2',
  '-msse4a',
  '-mssse3',
  '-mxop',
  '-mxsave',
  '-mxsavec',
  '-mxsaveopt',
  '-mxsaves',
]

gcc_target_help_compiler = find_program('c++', 'g++', required: false)
if cpp.get_id() == 'gcc' and gcc_target_help_compiler.found()
  target_help = run_command(
    gcc_target_help_compiler,
    native_args,
    '-Q',
    '--help=target',
    check: false
  )
  if target_help.returncode() == 0
    foreach line : target_help.stdout().split('\n')
      fields = line.strip().split()
      if fields.length() >= 2
        option = fields[0]
        value = fields[1]
        if option == '-march='
          target_cpu = value
        elif value == '[enabled]' and gcc_instruction_options.contains(option)
          allowed_instruction_set += [option.substring(2)]
        endif
      endif
    endforeach
  else
    warning('Could not query GCC target options with -Q --help=target')
  endif
endif

foreach native_arg : native_args
  if native_arg == '-march=native'
    if target_cpu == 'generic'
      target_cpu = 'native'
    endif
    if allowed_instruction_set.length() == 0
      if cpp.get_define('__AVX512F__', args: native_args) != ''
        allowed_instruction_set = ['avx512f']
      elif cpp.get_define('__AVX2__', args: native_args) != ''
        allowed_instruction_set = ['avx2']
      elif cpp.get_define('__AVX__', args: native_args) != ''
        allowed_instruction_set = ['avx']
      elif cpp.get_define('__SSE4_2__', args: native_args) != ''
        allowed_instruction_set = ['sse4.2']
      endif
    endif
  elif native_arg.startswith('/arch:')
    target_cpu = host_machine.cpu()
    allowed_instruction_set = [native_arg.split(':')[1]]
  endif
endforeach

allowed_instruction_set_literal = '['
foreach instruction : allowed_instruction_set
  if allowed_instruction_set_literal != '['
    allowed_instruction_set_literal += ', '
  endif
  allowed_instruction_set_literal += '"' + instruction + '"'
endforeach
allowed_instruction_set_literal += ']'

build_config.set(
  'NATIVE_OPTIMIZATION',
  get_option('native_optimization') ? 'True' : 'False'
)
build_config.set_quoted('TARGET_CPU', target_cpu)
build_config.set('ALLOWED_INSTRUCTION_SET', allowed_instruction_set_literal)
build_config.set_quoted('NATIVE_ARGUMENTS', ' '.join(native_args))
build_config.set_quoted('CPP_COMPILER_ID', cpp.get_id())
build_config.set_quoted('CPP_COMPILER_VERSION', cpp.version())
build_config.set_quoted('HOST_SYSTEM', host_machine.system())
build_config.set_quoted('HOST_CPU_FAMILY', host_machine.cpu_family())
build_config.set_quoted('HOST_CPU', host_machine.cpu())

configure_file(
  input: 'orgui/_build_config.py.in',
  output: '_build_config.py',
  configuration: build_config,
  install: true,
  install_dir: py.get_install_dir() / 'orgui'
)

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
  exclude_directories: ['__pycache__'],
  exclude_files: ['_build_config.py.in']
)
